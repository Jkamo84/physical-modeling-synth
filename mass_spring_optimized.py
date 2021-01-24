import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from initial_position_string import initial_position_string
from scipy.io.wavfile import write, read
from data import SEQUENCE
from datetime import datetime
from numba import njit, jit
from concurrent.futures import ThreadPoolExecutor
from model_settings import *
import roundrobin
from random import uniform, randrange
# sistema masa resorte
# completamente funcional ecuacion de onda en funcion de k
# u[i, j+1] = 2 * u[i, j] - u[i, j-1] + ((k * dt ** 2) / mass) * (u[i-1, j] - 2 * u[i, j] + u[i+1, j]) \
#     - 2 * ro_0 * dt * (u[i, j] - u[i, j-1]) \
#     + (ro_1 * dt / (h ** 2)) * (u[i+1,j] - 2*u[i,j] + u[i-1,j] - u[i+1,j-1] + 2*u[i,j-1] - u[i-1,j-1]) \
#     - (E * I * (dt ** 2) / (h ** 4)) * (u[i+2,j] - 4*u[i+1,j] + 6*u[i,j] - 4*u[i-1,j] + u[i-2,j])
#
# aleatoriedad:
# - fuerza de elongación
# - ubicación de punto de toque
# - filtro por amortiguamiento proporcional a punto de toque y fuerza
# - amortiguamiento normal segun fuerza
# - buenas y malas interpretaciones

class MassSpringSynth:
    def __init__(self, randomness=False, show_plots=True, ir_name='guitar_body', filename='nonlinear_string_song'):
        '''
        obtaining information from all strings
        '''
        self.randomness = randomness
        self.show_plots = show_plots
        self.ir_name = ir_name
        self.filename = filename
        self.string_length = L
        self.I = I
        self.string_masses = STRING_MASS
        self.amplitudes = AMPLITUDE
        self.ro_0 = RO_0
        self.ro_1 = RO_1
        self.local_ro_1 = LOCAL_RO_1
        self.thumb_decay = THUMB_COMPRESSION
        self.string_names = STRING_NAMES
        self.tempo = TEMPO
        self.young_modulus = E
        self.local_damping = DAMPED_SAMPLES_LOCALLY
        self.nx = NX
        self.frequency_factor = FREQUENCY_FACTOR

    def run(self):
        start = datetime.now().time()
        self.read_ir()
        self.generate_song()
        self.convolve_song()
        self.write_song()
        if self.show_plots:
            self.plot_graphs()
        print('---'*40)
        print(start, '--->', datetime.now().time())
        print('---'*40)

    def read_ir(self):
        self.fs, data_read = read(self.ir_name + '.wav')
        self.data_read = data_read / np.max(np.abs(data_read))

    def convolve_song(self):
        convolved = np.convolve(self.song, self.data_read[:,1], mode='full')#[:,1]
        self.convolved = convolved / np.max(np.abs(convolved))
        
    def write_song(self):
        write(self.filename + '.wav', self.fs, np.int16(self.convolved * 22000))
        write(self.filename + '_no_body.wav', self.fs, np.int16(self.song * 22000))

    def get_empty_song_matrix(self):
        total_duration = 0
        for note in SEQUENCE:
            if not note[4]:
                total_duration += note[2] * 60 / self.tempo
        total_duration += 2 * SEQUENCE[-1][3]
        cancion = np.zeros([6, int(total_duration * self.fs)])

        return cancion
    
    def get_random_parameters(self):
        pass

    def generate_song(self):

        cancion = self.get_empty_song_matrix()
        total_duration = 0
        length_sequence = len(SEQUENCE)
        #print([(0,int(length_sequence * 0.2)),(1,int(length_sequence * 0.6)),(2,int(length_sequence * 0.2))])

        robin_smooth_e = roundrobin.smooth([(0,int(length_sequence * 0.2)),(1,int(length_sequence * 0.6)),(2,int(length_sequence * 0.2))])
        performance_elongation = [robin_smooth_e() for _ in range(length_sequence)]

        robin_smooth_pp = roundrobin.smooth([(0,int(length_sequence * 0.1)),(1,int(length_sequence * 0.8)),(2,int(length_sequence * 0.1))])
        performance_pluck = [robin_smooth_pp() for _ in range(length_sequence)]

        robin_smooth_ph = roundrobin.smooth([(0,int(length_sequence * 0.25)),(1,int(length_sequence * 0.5)),(2,int(length_sequence * 0.25))])
        performance_phase = [robin_smooth_ph() for _ in range(length_sequence)]
        
        print('performance_elongation', performance_elongation)
        print('performance_pluck', performance_pluck)
        print('performance_phase', performance_phase)
        for index, note in enumerate(SEQUENCE):
            if self.randomness:
                elongation = uniform(*RANDOM_DEVIATION['elongation'][performance_elongation[index]])
                pluck_point = uniform(*RANDOM_DEVIATION['pluck_point'][performance_pluck[index]])
                human_phase = randrange(*RANDOM_DEVIATION['human_phase'][performance_phase[index]])
                print(elongation, pluck_point, human_phase)
            else:
                elongation = 0.003
                pluck_point = 0.43
                human_phase = 0
            string = note[0]
            print('------------------------- creating string ' + self.string_names[string])
            f_proportion = (2**(1.9*note[1]/12))
            f_factor_fret = self.frequency_factor[string]*f_proportion # factor que relaciona Young y longitud de la cuerda
            T = note[3]
            nt = int(self.fs * T)
            dt = T / nt
            nx_f = int(self.nx[string] / (2**(note[1]/12)))
            U = np.zeros([1, nt])[0]
            t = np.linspace(0, T, num=nt)
            x = np.linspace(0,  L, num=nx_f)
            masses = np.zeros([1, nx_f, 2])[0]
            self.initial_position = initial_position_string(self.string_length, nx_f-3, pe=elongation, pp=pluck_point)
            k = nx_f * self.young_modulus * self.I[string] * f_factor_fret/ self.string_length
            masses[1:-1, :] = self.initial_position
            masses[:, 0] = masses[:, 0] + x
            
            u_old = masses
            u = masses * 0.99998 # mejora distribución - suena resorte real
            #u[:,0] = u[:,0] * 0.99998
            u_new = masses * 0
            u_new[0] = masses[0]
            u_new[1] = masses[1]
            u_new[-1] = masses[-1]
            u_new[-2] = masses[-2]
            U[0] = np.linalg.norm(u[2, 0])
            U[1] = np.linalg.norm(u[2, 0])
            #linear density m/L
            
            simultaneous = note[4]
            print('nota ' + str(index + 1) + ' de ' + str(length_sequence) )
            if not simultaneous:
                total_duration += note[2]

            signal = generate_signal(nx_f, nt, k,
                                        self.string_masses[string]/(nx_f+1), 
                                        self.ro_0[string], 
                                        f_proportion * self.ro_1[string], 
                                        self.I[string], 
                                        self.amplitudes[string], 
                                        masses,
                                        u_old, u, u_new,
                                        U,
                                        self.local_ro_1[string],
                                        self.thumb_decay[string], T, self.fs, dt, t,
                                        self.local_damping,
                                        self.young_modulus,
                                        self.string_length)

            start = int((total_duration - note[2]) * self.fs * 60 / self.tempo)
            end = start + len(signal)
            human_phase += 960
            cancion[note[0],start + human_phase:end + human_phase] = signal
            pluck_point -= 0.04
        
        self.song = cancion.sum(axis=0)

    def plot_graphs(self):
        fig, axs = plt.subplots(4,2)

        axs[0,0].plot(self.song)
        axs[0,0].grid(True)

        X = np.abs(np.fft.fft(self.song))[:int(len(self.song)/2)]
        f = np.linspace(0, self.fs/2, num=len(X))
        axs[0,1].semilogx(f, X/np.max(X))
        axs[0,1].grid(True)

        axs[1,0].plot(self.data_read[:,1])
        axs[1,0].grid(True)

        X = np.abs(np.fft.fft(self.data_read[:,1]))[:int(len(self.data_read[:,1])/2)]
        f = np.linspace(0, self.fs/2, num=len(X))
        axs[1,1].semilogx(f, X/np.max(X))
        axs[1,1].grid(True)

        axs[2,0].plot(self.convolved)
        axs[2,0].grid(True)

        X = np.abs(np.fft.fft(self.convolved))[:int(len(self.convolved)/2)]
        f = np.linspace(0, self.fs/2, num=len(X))
        axs[2,1].semilogx(f, X/np.max(X))
        axs[2,1].grid(True)

        axs[3,0].plot(self.initial_position[:,0])
        axs[3,0].grid(True)

        axs[3,1].plot(self.initial_position[:,1])
        axs[3,1].grid(True)

        fig.tight_layout()
        plt.show()


@njit(nogil=True, fastmath=True)
def generate_signal(nx, nt, k, mass, ro_0, ro_1, I, amplitude, masses, u_old, u, u_new, U, local_ro, comp, T, fs, dt, t, local_limit, E, L):  
    h = L / nx
    ro = ro_1 * local_ro
    
    for j in range(1, nt - 1):
        if j == local_limit:
            ro = ro_1
        for i in range(2, nx - 2):
            alpha_minus = (k / mass) * (1 - (h / np.linalg.norm(u[i] - u[i-1])))
            alpha = (k / mass) * (1 - (h / np.linalg.norm(u[i+1] - u[i])))
            u_new[i] = (2*u[i] + (dt**2)*(alpha_minus*u[i-1] - (alpha_minus + alpha)*u[i] + alpha*u[i+1]) \
                - (1 - ro_0*dt)*u_old[i] \
                + (ro * dt / (h**2)) * (u[i+1] - 2*u[i] + u[i-1] - u_old[i+1] + 2*u_old[i] - u_old[i-1]) \
                - (E*I * (dt**2)/(h**4)) * (u[i+2] - 4*u[i+1] + 6*u[i] - 4*u[i-1] + u[i-2])) / (1 + ro_0*dt)

        U[j+1] = np.linalg.norm(u_new[2])
        for i in range(2, nx - 2):
            u_old[i] = u[i]
            u[i] = u_new[i]
            
    print('done')
    signal = U
    F = np.diff(signal) / np.diff(t)
    data = mass * np.diff(F) / np.diff(t[0:-1])
    signal = data
    signal[:local_limit] = comp * signal[:local_limit]

    return signal


if __name__ == "__main__":
    synth = MassSpringSynth(filename='definitive_normal')
    synth.run()
    synth = MassSpringSynth(randomness=True, filename='definitive_random')
    synth.run()