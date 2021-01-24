import math
TEMPO = 160
L = 0.65
NX = [65, 92, 119, 146, 173, 200] # aumenta la frecuencia y resolucion en fft
FREQUENCY_FACTOR = [1668000, 592000, 210000, 44600, 13100, 4200]
E = 195000000000
DAMPED_SAMPLES_LOCALLY = 600
STRING_DIAMETER = [0.000254, 0.0003302, 0.0004318, 0.0006604, 0.0009144, 0.0011684]
I = [(math.pi/4)*(x/2)**4 for x in STRING_DIAMETER]
A = [math.pi * (x / 2) ** 2 for x in STRING_DIAMETER]
STRING_MASS = [0.000401 * L, 0.000708 * L, 0.001140 * L, 0.002333 * L, 0.004466 * L, 0.00679 * L]
AMPLITUDE = [1, 1, 1, 0.8, 0.7, 0.4]
RO_0 = [1.5, 0.5, 0.1, 0.1, 0.1, 0.1]
LOCAL_RO_1 = [20, 20, 10, 10, 10, 10]
THUMB_COMPRESSION = [0.025, 0.05, 0.05, 0.05, 0.05, 0.05]
RO_1 =  [0.005, 0.005, 0.005, 0.003, 0.002, 0.002]
######################## nickel strings
# RO_0 = [1.5, 1.5, 0.1, 0.1, 0.1, 0.1]
# LOCAL_RO_1 = [60, 60, 40, 20, 10, 10]
# THUMB_COMPRESSION = [0.2, 0.2, 0.2, 0.2, 0.2, 0.5]
# RO_1 = [0.0008, 0.0005, 0.0008, 0.0004, 0.0004, 0.0005]# depende de la frecuecia
STRING_NAMES = ['E1', 'B', 'G', 'D', 'A', 'E2']
RANDOM_DEVIATION = {
    'elongation': [
            [0.002, 0.0025],
            [0.0026, 0.0035], # good
            [0.0036, 0.004]
        ],
    'pluck_point': [
            [0.33, 0.41],
            [0.42, 0.49],
            [0.5, 0.57]
        ],
    'human_phase': [
            [-600, -200],
            [-201, 200], # good
            [201, 600]
        ], # para 48k debe ser menos de 20ms
}