import numpy as np
import math

def initial_position_string(b, nx, pe=0.005, pp=0.5):
    pluck_elongation = pe
    pluck_point = pp
    m1 = pluck_elongation / pluck_point
    m2 = pluck_elongation / (b - pluck_point)
    xa = np.linspace(0, b, num=nx + 1)
    f = np.zeros([1, len(xa), 2])[0]
    l1 = np.linalg.norm([pluck_point, pluck_elongation])
    l2 = np.linalg.norm([b - pluck_point, pluck_elongation])
    down = False
    down_point = 0
    for i, x in enumerate(xa):
        if x <= pluck_point:
            if abs(x - pluck_point) < 0.01:
                fx = f[i-1][1]
            else:
                fx = x * m1
            f[i][1] = fx
        else:
            if not down:
                down = True
                down_point = i
            if abs(x - pluck_point) < 0.01:
                fx = f[i-1][1]
            else:
                fx = (b - x) * m2
            f[i][1] = fx
   
    for i, x in enumerate(xa):
        if i > 0 and i < len(xa) - 1:
            if i < down_point:
                if abs(x - pluck_point) < 0.01:
                    xx = 0
                else:
                    xx = (l1 / (down_point - 1)) - (b / nx)
                f[i][0] = xx * x**4
            elif i == down_point:
                f[i][0] = 0
            else:
                #xx = (b / nx) - (l2 / (len(xa) - down_point))
                if abs(x - pluck_point) < 0.01:
                    xx = 0
                else:
                    xx = (l2 / (len(xa) - down_point)) - (b / nx)
                f[i][0] = xx * (b - x)**4
    f[:,0] = (pe/4) * f[:,0] / np.max(np.abs(f[:,0]))
    
    return f