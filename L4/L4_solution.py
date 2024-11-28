import numpy as np
import scipy.io as sio

mat = sio.loadmat('dataset0.mat')
x_A = mat['xv']
y_A = mat['yv']
x_B = mat['xt']
y_B = mat['yt']

K = x_A.shape[0]
Q = x_B.shape[0]

f_A = np.reshape(mat['fv'], shape=(K))
f_B = np.reshape(mat['ft'], shape=(Q))

beta = 1.5

def gaussian_phi_fun(r):
    return np.exp(-(beta * r)**2)

def multiquadric_phi_fun(r):
    return np.sqrt(1+(beta * r)**2)

phi_function = gaussian_phi_fun

def phi(r):
    return phi_function(r)

def norm(p1, p2):
    return np.linalg.norm(np.subtract(p1, p2))

def create_base_fun():
    return lambda x,y: np.array([phi(norm([x_A[i], y_A[i]], [x,y])) for i in range(K)], )


def create_phi_matrix(base_function):
    return np.array([base_function(x_A[i], y_A[i]) for i in range(K)])

def calculate_f_b_i(calculated_a, base_function, i):
    return np.sum(np.multiply(calculated_a, base_function(x_B[i], y_B[i])))

def calculate_f_b(calculated_a, base_function):
    return np.array([calculate_f_b_i(calculated_a, base_function, i) for i in range(Q)])


base_fun = create_base_fun()
phi_matrix = create_phi_matrix(base_fun)
a = np.reshape(np.linalg.solve(phi_matrix, f_A), shape=(K))
calculated_f_B = calculate_f_b(a, base_fun)
errors = np.absolute(calculated_f_B - f_B)
max_absolute_error = np.max(errors)
mean_absolute_error = np.divide(np.sum(errors),Q)
print(max_absolute_error)
print(mean_absolute_error)
