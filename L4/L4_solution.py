import scipy.io as sio
import numpy as np

mat = sio.loadmat('dataset63332.mat')
x_A = mat['xt']
y_A = mat['yt']
x_B = mat['xv']
y_B = mat['yv']

K = x_A.shape[0]
Q = x_B.shape[0]

f_A = np.reshape(mat['ft'], shape=(K))
f_B = np.reshape(mat['fv'], shape=(Q))

beta = 1.45


def gaussian_phi_fun(r):
    return np.exp(-(beta * r) ** 2)


def multiquadric_phi_fun(r):
    return np.sqrt(1 + (beta * r) ** 2)


def polynomial_phi_fun(r):
    return r ** 2 + beta ** 2


phi_function = multiquadric_phi_fun


def phi(r):
    return phi_function(r)


def norm(p1, p2):
    return np.linalg.norm(np.subtract(p1, p2))


def create_base_fun():
    return lambda x, y: np.array([phi(norm([x_A[i], y_A[i]], [x, y])) for i in range(K)], )


def create_phi_matrix(base_function, x_vector, y_vector):
    return np.array([base_function(x_vector[i], y_vector[i]) for i in range(x_vector.shape[0])])


def calculate_f_b_2(calculated_a, base_function):
    return np.matmul(create_phi_matrix(base_function, x_B, y_B), calculated_a)


base_fun = create_base_fun()
phi_matrix = create_phi_matrix(base_fun, x_A, y_A)
a = np.reshape(np.linalg.solve(phi_matrix, f_A), shape=(K))
calculated_f_B = calculate_f_b_2(a, base_fun)
errors = np.absolute(calculated_f_B - f_B)
max_absolute_error = np.max(errors)
mean_absolute_error = np.divide(np.sum(errors), Q)
print(max_absolute_error)
print(mean_absolute_error)
