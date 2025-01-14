import matplotlib as mpl
import scipy.io as sio
import numpy as np
from matplotlib import pyplot as plt, cm

mpl.use('qt5agg')

x_min = 1
x_max = 6
y_min = 1
y_max = 5

plot_density = 100

mat = sio.loadmat('dataset63332.mat')
x_A = mat['xt']
y_A = mat['yt']
x_B = mat['xv']
y_B = mat['yv']

K = x_A.shape[0]
Q = x_B.shape[0]

x_A = np.reshape(x_A, shape=(K))
y_A = np.reshape(y_A, shape=(K))
x_B = np.reshape(x_B, shape=(Q))
y_B = np.reshape(y_B, shape=(Q))

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

np.array([x_val]*plot_density)
def phi(r):
    return phi_function(r)


def norm(p1, p2):
    return np.linalg.norm(np.subtract(p1, p2))


def create_base_fun():
    return lambda x, y: np.array([phi(norm([x_A[i], y_A[i]], [x, y])) for i in range(K)], )


def create_phi_matrix(base_function, x_vector, y_vector):
    return np.array([base_function(x_vector[i], y_vector[i]) for i in range(x_vector.shape[0])])


def calculate_f_b_2(calculated_a, base_function, x_vector, y_vector):
    return np.matmul(create_phi_matrix(base_function, x_vector, y_vector), calculated_a)


def plot_interpolated_surface(a_vector, base_function):
    ax = plt.figure().add_subplot(projection='3d')
    ax.scatter(x_A, y_A, f_A, color='r', marker='o', s=2)
    ax.scatter(x_B, y_B, f_B, color='g', marker='o', s=2)

    x = np.linspace(x_min, x_max, plot_density)
    y = np.linspace(y_min, y_max, plot_density)
    x_meshed, y_meshed = np.meshgrid(x, y)
    z = np.array([calculate_f_b_2(a_vector, base_function, x, np.array([y_val]*plot_density)) for y_val in y])
    ax.plot_surface(x_meshed, y_meshed, z, cmap=cm.coolwarm, rstride=1, cstride=1, alpha=0.6)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    plt.show(block=True)


base_fun = create_base_fun()
phi_matrix = create_phi_matrix(base_fun, x_A, y_A)
a = np.reshape(np.linalg.solve(phi_matrix, f_A), shape=(K))
calculated_f_B = calculate_f_b_2(a, base_fun, x_B, y_B)
errors = np.absolute(calculated_f_B - f_B)
max_absolute_error = np.max(errors)
mean_absolute_error = np.divide(np.sum(errors), Q)
print(max_absolute_error)
print(mean_absolute_error)

plot_interpolated_surface(a, base_fun)