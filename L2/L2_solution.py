import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt, cm

from L2.Bounce import Bounce
from L2.vector_methods import *

mpl.use('qt5agg')

GRAVITY = 10

A = 'A'
B = 'B'
C = 'C'
D = 'D'
E = 'E'
F = 'F'

v_x0 = 2
v_y0 = 4
v_z0 = 0
V_0 = create_vector(v_x0, v_y0, v_z0)

x_0 = -5
y_0 = 0
z_0 = 200
P_0 = create_vector(x_0, y_0, z_0)

surface_equation = lambda x, y: (x ** 2 + y ** 2)
surface_equation_parameters = {A: 1, B: 0, C: 1, D: 0, E: 0, F: 0}  # Ax^2 + Bx + Cy^2 + Dy + Exy + F

k = 1
plot_size = 15
plot_density = 100
bounce_resolution = 1000
number_of_bounces = 5


def calculate_x_in_time(x, v_x, t):
    return x + v_x * t


def calculate_y_in_time(y, v_y, t):
    return y + v_y * t


def calculate_z_in_time(z, v_z, t):
    return z + v_z * t - (GRAVITY / 2) * (t ** 2)


def quadratic_equation_results(a, b, delta):
    return (-b - (delta ** (1 / 2))) / (2 * a), (-b + (delta ** (1 / 2))) / (2 * a)


def calculate_time_of_impact(p_0, v_0, se_params):
    # at^2 + bt + c
    a = se_params[A] * (v_0[X] ** 2) + se_params[C] * (v_0[Y] ** 2) + se_params[E] * v_0[X] * v_0[Y] + GRAVITY / 2
    b = (2 * se_params[A] * p_0[X] * v_0[X] + se_params[B] * v_0[X] + 2 * se_params[C] * p_0[Y] * v_0[Y] + se_params[
        D] * v_0[Y]
         + se_params[E] * p_0[X] * v_0[Y] + se_params[E] * p_0[Y] * v_0[X] - v_0[Z])
    c = (se_params[A] * (p_0[X] ** 2) + se_params[B] * p_0[X] + se_params[C] * (p_0[Y] ** 2) + se_params[D] * p_0[Y]
         + se_params[E] * p_0[X] * p_0[Y] + se_params[F] - p_0[Z])
    delta = (b ** 2) - 4 * a * c
    if delta < 0:
        return None
    t_1, t_2 = quadratic_equation_results(a, b, delta)
    if t_1 > 0 and t_2 > 0:
        return t_1 if t_1 > t_2 else t_2
    elif t_1 > 0:
        return t_1
    elif t_2 > 0:
        return t_2
    else:
        return None


def print_bouncing_ball(surface, bounces):
    ax = plt.figure().add_subplot(projection='3d')

    for bounce in bounces:
        x = np.linspace(0, bounce.t_i, bounce_resolution)
        x = calculate_x_in_time(bounce.point_0[X], bounce.velocity_0[X], x)
        y = np.linspace(0, bounce.t_i, bounce_resolution)
        y = calculate_y_in_time(bounce.point_0[Y], bounce.velocity_0[Y], y)
        z = np.linspace(0, bounce.t_i, bounce_resolution)
        z = calculate_z_in_time(bounce.point_0[Z], bounce.velocity_0[Z], z)
        ax.plot(x, y, z, color='red', markersize=3)

    x = np.linspace(-plot_size, plot_size, plot_density)
    y = np.linspace(-plot_size, plot_size, plot_density)
    x, y = np.meshgrid(x, y)
    z = surface(x, y)
    ax.plot_surface(x, y, z, cmap=cm.coolwarm, rstride=1, cstride=1,
                    alpha=0.6)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    plt.show(block=True)


def calculate_bounce(p_0, v_0, se_params):
    t_1 = calculate_time_of_impact(p_0, v_0, se_params)
    x_prim = calculate_x_in_time(x_0, v_x0, t_1)
    y_prim = calculate_y_in_time(y_0, v_y0, t_1)
    z_prim = calculate_z_in_time(z_0, v_z0, t_1)
    print(x_prim, y_prim, z_prim, t_1, x_prim ** 2 + y_prim ** 2)
    return Bounce(p_0, create_vector(x_prim, y_prim, z_prim), v_0, t_1)


def calculate_new_velocity(p_i, v_0, se_params, t_i):
    # -Ax^2 - Bx - Cy^2 - Dy - Exy - F + z = 0
    v_z = v_0[Z] - GRAVITY * t_i
    v_i_before = create_vector(v_0[X], v_0[Y], v_z)
    N = create_vector(2 * se_params[A] * p_i[X] + se_params[B] + se_params[E] * p_i[Y],  # df/dx
                      2 * se_params[C] * p_i[Y] + se_params[D] + se_params[E] * p_i[X],  # df/dy
                      1)  # df/dz
    N_length = vector_module(N)
    n = multiply_vector(N, 1 / N_length)
    return add_vectors(v_i_before, multiply_vector(n, -2 * dot_product(v_0, n)))


def bounce_ball(bounces_number, p_0, v_0, se_params):
    bounces = []
    current_p = p_0
    current_v = v_0
    for i in range(bounces_number):
        print(current_p, current_v)
        bounce = calculate_bounce(current_p, current_v, se_params)
        current_p = bounce.point_i
        current_v = calculate_new_velocity(current_p, current_v, se_params, bounce.t_i)
        bounces.append(bounce)
    return bounces


calculated_bounces = bounce_ball(number_of_bounces,P_0, V_0, surface_equation_parameters)
print_bouncing_ball(surface_equation, calculated_bounces)
