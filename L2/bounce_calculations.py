from numpy import ndarray

from L2.Bounce import Bounce
from L2.vector_methods import *
from scipy.optimize import fsolve, root_scalar

GRAVITY = 10

A = 'A'
B = 'B'
C = 'C'
D = 'D'
E = 'E'
F = 'F'


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


def calculate_time_of_impact_v2(p_0, v_0, surface_equation, delta):
    result = root_scalar(surface_equation(p_0,v_0), bracket=(0 + delta,100))
    return result.root


def calculate_bounce(p_0, v_0, se_params):
    t_i = calculate_time_of_impact(p_0, v_0, se_params)
    x_prim = calculate_x_in_time(p_0[X], v_0[X], t_i)
    y_prim = calculate_y_in_time(p_0[Y], v_0[Y], t_i)
    z_prim = calculate_z_in_time(p_0[Z], v_0[Z], t_i)
    v_z_i = v_0[Z] - GRAVITY * t_i
    v_i = create_vector(v_0[X], v_0[Y], v_z_i)
    return Bounce(p_0, create_vector(x_prim, y_prim, z_prim), v_0, v_i, t_i)


def calculate_bounce_v2(p_0, v_0, surface_equation, delta):
    t_i = calculate_time_of_impact_v2(p_0, v_0, surface_equation, delta)
    x_prim = calculate_x_in_time(p_0[X], v_0[X], t_i)
    y_prim = calculate_y_in_time(p_0[Y], v_0[Y], t_i)
    z_prim = calculate_z_in_time(p_0[Z], v_0[Z], t_i)
    v_z_i = v_0[Z] - GRAVITY * t_i
    v_i = create_vector(v_0[X], v_0[Y], v_z_i)
    return Bounce(p_0, create_vector(x_prim, y_prim, z_prim), v_0, v_i, t_i)


def calculate_new_velocity(p_i, v_0, se_params, t_i, k):
    # -Ax^2 - Bx - Cy^2 - Dy - Exy - F + z = 0
    v_z = v_0[Z] - GRAVITY * t_i
    v_i_before = create_vector(v_0[X], v_0[Y], v_z)
    N = create_vector((2 * se_params[A] * p_i[X] + se_params[B] + se_params[E] * p_i[Y]),  # df/dx
                      (2 * se_params[C] * p_i[Y] + se_params[D] + se_params[E] * p_i[X]),  # df/dy
                      -1)  # df/dz
    N_length = vector_module(N)
    n = multiply_vector(N, 1 / N_length)
    return multiply_vector(add_vectors(v_i_before, multiply_vector(n, -2 * dot_product(v_i_before, n))), k**(1/2))


def calculate_new_velocity_v2(p_i, v_0, surface_equation, t_i, k, delta):
    v_z = v_0[Z] - GRAVITY * t_i
    v_i_before = create_vector(v_0[X], v_0[Y], v_z)
    N = create_vector((surface_equation(p_i[X] + delta, p_i[Y]) - surface_equation(p_i[X] - delta, p_i[Y]) )/(2*delta),  # df/dx
                      (surface_equation(p_i[X], p_i[Y] + delta) - surface_equation(p_i[X], p_i[Y] - delta) )/(2*delta),  # df/dy
                      -1)  # df/dz
    N_length = vector_module(N)
    n = multiply_vector(N, 1 / N_length)
    return multiply_vector(add_vectors(v_i_before, multiply_vector(n, -2 * dot_product(v_i_before, n))), k**(1/2))


def bounce_ball(bounces_number, p_0, v_0, se_params, k):
    bounces = []
    current_p = p_0
    current_v = v_0
    for i in range(bounces_number):
        bounce = calculate_bounce(current_p, current_v, se_params)
        current_p = bounce.point_i
        current_v = calculate_new_velocity(current_p, current_v, se_params, bounce.t_i, k)
        bounces.append(bounce)
    return bounces


def bounce_ball_v2(bounces_number, p_0, v_0, surface_equation, surface_equation_in_time, k, delta):
    bounces = []
    current_p = p_0
    current_v = v_0
    for i in range(bounces_number):
        bounce = calculate_bounce_v2(current_p, current_v, surface_equation_in_time, delta)
        current_p = bounce.point_i
        current_v = calculate_new_velocity_v2(current_p, current_v, surface_equation, bounce.t_i, k, delta)
        bounces.append(bounce)
    return bounces
