import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt, cm

from L2.Bounce import Bounce
from L2.bounce_calculations import calculate_x_in_time, calculate_y_in_time, calculate_z_in_time, bounce_ball, A, B, C, \
    D, E, F
from L2.vector_methods import *

mpl.use('qt5agg')

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


calculated_bounces = bounce_ball(number_of_bounces, P_0, V_0, surface_equation_parameters)
print_bouncing_ball(surface_equation, calculated_bounces)
