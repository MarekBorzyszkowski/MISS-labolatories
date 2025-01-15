import numpy as np
from matplotlib import pyplot as plt, cm
import matplotlib as mpl

# mpl.use('qt5agg')
np.random.seed(42)

OUTSIDE_CURVE = 0
INSIDE_CURVE_POSITIVE = 1
INSIDE_CURVE_NEGATIVE = 2

x_a, x_b = -1.0, 1.0
y_a, y_b = 0.0, 2.0
z_a, z_b = -3.5, 5.5
f_x = lambda x, y: x * x + x * y * y
F_x_integral = lambda x_a, x_b, y_a, y_b: ((x_b ** 3) * y_b / 3.0 + (x_b ** 2) * (y_b ** 3) / 6.0 -
                                           (((x_a ** 3) * y_b / 3.0) + (x_a ** 2) * (y_b ** 3) / 6.0) -
                                           (((x_b ** 3) * y_a / 3.0) + (x_b ** 2) * (y_a ** 3) / 6.0) +
                                           (((x_a ** 3) * y_a / 3.0) + (x_a ** 2) * (y_a ** 3) / 6.0))


def condition_checker(z, z_real):
    in_p, in_n, outside = [],[],[]
    for z_val, z_real_val in zip(z, z_real):
        if z_val >= 0:
            if z_val <= z_real_val:
                in_p.append(True)
                in_n.append(False)
                outside.append(False)
            else:
                in_p.append(False)
                in_n.append(False)
                outside.append(True)
        else:
            if z_val >= z_real_val:
                in_p.append(False)
                in_n.append(True)
                outside.append(False)
            else:
                in_p.append(False)
                in_n.append(False)
                outside.append(True)
    return in_p, in_n, outside

def plot_curve_with_points(points_outside_curve, points_inside_le_curve, points_inside_gt_curve):
    number_of_samples = 10 ** 1
    linspace_x = np.linspace(x_a, x_b, number_of_samples)
    linspace_y = np.linspace(y_a, y_b, number_of_samples)
    ax = plt.figure().add_subplot(projection='3d')
    ax.set_xlim([x_a - 1, x_b + 1])
    ax.set_ylim([y_a - 1, y_b + 1])
    ax.set_zlim([z_a - 1, z_b + 1])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    if len(points_inside_le_curve)>0:
        ax.scatter(*zip(*points_inside_le_curve), color='b', marker='o', s=2)
    if len(points_inside_gt_curve) > 0:
        ax.scatter(*zip(*points_inside_gt_curve), color='g', marker='o', s=2)
    if len(points_outside_curve) > 0:
        ax.scatter(*zip(*points_outside_curve), color='r', marker='o', s=2)
    x_meshed, y_meshed = np.meshgrid(linspace_x, linspace_y)
    z = np.array([f_x(linspace_x, np.array([y_val] * number_of_samples)) for y_val in linspace_y])
    ax.plot_surface(x_meshed, y_meshed, z, cmap=cm.coolwarm, rstride=1, cstride=1, alpha=0.6)
    plt.show(block=True)


def plot_error(errors_hm, error_mean, N_vals, title):
    fig2, ax2 = plt.subplots()
    ax2.set_title(title)
    ax2.set_xlabel('N')
    ax2.set_ylabel('error')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.plot(N_vals, errors_hm,  label='Hit and miss error')
    ax2.plot(N_vals, errors_mean,  label='Mean error')
    ax2.plot(N_vals, 1/np.sqrt(N_vals),  label='1/sqrt(N_vals)')
    ax2.legend()
    fig2.savefig(title.replace(" ", "_") + '.png')


def calc_integrals(N_min, N_max):
    N_vals, hit_mis_vals, mean_integral_vals = [], [], []
    for i in range(N_min, N_max + 1):
        print("calculating integrals for: ", i)
        N = 10 ** i
        N_vals.append(N)
        x_random_samples = (x_b - x_a) * np.random.random_sample(N) + x_a
        y_random_samples = (y_b - y_a) * np.random.random_sample(N) + y_a
        z_random_samples = (z_b - z_a) * np.random.random_sample(N) + z_a
        f_x_sampled = f_x(x_random_samples, y_random_samples)
        points = np.vstack((x_random_samples, y_random_samples, z_random_samples)).T
        in_p, in_n, outside = condition_checker(z_random_samples, f_x_sampled)
        points_outside_curve = points[outside]
        points_inside_le_curve = points[in_p]
        points_inside_gt_curve = points[in_n]

        N_b = len(points_inside_le_curve) - len(points_inside_gt_curve)
        hit_miss_integral = (x_b - x_a) * (y_b - y_a) * (z_b - z_a) * N_b / N
        mean_val = np.mean(f_x_sampled) * (x_b - x_a) * (y_b - y_a)
        print("hit miss: ", hit_miss_integral)
        print("mean value: ", mean_val)
        hit_mis_vals.append(hit_miss_integral)
        mean_integral_vals.append(mean_val)
        if i < 4:
            plot_curve_with_points(points_outside_curve, points_inside_le_curve, points_inside_gt_curve)
    return np.array(N_vals), np.array(hit_mis_vals), np.array(mean_integral_vals)


N_array, hit_mis_array, mean_integral_array = calc_integrals(1, 8)
integral_value = F_x_integral(x_a, x_b, y_a, y_b)
print("Analytical integral val", integral_value)
errors_hit_miss = abs(hit_mis_array - integral_value)
errors_mean = abs(mean_integral_array - integral_value)
plot_error(errors_hit_miss, errors_mean, N_array, "Algorithms absolute errors")
