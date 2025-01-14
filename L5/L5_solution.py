import numpy as np
import matplotlib.pyplot as plt

np.random.seed(42)

a = 0
b = 1
h = 3
f_x = lambda x: np.exp(x)
F_x = lambda x: np.exp(x)
condition_checker = lambda y, y_real: y <= y_real


def plot_curve_with_points(points_over_curve, points_le_curve):
    number_of_samples = 10 ** 4
    linspace_angle = np.linspace(a, b, number_of_samples)
    linspace_v_0 = np.linspace(0, h, number_of_samples)
    fig1, ax1 = plt.subplots()
    ax1.set_xlim([a - 1, b + 1])
    ax1.set_ylim([-1, h + 1])
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.plot(*zip(*points_le_curve), "b.")
    ax1.plot(*zip(*points_over_curve), "r.")
    ax1.plot(linspace_angle, f_x(linspace_angle), "g")
    ax1.plot(linspace_angle, [0] * len(linspace_angle), 'g')
    ax1.plot(linspace_angle, [h] * len(linspace_angle), 'g')
    ax1.plot([a] * len(linspace_v_0), linspace_v_0, 'g')
    ax1.plot([b] * len(linspace_v_0), linspace_v_0, 'g')
    fig1.show()


def plot_error(errors, N_vals, title):
    fig2, ax2 = plt.subplots()
    ax2.set_title(title)
    ax2.set_xlabel('N')
    ax2.set_ylabel('error')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.plot(N_vals, errors)
    fig2.show()


def calc_integrals(N_min, N_max):
    N_vals, hit_mis_vals, mean_integral_vals = [], [], []
    for i in range(N_min, N_max + 1):
        print("calculating integrals for: ", i)
        N = 10 ** i
        N_vals.append(N)
        x_random_samples = (b - a) * np.random.random_sample(N) + a
        y_random_samples = (h - 0) * np.random.random_sample(N) + 0
        f_x_sampled = f_x(x_random_samples)
        points = np.vstack((x_random_samples, y_random_samples)).T
        points_le_curve_condition = condition_checker(y_random_samples, f_x_sampled)
        points_over_curve = points[~points_le_curve_condition]
        points_le_curve = points[points_le_curve_condition]
        N_b = len(points_le_curve)
        hit_miss_integral = (b - a) * h * N_b / N
        mean_val = np.mean(f_x_sampled) * (b - a)
        print("hit miss: ", hit_miss_integral)
        print("mean value: ", mean_val)
        hit_mis_vals.append(hit_miss_integral)
        mean_integral_vals.append(mean_val)
        if i < 5:
            plot_curve_with_points(points_over_curve, points_le_curve)
    return np.array(N_vals), np.array(hit_mis_vals), np.array(mean_integral_vals)


N_array, hit_mis_array, mean_integral_array = calc_integrals(1, 8)
integral_value = F_x(b) - F_x(a)
errors_hit_miss = abs(hit_mis_array - integral_value)
errors_mean = abs(mean_integral_array - integral_value)
plot_error(errors_hit_miss, N_array, "Hit miss algorithm absolute errors")
plot_error(errors_mean, N_array, "Mean integral algorithm absolute errors")
