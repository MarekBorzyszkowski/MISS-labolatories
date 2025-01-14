import numpy as np

np.random.seed(42)

a = 0
b = 1
h = 3
f_x = lambda x : np.exp(x)
condition_checker = lambda y, y_real : y <= y_real
def calc_integrals(N_min, N_max):
    for i in range(N_min, N_max+1):
        print("calculating integrals for: ", i)
        N = 10 ** i
        x_random_samples = (b - a) * np.random.random_sample(N) + a
        y_random_samples = (h - 0) * np.random.random_sample(N) + 0
        f_x_sampled = f_x(x_random_samples)
        points = np.vstack((x_random_samples, y_random_samples)).T
        points_le_curve_condition = condition_checker(y_random_samples, f_x_sampled)
        # points_over_curve = points[~points_le_curve_condition]
        points_le_curve = points[points_le_curve_condition]
        N_b = len(points_le_curve)
        hit_miss_integral = (b - a) * h * N_b / N
        mean_val = np.mean(f_x_sampled) * (b - a)
        print("hit miss: ", hit_miss_integral)
        print("mean value: ", mean_val)

calc_integrals(4, 8)