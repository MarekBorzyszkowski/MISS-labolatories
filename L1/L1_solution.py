import random
import matplotlib.pyplot as plt
import numpy as np
from L1_methods import random_test_in_regions, middle_x_of_throw, quadratic_formula_from_throw

ALLOWED = "allowed"
FORBIDDEN = "forbidden"
GRAVITY = 10

regions = {ALLOWED: [(2, 1, 2, 4), (4, 1, 7, 1)], FORBIDDEN: [(3, 2, 3, 3)]}

v_0min_val = 1
v_0max_val = 15
alpha_min_val = 1
alpha_max_val = 89
number_of_tries_val = 100000

# v_0_test, alpha_test = 9.85, 68
# v_0_test, alpha_test = 12.013667928162864, 32.61695072706466


throws = random_test_in_regions(v_0min_val, v_0max_val, alpha_min_val, alpha_max_val, number_of_tries_val, regions)
fig1, ax1 = plt.subplots()
ax1.set_xlim([0, 90])
ax1.set_ylim([0, v_0max_val + 1])
ax1.plot(*zip(*throws), ".")
fig1.show()

fig2, ax2 = plt.subplots()
for throw in random.sample(throws, 5):
    alpha_val, v_0_val = throw
    x_max = middle_x_of_throw(v_0_val, alpha_val) * 2
    x_space = np.linspace(0, x_max, 10000)
    x_func = np.vectorize(quadratic_formula_from_throw(v_0_val, alpha_val))
    y_space = x_func(x_space)
    ax2.plot(x_space, y_space)
    # print(f"{x_max}, {v_0_val}, {alpha_val}")
for allowed in regions[ALLOWED]:
    region_x0, region_y0, region_x1, region_y1 = allowed
    ax2.plot([region_x0, region_x1], [region_y0, region_y1], "g")
for forbidden in regions[FORBIDDEN]:
    region_x0, region_y0, region_x1, region_y1 = forbidden
    ax2.plot([region_x0, region_x1], [region_y0, region_y1], "r")
fig2.show()

# print(is_throw_valid(v_0_test, alpha_test))
# for region_category in regions:
#     for region in regions[region_category]:
#         print(is_line_segment_touching_parabola(middle_x_of_throw(v_0_test, alpha_test),
#                                                 quadratic_formula_from_throw(v_0_test, alpha_test), region))
