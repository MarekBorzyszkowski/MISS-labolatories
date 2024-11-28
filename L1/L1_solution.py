import random
import matplotlib.pyplot as plt
import numpy as np
from L1_methods import random_test_in_regions, middle_x_of_throw, quadratic_formula_from_throw

ALLOWED = "allowed"
FORBIDDEN = "forbidden"
GRAVITY = 10

regions = {ALLOWED: [(2, 2, 2, 11), (7, 2, 7, 11) ], FORBIDDEN: [(3, 5, 6, 5), (3, 8, 6, 8)]}

v_0min_val = 1
v_0max_val = 30
alpha_min_val = 1
alpha_max_val = 89
number_of_tries_val = 100000

# v_0_test, alpha_test = 9.85, 68
# v_0_test, alpha_test = 12.013667928162864, 32.61695072706466

linspace_angle = np.linspace(alpha_min_val, alpha_max_val, number_of_tries_val)
linspace_v_0 = np.linspace(v_0min_val, v_0max_val, number_of_tries_val)

throws = random_test_in_regions(v_0min_val, v_0max_val, alpha_min_val, alpha_max_val, number_of_tries_val, regions)
fig1, ax1 = plt.subplots()
ax1.set_xlim([0, 90])
ax1.set_ylim([0, v_0max_val + 1])
ax1.set_xlabel('α')
ax1.set_ylabel('v_0')
ax1.plot(*zip(*throws), ".")
ax1.plot(linspace_angle, [v_0min_val] * len(linspace_angle), 'g')
ax1.plot(linspace_angle, [v_0max_val] * len(linspace_angle), 'g')
ax1.plot([alpha_min_val] * len(linspace_v_0), linspace_v_0, 'g')
ax1.plot([alpha_max_val] * len(linspace_v_0), linspace_v_0, 'g')
fig1.show()

fig2, ax2 = plt.subplots()
if len(throws) >=5:
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
ax2.set_xlabel('x')
ax2.set_ylabel('y')
fig2.show()

# results for exercise 3
print(len(throws) / number_of_tries_val)

# print(is_throw_valid(v_0_test, alpha_test))
# for region_category in regions:
#     for region in regions[region_category]:
#         print(is_line_segment_touching_parabola(middle_x_of_throw(v_0_test, alpha_test),
#                                                 quadratic_formula_from_throw(v_0_test, alpha_test), region))
