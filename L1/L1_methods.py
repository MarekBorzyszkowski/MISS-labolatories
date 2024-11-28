import math
import random

ALLOWED = "allowed"
FORBIDDEN = "forbidden"
GRAVITY = 10


def random_test_in_regions(v_0min, v_0max, alpha_min, alpha_max, number_of_tries, regions):
    successes = []
    for _ in range(number_of_tries):
        v_0 = random.uniform(v_0min, v_0max)
        alpha = random.uniform(alpha_min, alpha_max)
        if is_throw_valid(v_0, alpha, regions):
            successes.append((alpha, v_0))
    return successes


def is_throw_valid(v_0, alpha, regions):
    quadratic_func = quadratic_formula_from_throw(v_0, alpha)
    x_middle = middle_x_of_throw(v_0, alpha)
    for allowed_region in regions[ALLOWED]:
        if not is_line_segment_touching_parabola(x_middle, quadratic_func, allowed_region):
            return False
    for forbidden_region in regions[FORBIDDEN]:
        if is_line_segment_touching_parabola(x_middle, quadratic_func, forbidden_region):
            return False
    return True


def quadratic_formula_from_throw(v_0, alpha):
    sin_alpha = math.sin(math.radians(alpha))
    sin_2alpha = math.sin(math.radians(2 * alpha))
    x_0 = 0
    x_1 = x_0 + (v_0 ** 2) * sin_2alpha / GRAVITY
    x_middle = (x_0 + x_1) / 2
    h_max = (v_0 ** 2) * (sin_alpha ** 2) / (2 * GRAVITY)
    a = h_max / ((x_middle - x_0) * (x_middle - x_1))
    return lambda x: a * (x - x_0) * (x - x_1)


def middle_x_of_throw(v_0, alpha):
    sin_2alpha = math.sin(math.radians(2 * alpha))
    x_0 = 0
    x_1 = x_0 + (v_0 ** 2) * sin_2alpha / GRAVITY
    return (x_0 + x_1) / 2


def is_line_segment_touching_parabola(x_middle, quadratic_func, line_segment):
    line_x0, line_y0, line_x1, line_y1 = line_segment
    if (line_y0 - quadratic_func(line_x0)) * (line_y1 - quadratic_func(line_x1)) <= 0:
        return True
    if line_x0 < x_middle < line_x1:
        y_of_line_in_middle_of_parabola = (line_y1 - line_y0) / (line_x1 - line_x0) * (x_middle - line_x0) + line_y0
        if (line_y0 - quadratic_func(line_x0)) * (y_of_line_in_middle_of_parabola - quadratic_func(x_middle)) <= 0 or (
                line_y1 - quadratic_func(line_x1)) * (y_of_line_in_middle_of_parabola - quadratic_func(x_middle)) <= 0:
            return True
        else:
            return False
    else:
        return False
