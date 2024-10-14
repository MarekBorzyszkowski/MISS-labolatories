X = 'x'
Y = 'y'
Z = 'z'


def create_vector(x, y, z):
    return {X: x, Y: y, Z: z}


def dot_product(vector_1, vector_2):
    return vector_1[X] * vector_2[X] + vector_1[Y] * vector_2[Y] + vector_1[Z] * vector_2[Z]


def add_vectors(vector_1, vector_2):
    return create_vector(vector_1[X] + vector_2[X], vector_1[Y] + vector_2[Y], vector_1[Z] + vector_2[Z])


def multiply_vector(vector, scalar):
    return create_vector(vector[X] * scalar, vector[Y] * scalar, vector[Z] * scalar)


def vector_module(vector):
    return (vector[X] ** 2 + vector[Y] ** 2 + vector[Z] ** 2) ** 0.5
