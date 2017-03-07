from decimal import Decimal, getcontext
from copy import deepcopy

from Vectors import Vector
from Planes import Plane

getcontext().prec = 30


class LinearSystem(object):

    ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG = 'All planes in the system should live in the same dimension'
    NO_SOLUTIONS_MSG = 'No solutions'
    INF_SOLUTIONS_MSG = 'Infinitely many solutions'

    def __init__(self, planes):
        try:
            d = planes[0].dimension
            for p in planes:
                assert p.dimension == d

            self.planes = planes
            self.dimension = d

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    def swap_rows(self, row1, row2):
        self[row1], self[row2] = self[row2], self[row1]

    def multiply_coefficient_and_row(self, coefficient, row):
        new_coordinates = [x*coefficient for x in self[row].normal_vector.coordinates]
        self[row] = Plane(normal_vector=Vector(new_coordinates), constant_term=(self[row].constant_term * coefficient))

    def add_multiple_times_row_to_row(self, coefficient, row_to_add, row_to_be_added_to):
        new_coordinates = [x*coefficient+y for x,y in zip(self[row_to_add].normal_vector.coordinates, self[row_to_be_added_to].normal_vector.coordinates)]
        self[row_to_be_added_to] = Plane(normal_vector=Vector(new_coordinates), constant_term=(self[row_to_be_added_to].constant_term + (self[row_to_add].constant_term*coefficient)))

    def compute_triangular_form(self):
        system = deepcopy(self)

        for current_row in range(len(system.planes)):

            system.check_for_swaps(current_row)
            indices = system.indices_of_first_nonzero_terms_in_each_row()
            first_term_pos = indices[current_row]
            while first_term_pos < current_row and first_term_pos != -1:
                coef = -float(system[current_row].normal_vector.coordinates[first_term_pos] / system[first_term_pos].normal_vector.coordinates[first_term_pos])
                system.add_multiple_times_row_to_row(coef, first_term_pos, current_row)
                indices = system.indices_of_first_nonzero_terms_in_each_row()
                if first_term_pos == indices[current_row]:
                    break
                else:
                    first_term_pos = indices[current_row]

        return system

    def check_for_swaps(self, row):

        indices = self.indices_of_first_nonzero_terms_in_each_row()

        if indices[row] > row:
            for i in range(row, len(indices)):
                if indices[i] <= row:
                    self.swap_rows(row, i)
                    break


    def indices_of_first_nonzero_terms_in_each_row(self):
        num_equations = len(self)
        num_variables = self.dimension

        indices = [-1] * num_equations

        for i,p in enumerate(self.planes):
            try:
                indices[i] = p.first_nonzero_index(p.normal_vector)
            except Exception as e:
                if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                    continue
                else:
                    raise e

        return indices

    def compute_rref(self):
        tf = self.compute_triangular_form()

        loop_num = 0
        if (len(tf[0].normal_vector.coordinates) > len(tf.planes)):
            loop_num = len(tf.planes)
        else:
            loop_num = len(tf[0].normal_vector.coordinates)
        for current_row in range(loop_num):

            indices = tf.indices_of_first_nonzero_terms_in_each_row()
            first_term_pos = indices[current_row]

            # make sure the coefficient of the column term is 1
            if (tf[current_row].normal_vector[current_row] != 1) and (tf[current_row].normal_vector[current_row] != 0):
                tf.multiply_coefficient_and_row((1.0/tf[current_row].normal_vector.coordinates[current_row]), current_row)

            i = current_row + 1
            while i < len(tf[current_row].normal_vector.coordinates):
                if (tf[current_row].normal_vector.coordinates[i] != 0) and (len(tf.planes) > i) and (tf[i].normal_vector.coordinates[i] != 0):
                    coef = -float(tf[current_row].normal_vector.coordinates[i] / tf[i].normal_vector.coordinates[i])
                    tf.add_multiple_times_row_to_row(coef, i, current_row)
                i += 1   

        return tf


    def __len__(self):
        return len(self.planes)


    def __getitem__(self, i):
        return self.planes[i]


    def __setitem__(self, i, x):
        try:
            assert x.dimension == self.dimension
            self.planes[i] = x

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    def __str__(self):
        ret = 'Linear System:\n'
        temp = ['Equation {}: {}'.format(i+1,p) for i,p in enumerate(self.planes)]
        ret += '\n'.join(temp)
        return ret


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps

def main():


    p1 = Plane(normal_vector=Vector([5.862,1.178,-10.366]), constant_term=-8.15)
    p2 = Plane(normal_vector=Vector([-2.931,-.589,5.183]), constant_term=-4.075)
    # p3 = Plane(normal_vector=Vector([1,2,-5]), constant_term=3)
    s = LinearSystem([p1,p2])
    print s
    print s.compute_rref()

    p1 = Plane(normal_vector=Vector([8.631,5.112,-1.816]), constant_term=-5.113)
    p2 = Plane(normal_vector=Vector([4.315,11.132,-5.27]), constant_term=-6.775)
    p3 = Plane(normal_vector=Vector([-2.158,3.01,-1.727]), constant_term=-0.831)
    s = LinearSystem([p1,p2,p3])
    print s
    print s.compute_rref()

    p1 = Plane(normal_vector=Vector([5.262,2.739,-9.878]), constant_term=-3.441)
    p2 = Plane(normal_vector=Vector([5.111,6.358,7.638]), constant_term=-2.152)
    p3 = Plane(normal_vector=Vector([2.016,-9.924,-1.367]), constant_term=-9.278)
    p4 = Plane(normal_vector=Vector([2.167,-13.593,-18.883]), constant_term=-10.567)
    s = LinearSystem([p1,p2,p3,p4])
    print s
    print s.compute_rref()


if __name__=='__main__':
    main()