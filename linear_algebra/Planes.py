from decimal import Decimal, getcontext
import numpy as np
from Vectors import Vector, areParallel, areOrthogonal

getcontext().prec = 30

class Plane(object):

	NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'

	def __init__(self, normal_vector=None, constant_term=None):
		self.dimension = 3

		if not normal_vector:
			all_zeros = ['0']*self.dimension
			normal_vector = Vector(all_zeros)
		self.normal_vector = normal_vector

		if not constant_term:
			constant_term = float(0.0)
		self.constant_term = float(constant_term)

		self.set_basepoint()

	def set_basepoint(self):
		try:
			n = self.normal_vector
			c = self.constant_term
			basepoint_coords = [0.0]*self.dimension

			initial_index = Plane.first_nonzero_index(n)
			initial_coefficient = n[initial_index]

			basepoint_coords[initial_index] = c/float(initial_coefficient)
			self.basepoint = Vector(basepoint_coords)

		except Exception as e:
			if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
				self.basepoint = None
			else:
				raise e

	def __str__(self):

		num_decimal_places = 3

		def write_coefficient(coefficient, is_initial_term=False):
			coefficient = round(coefficient, num_decimal_places)
			if coefficient % 1 == 0:
				coefficient = int(coefficient)

			output = ''

			if coefficient < 0:
				output += '-'
			if coefficient > 0 and not is_initial_term:
				output += '+'

			if not is_initial_term:
				output += ''

			if abs(coefficient) != 1:
				output += '{}'.format(abs(coefficient))

			return output

		n = self.normal_vector

		try:
			initial_index = Plane.first_nonzero_index(n)
			terms = [write_coefficient(n[i], is_initial_term=(i==initial_index)) + 'x_{}'.format(i+1)
					for i in range(self.dimension) if round(n[i], num_decimal_places) != 0]
			output = ' '.join(terms)

		except Exception as e:
			if str(e) == self.NO_NONZERO_ELTS_FOUND_MSG:
				output = '0'
			else:
				raise e

		constant = round(self.constant_term, num_decimal_places)
		if constant % 1 == 0:
			constant = int(constant)
		output += ' = {}'.format(constant)

		return output

	@staticmethod
	def first_nonzero_index(iterable):
		for k, item in enumerate(iterable):
			if not MyDecimal(item).is_near_zero():
				return k
		raise Exception(Plane.NO_NONZERO_ELTS_FOUND_MSG)

class MyDecimal(Decimal):
	def is_near_zero(self, eps=1e-10):
		return abs(self) < eps

def planesAreParallel(plane1, plane2):
	return areParallel(plane1.normal_vector, plane2.normal_vector)

def planesAreCoincident(plane1, plane2):
	if (not planesAreParallel(plane1, plane2)):
		return False
	return areOrthogonal((plane1.basepoint - plane2.basepoint), plane2.normal_vector)

def main():

	plane1 = Plane(Vector([-0.412, 3.806, 0.728]), -3.46)
	plane2 = Plane(Vector([1.03, -9.515, -1.82]), 8.65)
	print planesAreParallel(plane1, plane2)
	print planesAreCoincident(plane1, plane2)

	plane1 = Plane(Vector([2.611, 5.528, 0.283]), 4.6)
	plane2 = Plane(Vector([7.715, 8.306, 5.342]), 3.76)
	print planesAreParallel(plane1, plane2)
	print planesAreCoincident(plane1, plane2)

	plane1 = Plane(Vector([-7.926, 8.625, -7.217]), -7.952)
	plane2 = Plane(Vector([-2.642, 2.875, -2.404]), -2.443)
	print planesAreParallel(plane1, plane2)
	print planesAreCoincident(plane1, plane2)



if __name__ == '__main__':
	main()