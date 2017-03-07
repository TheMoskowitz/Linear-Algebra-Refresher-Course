from decimal import Decimal, getcontext
import numpy as np
from Vectors import Vector, areParallel, areOrthogonal

getcontext().prec = 30

class Line(object):

	NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'

	def __init__(self, normal_vector=None, constant_term=None):
		self.dimension = 2

		if not normal_vector:
			all_zeros = ['0']*self.dimension
			normal_vector = Vector(all_zeros)
		self.normal_vector = normal_vector
		self.direction_vector = Vector([-normal_vector[1], normal_vector[0]])

		if not constant_term:
			constant_term = float(0.0)
		self.constant_term = float(constant_term)

		self.set_basepoint()

	def set_basepoint(self):
		try:
			n = self.normal_vector
			c = self.constant_term
			basepoint_coords = [0.0]*self.dimension

			initial_index = Line.first_nonzero_index(n)
			initial_coefficient = n[initial_index]

			basepoint_coords[initial_index] = c/float(initial_coefficient)
			self.basepoint = Vector(basepoint_coords)

		except Exception as e:
			if str(e) == Line.NO_NONZERO_ELTS_FOUND_MSG:
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
			initial_index = Line.first_nonzero_index(n)
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
		raise Exception(Line.NO_NONZERO_ELTS_FOUND_MSG)

class MyDecimal(Decimal):
	def is_near_zero(self, eps=1e-10):
		return abs(self) < eps

def linesAreParallel(line1, line2):
	return areParallel(line1.normal_vector, line2.normal_vector)

def linesAreCoincident(line1, line2):
	if (not linesAreParallel(line1, line2)):
		return False
	return areOrthogonal((line1.basepoint - line2.basepoint), line2.normal_vector)

def findIntersection(line1, line2):
	if (linesAreParallel(line1, line2)):
		if (linesAreCoincident(line1, line2))
			return line1
		else:
			raise Exception('Parallel lines have no intersection point')
	x = (line2.normal_vector[1]*line1.constant_term - line1.normal_vector[1]*line2.constant_term)/(line1.normal_vector[0]*line2.normal_vector[1] - line1.normal_vector[1]*line2.normal_vector[0])
	y = (-line2.normal_vector[0]*line1.constant_term + line1.normal_vector[0]*line2.constant_term)/(line1.normal_vector[0]*line2.normal_vector[1] - line1.normal_vector[1]*line2.normal_vector[0])
	return Vector((x,y))


def main():

	line1 = Line(normal_vector=Vector([4.046, 2.836]), constant_term=1.21)
	line2 = Line(normal_vector=Vector([10.115, 7.09]), constant_term=3.025)
	print linesAreParallel(line1, line2)
	print linesAreCoincident(line1, line2)
	if (not linesAreParallel(line1, line2)):
		print findIntersection(line1, line2)

	line1 = Line(normal_vector=Vector([7.204, 3.182]), constant_term=8.68)
	line2 = Line(normal_vector=Vector([8.172, 4.114]), constant_term=9.883)
	print linesAreParallel(line1, line2)
	print linesAreCoincident(line1, line2)
	if (not linesAreParallel(line1, line2)):	
		print findIntersection(line1, line2)

	line1 = Line(normal_vector=Vector([1.182, 5.562]), constant_term=6.744)
	line2 = Line(normal_vector=Vector([1.773, 8.343]), constant_term=9.525)
	print linesAreParallel(line1, line2)
	print linesAreCoincident(line1, line2)
	if (not linesAreParallel(line1, line2)):
		print findIntersection(line1, line2)




if __name__=='__main__':
	main()












