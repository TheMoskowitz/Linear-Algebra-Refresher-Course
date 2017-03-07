import numpy as np
from decimal import Decimal

class Vector(object):
	
	def __init__(self, coordinates):
		try:
			if not coordinates:
				raise ValueError
			self.coordinates = tuple(coordinates)
			self.dimension = len(coordinates)
		except ValueError:
			raise ValueError("The coordinates must be nonempty")
		except TypeError:
			raise TypeError('The coordinates must be an iterable')

	def __str__(self):
		return 'Vector: {}'.format(self.coordinates)

	def __iter__(self):
		return iter(self.coordinates)

	def __getitem__(self, index):
		return self.coordinates[index]

	def __len__(self):
		return self.dimension

	def __eq__(self, v):
		return self.coordinates == v.coordinates

	def __add__(self, v):
		new_coordinates = [x+y for x,y in zip(self.coordinates, v.coordinates)]
		return Vector(new_coordinates)

	def __sub__(self, v):
		new_coordinates = [x-y for x,y in zip(self.coordinates, v.coordinates)]
		return Vector(new_coordinates)

	def __mul__(self, num):
		new_coordinates = [x*num for x in self.coordinates]
		return Vector(new_coordinates)

	def getMagnitude(self):
		total = np.sum([x**2 for x in self.coordinates])
		return np.sqrt(total)

	def normalize(self):
		try:
			mag = self.getMagnitude()
			new_coordinates = [(x/mag) for x in self.coordinates]
			return Vector(new_coordinates)
		except ZeroDivisionError:
			raise Exception("Cannot normalize the zero vetor")

def dotProduct(v1, v2):
	try:
		return np.sum([x*y for x,y in zip(v1.coordinates, v2.coordinates)])
	except:
		raise Exception('The dimensions of the two vectors mut be compatible')

def getAngle(v1, v2, radians=True):
	dp = dotProduct(v1, v2)
	mags = v1.getMagnitude()*v2.getMagnitude()
	if np.isnan(np.arccos(dp/mags)):
		return 0
	if (radians):
		return np.arccos(dp/mags)
	else:
		return (np.arccos(dp/mags)*(180/np.pi))

def areParallel(v1, v2):
	if abs(np.round(dotProduct(v1, v2), 2)) == np.round(v1.getMagnitude()*v2.getMagnitude(), 2):
		return True
	else:
		return False

def areOrthogonal(v1, v2):
	if np.round(dotProduct(v1, v2), 5) == 0:
		return True
	else:
		return False

def getProjection(v1, v2):
	return v2.normalize()*(dotProduct(v1, v2.normalize()))

def getOrthoComponent(v1, v2):
	return v1 - getProjection(v1, v2)

def crossProduct(v1, v2):
	if v1.dimension != 3 or v2.dimension != 3:
		raise TypeError("Cross product only makes sense in three dimensions")
	dims = ()
	dims += (v1.coordinates[1]*v2.coordinates[2] - v1.coordinates[2]*v2.coordinates[1],)
	dims += (v2.coordinates[0]*v1.coordinates[2] - v1.coordinates[0]*v2.coordinates[2],)
	dims += (v1.coordinates[0]*v2.coordinates[1] - v2.coordinates[0]*v1.coordinates[1],)
	return Vector(dims)

def parallelogramArea(v1, v2):
	return crossProduct(v1, v2).getMagnitude()

def triangleArea(v1, v2):
	return 0.5*parallelogramArea(v1, v2)

def main():

	vec1 = Vector([8.462, 7.893, -8.187])
	vec2 = Vector([6.984, -5.975, 4.778])
	print crossProduct(vec1, vec2)

	vec1 = Vector([-8.987, -9.838, 5.031])
	vec2 = Vector([-4.268, -1.861, -8.866])
	print parallelogramArea(vec1, vec2)

	vec1 = Vector([1.5, 9.547, 3.691])
	vec2 = Vector([-6.007, 0.124, 5.772])
	print triangleArea(vec1, vec2)



if __name__ == '__main__':
	main()