'''Program for generation of random crystal structures.
by Scott Fredericks, Spring 2018
Given a space group number between 1 and 230,
and a number N of atoms in the primitive cell,
produces a crystal structure with random atomic coordinates.
Outputs a cif file with conventional setting'''

import pymatgen
from random import uniform as rand
from random import choice as choose
from math import sqrt
from math import pi
from math import sin
from math import cos
from math import acos
from math import fabs
import database.make_sitesym as make_sitesym
import database.hall as hall

#Define variables
#------------------------------
deg = 2.*pi/360. #radian conversion factor for sin and cos
tol = 1.0 #seperation tolerance in Angstroms
max1 = 10 #Attempts for generating lattices
max2 = 3 #Attempts for a given lattice
max3 = 10 #Attempts for a given Wyckoff position

#Define functions
#------------------------------
def displacement(xyz1, xyz2): #Displacement vector
	return [xyz2[0]-xyz1[0], xyz2[1]-xyz1[1], xyz2[2]-xyz1[2]]

def dsquared(xyz1, xyz2): #Distance squared
	v = displacement(xyz1, xyz2)
	return (v[0]**2 + v[1]**2 + v[2]**2)

def distance(xyz1, xyz2): #Euclidean distance
	return sqrt(dsquared(xyz1, xyz2))

def rand_coords():
	return [rand(0,1), rand(0,1), rand(0,1)]

def filter_site(v):
	w = v
	for i in range(len(w)):
		while w[i]<0: w[i] += 1
		while w[i]>=1: w[i] -= 1
	return w

#Choose random lattice parameters consistent with the lattice type
#Uses the conventional setting; pymatgen handles this with the Structure class
def choose_lattice(sg, N):
	v = 15. #volume per atom
	m = 2. #minimum lattice spacing
	#Triclinic
	if sg <= 2:
		gamma = rand(30., 150.)
		beta = rand(30., 150.)
		alpha = 0
		while (alpha < 30. or alpha > 150.):
			theta = rand(30., 150.) #alpha is restricted by beta and gamma, so an intermediary value is needed
			if fabs(cos(gamma*deg)*cos(beta*deg)+sin(gamma*deg)*cos(theta*deg)) <= 1.0:
				alpha = acos(cos(gamma*deg)*cos(beta*deg)+sin(gamma*deg)*cos(theta*deg)) / deg
		x = sqrt(1. - cos(alpha*deg)**2 - cos(beta*deg)**2 - cos(gamma*deg)**2 + 2.*cos(alpha*deg)*cos(beta*deg)*cos(gamma*deg))
		#P lattice
		a = rand(2., v*N/(m*m*x))
		b = rand(2., v*N/(a*m*x))
		c = v*N/(a*b*x)
	#Monoclinic
	elif sg <= 15:
		alpha, gamma = 90., 90.
		beta = rand(30., 150.)
		x = sin(beta*deg)
		#C lattice
		if sg in [5, 8, 9, 12, 15]:
			a = rand(2., 2.*v*N/(m*m*x))
			b = rand(2., 2.*v*N/(a*m*x))
			c = 2.*v*N/(a*b*x)
		#P lattice
		else:
			a = rand(2., v*N/(m*m*x))
			b = rand(2., v*N/(a*m*x))
			c = v*N/(a*b*x)
	#Orthorhombic
	elif sg <= 74:
		alpha, beta, gamma = 90., 90., 90.
		#C lattice
		if sg in [20, 21, 35, 36, 37, 63, 64, 65, 66, 67, 68]:
			a = rand(2., 2.*v*N/(m*m))
			b = rand(2., 2.*v*N/(a*m))
			c = 2.*v*N/(a*b)
		#F lattice
		elif sg in [22, 42, 43, 69, 70]:
			a = rand(2., 4.*v*N/(m*m))
			b = rand(2., 4.*v*N/(a*m))
			c = 4.*v*N/(a*b)
		#I lattice
		elif sg in [23, 24, 44, 45, 46, 71, 72, 73, 74]:
			a = rand(2., 2.*v*N/(m*m))
			b = rand(2., 2.*v*N/(a*m))
			c = 2.*v*N/(a*b)
		#A lattice
		elif sg in [38, 39, 40, 41]:
			a = rand(2., 2.*v*N/(m*m))
			b = rand(2., 2.*v*N/(a*m))
			c = 2.*v*N/(a*b)
		#P lattice
		else:
			a = rand(2., v*N/(m*m))
			b = rand(2., v*N/(a*m))
			c = v*N/(a*b)
	#Tetragonal
	elif sg <= 142:
		alpha, beta, gamma = 90., 90., 90.
		#I lattice
		if sg in [79, 80, 82, 87, 88, 97, 98, 107, 108, 109, 110, 119, 120, 121, 122, 139, 140, 141, 142]:
			c = rand(2., 2*v*N/(m*m))
			a = sqrt(2*v*N/c)
		#P lattice
		else:
			c = rand(2., v*N/(m*m))
			a = sqrt(v*N/c)
	#Trigonal/Rhombohedral/Hexagonal
	elif sg <= 194:
		alpha, beta, gamma = 90., 90., 120.
		x = sqrt(3.)/2.
		#Rhombohedral (R) lattice
		if sg in [146, 148, 155, 160, 161, 166, 167]:
			c = rand(2., 3.*v*N/(m*m*x))
			a = sqrt(3.*v*N/(c*x))
		#P lattice
		else:
			c = rand(2., v*N/(m*m*x))
			a = sqrt(v*N/(c*x))
			a, b, c = a, a, c
	#Cubic
	else:
		alpha, beta, gamma = 90., 90., 90.
		#F lattice
		if sg in [196, 202, 203, 209, 210, 216, 219, 225, 226, 227, 228]:
			s = (4*N*v) ** (1./3.)
		#I lattice
		elif sg in [197, 199, 204, 206, 211, 214, 217, 220, 229, 230]:
			s = (2*N*v) ** (1./3.)
		#P lattice
		else:
			s = (N*v) ** (1./3.)
		a, b, c = s, s, s
	return [[a, b, c], [alpha, beta, gamma]]

def get_wyckoff_positions(international_number):
	array = []
	hall_number = hall.hall_from_hm(international_number)
	wyckoff_positions = make_sitesym.get_wyckoff_position_operators('database/Wyckoff.csv', hall_number)
	for x in wyckoff_positions:
		temp = []
		for y in x:
			temp.append(pymatgen.core.operations.SymmOp.from_rotation_and_translation(list(y[0]), filter_site(y[1]/24)))
		array.append(temp)
	return array

def connected_components(graph): #Return a set of connected components for an undirected graph
	def add_neighbors(el, seen=[]):
		if seen == []: seen = [el]
		for x in graph[el]:
			if x not in seen:
				seen.append(x)
				add_neighbors(x, seen)
		return seen
	unseen = list(range(len(graph)))
	sets = []
	i = 0
	while (unseen != []):
		x = unseen.pop()
		sets.append([])
		for y in add_neighbors(x):
			sets[i].append(y)
			if y in unseen: unseen.remove(y)
		i += 1
	return sets

def generate_structure(space_group_num, species, num_atoms):
	wyckoffs = get_wyckoff_positions(space_group_num) #2d array of Symm_Ops
	i = 0
	j = 0
	wyckoffs_organized = [[]] #2D Array of Wyckoff positions organized by multiplicity
	old = len(wyckoffs[0])
	for x in wyckoffs:
		mult = len(wyckoffs[i])
		if mult != old:
			i += 1
			old = mult
		wyckoffs_organized[i].append(x)
	y = 0 #total number of atoms in primitive cell
	for x in num_atoms:
		y += x
	attempts1 = 0 #Num of attempted lattices
	success = False #True once successful structure has been generated
	while (attempts1 < max1 and not success):
		#Generate a lattice and empty structure object
		lat = choose_lattice(space_group_num, y)
		conventional_lattice = pymatgen.core.lattice.Lattice.from_lengths_and_angles(lat[0], lat[1])
		attempts2 = 0 #Num of attempts for a given lattice
		while (attempts2 < max2 and not success):
			struct1 = pymatgen.core.structure.Structure.from_spacegroup(space_group_num, conventional_lattice, [], [])
			#Add atoms one species at a time
			i = 0
			success_atom = True
			while (i < len(species) and success_atom):
				attempts3 = 0 #Num of attempts for a given species
				#Strangely, just using struct2 = struct1 causes methods on struct2 to affect struct1
				success_atom = False
				element = species[i]
				needed = num_atoms[i]
				total = 0
				while (not success_atom and attempts3 < max3):
					loop = True
					j = 0
					while loop:
						mult = len(wyckoffs_organized[j])
						if (needed - total <= mult):
							#Choose a random Wyckoff position for given multiplicity
							ops = choose(wyckoffs_organized[j])
							#Generate a list of coords from ops
							point = rand_coords()
							new_coords = []
							for x in ops:
								new_coords.append(x.operate(point))
							#Create a graph from new_coords
							graph = []
							struct2 = pymatgen.core.structure.Structure.from_spacegroup(1, conventional_lattice, ['H'], new_coords)
							for x in struct2.sites:
								pass
							loop = False
						else:
							j += 1
							if j >= len(wyckoffs_organized): loop = False
					attempts3 += 1
				i += 1
			if not success_atom: attempts2 += 1
		if not success: attempts1 += 1
	return struct1

#-----------------Test Functionality------------------
'''space_group_num = int(input("International number of space group: "))
if space_group_num < 1 or space_group_num > 230:
	print("Error: Invalid space group number.")
	quit()
species = input("Atomic species (separated by commas): ").split(',')
y = input("Number of atoms in the primitive cell (separated by commas): ").split(',')
num_atoms = []
for x in y:
	num_atoms.append(int(x))
if len(species) != len(num_atoms):
	print("Error: Number of atoms does not match number of species.")
	quit()'''
print(generate_structure(1, ['H'], [1]))

