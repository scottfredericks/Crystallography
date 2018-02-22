'''Program for generation of random crystal structures.
by Scott Fredericks, Spring 2018
Given a space group number between 1 and 230,
and a number N of atoms in the primitive cell,
produces a crystal structure with random atomic coordinates.
Outputs a cif file with conventional setting'''

import pymatgen
from random import uniform as rand
from math import sqrt
from math import pi
from math import sin
from math import cos
import database.make_sitesym as make_sitesym
import database.hall as hall

#Define variables
deg = 2.*pi/360. #radian conversion factor for sin and cos
min_separation = 1.0 #seperation tolerance in Angstroms
num_attempts1 = 3
num_attempts2 = 10
space_group_num = int(input("International number of space group: "))
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
	quit()

#Choose random lattice parameters consistent with the lattice type
#Uses the conventional setting; pymatgen handles this with the Structure class
def choose_lattice(sg, N):
	v = 15. #volume per atom
	m = 2. #minimum lattice spacing
	#Triclinic
	if sg <= 2:
		alpha = rand(30., 150.)
		beta = rand(30., 150.)
		gamma = rand(30., 150.)
		x = sqrt(1. - cos(alpha*deg)**2 - cos(beta*deg)**2 - cos(gamma*deg)**2 + 2*cos(alpha*deg)*cos(beta*deg)*cos(gamma*deg))
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

def rand_coords():
	return [rand(0,1), rand(0,1), rand(0,1)]

def get_wyckoff_positions(hall_number):
	array = []
	wyckoff_positions = make_sitesym.get_wyckoff_position_operators('database/Wyckoff.csv', hall_number)
	for x in wyckoff_positions:
		temp = []
		for y in x:
			temp.append(pymatgen.core.operations.SymmOp.from_rotation_and_translation(list(y[0]), filter_site(y[1]/24)))
		array.append(temp)
	return array

#Create 1-atom structure with pymatgen
y = 0
for x in num_atoms:
	y += x
params = choose_lattice(space_group_num, y)
conventional_lattice = pymatgen.core.lattice.Lattice.from_lengths_and_angles( params[0], params[1] )

#-----------------Test Functionality------------------
struct1 = pymatgen.core.structure.Structure.from_spacegroup(space_group_num, conventional_lattice, [species[0]], [rand_coords()])
	#Add largest Wyckoff positions first (n<num_atoms)
		#Check distance with periodic_site, merge if needed
	#If n == num_atoms, success
	#else, try again up to num_attempts_2 times
	#If still failed, start over up to num_attempts_1 times

#If unsuccessful, output error message
#If successful, output Cif file (convert from primitive cell if needed)
