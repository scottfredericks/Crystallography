'''Program for generation of random crystal structures.
by Scott Fredericks, Spring 2018
Given a space group number between 1 and 230,
and a number N of atoms in the primitive cell,
produces a crystal structure with random atomic coordinates.
Outputs a cif file with conventional setting'''

import pymatgen as pmg

#Define variables
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
tol = 1.0 #seperation tolerance in Angstroms
num_attempts1 = 3
num_attempts2 = 10

#Choose random lattice parameters consistent with the lattice type
#Uses the conventional setting; pymatgen handles this with the Structure class
def choose_lattice(space_group_num, N):
	sg = space_group_num
	#Triclinic
	#P lattice
	if sg <= 2:
		pass 
	#Monoclinic
	elif sg <= 15:
		#C lattice
		if sg in [5, 8, 9, 12, 15]:
			pass
		#P lattice
		else:
			pass
	#Orthorhombic
	elif sg <= 74:
		#C lattice
		if sg in [20, 21, 35, 36, 37, 63, 64, 65, 66, 67, 68]:
			pass
		#F lattice
		elif sg in [22, 42, 43, 69, 70]:
			pass
		#I lattice
		elif sg in [23, 24, 44, 45, 46, 71, 72, 73, 74]:
			pass
		#A lattice
		elif sg in [38, 39, 40, 41]:
			pass
		#P lattice
		else:
			pass
	#Tetragonal
	elif sg <= 142:
		#I lattice
		if sg in [79, 80, 82, 87, 88, 97, 98, 107, 108, 109, 110, 119, 120, 121, 122, 139, 140, 141, 142]:
			pass
		#P lattice
		else:
			pass
	#Trigonal/Rhombohedral/Hexagonal
	elif sg <= 194:
		#Rhombohedral (R) lattice
		if sg in [146, 148, 155, 160, 161, 166, 167]:
			pass
		#P lattice
		else:
			pass
	#Cubic
	else:
		#F lattice
		if sg in [196, 202, 203, 209, 210, 216, 219, 225, 226, 227, 228]:
			s = (4*N) ** (1./3.)
		#I lattice
		elif sg in [197, 199, 204, 206, 211, 214, 217, 220, 229, 230]:
			s = (2*N) ** (1./3.)
		#P lattice
		else:
			s = N ** (1./3.)
		a, b, c = s, s, s
		alpha, beta, gamma = 90., 90., 90.
			

#Create structure with pymatgen
#primitive_lattice = pmg.core.lattice.from_lengths_and_angles()
#pmg.core.structure.Structure.from_spacegroup(sg, lattice, species, coords)
	#Add largest Wyckoff positions first (n<num_atoms)
		#Check distance with periodic_site, merge if needed
	#If n == num_atoms, success
	#else, try again up to num_attempts_2 times
	#If still failed, start over up to num_attempts_1 times

#If unsuccessful, output error message
#If successful, output Cif file (convert from primitive cell if needed)
