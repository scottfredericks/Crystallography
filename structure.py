'''Program for generation of random crystal structures.
by Scott Fredericks, Spring 2018
Given a space group number between 1 and 230,
and a number N of atoms in the primitive cell,
produces a crystal structure with random atomic coordinates.'''

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

def choose_lattice(space_group_num, N):
	sg = space_group_num
	#Triclinic
	if <= 2:
		pass 
	#Monoclinic
	elif sg <= 15:
		pass
	#Orthorhombic
	elif sg <= 74:
		pass
	#Tetragonal
	elif sg <= 142:
		pass
	#Trigonal/Rhombohedral
	elif sg <= 167:
		#Rhombohedral
		if sg in [146, 148, 155, 160, 161, 166, 167]:
			pass
		#Trigonal
		else:
			pass
	#Hexagonal
	elif sg <= 194:
		pass
	#Cubic
	else:
		pass

#Create structure with pymatgen

	#Add largest Wyckoff positions first (n<num_atoms)
		#Check distance with periodic_site, merge if needed
	#If n == num_atoms, success
	#else, try again up to num_attempts_2 times
	#If still failed, start over up to num_attempts_1 times

#If unsuccessful, output error message
#If successful, output Cif file (convert from primitive cell if needed)
