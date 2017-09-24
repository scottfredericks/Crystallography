'''
Pseudosymmetry Search Program for Crystal Structures
By: Scott Fredericks 2017, UNLV Physics Department
Advisor: Dr. Qiang Zhu
Provided under MIT license
License details can be found in "Licenses" folder
Version 0.1, under continuous development
Based on "PSEUDO" from the Bilbao Crystallographic Server:
http://www.cryst.ehu.es/cryst/pseudosymmetry.html
C. Capillas, E.S. Tasci, G. de la Flor, D. Orobengoa, J.M. Perez-Mato and M.I. Aroyo. "A new computer tool at the Bilbao Crystallographic Server to detect and characterize pseudosymmetry". Z. Krist. (2011), 226(2), 186-196 DOI:10.1524/zkri.2011.1321.
---Dependencies---
Python 3.6
PYCIFRW version 4.3
	https://pypi.python.org/pypi/PyCifRW/4.3
Numpy version 1.12.1
	http://www.numpy.org/
---Current objectives---
	Class tr: Define transformation matrix-column pair from (x, y, z) format
	Convert H-M symbol into list of irreducible symm opp's
	Find the irreducible symmetry operations for a group/supergroup
-----------------------------------------------------'''


import sys
import os.path
import numpy as np
import CifFile as cf

#shorthand and constants
pi = np.pi
def cos(x): return np.cos(x)
def sin(x): return np.sin(x)
def sqrt(x): return np.sqrt(x)

'''Parameters needed to define a structure'''
expected = ["space group name", "cell length x", "cell length y",
	"cell length z", "cell angle a", "cell angle b", "cell angle c",
	"number of atoms in unit", "atomic info"]

'''corresponding variables in crystal class'''
exdic = {'space group name': "gname",
	'cell length x': "cx",
	'cell length y': "cy",
	'cell length z': "cz",
	'cell angle a': "ca",
	'cell angle b': "cb",
	'cell angle c': "cc",
	'number of atoms in unit': "nunit",
	'atomic info': "info"}

'''Store structure data for a given material'''
class crystal:
	data = {
	#Index of space group
	'gname': 0.0,
	#cell parameters: lengths x,y,z, angles a,b,c
	'cx': 0.0,
	'cy': 0.0,
	'cz': 0.0,
	'ca': 0.0,
	'cb': 0.0,
	'cc': 0.0,
	#number of atoms in unit
	'nunit': 0.0}
	label = []
	pos = [[]]
	#Add list of symmetry operations... irreducible?

'''Get user input for needed values. (x) is a crystal to be updated by user
input. Using "new" in place of x generates a new crystal structure. (needed)
is a list of strings which we will ask for, usually (expected). Outputs a new
crystal object with the provided info.'''
def askdata(x, needed):
	if x == "new": y = crystal
	else: y = x
	print("~Please enter the needed information:")
	i = 0
	while i < len(needed):
		#space group name is not a float
		if needed[i] == "space group name":
			y.data[exdic[needed[i]]] = input(needed[i] + ": ")
		#For atomic labels and positions, we must parse the input
		elif needed[i] == "atomic info":
			y.label = []
			y.pos = []
			print("~Please enter the atomic symbols, labels, and positions in the form:")
			print("Symbol # x y z")
			text = []
			textlist = []
			j = 0
			while j < y.data['nunit']:
				text.append(input(str(j+1) + ": "))
				textlist.append(text[j].split())
				j += 1
			j = 0
			while j < y.data['nunit']:
				k = 0
				while k < 5:
					if k == 0: y.label.append(str(textlist[j][k]))
					elif k == 1: y.label[j] += str(textlist[j][k])
					else: y.pos[j].append(float(textlist[j][k]))
					k += 1
				if j != (y.data['nunit']-1): y.pos.append([])
				j += 1
		#all other crystal data can be stored as floats
		elif needed[i] != "atomic info":
			y.data[exdic[needed[i]]] = float(input(needed[i] + ": "))
		i += 1
	return y

names = {'space group name': "_symmetry_space_group_name_H-M",
	'cell length x': "_cell_length_a",
	'cell length y': "_cell_length_b",
	'cell length z': "_cell_length_c",
	'cell angle a': "_cell_angle_alpha",
	'cell angle b': "_cell_angle_beta",
	'cell angle c': "_cell_angle_gamma"}

'''Create a crystal object by parsing a cif file from (path).'''
def loadcif(path):
	print("Loading file: " + path)
	x = crystal
	missing = []
	myfile = cf.ReadCif(path)
	info = myfile.first_block()
	
	for s in expected:
		if (s != "atomic info") and (s !="number of atoms in unit"):
			if names[s] in info:
				x.data[exdic[s]] = info[names[s]]
			else: missing.append(s)			
		elif s == "number of atoms in unit":
			if '_atom_site_label' in info:
				x.data['gunit'] = len(info.GetItemValue('_atom_site_label'))
		elif s == "atomic info":
			if '_atom_site_label' in info\
			and '_atom_site_fract_x' in info\
			and '_atom_site_fract_y' in info\
			and '_atom_site_fract_z' in info:
				i = 0
				while i < len(info.GetItemValue('_atom_site_fract_x')):
					x.label.append(info.GetItemValue('_atom_site_label')[i])
					x.pos[i].append(float(info.GetItemValue('_atom_site_fract_x')[i]))
					x.pos[i].append(float(info.GetItemValue('_atom_site_fract_y')[i]))
					x.pos[i].append(float(info.GetItemValue('_atom_site_fract_z')[i]))
					if i < len(info.GetItemValue('_atom_site_fract_x')) - 1:
						x.pos.append([])
					i += 1
	#If any parameters are missing, get directly from user
	if missing != []:
		print("Could not find the following: " + str(missing))
		x = askdata(x, missing)
	'''print("Entered info:") #Check that the file loaded correctly
	print(x.data)
	print(x.label)
	print(x.pos)'''
	print("Finished loading " + path)
	return x

'''Class for storing a set of symmetry group operations...?'''
#class sgroup:

'''A class for storing matrix-column transformations (W,w)
The first 3 vectors are the matrix W, and the 4th vector
represents the translation w. Can be initialized using
text string in "a+bx,c+dy,e+fz" format. We also define
some basic transformations for testing.'''
class tr:
	#Basic transformations for testing
	identity = [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0], [0.0,0.0,0.0]]
	inversion = [[-1.0,0.0,0.0],[0.0,-1.0,0.0],[0.0,0.0,-1.0], [0.0,0.0,0.0]]
	reflectx = [[-1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0], [0.0,0.0,0.0]]
	reflecty = [[1.0,0.0,0.0],[0.0,-1.0,0.0],[0.0,0.0,1.0], [0.0,0.0,0.0]]
	reflectz = [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,-1.0], [0.0,0.0,0.0]]
	#rotation by 90 degrees clockwise
	rotatex = [[1.0,0.0,0.0],[0.0,0.0,1.0],[0.0,-1.0,0.0], [0.0,0.0,0.0]]
	rotatey = [[0.0,0.0,-1.0],[0.0,1.0,0.0],[1.0,0.0,0.0], [0.0,0.0,0.0]]
	rotatez = [[0.0,1.0,0.0],[-1.0,0.0,0.0],[0.0,0.0,1.0], [0.0,0.0,0.0]]
	#rotation by 90 degrees counter-clockwise
	rotatexc = [[1.0,0.0,0.0],[0.0,0.0,-1.0],[0.0,1.0,0.0], [0.0,0.0,0.0]]
	rotateyc = [[0.0,0.0,1.0],[0.0,1.0,0.0],[-1.0,0.0,0.0], [0.0,0.0,0.0]]
	rotatezc = [[0.0,-1.0,0.0],[1.0,0.0,0.0],[0.0,0.0,1.0], [0.0,0.0,0.0]]
	#translation by a vector v
	def translate(v):
		return [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0], [float(v[0]), float(v[1]), float(v[2])]]
	#reflection about a plane defined by v
	#rotation by 2pi/n about an axis v	
	def rotate(n, v):
		if n == 0: return identity
		#rotation angle a
		a = 2.0 * pi / n
		#normalize v
		v = [float(v[0]), float(v[1]), float(v[2])]
		absv = sqrt((v[0]**2 + v[1]**2 + v[2]**2))
		v = [v[0]*absv, v[1]*absv, v[2]*absv]
		return [[cos(a)+(1-cos(a))*v[0]*v[0], (1-cos(a))*v[0]*v[1] - sin(a)*v[2], (1-cos(a))*v[0]*v[2] + sin(a)*v[1]],
			[(1-cos(a))*v[1]*v[0] + sin(a)*v[2], cos(a) + (1-cos(a))*v[1]*v[1], (1-cos(a))*v[1]*v[2] - sin(a)*v[0]],
			[(1-cos(a))*v[2]*v[0] - sin(a)*v[1], (1-cos(a))*v[2]*v[1] + sin(a)*v[0], cos(a) + (1-cos(a))*v[2]*v[2]],
			[0,0,0]]
	
	#Allow initialization from a string of the form "a+bx,c+dy,e+fz"

'''Takes in a float vector (a), and a transformation matrix-column
pair (b). Outputs a new coordinate triplet. 'l' means left, as we are
applying a transormation 'g' to a vector 'a': ga, as opposed to ag'''
#convert to numpy
def ltransform(b, a):
	for i in range(4):
		for j in range(3):
			if b[i][j] is not float: b[i][j] = float(b[i][j])
	for i in range(3):
		if a[i] is not float: a[i] = float(a[i])
	c = [float(0), float(0), float(0)]
	i = 0
	while i < 3:
		j = 0
		while j < 3:
			c[i] += b[i][j]*a[j]
			j += 1
		c[i] += b[3][i]
		i += 1
	return c

'''rtransform'''
	#convert to numpy

'''Outputs the set of minimal supergroups for a given symmetry group
Takes in a data object (a), outputs a potential symmetry operation.
Needs to use irreducible representations to avoid redundancy.'''
#def super(a):

'''Checks the relative positions of atoms for 2 data sets. Outputs a real number
representing the maximum difference (du) in positions.'''
#def compare(a, b)


#-----------------
#Main program loop
#-----------------


print("---Pseudosymmetry search program---")

'''The branching set of pseudosymmetric groups, including the user-defined input group'''
tree = [[]]

'''Obtain a data object, then set it as the first element of tree.
If CIF file passed as the first argument, set it as 1st tree element'''
if len(sys.argv) > 1:
	if os.path.isfile(sys.argv[1]):
		tree = [[loadcif(sys.argv[1])]]
	else:
		print("Error: Could not read file.\n---Closing program.---")
		sys.exit()
else:
	#Let the user choose a CIF file in-program
	choice = input("Type '1' to choose a CIF file, or type '2' to enter structure data manually:\n")
	if choice == "1":
		usrpath = input("Enter the file's path:\n")
		if os.path.isfile(usrpath):
			tree = [[loadcif(usrpath)]]
		else:
			print("Error: Could not read file.\n---Closing program.---")
			sys.exit()
	elif choice == "2":
		tree = [[askdata("new", expected)]]
#Let the user choose a maximum tolerance tol (default 1 Angstrom)
#Loop over (a) in tree:
#Loop over the supergroup list super(a):
	#dumax = 0
	#Loop over the operations gi for all i=index(G)
		#Apply the symmetry operation to the structure, output b
		#Loop over each point in b:
			#set dumax to 1st difference calculated
			#update dumax if difference > dumax
		#If dumax > tol: go to next supergroup in a
	#If dumax <= tol: add b to tree
#Once no more pseudosymmetries are found, output the elements of the tree
	#Display symmetry group, positions, dumax, index
print("Before:")
print(tree[0][0].pos)
print("After:")
tree[0].append(crystal)
for i in range(len(tree[0][0].pos)):
	tree[0][1].pos[i] = ltransform(tr.rotate(4, [1,0,0]), tree[0][0].pos[i])
print(tree[0][1].pos)

print("---Closing program---")
