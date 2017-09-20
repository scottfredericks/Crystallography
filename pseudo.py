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
	Define a class "data" for structure information
	Define __init__ for data using cif file as input
	Define __init__ for data using manual user input
-----------------------------------------------------'''


import sys
import os.path
import numpy as np
import CifFile as cf

'''Class for storing a set of symmetry group operations.'''
#class sgroup:

'''Function for transforming a position to a new position.
Takes in a coordinate triplet (a), and a transformation matrix-column
pair (b). Outputs a new coordinate triplet.'''
#def transform(a, b):
	
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

#Create a 2d tree (tree) of CifDic objects for our set of pseudosymmetric groups.
tree = [[]]

#Obtain a data object, then set it as the first element of tree
#If CIF file passed as the first argument, set it as 1st tree element
if len(sys.argv) > 1:
	if os.path.isfile(sys.argv[1]):
		file = cf.ReadCif(sys.argv[1])
		#tree[0].append(...) #We need to convert the file into data
	else:
		print("Error: Could not read file.\n---Closing program.---")
		sys.exit()
else:
	#Let the user choose a CIF file in-program
	choice = input("Type '1' to choose a CIF file, or type '2' to enter structure data manually:\n")
	if choice == "1":
		usrpath = input("Enter the file's path:\n")
		if os.path.isfile(usrpath):
			file = cf.ReadCif(usrpath)
			#tree[0].append(...
		else:
			print("Error: Could not read file.\n---Closing program.---")
			sys.exit()
	elif choice == "2":
		#let the user input structure data directly
	#OR
	#Allow input from the user

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
print("---Closing program---")
