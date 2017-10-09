'''
by Scott Fredericks 2017
Module for importing structure data from various input types.
Using pymatgen (pymatgen.org)
'''
import pymatgen as mg

def fromFile(path):
	return mg.core.structure.IStructure.from_file(path, primitive=False, sort=False, merge_tol=0.0)
