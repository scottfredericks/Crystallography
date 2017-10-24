'''
by Scott Fredericks 2017
Module for importing structure data from various input types.
Using pymatgen (pymatgen.org)
'''
import pymatgen

def from_file(path):
	return pymatgen.core.structure.Structure.from_file(path, primitive=False, sort=False, merge_tol=0.01)
