'''
Replacement file for pseuso. Now using pymatgen.
'''

import pymatgen as mg
import getdata

struct1 = getdata.fromFile("test.cif")
print(struct1)
