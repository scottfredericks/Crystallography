import pymatgen
import spglib
import spglib._spglib as spg

mystruct = pymatgen.core.structure.Structure.from_file('test.cif')
ds = mystruct.as_dict()
for x in ds:
	print(x)
#data = spg.dataset(ds['lattice'], ds['positions'], ds['numbers'], ds['hall_number'], .001, 5)
