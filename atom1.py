import numpy as np
from pyscf import gto, scf, mcscf,dft
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, PyMol



m = Chem.MolFromSmiles('CCC')
m = Chem.AddHs(m)
#print (Chem.MolToMolBlock(m))
AllChem.EmbedMolecule(m)
AllChem.MMFFOptimizeMolecule(m)

bloc = Chem.MolToMolBlock(m).split("\n")


init_coords = ""
for ix, line in enumerate(bloc[4:m.GetNumAtoms()+4]):
    spl = line.split()
    init_coords += spl[3] + " " + spl[0]+ " " + spl[1]+ " "+ spl[2] + "; \n"

#v = PyMol.MolViewer()
#v.ShowMol(ibuH)
#print (init_coords)

mol = gto.M(atom= init_coords, basis = 'sto-3g')
mf = scf.RHF(mol).run()
#mc = mcscf.CASSCF(mf, 2 ,2)
#mc.kernel()


