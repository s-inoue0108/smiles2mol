import os
import sys
import site
import subprocess

def ensure_openbabel_wheel():
    if os.environ.get("LD_LIBRARY_PATH_INITIALIZED"):
        return

    paths = []
    paths += site.getsitepackages()
    paths += [site.getusersitepackages()]

    for p in paths:
        wheel_lib = os.path.join(p, "openbabel")
        if os.path.isdir(wheel_lib):
            os.environ["LD_LIBRARY_PATH"] = wheel_lib
            os.environ["LD_LIBRARY_PATH_INITIALIZED"] = "1"
            os.execv(sys.executable, [sys.executable] + sys.argv)
            break

ensure_openbabel_wheel()

try:
    from openbabel import pybel
except ImportError:
    import pybel
    
from rdkit import Chem

class SmilesIO:
    def __init__(self, smiles):
        self.smiles = smiles
        self.rdmol = self.gen_rdmol()
        
    def gen_rdmol(self):
        rdmol = Chem.MolFromSmiles(self.smiles)
        rdmol = Chem.AddHs(rdmol)
        return rdmol
    
    def get_rdmol(self):
        return self.rdmol

class ChemicalFormatIO:
    def __init__(self, xyzstr, output):
        self.xyzstr = xyzstr
        self.output = output
        self.filename = os.path.basename(self.output)
        self.basename = os.path.splitext(self.filename)[0]
        self.format = os.path.splitext(self.filename)[1][1:]
        
    def get_filename(self):
        return self.filename
    
    def get_basename(self):
        return self.basename
        
    def get_format(self):
        return self.format
        
    def get_rdmol(self):
        mol = pybel.readstring("xyz", self.xyzstr)
        mol_block = mol.write("mol").rstrip()
        rdmol = Chem.MolFromMolBlock(mol_block, removeHs=False)
        return rdmol
        
    def export(self):
        mol = pybel.readstring("xyz", self.xyzstr)
        mol.write(self.format, self.filename, overwrite=True)
        print(f"Generate chemical structure: {self.filename}")
        return self
