# ---------------------------------------------------------------------------
# @author Shota Inoue
# @date   2025-11-02
#
# Generate 3D Molecule structure by SMILES strings.
#
# Package dependencies:
#    - NumPy
#    - RDKit
#    - OpenBabel
#    - NetworkX
# ---------------------------------------------------------------------------

from .lib.io import SmilesIO, ChemicalFormatIO
from .lib.structure import Structure
from .lib.render import Renderer

# main function
def process(args):
    smio = SmilesIO(args.smiles)
    rdmol = smio.get_rdmol()
    sc = Structure(rdmol)
    
    best_mol = sc.generate3Drdmol(nconfs=args.iter)
    xyzstr = sc.get_xyzstr(best_mol)
    
    cfio = ChemicalFormatIO(xyzstr, output=args.output)
    
    if args.render:
        rdmol = cfio.get_rdmol()
        renderer = Renderer(rdmol)
        renderer.run()
    else:
        cfio.export()
