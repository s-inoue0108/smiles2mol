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

import argparse as ap
import fnmatch
from smiles2mol import __version__
from .process import process

# extension validator
def ext_validator(ptns):
    def extension(val):
        for ptn in ptns:
            if fnmatch.fnmatch(val, ptn):
                return val
        raise ap.ArgumentTypeError(f"{val} is an invalid extension file.")
    return extension
    
# negative values validator
def neg_validator(val):
    fval = int(val)
    if fval >= 0:
        return fval
    raise ap.ArgumentTypeError(f"{val!r} is not positive.")

# main function (command line parser)
def main():
    ptns_out = ["*.mol", "*.mol2", "*.pdb", "*.xyz", "*.com", "*.gjf", "*.orcainp", "*.svg"]
    parser = ap.ArgumentParser(
        prog="smiles2mol",
        description=(
            f"Chemical structure generator from SMILES strings.\n\n"
            "----------------------------------------------------\n"
            "author  : Shota Inoue\n"
            "license : MIT\n"
            "guide   : https://github.com/s-inoue0108/smiles2mol\n"
            "----------------------------------------------------"
        ),
        formatter_class=ap.RawTextHelpFormatter,
    )
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("-s", "--smiles", required=True, help="Input SMILES strings.")
    parser.add_argument("-o", "--output", default="molecule.mol", type=ext_validator(ptns_out), help="Output file name.")
    parser.add_argument("--iter", type=neg_validator, default=1000, help="Number of optimize iteration.")
    parser.add_argument("--render", action="store_true", help="Render chemical structure instead of export.")

    # execute
    args = parser.parse_args()
    process(args)

if __name__ == "__main__":
    main()
