# smiles2mol

## Installation

```bash
pip install smiles2mol
```

## Usage

### Generate Chemical Structure

Specify a SMILES string with the `-s` or `--smiles` option.

```bash
# e.g.) Ethanol
smiles2mol -s "CCO" -o ethanol.mol
```

The following formats can be output:

| Format | Explanation |
| :----- | :---------- |
| `mol`  | MDL mol file. |
| `mol2` | Sybyl mol2 file. |
| `pdb`  | Protein data bank file. |
| `gjf` or `com` | Gaussian input file. |
| `orcainp` | ORCA input file. |
| `svg` | 2D vector graphics. |

### Render Chemical Structure

Specifying `--render` renders the 3D chemical structure using Matplotlib.

```bash
smiles2mol -s "CCO" --render
```

### Other options

```bash
- `--iter`: Specify number of iterations for structural optimization.
- `-v` or `--version`: Display version.
- `-h` or `--help`: Display help.
```
