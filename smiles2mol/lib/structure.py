import os
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import networkx as nx

class Structure:
    def __init__(self, rdmol):
        self.rdmol = rdmol
        
    def generate3Drdmol(self, nconfs=10, seed=1):
        params = AllChem.ETKDGv3()
        params.randomSeed = seed
        rdmol = self.rdmol

        confs = AllChem.EmbedMultipleConfs(rdmol, numConfs=nconfs, params=params)

        # Optimization with MMFF
        for cid in confs:
            AllChem.MMFFOptimizeMolecule(rdmol, confId=cid)

        # sort with energy
        props = []
        for cid in confs:
            mmff = AllChem.MMFFGetMoleculeForceField(rdmol, AllChem.MMFFGetMoleculeProperties(rdmol), confId=cid)
            props.append((cid, mmff.CalcEnergy()))
        props.sort(key=lambda x: x[1])
        best_cid = props[0][0]
        
        # generate rdmol object
        mol_cid = Chem.Mol(rdmol)
        best_conf = rdmol.GetConformer(best_cid)
        mol_cid.RemoveAllConformers()
        mol_cid.AddConformer(best_conf, assignId=True)
        return mol_cid

    def get_coords(self, rdmol):
        coords = []
        conf = rdmol.GetConformer()
        for atom in rdmol.GetAtoms():
            idx = atom.GetIdx()
            symbol = atom.GetSymbol()
            pos = conf.GetAtomPosition(idx)
            coords.append([symbol, pos.x, pos.y, pos.z])
        return coords
        
    def get_bonds(self, rdmol):
        bonds = []
        for bond in rdmol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            order = bond.GetBondTypeAsDouble()
            bonds.append([order, (a1, a2)])
        return bonds
    
    def get_xyzstr(self, rdmol):
        coords = self.get_coords(rdmol)
        strs = ["\t".join([c[0], f"{c[1]:.6f}", f"{c[2]:.6f}", f"{c[3]:.6f}"]) for c in coords]
        return "\n".join([f"{len(coords)}", ""] + strs)
    
    def get_graph(self, rdmol):
        graph = nx.Graph()
        coords = self.get_coords(rdmol)
        bonds = self.get_bonds(rdmol)
        
        pos_list = np.array([[c[1], c[2], c[3]] for c in coords])

        # add node
        for i, (s, x, y, z) in enumerate(coords):
            symbol = rdmol.GetAtomWithIdx(i).GetSymbol()
            graph.add_node(i, symbol=s, pos=np.array([x, y, z]))

        # add edge
        nodes = graph.nodes(data=True)
        for b in bonds:
            i = b[1][0]
            j = b[1][1]
            dist = np.linalg.norm(graph.nodes[i]["pos"] - graph.nodes[j]["pos"])
            graph.add_edge(i, j, distance=dist, order=b[0])
        return graph
        
