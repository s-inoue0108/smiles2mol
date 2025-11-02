import os
import sys
from importlib.resources import files
import numpy as np
import networkx as nx
from matplotlib import font_manager as fm
import matplotlib.pyplot as plt
import matplotlib.colors as mcl

from .structure import Structure

class Renderer:
    def __init__(self, rdmol):
        self.rdmol = rdmol
        
    # visualize
    def run(self):
        atom_colors = {
            "H": "lightgray",
            "B": "pink",
            "C": "gray",
            "N": "blue",
            "O": "red",
            "F": "cyan",
            "Si": "lightblue",
            "P": "darkorange",
            "S": "yellow",
            "Cl": "green",
            "Se": "sandybrown",
            "Br": "salmon",
            "I": "blueviolet",
            "Sc": "#E6E6F2",
            "Ti": "#B3B3CC",
            "V": "#A680B3",
            "Cr": "#8A99C7",
            "Mn": "#9C4FA3",
            "Fe": "orange",
            "Co": "#4D80FF",
            "Ni": "#4FB34F",
            "Cu": "#C78033",
            "Zn": "#7D7FB0",
            "Y": "#94FFFF",
            "Zr": "#94E0E0",
            "Nb": "#73C2C2",
            "Mo": "#619999",
            "Tc": "#4D8080",
            "Ru": "#4069E0",
            "Rh": "#4069E0",
            "Pd": "silver",
            "Ag": "#BFBFBF",
            "Cd": "#57A6C2",
            "La": "#70B0F2",
            "Hf": "#4D8080",
            "Ta": "#4D8080",
            "W": "#1F1F1F",
            "Re": "#404040",
            "Os": "#4D4D4D",
            "Ir": "#808080",
            "Pt": "silver",
            "Au": "gold",
            "Hg": "#B3B3B3",
            "Tl": "#A8534D",
            "Pb": "#565656",
            "Bi": "#9E4F48",
            "Po": "#993333",
            "At": "#760076",
            "Rn": "#3FDFFF",
            "*": "purple",
        }
        
        background = "#ffffff"
        foreground = "#000000"
        negative = "#9900fa"
        positive = "#152eff"

        arial_path = files("smiles2mol.assets") / "Arial.ttf"
        arial = fm.FontProperties(fname=arial_path)
        
        plt.rcParams["font.family"] = arial.get_name()
        plt.rcParams["mathtext.fontset"] = "cm"

        plt.rcParams["figure.facecolor"] = background
        plt.rcParams["axes.facecolor"]   = background
        plt.rcParams["text.color"] = foreground
        plt.rcParams["axes.labelcolor"] = foreground
        plt.rcParams["xtick.color"] = foreground
        plt.rcParams["ytick.color"] = foreground
        
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection="3d")
        fig.canvas.manager.set_window_title(f"smiles2mol | renderer")
        
        # get graph
        sc = Structure(self.rdmol)
        graph = sc.get_graph(self.rdmol)

        # get nodes coordinate
        pos = nx.get_node_attributes(graph, "pos")
        symbols = nx.get_node_attributes(graph, "symbol")

        # draw nodes
        for i, (x, y, z) in pos.items():
            size = 100 if symbols[i] == "H" else 300
            fsize = 16 if symbols[i] == "H" else 24
            color = atom_colors[symbols[i]]
            
            # display node
            ax.scatter(x, y, z, s=size, c=color)

        # draw edges
        nodes = graph.nodes(data=True)
        for i, j, data in graph.edges(data=True):
            symbol_i = nodes[i]["symbol"]
            symbol_j = nodes[j]["symbol"]
            
            x = [pos[i][0], pos[j][0]]
            y = [pos[i][1], pos[j][1]]
            z = [pos[i][2], pos[j][2]]
            
            if data["order"] == 3:
                ax.plot(x, y, z, color="gray", linewidth=7.5)
                bv = np.array([x[1] - x[0], y[1] - y[0], z[1] - z[0]])
                
                # normal vec
                tmp = np.array([0, 0, 1]) if not np.allclose(bv[:2], 0) else np.array([0, 1, 0])
                u = np.cross(bv, tmp)
                u /= np.linalg.norm(u)
                offset = 0.05
                for sft in [-offset, offset]:
                    s = np.array((x[0], y[0], z[0])) + sft * u
                    e = np.array((x[1], y[1], z[1])) + sft * u
                    ax.plot([s[0], e[0]], [s[1], e[1]], [s[2], e[2]], color=background, linewidth=1.5)
            elif data["order"] == 2:
                ax.plot(x, y, z, color="gray", linewidth=4.5)
                ax.plot(x, y, z, color=background, linewidth=1.5)
            elif data["order"] == 1.5:
                ax.plot(x, y, z, color="gray", linestyle=(0, (1, 1)), linewidth=4.5)
                ax.plot(x, y, z, color=background, linewidth=1.5)
            else:
                ax.plot(x, y, z, color="gray", linewidth=2)
        
        x_max = np.max([item[1][0] for item in pos.items()])
        x_min = np.min([item[1][0] for item in pos.items()])
        y_max = np.max([item[1][1] for item in pos.items()])
        y_min = np.min([item[1][1] for item in pos.items()])
        z_max = np.max([item[1][2] for item in pos.items()])
        z_min = np.min([item[1][2] for item in pos.items()])
            
        ax.set_axis_off()
        ax.set_box_aspect([x_max - x_min, y_max - y_min, z_max - z_min])
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_zlim(z_min, z_max)
        plt.tight_layout()
        
        plt.show()
