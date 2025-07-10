"""
burbuja.py

Manage all stages of Burbuja, including use as a command-line tool or
as the API.
"""

import time
import pathlib
import argparse

import mdtraj

import Burbuja.modules.base as base
import Burbuja.modules.structures as structures

def burbuja(
        structure: mdtraj.Trajectory,
        grid_resolution: float = 0.1,
        use_cupy: bool = False
        ) -> structures.Bubble:
    """
    Perform bubble detection and analysis on the structure.
    """

    # TODO: make structure able to take a list or a structure object from mdtraj
    #  the list being for either multiple frames or pieces of a very large structure.
    lengths = base.reshape_atoms_to_orthorombic(structure)
    box_grid = structures.Grid(
        approx_grid_space=grid_resolution,
        boundaries=lengths)
    box_grid.initialize_cells(use_cupy=use_cupy)
    box_grid.calculate_cell_masses(structure)
    box_grid.calculate_densities_cuda()
    bubble = box_grid.generate_bubble_object()
    return bubble

def has_bubble(
        structure: mdtraj.Trajectory,
        grid_resolution: float = 0.1,
    ) -> bool:
    """
    Check if the structure has a bubble based on density threshold.
    """
    bubble = burbuja(structure, grid_resolution)
    if bubble.total_bubble_volume > base.MINIMUM_BUBBLE_VOLUME:
        return True
    else:
        return False

def main():
    argparser = argparse.ArgumentParser(
        description="Automatically detect bubbles and vapor pockets and local "
            "voids within molecular dynamics simulation structures making use "
            "of explicit solvent.")
    argparser.add_argument(
        "structure_file", metavar="STRUCTURE_FILE", type=str, 
        help="Path to the input structure or trajectory file - most likely a "
        "single-frame PDB file.")
    argparser.add_argument(
        "-t", "--topology", dest="topology",
        metavar="TOPOLOGY", type=str, default=None,
        help="Path to the topology file used by mdtraj to characterize atoms. "
        "Default: None.")
    argparser.add_argument(
        "-r", "--grid_resolution", dest="grid_resolution",
        metavar="GRID_RESOLUTION", type=float, default=0.1,
        help="Grid resolution in nanometers for the bubble detection. "
        "Default: 0.1 nm.")
    args = argparser.parse_args()
    args = vars(args)
    structure_file = pathlib.Path(args["structure_file"])
    topology = pathlib.Path(args["topology"]) if args["topology"] else None
    grid_resolution = args["grid_resolution"]
    # TODO: replace with a function to detect file size and load into multiple mdtrajs if necessary.
    if topology is None:
        structure = mdtraj.load(structure_file)
    else:
        structure = mdtraj.load(structure_file, top=topology)
    #bubble = burbuja(structure)
    has_bubble_result = has_bubble(structure, grid_resolution)
    if has_bubble_result:
        print("The structure has a bubble.")
    else:
        print("No bubble detected in structure.")

if __name__ == "__main__":
    main()