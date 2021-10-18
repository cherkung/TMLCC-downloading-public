import sys
import re
import copy
import string
import numpy as np
import pandas as pd
import math
from ase.io import read
from operator import itemgetter
from itertools import groupby
from icecream import ic
import os
from pathlib import Path

# default global numbers
E = 1e-3
PI = np.pi
angstrom = 1e-10

vdw_radii = {'H': 1.1, 'He': 1.4, 'Li': 1.82, 'Be': 1.53, 'B': 1.92, 'C': 1.7, 'N': 1.55,
             'O': 1.52, 'F': 1.47, 'Ne': 1.54, 'Na': 2.27, 'Mg': 1.73, 'Al': 1.84, 'Si': 2.1,
             'P': 1.8, 'S': 1.8, 'Cl': 1.75, 'Ar': 1.88, 'K': 2.75, 'Ca': 2.31, 'Sc': 2.15, 'Ti': 2.11,
             'V': 2.07, 'Cr': 2.06, 'Mn': 2.05, 'Fe': 2.04, 'Co': 2.0, 'Ni': 1.97, 'Cu': 1.96, 'Zn': 2.01,
             'Ga': 1.87, 'Ge': 2.11, 'As': 1.85, 'Se': 1.9, 'Br': 1.85, 'Kr': 2.02, 'Rb': 3.03, 'Sr': 2.49,
             'Y': 2.32, 'Zr': 2.23, 'Nb': 2.18, 'Mo': 2.17, 'Tc': 2.16, 'Ru': 2.13, 'Rh': 2.1, 'Pd': 2.1,
             'Ag': 2.11, 'Cd': 2.18, 'In': 1.93, 'Sn': 2.17, 'Sb': 2.06, 'Te': 2.06, 'I': 1.98, 'Xe': 2.16,
             'Cs': 3.43, 'Ba': 2.68, 'La': 2.43, 'Ce': 2.42, 'Pr': 2.4, 'Nd': 2.39, 'Pm': 2.38, 'Sm': 2.36,
             'Eu': 2.35, 'Gd': 2.34, 'Tb': 2.33, 'Dy': 2.31, 'Ho': 2.3, 'Er': 2.29, 'Tm': 2.27, 'Yb': 2.26,
             'Lu': 2.24, 'Hf': 2.23, 'Ta': 2.22, 'W': 2.18, 'Re': 2.16, 'Os': 2.16, 'Ir': 2.13, 'Pt': 2.13,
             'Au': 2.14, 'Hg': 2.23, 'Tl': 1.96, 'Pb': 2.02, 'Bi': 2.07, 'Po': 1.97, 'At': 2.02, 'Rn': 2.2,
             'Fr': 3.48, 'Ra': 2.83, 'Ac': 2.47, 'Th': 2.45, 'Pa': 2.43, 'U': 2.41, 'Np': 2.39, 'Pu': 2.43,
             'Am': 2.44, 'Cm': 2.45, 'Bk': 2.44, 'Cf': 2.45, 'Es': 2.45, 'Fm': 2.45, 'Md': 2.46, 'No': 2.46,
             'Lr': 2.46}


def get_box_parameters(cifname):
    content = open(cifname, 'r').readlines()
    for line in content:
        if "_cell_length_a" in line:
            a = float(line.split()[1])
        if "_cell_length_b" in line:
            b = float(line.split()[1])
        if "_cell_length_c" in line:
            c = float(line.split()[1])
        if "_cell_angle_alpha" in line:
            alpha_r = float(line.split()[1]) * PI / 180.0
        if "_cell_angle_beta" in line:
            beta_r = float(line.split()[1]) * PI / 180.0
        if "_cell_angle_gamma" in line:
            gamma_r = float(line.split()[1]) * PI / 180.0
    return a, b, c, alpha_r, beta_r, gamma_r


def get_box_matrix(cifname):
    a, b, c, alpha_r, beta_r, gamma_r = get_box_parameters(cifname)
    a_x = a
    a_y = 0.0
    a_z = 0.0
    b_x = b * np.cos(gamma_r)
    b_y = b * np.sin(gamma_r)
    b_z = 0.0
    c_x = c * np.cos(beta_r)
    c_y = (b * c * np.cos(alpha_r) - b_x * c_x) / b_y
    c_z = np.sqrt(c ** 2 - c_x ** 2 - c_y ** 2)
    cell_matrix = np.matrix([
        [a_x, a_y, a_z],
        [b_x, b_y, b_z],
        [c_x, c_y, c_z]],
        dtype=float)
    crs = np.cross([b_x, b_y, b_z], [c_x, c_y, c_z])
    inverse_matrix = cell_matrix.I
    unit_cell_volume = a_x * crs[0] + a_y * crs[1] + a_z * crs[2]
    return cell_matrix, inverse_matrix, unit_cell_volume


def get_atomic_coordinates(cifname):
    ase_atoms = read(cifname)
    return ase_atoms.get_chemical_symbols(), ase_atoms.get_positions()


def RandomNumberOnUnitSphere():
    thetha = 0.0
    phi = 0.0
    theta = 2 * PI * np.random.random_sample()
    phi = np.arccos(2 * np.random.random_sample() - 1.0)
    x = np.cos(theta) * np.sin(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(phi)
    return x, y, z


def _dot_product(unit_cell_matrix, atom_positions):
    s_x = unit_cell_matrix.item(0) * atom_positions[0] + unit_cell_matrix.item(3) * atom_positions[
        1] + unit_cell_matrix.item(6) * atom_positions[2]
    s_y = unit_cell_matrix.item(1) * atom_positions[0] + unit_cell_matrix.item(4) * atom_positions[
        1] + unit_cell_matrix.item(7) * atom_positions[2]
    s_z = unit_cell_matrix.item(2) * atom_positions[0] + unit_cell_matrix.item(5) * atom_positions[
        1] + unit_cell_matrix.item(8) * atom_positions[2]
    new_coord = np.array([s_x, s_y, s_z], float)
    return new_coord


def ApplyBoundaryConditions(vec2, pos, cell_matrix, inverse_matrix):
    w = [0, 0, 0]
    x = [0, 0, 0]
    vec2 = copy.deepcopy(np.asarray(vec2))
    inverse_matrix2 = copy.deepcopy(np.asarray(inverse_matrix))
    fractional = _dot_product(inverse_matrix, vec2)
    # apply boundary conditions
    x[0] = fractional[0] - np.rint(fractional[0])
    x[1] = fractional[1] - np.rint(fractional[1])
    x[2] = fractional[2] - np.rint(fractional[2])
    cartesian = _dot_product(cell_matrix, x)
    return cartesian


def CheckSurfaceAreaOverlap(pAtom, pos, cell_matrix, inverse_matrix, atom_type, vdW_pAtom, fAtom_b):
    well_depth_factor = 1.0
    fAtom_o = [0, 0, 0]
    distance = [0, 0, 0]
    # start enumerating all atoms in the object
    for i, elem in enumerate(pos):
        fAtom_o = copy.deepcopy(elem)
        if not np.array_equal(fAtom_o, fAtom_b):
            vdW_fAtom = vdw_radii.get(atom_type[i])
            equilibrium_distance = well_depth_factor * 0.5 * (vdW_pAtom + vdW_fAtom)
            fAtom_o = _dot_product(cell_matrix, fAtom_o)
            distance[0] = pAtom[0] - fAtom_o[0]
            distance[1] = pAtom[1] - fAtom_o[1]
            distance[2] = pAtom[2] - fAtom_o[2]
            dr = ApplyBoundaryConditions(distance, pos, cell_matrix, inverse_matrix)
            rr = (dr[0] * dr[0]) + (dr[1] * dr[1]) + (dr[2] * dr[2])
            if rr < (equilibrium_distance * equilibrium_distance):
                return True

    return False


def vsa(atom_type, pos, cell_matrix, inverse_matrix, unit_cell_volume, probe_diameter=1.0, nSample=20, total=0,
        equilibrium_distance=0.0):
    # Computes the volumetric surface area (VSA) of atoms. Return VSA in m^2/cm3.
    #
    # atom_type : type of atom
    # pos : coordinates
    # cell_matrix : Crystall cell matrix (A, B, C)
    # probe_diameter : Diameter of a probe
    # nSample : Sample probes
    # counted = 0.0
    vec = np.zeros((1, 3), float)
    SurfaceAreaAverage = 0.0
    well_depth_factor = 1.0
    fractional_list = []

    # start enumerating all atoms in the framework
    for i, elem in enumerate(atom_type):
        total = 0.0
        counted = 0.0
        vdW_fAtom = vdw_radii.get(elem)
        # Lorentz rule of mixing hard spheres
        equilibrium_distance = well_depth_factor * 0.5 * (probe_diameter + vdW_fAtom)

        # start MC sampling for a given #
        for attempt in range(nSample):
            pAtom = np.array([0.0, 0.0, 0.0], float)
            total += 1
            # find a random number around the sphere
            vec = RandomNumberOnUnitSphere()
            # calculate the coordinates of the center of probe using vec and equilibrium distance
            fAtom_xyz = _dot_product(cell_matrix, pos[i])
            pAtom[0] = fAtom_xyz[0] + vec[0] * equilibrium_distance
            pAtom[1] = fAtom_xyz[1] + vec[1] * equilibrium_distance
            pAtom[2] = fAtom_xyz[2] + vec[2] * equilibrium_distance

            # check for the overlap between probe and the other framework atoms
            overlap = CheckSurfaceAreaOverlap(pAtom, pos, cell_matrix, inverse_matrix, atom_type, probe_diameter, pos[i])
            if not overlap:
                counted += 1

        tmp = (counted / total) * 4.0 * PI * (equilibrium_distance * equilibrium_distance)

        SurfaceAreaAverage += tmp
        fractional_list.append([elem[0], counted / total, tmp, 1.0e4 * tmp / unit_cell_volume])

    return 1.0e4 * SurfaceAreaAverage / unit_cell_volume, SurfaceAreaAverage, atom_type


def get_surface(cifname):
    cell_matrix, inverse_matrix, unit_cell_volume = get_box_matrix(cifname)
    atom_type, pos = get_atomic_coordinates(cifname)
    vol, surface_area, a = vsa(atom_type, pos, cell_matrix, inverse_matrix, unit_cell_volume, probe_diameter=3.72, nSample=200)
    return vol, surface_area, a


def main():
    # Split the whole calculation in to a smaller part, defined by round_number and max_round
    # Each round except the last round will have equal number of samples to work with, with the last round is the remainder
    repo_folder = Path(os.path.realpath(__file__)).parents[2]
    data = pd.read_csv(os.path.join(repo_folder, "data", "train.csv"))
    try:
        if int(sys.argv[2]) == 1:
            raise IndexError
        round_number = int(sys.argv[1])
        max_round = int(sys.argv[2])
        if(round_number != max_round):
            round_start_idx = len(data.index) // (max_round-1) * (round_number-1)
            round_end_idx = round_start_idx + len(data.index) // (max_round-1) - 1
        else:
            round_start_idx = len(data.index) - len(data.index) // (max_round-1)
            round_end_idx = len(data.index) - 1
        pass
    except IndexError:
        round_number = 0
        round_start_idx = 0
        round_end_idx = len(data.index) - 1
    surface_data_presort = []

    for i in range(round_start_idx, round_end_idx):
        mof_name = "mof_unit_{}.cif".format(i+1)
        ic(mof_name)
        mof_cif = os.path.join(repo_folder, "data/mof_cif_train/", mof_name)
        vol, surface_area, a = get_surface(mof_cif)

        # add to dictionary
        surface_data_presort.append({"MOFname": mof_name, "vsa": [vol], "surface_area": [surface_area], "atom": [a]})
    surface_data = pd.DataFrame(surface_data_presort).sort_values(by='MOFname')
    surface_data.to_csv(os.path.join(repo_folder, "data", "clean_{}.csv".format(round_number)))


# if __name__ == '__main__':
#     main()
