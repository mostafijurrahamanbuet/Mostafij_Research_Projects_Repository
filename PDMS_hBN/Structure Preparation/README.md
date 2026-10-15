
# Molecular Structure Preparation

This directory contains Python scripts used to prepare molecular structures for LAMMPS simulations of PDMS-based nanocomposites containing two-dimensional fillers.

## PDB-to-LAMMPS Data Conversion

The Python conversion script processes molecular structures in PDB format and generates LAMMPS-compatible data files.

The processing workflow includes:

1. Reading and processing atomic coordinates from PDB files.
2. Identifying PDMS and two-dimensional filler atoms.
3. Assigning atomic types, masses, and partial charges.
4. Constructing molecular topology, including bonds, angles, and dihedrals.
5. Defining simulation box dimensions.
6. Generating LAMMPS data files for subsequent molecular dynamics simulations.

**Programming Language:** Python

**Libraries:** MDAnalysis, NumPy

The scripts are intended to support molecular structure preparation prior to energy minimization and equilibration.
