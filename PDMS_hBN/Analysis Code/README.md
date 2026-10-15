
# Molecular Dynamics Data Analysis

This folder contains Python scripts developed for post-processing LAMMPS simulation outputs obtained from non-equilibrium molecular dynamics (NEMD) simulations of PDMS-based nanocomposites containing two-dimensional fillers.

The scripts support the extraction and processing of temperature and energy-exchange data used in thermal conductivity calculations.

## 1. Temperature Profile Analysis

**Script:** `Temperature_Profile.py`

This Python script processes spatially resolved temperature data generated during NEMD simulations.

The analysis includes:

- Reading LAMMPS temperature-profile output files.
- Extracting temperature data within a selected timestep interval.
- Calculating average temperatures for individual spatial bins.
- Generating processed temperature data for subsequent temperature-gradient analysis.

The processed results can be used to construct steady-state temperature profiles and determine temperature gradients through linear fitting.

## 2. Heat Flux Data Processing

**Script:** `Flux_Profile.py`

This Python script extracts and processes thermostat energy-exchange data generated during NEMD simulations.

The analysis includes:

- Reading LAMMPS energy-exchange output files.
- Filtering data within a specified timestep range.
- Extracting energy-exchange information associated with hot and cold regions.
- Exporting the filtered data to Excel format for further analysis.

The processed energy-exchange data can be used to evaluate the steady-state energy transfer rate and heat flux.

## 3. Thermal Conductivity Evaluation

The effective thermal conductivity is calculated using Fourier's law by combining the steady-state heat flux and the temperature gradient.

The Python scripts in this directory provide the processed data needed for subsequent thermal conductivity calculations.

## Programming Tools

- Python
- NumPy
- Pandas
- OpenPyXL

## Research Application

These analysis scripts support investigations of thermal transport in PDMS nanocomposites containing h-BN, graphene, graphene/h-BN hybrid structures (BNC), and carbon-doped h-BN fillers.
