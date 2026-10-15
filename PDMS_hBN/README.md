
# Thermal Conductivity Enhancement in PDMS-Based Nanocomposites with 2D Fillers

**Undergraduate Thesis | Mechanical Engineering | BUET**

**Research Area:** Nanoscale Thermal Transport | Polymer Nanocomposites | Thermal Interface Materials

**Simulation Software:** LAMMPS  
**Method:** Non-Equilibrium Molecular Dynamics (NEMD)

---

## Research Overview

This research investigates the effects of filler composition, volume fraction, and morphology on the thermal conductivity of polydimethylsiloxane (PDMS)-based nanocomposites using molecular dynamics simulations in LAMMPS.

Four types of two-dimensional fillers are considered:

- Hexagonal boron nitride (h-BN)
- Graphene
- Graphene/h-BN in-plane hybrid heterostructures (BNC)
- Carbon-doped h-BN

Both united-atom (UA) and all-atom (OPLS-AA) representations of PDMS were considered during the preliminary investigation. The hybrid UA model was selected for subsequent investigations because of its reduced computational cost.

The study examines how filler loading, internal composition, and geometric configuration influence thermal transport, with the objective of identifying promising nanocomposite designs for thermal interface material (TIM) applications.

---

## Graphical Abstract

![PDMS Nanocomposite Graphical Abstract](graphical_abstract.png)

*Graphical overview of the molecular dynamics workflow for investigating thermal transport in PDMS nanocomposites containing two-dimensional fillers.*

---

## Representative Molecular Structure

<p align="center">
  <img src="Circular_BNC_PDMS_composite.png" alt="Circular graphene/h-BN BNC–PDMS nanocomposite" width="420">
</p>

**Figure:** Representative atomistic configuration of a PDMS-based nanocomposite containing graphene/h-BN hybrid filler sheets with circular graphene domains.

The composite consists of three planar filler sheets embedded within a PDMS matrix, with polymer chains distributed between and around the filler layers.

---

## Simulation Methodology

### 1. Model Construction

The initial molecular structures were generated using molecular modeling tools and packed using PACKMOL.

The investigated systems consist of PDMS polymer chains and aligned two-dimensional filler sheets.

For the comparative study at 9% filler volume fraction, three filler sheets were placed parallel to the XY plane and separated along the z-direction.

### 2. Interatomic Potentials

The molecular interactions were modeled using:

- **Tersoff potential:** Covalent interactions within h-BN, graphene, and BNC fillers.
- **Hybrid UA / OPLS-AA force fields:** Molecular interactions within PDMS.
- **Lennard-Jones potential:** Nonbonded interactions between PDMS and filler sheets.
- **PPPM method:** Long-range Coulombic interactions.

### 3. Energy Minimization and Equilibration

The initial composite structures were relaxed using conjugate-gradient energy minimization to remove unfavorable atomic contacts.

The systems were subsequently equilibrated using temperature- and pressure-controlled molecular dynamics simulations.

The equilibration procedure includes NPT densification followed by NVE relaxation.


### 4. Thermal Conductivity Calculation

The effective thermal conductivity of the PDMS-based nanocomposites was calculated using non-equilibrium molecular dynamics (NEMD) simulations in LAMMPS.

A hot region (300 K) was placed at the center of the simulation domain, while two cold regions (260 K) were positioned at opposite ends along the x-direction.


The steady-state thermal conductivity was calculated using Fourier's law:

```math
\kappa = \frac{J}{\left|\frac{dT}{dx}\right|}
```

where:

- **κ** = Effective thermal conductivity (W m⁻¹ K⁻¹).
- **J** = Steady-state heat flux (W m⁻²).
- **dT/dx** = Temperature gradient (K m⁻¹).

The heat flux was calculated from the energy transfer rate:

```math
J = \frac{1}{2A}\frac{dE}{dt}
```

The cross-sectional area perpendicular to the heat-flow direction was calculated as:

```math
A = L_y L_z
```

The factor of two accounts for heat flow in two opposite directions from the central hot region under periodic boundary conditions.


The temperature gradient was obtained by linear fitting of the steady-state temperature profiles on both sides of the hot region. The average magnitude of the two fitted slopes was used to calculate the thermal conductivity.


## Research Files

The following folders contain the simulation inputs, molecular structures, and force-field data used in this research.

| Category | Description | Files |
|---|---|---|
| Simulation Code | LAMMPS inputs for energy minimization, equilibration, and NEMD simulations | [View Code](MD%20code/) |
| Simulation Data | PDB molecular structures and LAMMPS data files | [View Data](Data%20file/) |
| Force Field | B–N–C Tersoff potential parameters | [View Force Field](Force%20field/) |
| Structure Preparation | Python script for converting PDB structures into LAMMPS data files | Upload pending |
| Data Analysis | Python scripts for processing temperature profiles and thermostat energy-exchange data | Upload pending |

---

## Data Analysis

### Temperature Profile Analysis

Python scripts are used to extract and average spatial temperature data from LAMMPS outputs.

The processed data can be used to construct temperature profiles and determine the temperature gradient through subsequent linear fitting.

### Heat Flux Analysis

Python scripts are used to extract and filter thermostat energy-exchange data from LAMMPS outputs.

The processed data can be used to evaluate the steady-state energy transfer rate and calculate heat flux.

The thermal conductivity is then obtained from the heat flux and temperature gradient using Fourier's law.

*Analysis scripts will be added to the repository.*

---

## Conference Presentation

**Title:** Molecular Dynamics Study of Thermal Conductivity Enhancement in PDMS Nanocomposites Using H-BN Based Fillers

**Conference:** ASME International Mechanical Engineering Congress and Exposition (IMECE 2026)

**Paper Number:** 198871

**Status:** Technical presentation accepted.

[**View Official Conference Abstract →**](https://imece.secure-platform.com/a/solicitations/281/sessiongallery/24922/application/198871)

**Presentation Slides:** To be uploaded after the conference presentation.

---

## Research Status

**Manuscript:** In preparation.

Simulation codes, supporting data, and analysis scripts are being organized for research documentation and reproducibility.

---

## Authors

**1. Mostafijur Rahaman**  
Department of Mechanical Engineering  
Bangladesh University of Engineering and Technology (BUET), Dhaka, Bangladesh  
*First Author*

**2. Arman Hassan**  
Department of Mechanical Engineering  
Bangladesh University of Engineering and Technology (BUET), Dhaka, Bangladesh  
*Co-Author*

**3. Prof. Dr. Md. Ashiqur Rahman**  
Department of Mechanical Engineering  
Bangladesh University of Engineering and Technology (BUET), Dhaka, Bangladesh  
*Corresponding Author*

[Email](mailto:mostafijur.rahaman.buet@gmail.com)
