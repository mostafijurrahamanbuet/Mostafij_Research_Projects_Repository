
# CO₂ Self-Diffusion in Metal–Organic Frameworks under N₂ Co-Loading and Varying Pore Sizes: A Molecular Dynamics Study

**Ongoing Research | Computational Materials Science | Molecular Diffusion in Porous Materials**

**Research Areas:** Metal–Organic Frameworks (MOFs) | Gas Transport | CO₂ Self-Diffusion | Pore-Size Effects | Molecular Dynamics

**Simulation Software:** LAMMPS  
**Method:** Molecular Dynamics (MD) | Mean-Square Displacement (MSD) Analysis

---

## Research Overview

This ongoing research investigates CO₂ self-diffusion in metal–organic frameworks (MOFs) using molecular dynamics simulations in LAMMPS.

The study aims to understand how **N₂ co-loading, framework composition, and pore-size variations** influence the molecular mobility and self-diffusion behavior of CO₂ within porous materials.

The investigated and planned framework systems include:

- **MOF-5:** Zn-based metal–organic framework.
- **Ni-based MOF:** Nickel-containing porous framework.
- **Cu-based MOF:** Additional framework considered for future pore-size comparisons.

The broader objective is to characterize gas transport in MOFs and explore its relevance to CO₂ capture and gas-separation applications.

---

## Graphical Abstract

![CO2 Self-Diffusion in Metal-Organic Frameworks](graphical%20abstruct.png)

*Graphical overview of the molecular dynamics workflow for investigating CO₂ self-diffusion in MOFs under varying N₂ co-loading conditions and framework pore characteristics. The illustration includes planned research extensions.*

---

## Simulation Methodology

### 1. Molecular Model Construction

Atomistic models of MOF-5 and a Ni-based MOF containing CO₂ molecules have been prepared for molecular dynamics simulations.

The MOF frameworks are treated as **rigid structures**, while CO₂ molecules are allowed to move within the framework pores.

The research will examine the influence of different framework structures and gas-loading conditions on molecular diffusion.

### 2. Interatomic Interactions

The molecular interactions are modeled using framework-specific force-field parameters and nonbonded interactions.

The computational approach incorporates:

- **TraPPE model:** Molecular representation of CO₂.
- **Lennard-Jones potential:** Short-range intermolecular interactions.
- **Coulombic interactions:** Electrostatic interactions between atomic sites.
- **PPPM method:** Long-range electrostatic calculations.
- **Lorentz–Berthelot mixing rules:** Cross-interaction parameters.

Framework-specific parameters are defined in the corresponding LAMMPS simulation inputs.

### 3. Molecular Dynamics Simulation

The computational workflow includes:

1. Initial molecular relaxation.
2. NVT equilibration at the specified simulation temperature.
3. NVE equilibration.
4. Molecular dynamics production runs for diffusion analysis.

The MOF framework is maintained as rigid during the simulations.

### 4. CO₂ Self-Diffusivity Calculation

CO₂ self-diffusion is investigated through mean-square displacement (MSD) analysis.

The self-diffusion coefficient is determined using the Einstein relation for three-dimensional diffusion:

```math
D_s = \frac{1}{6}\lim_{t\to\infty}\frac{d}{dt}\left\langle |\mathbf{r}(t)-\mathbf{r}(0)|^2 \right\rangle
```

where:

- **Dₛ:** CO₂ self-diffusion coefficient.
- **r(t):** Molecular position at time t.
- **r(0):** Initial molecular position.

The diffusion coefficient is estimated from the slope of the MSD curve within the linear diffusive regime.

Directional MSD components may also be examined to characterize anisotropic molecular transport.

---

## Research Files

The following folders contain the available computational models and molecular dynamics simulation scripts.

| Research Materials | Description | Access |
|---|---|---|
| **LAMMPS Simulation Code and Force Fields** | Molecular relaxation, equilibration, MD simulations, force-field parameters, and MSD calculations | [Explore Simulation Code →](MD%20code%20with%20force%20field/) |
| **MOF Structures and Simulation Data** | Molecular configurations of MOF frameworks containing CO₂ molecules | [Explore Simulation Data →](Data%20file/) |

---

## Research Status

**Status: Ongoing Research**

The repository currently documents the initial molecular models and computational simulation methodology.

The investigation of N₂ co-loading, additional MOF structures, and systematic pore-size variations is under development.

Numerical diffusivity results and comparative research findings are not included at this stage.

Additional simulation files and supporting documentation will be uploaded as the research progresses.

---

## Contact

**Mostafijur Rahaman**  
Department of Mechanical Engineering  
Bangladesh University of Engineering and Technology (BUET)  
Dhaka, Bangladesh

**Email:** [mostafijur.rahaman.buet@gmail.com](mailto:mostafijur.rahaman.buet@gmail.com)
