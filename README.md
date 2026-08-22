 # Mostafijur Rahaman — Research Projects Repository

Molecular dynamics and materials simulation projects from my undergraduate research at Bangladesh University of Engineering and Technology (BUET), Department of Mechanical Engineering. This repository collects the code, input files, and analysis scripts for three ongoing/completed projects, all centered on thermal transport and diffusion in nanoscale and porous materials.

**Tools used:** LAMMPS, OVITO, Atomsk, Materials Studio, VESTA, Python (NumPy, Matplotlib, MDAnalysis)

---

## 1. PDMS/2D-Filler Polymer Nanocomposites — Thermal Conductivity Enhancement
`/PDMS_hBN`

Undergraduate thesis project. Uses non-equilibrium molecular dynamics (NEMD) to study thermal conductivity enhancement in PDMS-based polymer nanocomposites reinforced with four 2D fillers — h-BN, graphene, BNC, and carbon-doped h-BN — at 9% volume fraction. Systems built with 56 PDMS chains × 50 monomers using a united-atom force field, sandwich-geometry filler placement, multi-stage NPT densification, and two-temperature NEMD (300 K / 260 K) to extract thermal conductivity. Currently being written up as a journal manuscript.

## 2. Defect Engineering in h-BN/Black Phosphorus Heterobilayers
`/2D material Heterobylayer new`

Investigates how point defects (nitrogen and boron vacancies at varying density/clustering) affect thermal conductivity in h-BN/BP heterobilayers, using reverse non-equilibrium MD (Müller-Plathe method). Notable finding: nitrogen vacancies increase thermal conductivity relative to boron vacancies — the opposite trend from isolated monolayer literature — with non-monotonic behavior at certain defect densities under further investigation. LAMMPS scripts use a hybrid/overlay potential combining Stillinger-Weber, Tersoff (BNC), and LJ/cut pair styles.

## 3. CO₂/N₂ Co-Diffusion in Metal-Organic Frameworks
`/MOF_CO2_Difussion_New`

Studies competitive diffusion of CO₂ and N₂ in Ni₃(HITP)₂ MOFs using LAMMPS with UFF4MOF and TraPPE force fields. Includes rigid-body CO₂ treatment (`fix rigid/nve/small`), DDEC6 partial charges, and MSD-based diffusivity extraction with per-axis anisotropy analysis.

---

## Repository Structure

Each project folder contains its own README with build details, LAMMPS input scripts, analysis code, and (where applicable) representative figures.

## Contact

Mostafijur Rahaman
Department of Mechanical Engineering, BUET
[LinkedIn](https://www.linkedin.com/in/mostafijur-rahaman-2b1835251/) · [Email](mailto:mostafijur.rahaman.buet@gmail.com)
