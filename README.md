# CFD-Analysis-of-Minimum-Length-Nozzle-Designed-with-MOC

Design a minimum-length converging–diverging nozzle with the **Method of Characteristics (MoC)**, generate a block-structured quadrilateral mesh in **Gmsh**, and solve the compressible flow using **SU2 (RANS–SST)**.  
This repo contains the MATLAB design script, a Python mesher, an SU2 config, and example figures/results.

<p align="center"><em>MoC characteristics → nozzle contour → conformal blocks → SU2 solution</em></p>

---

## Table of Contents

- [Features](#features)
- [Repository Layout](#repository-layout)
- [Requirements](#requirements)
- [Quick Start](#quick-start)
- [Workflow](#workflow)
- [Design Script (MATLAB)](#design-script-matlab)
- [Meshing Script (Python--Gmsh)](#meshing-script-python--gmsh)
- [SU2 Setup](#su2-setup)
- [Results](#results)
- [How to Reproduce](#how-to-reproduce)
- [Troubleshooting](#troubleshooting)
- [Cite / Acknowledgments](#cite--acknowledgments)
- [License](#license)

---

## Features

- **MoC nozzle design** with selectable converging profiles:
  - Circular Arc, Witoszyński (adjustable exponents), Bicubic, Conical, Cosine-bell, Quintic (“NASA smooth”), or “No converging section”.
- **Automatic characteristic net** and **minimum-length exit** (parallel, shock-free by construction).
- **Complete data export**: wall & full coordinates, *all* node properties, wall node properties, geometry features, and performance metrics (CSV + TXT).
- **Block-structured, transfinite quad mesh** in Gmsh with re-used vertical interfaces and symmetry axis.
- **SU2 RANS–SST** configuration for air with total-conditions inlet and fixed static outlet pressure.

---

## Repository Layout

