# Project on Computational Fluid Dynamics (CFD)

## 📌 Overview
This project focuses on the numerical simulation and analysis of fluid flow using **Computational Fluid Dynamics (CFD)** techniques.  
The current phase of the project is implemented using **Force 2.0**, and the project is designed to evolve further with advanced physical models, numerical schemes, and post-processing capabilities.

The repository contains simulation cases, parameter studies, and result datasets related to different flow regimes and conditions.

---

## ⚙️ Current Status
- **CFD Engine:** Force 2.0  
- **Project Phase:** Active development (Initial & Intermediate stages)
- **Simulation Type:** Numerical flow analysis
- **Platform:** Desktop-based computational environment

---

## 🧪 Features Implemented (Force 2.0)
- Grid-based numerical simulation
- Parametric studies for:
  - Reynolds number variations
  - Time-dependent flow behavior
  - Initial and boundary condition analysis
- Organized simulation cases for:
  - Different Reynolds numbers
  - Different time steps
  - Different inflow/outflow conditions
- Output data stored in structured `.dat` files

---

## 📂 Project Structure
```text
CFD_PROJECT/
│
├── Basic_time_0/           # Baseline simulation case <br>
├── nx_ny_lx_ly_*/          # Grid and domain configurations <br>
├── re_50_time_*            # Reynolds number = 50 cases <br>
├── re_100_time_*           # Reynolds number = 100 cases <br>
├── re_500_time_*           # Reynolds number = 500 cases <br>
├── ri_*_time_*             # Richardson number / related parameter studies <br>
├── .vscode/                # Editor configuration <br>
└── README.md               # Project documentation <br>
📊 Output Data
Velocity fields

Flow variables stored as .dat files

Time-step based numerical results

Structured for post-processing and visualization

🚀 Future Scope (Next Phases)
The project will be extended beyond Force 2.0 with the following goals:

Implementation of advanced numerical solvers

Improved stability and accuracy schemes

Turbulence modeling

Multi-physics coupling

Optimization for larger grids

Enhanced visualization and post-processing

Validation against benchmark CFD problems

🎯 Objectives
Understand fluid flow behavior numerically

Analyze the effect of physical parameters on flow

Develop a scalable and extensible CFD framework

Prepare the foundation for research and advanced simulations

🛠 Tools & Technologies
Numerical methods for CFD

Structured grids

Data-driven post-processing

Version control using Git & GitHub

📌 Notes
This repository is under continuous development.

Output files may be large due to numerical simulations.

Folder naming reflects simulation parameters for clarity.

📜 License
This project is intended for academic and research purposes.

✨ Author
Jeet Saha
