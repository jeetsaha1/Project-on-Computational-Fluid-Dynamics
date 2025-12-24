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
├── Basic_time_0/           # Baseline simulation case 
├── nx_ny_lx_ly_*/          # Grid and domain configurations 
├── re_50_time_*            # Reynolds number = 50 cases 
├── re_100_time_*           # Reynolds number = 100 cases 
├── re_500_time_*           # Reynolds number = 500 cases 
├── ri_*_time_*             # Richardson number / related parameter studies 
├── .vscode/                # Editor configuration 
└── README.md               # Project documentation 
📊 Output Data
Velocity fields
```
📊 Simulation Outputs & Data Organization

Flow variables are stored in structured .dat files

Results are generated in a time-step–based numerical framework

Data layout is optimized for post-processing, visualization, and analysis

Output structure supports reproducibility and parametric comparison across cases

🚀 Future Scope & Project Roadmap

This project is designed as a continuously evolving CFD framework.
Future development beyond Force 2.0 will focus on:

Integration of advanced numerical solvers

Development of higher-order stability and accuracy schemes

Implementation of turbulence models

Multi-physics coupling (heat transfer, buoyancy-driven flows, etc.)

Performance optimization for large-scale computational grids

Advanced flow visualization and post-processing pipelines

Validation and benchmarking against standard CFD test cases

🎯 Project Objectives

Numerically investigate and understand fluid flow behavior

Analyze the influence of physical and numerical parameters on flow characteristics

Design a scalable, modular, and extensible CFD framework

Establish a strong foundation for advanced simulations and research-oriented work

🛠 Tools & Technologies

Numerical methods applied to Computational Fluid Dynamics

Structured grid–based discretization techniques

Data-driven post-processing methodologies

Version control and collaboration using Git & GitHub

📌 Important Notes

This repository is under active and continuous development

Simulation output files may be large due to numerical resolution and time-stepping

Folder naming conventions clearly reflect simulation parameters for traceability

📜 License

This project is intended solely for academic and research purposes.

✨ Author & Supervision

Author: Jeet Saha
Supervisor: Prof. (Dr.) Subhasree Dutta
