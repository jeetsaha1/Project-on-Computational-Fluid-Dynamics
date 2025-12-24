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
## 📊 Simulation Output & Data Handling
- Flow variables are stored in structured **`.dat` files**
- Results are generated in a **time-step–based numerical framework**
- Data organization is optimized for:
  - Efficient post-processing
  - Scientific visualization
  - Comparative analysis across cases

---

## 🚀 Future Scope (Next Phases)
This project is designed as a **progressive CFD framework** and will be extended beyond **Force 2.0** with the following planned enhancements:

- Integration of **advanced numerical solvers**
- Improved **stability and accuracy schemes**
- Introduction of **turbulence models**
- **Multi-physics coupling** (e.g., heat transfer, buoyancy effects)
- Optimization for **large-scale computational grids**
- Enhanced **visualization and post-processing pipelines**
- Validation against **standard benchmark CFD problems**

---

## 🎯 Project Objectives
- To numerically understand and analyze **fluid flow behavior**
- To study the impact of **physical and non-dimensional parameters** on flow
- To develop a **scalable, modular, and extensible CFD framework**
- To establish a strong foundation for **advanced simulations and research work**

---

## 🛠 Tools & Technologies
- Numerical methods for **Computational Fluid Dynamics**
- **Structured grid**–based discretization
- Data-driven **post-processing and visualization**
- Version control and collaboration using **Git & GitHub**

---

## 📌 Important Notes
- This repository is under **continuous development**
- Output files may be **large in size** due to numerical simulations
- Folder naming conventions reflect **simulation parameters** for clarity and reproducibility

---

## 📜 License
This project is intended strictly for **academic and research purposes**.

---

## ✨ Author & Supervision
**Author:** Jeet Saha  
**Supervisor:** Prof. (Dr.) Subhasree Dutta
