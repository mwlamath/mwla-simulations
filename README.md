# MWLA Simulations

This repository contains reproducible **numerical experiments** validating the  
**Mathematical Work of Logical Algebra (MWLA)** — a unified ledger framework for  
conservation, regularity, and computational balance across physical and informational systems.

Each simulation verifies that the discrete MWLA balance  


d/dt(|Vᵢ| qᵢ) + Σⱼ Jᵢⱼ = Sᵢ + Dᵢ

remains numerically closed (Σ Dᵢ ≈ 0) and stable across refinement levels.

---

## 🌐 Project Website

🔗 **[mwlamath.com](https://mwlamath.com)**  
Visit the MWLA research site for the full paper, formal definitions, and derivations.

---

## 🧠 Overview

This repository is organized by simulation type.  
Each folder contains:
- A self-contained simulation script (`code/`)
- Output data or zipped results (`results/`)
- Figures or images (`media/`)
- A localized `INSTRUCTIONS.md` for easy replication

---

## 🌀 Completed Simulation

### **Sim 1 — 2-D Lid-Driven Cavity (Ledger Audit)**

**Goal:** Verify that the discrete MWLA update preserves global conservation  
(Σ Dᵢ ≈ 0) and maintains a small normalized defect metric (Δ̄ₕ ≈ 10⁻⁹).

**Model**


∂t q + div(q u − D ∇q) = S, Ω = [0, 1] × [0, 1]
u = (1, 0) on lid; u = 0 elsewhere
D = ν I; S = 0


**Artifacts**
- 🧩 [Source code](./sims/cavity/code/mwla_cavity_v2_3_1.py)
- 📦 [Results archive](./sims/cavity/results/mwla_outputssuccess1.zip)
- 🖼️ [Visualization](./sims/cavity/media/cavity_v2_3_1.png)
- 📘 [Instructions](./sims/cavity/INSTRUCTIONS.md)

---

## ⚙️ Requirements

- Python 3.10 or newer  
- Packages: `numpy`, `matplotlib`  
(Optional) Git LFS for large `.zip` or `.png` assets

---

## 🧩 Folder Layout



MWLA-sims/
├─ README.md
├─ sims/
│ └─ cavity/
│ ├─ code/
│ │ └─ mwla_cavity_v2_3_1.py
│ ├─ results/
│ │ └─ mwla_outputssuccess1.zip
│ ├─ media/
│ │ └─ mwla_cavity.png
│ └─ INSTRUCTIONS.md


Future simulations (e.g., diffusion, quantum, interface) can be added under  
`sims/` using the same folder structure.

---

## 📜 License

MIT License — open for collaboration and derivative research.  
Please cite: *Valdez, J. (2025). “Mathematical Work of Logical Algebra: Unified Conservation Ledger Framework.”*