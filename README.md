# MWLA Simulations

This repository contains reproducible **numerical experiments** validating the  
**Mathematical Work of Logical Algebra (MWLA)** — a unified ledger framework for conservation,  
regularity, and computational balance across physical and informational systems.

Each simulation verifies that the discrete MWLA balance  
d/dt(|Vᵢ| qᵢ) + Σⱼ Jᵢⱼ = Sᵢ + Dᵢ

yaml
Copy code
remains numerically closed (Σ Dᵢ ≈ 0) and stable across refinement levels.

---

## 🌐 Project Website

🔗 **[mwlamath.com](https://mwlamath.com)**  
Visit the official MWLA site for the research paper, derivations, and theoretical background.

---

## 🌀 Simulation 1 — Lid-Driven Cavity (Ledger Audit)

**Goal:** Confirm that the discrete MWLA update preserves global conservation  
and keeps the normalized defect metric small and stable.

**Model**
∂t q + div(q u − D ∇q) = S, Ω = [0, 1] × [0, 1]
u = (1, 0) on the lid, u = 0 elsewhere
D = ν I, S = 0

yaml
Copy code

**Artifacts**
- 🧩 [Source code](./sims/cavity/code/mwla_cavity_v2_3_1.py)  
- 📦 [Results archive](./sims/cavity/results/mwla_outputssuccess1.zip)  
- 🖼️ [Visualization](./sims/cavity/media/mwla_cavity.png)

---

## ⚙️ Requirements

- Python 3.10 or newer  
- Packages: `numpy`, `matplotlib`  
(Optional) Git LFS for large `.zip` or `.png` files

---

## 🧠 How to Reproduce

See [`INSTRUCTIONS.md`](./INSTRUCTIONS.md) for a full, step-by-step guide.

---

## 📁 Folder Layout

MWLA-sims/
├─ README.md
├─ INSTRUCTIONS.md
├─ sims/
│ └─ cavity/
│ ├─ code/
│ │ └─ mwla_cavity_v2_3_1.py
│ ├─ results/
│ │ └─ mwla_outputssuccess1.zip
│ └─ media/
│ └─ mwla_cavity.png

yaml
Copy code

Future simulations (e.g., diffusion, quantum, interface) can be added under `sims/`
following the same three-folder structure.

---

## 📜 License

MIT License — open for collaboration and derivative research.  
Please cite: *Valdez, J. (2025). “Mathematical Work of Logical Algebra: Unified Conservation Ledger Framework.”*