# MWLA Simulations

This repository hosts reproducible **numerical experiments** validating the  
**Mathematical Work of Logical Algebra (MWLA)** — a unified framework that treats  
physics, probability, and computation as one system of balance equations.

Every experiment checks that a simple accounting rule holds true:

d/dt(|Vᵢ| qᵢ) + Σⱼ Jᵢⱼ = Sᵢ + Dᵢ

markdown
Copy code

Here:
- **q** = quantity (mass, energy, probability, etc.)
- **J** = flux (flow of q)
- **S** = source or sink  
- **D** = “defect” — imbalance or conservation error

When `D → 0`, the system is balanced and smooth; when `D` grows, irregularities appear.  
MWLA studies how keeping this ledger balanced guarantees stability in physics and computation.

---

## 🌐 Project Website

🔗 **[mwlamath.com](https://mwlamath.com)**  
Visit the MWLA research site for formal definitions, proofs, and visual demonstrations.

---

## 🧠 Repository Overview

Each simulation is stored in its own folder under `sims/` and includes:
- `code/` — Python scripts and configuration files  
- `results/` — numerical outputs or compressed data  
- `media/` — figures, plots, and visualizations  
- `INSTRUCTIONS.md` — short guide to reproduce the experiment  

---

## 🌀 Completed Simulation

### **Sim 1 — 2-D Lid-Driven Cavity (Ledger Audit)**

A classic fluid-dynamics test adapted for MWLA.  
Instead of only measuring flow velocity or pressure, this run tracks how well the **numerical ledger** closes at every timestep.

**Plain-language summary:**
> Imagine a box of fluid with a moving lid on top.  
> MWLA checks whether every bit of fluid that moves out somewhere also moves in somewhere else —  
> just like balancing a checkbook for motion and energy.  
> The result: near-perfect accounting accuracy, showing that the discrete ledger works.

**Model**
∂t q + div(q u − D ∇q) = S, Ω = [0, 1] × [0, 1]
u = (1, 0) on the top lid; u = 0 elsewhere
D = ν I; S = 0

markdown
Copy code

**Results**
- Total defect Σ Dᵢ stayed within machine epsilon (~10⁻¹⁵)  
- Normalized defect Δ̄ₕ ≈ 10⁻⁹ — excellent numerical closure  
- Confirms theoretical predictions from MWLA Section 6

**Artifacts**
- 🧩 [Source code](./sims/cavity/code/mwla_cavity_v2_3_1.py)  
- 📦 [Results archive](./sims/cavity/results/mwla_outputssuccess1.zip)  
- 🖼️ [Visualization](./sims/cavity/media/mwla_cavity.png)  
- 📘 [Instructions](./sims/cavity/INSTRUCTIONS.md)

---

## ⚙️ Requirements

- Python 3.10 or newer  
- Packages: `numpy`, `matplotlib`  
- (Optional) Git LFS for large `.zip` and `.png` files  

---

## 🧩 Folder Layout

```text
MWLA-sims/
├─ README.md
├─ .gitignore
├─ .gitattributes
└─ sims/
   └─ cavity/
      ├─ code/
      │   └─ mwla_cavity_v2_3_1.py
      ├─ results/
      │   └─ mwla_outputssuccess1.zip
      ├─ media/
      │   └─ mwla_cavity.png
      └─ INSTRUCTIONS.md
Future experiments — diffusion, quantum continuity, coupled PDE-ODE systems —
can be added in the same structure under sims/.