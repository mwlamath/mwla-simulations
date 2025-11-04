# 🌀 Simulation 1 — Lid-Driven Cavity (MWLA Ledger Audit)

This folder reproduces the **2-D lid-driven cavity test** used to verify the discrete  
**Mathematical Work of Logical Algebra (MWLA)** ledger update and defect stability.

---

## 📘 1. Overview

The simulation integrates the MWLA balance equation

∂t q + div(q u − D ∇q) = S

pgsql
Copy code

over a unit square cavity domain.  
It checks that the discrete defect term

Dᵢⁿ⁺½ = |Vᵢ|(qᵢⁿ⁺¹ − qᵢⁿ) + Δt(ΣⱼJᵢⱼⁿ⁺½ − Sᵢⁿ⁺½)

yaml
Copy code

satisfies global closure (`Σ Dᵢ ≈ 0`) and that the normalized defect  
`Δ̄ₕ` remains small (≈10⁻⁹ to 10⁻⁸).

---

## ⚙️ 2. Requirements

- **Python:** 3.10 or newer  
- **Packages:**  
  ```bash
  pip install numpy matplotlib
No GPU libraries are required — this version runs fully on CPU.

▶️ 3. Run the Simulation
From the repository root, execute:

bash
Copy code
python sims/cavity/code/mwla_cavity_v2_3_1.py \
  --N 160 --Re 400 --T_end 20 \
  --ramp_steps 300 \
  --omega 1.9 \
  --p_abs 1e-7 --p_rel 1e-6 --p_iters_max 8000 \
  --CFL_conv 0.6 --CFL_diff 0.45 \
  --print_every 50 \
  --outdir sims/cavity/results \
  --name cavity_v2_3_1
This will:

Print per-step solver diagnostics (time, dt, residuals)

Save logs and NumPy arrays under sims/cavity/results/

Generate a visualization similar to media/mwla_cavity.png

📊 4. Expected Results
Quantity	Expected Behavior
Global defect Σ Dᵢ	≈ 0 (within 1e-15)
Normalized defect Δ̄ₕ	1e-9 to 1e-8
Flow field	Stable recirculation pattern

A successful run confirms ledger closure and numerical smoothness predicted by MWLA Section 6.

📦 5. Files
bash
Copy code
code/mwla_cavity_v2_3_1.py       → simulation script
results/mwla_outputssuccess1.zip → archived output data
media/mwla_cavity.png            → sample visualization
INSTRUCTIONS.md                  → this file
🌐 6. References
More about the theory and formalism:
🔗 MWLA Research Website

