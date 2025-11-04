# mwla_cavity_v2_3_1.py — SOR Poisson with abs+rel tolerance, iter logging, adaptive omega
import sys, os, csv, argparse
import numpy as np
import matplotlib
try:
    if sys.platform.startswith("win"):
        matplotlib.use("TkAgg", force=True)
except Exception:
    pass
import matplotlib.pyplot as plt
from matplotlib import cm

def ramp(v_final, k, k_ramp):
    if k_ramp <= 0: return v_final
    a = min(1.0, max(0.0, k/float(k_ramp)))
    return a*v_final

def upwind_x(phi_c, phi_w, phi_e, u_c):
    pos = (u_c >= 0.0)
    return np.where(pos, phi_c - phi_w, phi_e - phi_c)

def upwind_y(phi_c, phi_s, phi_n, v_c):
    pos = (v_c >= 0.0)
    return np.where(pos, phi_c - phi_s, phi_n - phi_c)

def rms(a):
    a = np.nan_to_num(a, nan=0.0, posinf=0.0, neginf=0.0)
    return float(np.sqrt(np.mean(a*a)))

def poisson_sor_pressure(p, rhs, dx, iters_max=8000, omega=1.9, p_abs=1e-7, p_rel=1e-6,
                         adapt_omega=True, omega_min=1.5, stall_window=50):
    """
    Solve ∇²p = rhs on interior (Neumann dp/dn=0) using red–black SOR.
    Stop when res_rms <= p_abs + p_rel * rhs_rms.
    Returns: (iters_used, res_rms)
    """
    N = p.shape[0]
    h2 = dx*dx
    pc = p[1:-1, 1:-1]
    iy = np.arange(1, N-1)[:, None]
    ix = np.arange(1, N-1)[None, :]
    red_mask = ((iy + ix) % 2) == 0
    black_mask = ~red_mask

    rhs_rms = rms(rhs)

    best_res = np.inf
    worsened = 0

    for it in range(iters_max):
        # neighbors
        pE = p[1:-1, 2:]; pW = p[1:-1, 0:-2]
        pN = p[2:, 1:-1]; pS = p[0:-2, 1:-1]

        # red
        new_red = 0.25*(pE + pW + pN + pS - h2*rhs)
        pc[red_mask] = (1.0 - omega)*pc[red_mask] + omega*new_red[red_mask]

        p[:, 0]  = p[:, 1]; p[:, -1] = p[:, -2]
        p[0, :]  = p[1, :]; p[-1, :] = p[-2, :]

        # neighbors again
        pE = p[1:-1, 2:]; pW = p[1:-1, 0:-2]
        pN = p[2:, 1:-1]; pS = p[0:-2, 1:-1]

        # black
        new_blk = 0.25*(pE + pW + pN + pS - h2*rhs)
        pc[black_mask] = (1.0 - omega)*pc[black_mask] + omega*new_blk[black_mask]

        p[:, 0]  = p[:, 1]; p[:, -1] = p[:, -2]
        p[0, :]  = p[1, :]; p[-1, :] = p[-2, :]

        # residual
        pE = p[1:-1, 2:]; pW = p[1:-1, 0:-2]
        pN = p[2:, 1:-1]; pS = p[0:-2, 1:-1]
        lap = (pE + pW + pN + pS - 4.0*pc) / h2
        res_rms = rms(lap - rhs)

        # abs + rel stopping
        if res_rms <= (p_abs + p_rel * rhs_rms):
            return it + 1, res_rms

        # simple adaptive omega if stalled
        if adapt_omega:
            if res_rms < best_res:
                best_res = res_rms
                worsened = 0
            else:
                worsened += 1
                if worsened >= stall_window and omega > omega_min:
                    omega = max(omega_min, omega - 0.05)
                    worsened = 0  # reset and try again

    return iters_max, res_rms

def run_cavity(N=96, Re=100.0, U_lid=1.0, T_end=10.0, ramp_steps=100,
               CFL_conv=0.4, CFL_diff=0.45, p_iters_max=8000, omega=1.9,
               p_abs=1e-7, p_rel=1e-6, print_every=100):
    L = 1.0
    dx = dy = L/(N-1)
    nu = U_lid*L/Re

    u = np.zeros((N,N)); v = np.zeros((N,N)); p = np.zeros((N,N))

    denom = (1.0/(dx*dx) + 1.0/(dy*dy))
    t = 0.0; k = 0

    print(f"Grid {N}x{N} | Re={Re} nu={nu:.5g} | ramp={ramp_steps} | SOR ω={omega} "
          f"| p_abs={p_abs} p_rel={p_rel} | p_iters_max={p_iters_max} | CFLc={CFL_conv} CFLd={CFL_diff}")

    R_u_rms_hist, div_star_hist, dt_hist, time_hist, sor_iters_hist, pres_res_hist = [], [], [], [], [], []

    while t < T_end:
        u_old = u.copy(); v_old = v.copy()

        umax = float(np.max(np.abs(u_old))); vmax = float(np.max(np.abs(v_old)))
        Umax = max(umax, vmax, 1e-12)
        dt_conv = CFL_conv*min(dx,dy)/Umax
        dt_diff = CFL_diff/(nu*denom + 1e-16)
        dt = min(dt_conv, dt_diff, T_end - t)

        # interior aliases
        uc = u[1:-1,1:-1]; vc = v[1:-1,1:-1]
        uw = u[1:-1,:-2];  ue = u[1:-1,2:]
        us = u[:-2,1:-1];  un = u[2:,1:-1]
        vw = v[1:-1,:-2];  ve = v[1:-1,2:]
        vs = v[:-2,1:-1];  vn = v[2:,1:-1]

        # upwind advection
        du_dx = upwind_x(uc, uw, ue, uc)/dx
        du_dy = upwind_y(uc, us, un, vc)/dy
        dv_dx = upwind_x(vc, vw, ve, uc)/dx
        dv_dy = upwind_y(vc, vs, vn, vc)/dy
        u_conv = uc*du_dx + vc*du_dy
        v_conv = uc*dv_dx + vc*dv_dy

        # diffusion
        u_diff = nu*((ue - 2*uc + uw)/(dx*dx) + (un - 2*uc + us)/(dy*dy))
        v_diff = nu*((ve - 2*vc + vw)/(dx*dx) + (vn - 2*vc + vs)/(dy*dy))

        # predictor
        F_u = uc + dt*(-u_conv + u_diff)
        F_v = vc + dt*(-v_conv + v_diff)

        # Poisson RHS
        dFu_dx = (F_u[:,1:] - F_u[:,:-1])/dx
        dFv_dy = (F_v[1:,:] - F_v[:-1,:])/dy
        rhs_core = dFu_dx[1:,:] + dFv_dy[:,1:]
        rhs = np.zeros_like(F_u); rhs[1:,1:] = rhs_core

        iters_used, res_rms = poisson_sor_pressure(
            p, rhs, dx, iters_max=p_iters_max, omega=omega, p_abs=p_abs, p_rel=p_rel,
            adapt_omega=True
        )

        # corrector
        u[1:-1,1:-1] = F_u - dt*(p[1:-1,2:] - p[1:-1,:-2])/(2*dx)
        v[1:-1,1:-1] = F_v - dt*(p[2:,1:-1] - p[:-2,1:-1])/(2*dy)

        # boundaries
        Utop = ramp(U_lid, k, ramp_steps)
        u[0,:]=0.0; u[-1,:]=Utop; u[:,0]=0.0; u[:,-1]=0.0
        v[0,:]=0.0; v[-1,:]=0.0; v[:,0]=0.0; v[:,-1]=0.0

        # audits
        J_net = -u_conv + u_diff
        S_src = -(p[1:-1,2:] - p[1:-1,:-2])/(2*dx)
        R_u_rms = rms((u[1:-1,1:-1] - u_old[1:-1,1:-1])/dt + J_net - S_src)

        du_dx_c = (u[1:-1,2:] - u[1:-1,0:-2])/(2*dx)
        dv_dy_c = (v[2:,1:-1] - v[0:-2,1:-1])/(2*dy)
        div_star = rms(du_dx_c + dv_dy_c) * (1.0 / (U_lid / 1.0))

        R_u_rms_hist.append(R_u_rms)
        div_star_hist.append(div_star)
        dt_hist.append(dt)
        time_hist.append(t)
        sor_iters_hist.append(iters_used)
        pres_res_hist.append(res_rms)

        if (k % print_every) == 0:
            print(f"step {k:5d}  t={t:7.3f}  dt={dt:.3e}  Utop={Utop:.3f}  "
                  f"SOR iters={iters_used:4d}  p_res_rms={res_rms:.2e}  "
                  f"R_u_rms={R_u_rms:.2e}  div*={div_star:.2e}")

        k += 1
        t += dt

    return u, v, p, np.array(R_u_rms_hist), np.array(div_star_hist), np.array(dt_hist), np.array(time_hist), \
           np.array(sor_iters_hist), np.array(pres_res_hist), dx, dy

def save_outputs(outdir, name, u, v, p, R_u_rms, div_star, dt_hist, time_hist, sor_iters, pres_res, dx, dy):
    os.makedirs(outdir, exist_ok=True)
    csv_path = os.path.join(outdir, f"{name}_audit.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["step","time","dt","R_u_rms","div_star","sor_iters","p_res_rms"])
        for i in range(len(R_u_rms)):
            w.writerow([i, float(time_hist[i]), float(dt_hist[i]), float(R_u_rms[i]),
                        float(div_star[i]), int(sor_iters[i]), float(pres_res[i])])
    print(f"[Saved] {csv_path}")
    npz_path = os.path.join(outdir, f"{name}_fields.npz")
    np.savez_compressed(npz_path, u=u, v=v, p=p, R_u_rms=R_u_rms, div_star=div_star,
                        dt_hist=dt_hist, time_hist=time_hist,
                        sor_iters=sor_iters, p_res_rms=pres_res, dx=dx, dy=dy)
    print(f"[Saved] {npz_path}")

def plot_and_save(outdir, name, u, v, R_u_rms, div_star):
    N = u.shape[0]
    X, Y = np.meshgrid(np.linspace(0,1,N), np.linspace(0,1,N))
    vel = np.sqrt(u*u + v*v)
    fig, axes = plt.subplots(1,2, figsize=(14,6))

    ax = axes[0]
    cf = ax.contourf(X, Y, vel, 15, cmap=cm.viridis)
    skip = max(1, N//15)
    ax.quiver(X[::skip,::skip], Y[::skip,::skip], u[::skip,::skip], v[::skip,::skip], color='white', scale=5)
    ax.set_title("2D Lid-Driven Cavity (|u|)"); ax.set_xlabel("X"); ax.set_ylabel("Y"); ax.axis('equal')
    fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.04)

    ax = axes[1]
    ax.plot(R_u_rms, label="u-momentum residual RMS")
    ax.plot(div_star, label="divergence* (normalized)")
    ax.set_yscale("log"); ax.grid(True, which="both", ls="--", alpha=0.5)
    ax.set_title("Residual & Normalized Divergence"); ax.set_xlabel("Timestep"); ax.set_ylabel("RMS (log)")
    ax.legend()

    plt.tight_layout()
    png = os.path.join(outdir, f"{name}.png")
    fig.savefig(png, dpi=180)
    print(f"[Saved] {png}")
    try:
        gui_backends = {'Qt5Agg','QtAgg','TkAgg','MacOSX','WXAgg'}
        if matplotlib.get_backend() in gui_backends:
            plt.show(block=True)
        elif sys.platform.startswith("win"):
            os.startfile(os.path.abspath(png))  # type: ignore[attr-defined]
    except Exception as e:
        print(f"[INFO] Could not show image: {e}")

def main():
    ap = argparse.ArgumentParser(description="MWLA cavity v2.3.1 (abs+rel Poisson tol + iter logging)")
    ap.add_argument("--N", type=int, default=96)
    ap.add_argument("--Re", type=float, default=100.0)
    ap.add_argument("--T_end", type=float, default=10.0)
    ap.add_argument("--ramp_steps", type=int, default=100)
    ap.add_argument("--CFL_conv", type=float, default=0.4)
    ap.add_argument("--CFL_diff", type=float, default=0.45)
    ap.add_argument("--omega", type=float, default=1.9)
    ap.add_argument("--p_abs", type=float, default=1e-7)
    ap.add_argument("--p_rel", type=float, default=1e-6)
    ap.add_argument("--p_iters_max", type=int, default=8000)
    ap.add_argument("--print_every", type=int, default=100)
    ap.add_argument("--outdir", default="mwla_outputs")
    ap.add_argument("--name", default="cavity_v2_3_1")
    args = ap.parse_args()

    u, v, p, R_u_rms, div_star, dt_hist, time_hist, sor_iters, pres_res, dx, dy = run_cavity(
        N=args.N, Re=args.Re, U_lid=1.0, T_end=args.T_end, ramp_steps=args.ramp_steps,
        CFL_conv=args.CFL_conv, CFL_diff=args.CFL_diff,
        p_iters_max=args.p_iters_max, omega=args.omega,
        p_abs=args.p_abs, p_rel=args.p_rel, print_every=args.print_every
    )
    save_outputs(args.outdir, args.name, u, v, p, R_u_rms, div_star, dt_hist, time_hist, sor_iters, pres_res, dx, dy)
    plot_and_save(args.outdir, args.name, u, v, R_u_rms, div_star)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print("\n[ERROR]", e)
        import traceback; traceback.print_exc()
        try:
            if not sys.stdin or not sys.stdin.isatty():
                input("\nPress Enter to exit...")
        except Exception:
            pass
        raise
