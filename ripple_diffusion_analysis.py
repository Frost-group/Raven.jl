import os
import argparse
from pathlib import Path
import numpy as np
import h5py
from ase.io import read, iread

from ripple.prepare_diffusion import prepare_diffusion

# Fixed analysis defaults
FIT_LAG_FRAC = (0.05, 0.95)   # your requested default
FIT_N_LAGS_OLS = 200
FIT_N_LAGS_MCMC = 100
PSD_PROJECT  = False          # keep False unless you explicitly want PSD projection for OLS

# MCMC defaults (edit here only if needed)
MCMC_CFG = dict(
    n_blocks=50,
    min_block_len=10,
    n_samples=1500,
    n_walkers=32,
    n_burn=400,
    n_thin=10,
    seed=42,
    progress=False,
    cond_max=1e16,
)


def _set_env_single_thread():
    # Avoid oversubscription on HPC
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")


def _write_fit_group(parent: h5py.Group, name: str, fit):
    g = parent.require_group(name)
    g.create_dataset("lag_min", data=np.int32(fit.lag_min))
    g.create_dataset("lag_max", data=np.int32(fit.lag_max))
    g.create_dataset("time_min_ps", data=np.float64(fit.time_min_ps))
    g.create_dataset("time_max_ps", data=np.float64(fit.time_max_ps))
    g.create_dataset(f"D_tensor_A2_per_ps", data=np.asarray(fit.D_tensor, dtype=np.float64))
    g.create_dataset(f"intercept_tensor_A2", data=np.asarray(fit.intercept_tensor, dtype=np.float64))
    g.create_dataset(f"r2_fro", data=float(fit.r2_fro) if fit.r2_fro is not None else np.nan)
    if getattr(fit, "D_std", None) is not None:
        g.create_dataset(f"D_std", data=np.asarray(fit.D_std, dtype=np.float64))
    if getattr(fit, "D_ci", None) is not None:
        g.create_dataset(f"D_ci_16_84", data=np.asarray(fit.D_ci, dtype=np.float64))
    if getattr(fit, "intercept_std", None) is not None:
        g.create_dataset(f"intercept_std", data=np.asarray(fit.intercept_std, dtype=np.float64))
    if getattr(fit, "intercept_ci", None) is not None:
        g.create_dataset(f"intercept_ci_16_84", data=np.asarray(fit.intercept_ci, dtype=np.float64))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", type=str, required=True, help="--path /gpfs/home/ll2121/Lammps_IAP/Project_LiPS/LiPS_0.2v_750K")
    parser.add_argument("--temps", type=int, required=True, help="--temps 700")
    parser.add_argument("--vacancies", type=str, required=True, help="--vacancies 0 or 0.05")
    parser.add_argument("--specie", type=str, required=True, help="--specie Li")
    parser.add_argument("--terminal", type=int, required=True, help="--terminal 24000")
    parser.add_argument("--frame_dt_fs", type=float, default=100.0, help="MD timestep * dump_every, in fs")
    parser.add_argument("--workers", type=int, default=2, help="--workers 64")
    args = parser.parse_args()
    
    T = int(args.temps)
    v_str = args.vacancies.strip()
    specie = args.specie.strip()
    terminal = int(args.terminal)
    frame_dt_ps = float(args.frame_dt_fs) / 1000.0
    n_workers = int(args.workers)
    cwd = Path(args.path)
    traj_name = f"traj_LiPS_844_{v_str}v_{T}K.dump"
    traj_path = cwd / traj_name
    if not traj_path.exists():
        raise FileNotFoundError(f"Trajectory not found: {traj_path}")
    out_path = cwd / f"RippleAnalysis_{traj_name}_f{terminal}.h5"
    print(f"[info] CWD: {cwd}")
    print(f"[info] traj: {traj_path.name}")
    print(f"[info] read first {terminal} frames, dt = {frame_dt_ps} ps, workers = {n_workers}")
    # ----- load trajectory -----
    traj_factory = lambda: iread(str(traj_path), format="lammps-dump-text", units="metal", index=f":{terminal}")
    # iread("traj_LiPS_844_0.2v_750K.dump", format="lammps-dump-text", units="metal", index=f":2000")
    
    # ----- prepare memmaps + compute curves -----
    traj_tag = f"{v_str}v_{T}K_f{terminal}"
    traj_diff = prepare_diffusion(
        traj_factory,
        target_atoms=specie,
        timestep_ps=frame_dt_ps,
        n_frames=terminal,
        save_dir=str(cwd),
        trajectory_tag=traj_tag,
        pos_dtype="float64",
        reuse_memmap=True,
    )
    msd_res = traj_diff.msd_cal(n_workers=n_workers)
    cmsd_res = traj_diff.cmsd_cal(n_workers=n_workers)
    TIME = np.asarray(msd_res.times_ps, dtype=np.float64)           # (T,)
    MSD  = np.asarray(msd_res.msd_tensor, dtype=np.float64)         # (T,3,3)
    cMSD = np.asarray(cmsd_res.cmsd_tensor, dtype=np.float64)       # (T,3,3)
    
    # ----- fit OLS + MCMC (self) -----
    fit_self_ols = msd_res.diffusion_tensor(
        method="ols",
        fit_lag_frac=FIT_LAG_FRAC,
        fit_n_lags=FIT_N_LAGS_OLS,
        psd=PSD_PROJECT,
    )
    fit_self_mcmc = msd_res.diffusion_tensor(
        method="mcmc",
        fit_lag_frac=FIT_LAG_FRAC,
        fit_n_lags=FIT_N_LAGS_MCMC,
        psd=PSD_PROJECT,
        return_chain=False,
        **MCMC_CFG,
    )
    # ----- fit OLS + MCMC (sigma) -----
    fit_sigma_ols = cmsd_res.diffusion_tensor(
        method="ols",
        fit_lag_frac=FIT_LAG_FRAC,
        fit_n_lags=FIT_N_LAGS_OLS,
        psd=PSD_PROJECT,
    )
    fit_sigma_mcmc = cmsd_res.diffusion_tensor(
        method="mcmc",
        fit_lag_frac=FIT_LAG_FRAC,
        fit_n_lags=FIT_N_LAGS_MCMC,
        psd=PSD_PROJECT,
        return_chain=False,
        **MCMC_CFG,
    )
    
    # ----- save HDF5 -----
    with h5py.File(out_path, "w") as h5f:
        # metadata
        s = h5py.string_dtype(encoding="utf-8")
        meta = h5f.require_group("meta")
        meta.create_dataset("trajectory_file", data=traj_path.name, dtype=s)
        meta.create_dataset("temperature_K", data=np.int32(T))
        meta.create_dataset("vacancy", data=np.float64(v_str))   # 或 v_val
        meta.create_dataset("specie", data=specie, dtype=s)
        meta.create_dataset("terminal_frames", data=np.int32(terminal))
        meta.create_dataset("frame_dt_fs", data=np.float64(args.frame_dt_fs))
        meta.create_dataset("n_sel", data=np.int32(msd_res.n_sel))
        meta.create_dataset("workers", data=np.int32(n_workers))
        # curves
        grp = h5f.require_group("curves")
        grp.create_dataset("time_ps", data=TIME, compression="gzip", compression_opts=4, shuffle=True)
        grp.create_dataset("msd_tensor_A2", data=MSD, compression="gzip", compression_opts=4, shuffle=True)
        grp.create_dataset("cmsd_tensor_A2", data=cMSD, compression="gzip", compression_opts=4, shuffle=True)
        # fits
        fitgrp = h5f.require_group("fits")
        fitgrp.create_dataset("fit_lag_frac", data=np.asarray(FIT_LAG_FRAC, dtype=np.float64))
        fitgrp.create_dataset("fit_n_lags_ols", data=np.int32(FIT_N_LAGS_OLS))
        fitgrp.create_dataset("fit_n_lags_mcmc", data=np.int32(FIT_N_LAGS_MCMC))
        fitgrp.create_dataset("psd_project", data=np.int8(int(PSD_PROJECT)))
        _write_fit_group(fitgrp, "D_self_ols", fit_self_ols)
        _write_fit_group(fitgrp, "D_self_mcmc", fit_self_mcmc)
        _write_fit_group(fitgrp, "D_sigma_ols", fit_sigma_ols)
        _write_fit_group(fitgrp, "D_sigma_mcmc", fit_sigma_mcmc)

        h5f.flush()

    print(f"[done] wrote: {out_path}")


if __name__ == "__main__":
    main()
