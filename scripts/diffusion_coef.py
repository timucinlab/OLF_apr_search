#!/usr/bin/env python3
"""
Diffusion coefficient calculator for urea and water.

python ../diffusion_coef.py \
    --psf input.psf \
    --dcd unwrapped.dcd \
    --urea_sel "resname URE" \
    --water_sel "resname TIP3" \
    --tau-min 500 \
    --tau-max 2000 \
    --dt-ps 100 \
    --sample-urea 500 \
    --sample-water 1500 \
    --boot 300
"""

import argparse
import sys
import numpy as np
import MDAnalysis as mda

# ------------------------ Arguments ------------------------
def parse_args():
    ap = argparse.ArgumentParser(
        description="diffusion coefficients for urea and water",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    ap.add_argument("--psf", required=True, help="PSF topology file")
    ap.add_argument("--dcd", required=True, help="DCD trajectory (unwrapped)")
    ap.add_argument("--urea_sel", default="resname URE", help="Urea selection string")
    ap.add_argument("--water_sel", default="resname TIP3 TIP3P SOL HOH", 
                   help="Water selection string")
    ap.add_argument("--tau-min", type=float, default=500.0, 
                   help="Fit window start (ps)")
    ap.add_argument("--tau-max", type=float, default=2000.0, 
                   help="Fit window end (ps)")
    ap.add_argument("--max-lag", type=float, default=3000.0, 
                   help="Max tau considered (ps)")
    ap.add_argument("--stride", type=int, default=1, 
                   help="Analyze every Nth frame")
    ap.add_argument("--start", type=int, default=None, 
                   help="Starting frame")
    ap.add_argument("--stop", type=int, default=None, 
                   help="Stopping frame")
    ap.add_argument("--dt-ps", type=float, default=None, 
                   help="Frame spacing (ps); overrides DCD if given")
    ap.add_argument("--sample-urea", type=int, default=None, 
                   help="Number of urea molecules to sample (None=all)")
    ap.add_argument("--sample-water", type=int, default=None, 
                   help="Number of water molecules to sample (None=all)")
    ap.add_argument("--boot", type=int, default=300, 
                   help="Bootstrap replicates for error estimation")
    ap.add_argument("--seed", type=int, default=42, 
                   help="Random seed")
    ap.add_argument("--output", type=str, default=None,
                   help="Output file for results (default: print to screen)")
    return ap.parse_args()

# ------------------------ Helper Functions ------------------------
def fit_slope(x, y):
    """
    Linear fit y ~ a + b*x ignoring NaNs.
    Returns: slope b, R², number of points
    """
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return np.nan, np.nan, int(mask.sum())
    
    X = np.vstack([np.ones(mask.sum()), x[mask]]).T
    beta, *_ = np.linalg.lstsq(X, y[mask], rcond=None)
    
    y_pred = X @ beta
    residuals = y[mask] - y_pred
    ss_res = float(np.sum(residuals**2))
    ss_tot = float(np.sum((y[mask] - np.mean(y[mask]))**2))
    r2 = 1.0 - ss_res/ss_tot if ss_tot > 0 else np.nan
    
    return float(beta[1]), r2, int(mask.sum())

def compute_msd_per_molecule(positions, max_lag_frames):
    """
    Compute MSD for each molecule separately.
    
    Parameters
    ----------
    positions : ndarray, shape (n_frames, n_molecules, 3)
        Unwrapped COM positions in Angstroms
    max_lag_frames : int
        Maximum lag time in frames
    
    Returns
    -------
    msd_per_mol : ndarray, shape (n_molecules, max_lag_frames)
        MSD curve for each molecule
    lags : ndarray, shape (max_lag_frames,)
        Lag times in frames
    """
    n_frames, n_mols, _ = positions.shape
    max_lag = min(max_lag_frames, n_frames - 1)
    lags = np.arange(1, max_lag + 1, dtype=int)
    
    msd_per_mol = np.full((n_mols, max_lag), np.nan, dtype=float)
    
    for mol_idx in range(n_mols):
        pos_mol = positions[:, mol_idx, :]  # (n_frames, 3)
        
        for lag_idx, lag in enumerate(lags):
            # Compute displacements for this lag
            displacements = pos_mol[lag:] - pos_mol[:-lag]  # (n_frames-lag, 3)
            squared_disp = np.sum(displacements**2, axis=1)  # (n_frames-lag,)
            msd_per_mol[mol_idx, lag_idx] = np.mean(squared_disp)
    
    return msd_per_mol, lags

def bootstrap_diffusion(msd_per_mol, lags_ps, tau_min, tau_max, n_boot, rng):
    """
    Compute diffusion coefficient with bootstrap error estimation.
    
    Parameters
    ----------
    msd_per_mol : ndarray, shape (n_molecules, n_lags)
        MSD for each molecule at each lag
    lags_ps : ndarray
        Lag times in picoseconds
    tau_min, tau_max : float
        Fit window in picoseconds
    n_boot : int
        Number of bootstrap replicates
    rng : numpy.random.Generator
        Random number generator
    
    Returns
    -------
    results : dict
        D (diffusion coefficient in m²/s), 
        D_std (standard error),
        ci (95% confidence interval),
        R2 (coefficient of determination),
        n_fit (number of points in fit)
    """
    n_mols, n_lags = msd_per_mol.shape
    
    # Define fit window
    fit_mask = (lags_ps >= tau_min) & (lags_ps <= tau_max)
    if fit_mask.sum() < 3:
        return {
            'D': np.nan, 'D_std': np.nan, 'ci': (np.nan, np.nan),
            'R2': np.nan, 'n_fit': 0
        }
    
    lags_fit = lags_ps[fit_mask]
    
    # Point estimate using all molecules
    msd_mean = np.nanmean(msd_per_mol[:, fit_mask], axis=0)
    slope, r2, n_pts = fit_slope(lags_fit, msd_mean)
    D_point = (slope / 6.0) * 1e-8  # Å²/ps to m²/s: 1 Å²/ps = 1e-8 m²/s
    
    # Bootstrap
    D_boot = []
    for _ in range(n_boot):
        # Resample molecules with replacement
        idx = rng.integers(0, n_mols, size=n_mols)
        msd_boot = np.nanmean(msd_per_mol[idx][:, fit_mask], axis=0)
        slope_boot, _, _ = fit_slope(lags_fit, msd_boot)
        D_boot.append((slope_boot / 6.0) * 1e-8)
    
    D_boot = np.array(D_boot)
    D_boot = D_boot[np.isfinite(D_boot)]  # Remove NaNs
    
    if len(D_boot) > 0:
        D_mean = float(np.mean(D_boot))
        D_std = float(np.std(D_boot))
        ci_low = float(np.percentile(D_boot, 2.5))
        ci_high = float(np.percentile(D_boot, 97.5))
    else:
        D_mean = D_point
        D_std = np.nan
        ci_low = np.nan
        ci_high = np.nan
    
    return {
        'D': D_point,
        'D_mean': D_mean,
        'D_std': D_std,
        'ci': (ci_low, ci_high),
        'R2': float(r2),
        'n_fit': int(n_pts)
    }

def analyze_species(u, mol_sel, sample_n, dt_ps, stride, start, stop, 
                   max_lag_ps, tau_min, tau_max, n_boot, rng, species_name):
    """
    Analyze diffusion for species.
    
    Returns
    -------
    results : dict
        Contains diffusion coefficient and statistics
    """
    print(f"\n{'='*70}")
    print(f"Analyzing {species_name}...")
    print(f"{'='*70}")
    
    # Select molecules
    mol = u.select_atoms(mol_sel)
    if mol.n_residues == 0:
        print(f"ERROR: No molecules found with selection: {mol_sel}")
        return None
    
    print(f"Total {species_name} molecules: {mol.n_residues}")
    
    # Subsample if requested
    residues = list(mol.residues)
    if sample_n is not None and sample_n < len(residues):
        idx = rng.choice(len(residues), size=sample_n, replace=False)
        residues = [residues[i] for i in sorted(idx)]
        print(f"Subsampled to: {len(residues)} molecules")
    
    n_mols = len(residues)
    
    # Extract trajectory
    traj = u.trajectory[start:stop:stride]
    print(f"Frames to analyze: {len(traj)}")
    
    # Collect COM positions (unwrapped)
    positions = []
    for ts in traj:
        coms = np.array([res.atoms.center_of_mass() for res in residues])
        positions.append(coms)
    
    positions = np.array(positions)  # (n_frames, n_mols, 3)
    n_frames = positions.shape[0]
    
    print(f"Collected positions: {n_frames} frames × {n_mols} molecules")
    
    # Compute MSD
    max_lag_frames = int(max_lag_ps / dt_ps)
    print(f"Computing MSD up to {max_lag_frames} frames ({max_lag_ps:.0f} ps)...")
    
    msd_per_mol, lags = compute_msd_per_molecule(positions, max_lag_frames)
    lags_ps = lags * dt_ps
    
    # Fit and bootstrap
    print(f"Fitting diffusion coefficient (window: {tau_min:.0f}-{tau_max:.0f} ps)...")
    results = bootstrap_diffusion(msd_per_mol, lags_ps, tau_min, tau_max, 
                                  n_boot, rng)
    
    # Print results
    print(f"\nResults for {species_name}:")
    print(f"  D = {results['D']:.3e} m²/s")
    print(f"  D (bootstrap mean) = {results['D_mean']:.3e} m²/s")
    print(f"  Standard error = {results['D_std']:.3e} m²/s")
    print(f"  95% CI = [{results['ci'][0]:.3e}, {results['ci'][1]:.3e}] m²/s")
    print(f"  R² = {results['R2']:.4f}")
    print(f"  Fit points = {results['n_fit']}")
    
    return results

# ------------------------ Main ------------------------
def main():
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    
    print("="*70)
    print("Diffusion Coefficient Calculator")
    print("="*70)
    print(f"\nInput files:")
    print(f"  PSF: {args.psf}")
    print(f"  DCD: {args.dcd}")
    
    # Load trajectory
    u = mda.Universe(args.psf, args.dcd)
    
    # Get time step
    if args.dt_ps is not None:
        dt_ps = args.dt_ps
    else:
        dt_ps = getattr(u.trajectory, 'dt', None)
        if dt_ps is None:
            print("\nERROR: Cannot determine time step. Use --dt-ps")
            return 1
    
    dt_ps = float(dt_ps) * args.stride
    
    print(f"\nTrajectory info:")
    print(f"  Total frames: {len(u.trajectory)}")
    print(f"  Stride: {args.stride}")
    print(f"  Effective dt: {dt_ps:.3f} ps")
    print(f"  Analysis window: frames {args.start or 0} to {args.stop or len(u.trajectory)}")
    
    print(f"\nAnalysis parameters:")
    print(f"  Fit window: {args.tau_min:.0f} - {args.tau_max:.0f} ps")
    print(f"  Max lag: {args.max_lag:.0f} ps")
    print(f"  Bootstrap replicates: {args.boot}")
    print(f"  Random seed: {args.seed}")
    
    # Analyze urea
    results_urea = analyze_species(
        u, args.urea_sel, args.sample_urea, dt_ps, args.stride,
        args.start, args.stop, args.max_lag, args.tau_min, args.tau_max,
        args.boot, rng, "UREA"
    )
    
    # Analyze water (reload universe to be safe)
    u2 = mda.Universe(args.psf, args.dcd)
    results_water = analyze_species(
        u2, args.water_sel, args.sample_water, dt_ps, args.stride,
        args.start, args.stop, args.max_lag, args.tau_min, args.tau_max,
        args.boot, rng, "WATER"
    )
    
    # Summary
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    
    if results_urea:
        print(f"\nUREA:")
        print(f"  D = {results_urea['D']:.3e} ± {results_urea['D_std']:.3e} m²/s")
        print(f"  95% CI = [{results_urea['ci'][0]:.3e}, {results_urea['ci'][1]:.3e}]")
        print(f"  R² = {results_urea['R2']:.4f}")
    
    if results_water:
        print(f"\nWATER:")
        print(f"  D = {results_water['D']:.3e} ± {results_water['D_std']:.3e} m²/s")
        print(f"  95% CI = [{results_water['ci'][0]:.3e}, {results_water['ci'][1]:.3e}]")
        print(f"  R² = {results_water['R2']:.4f}")
    
    # Write output file if requested
    if args.output:
        with open(args.output, 'w') as f:
            f.write("# Diffusion Coefficient Results\n")
            f.write(f"# PSF: {args.psf}\n")
            f.write(f"# DCD: {args.dcd}\n")
            f.write(f"# dt: {dt_ps:.3f} ps\n")
            f.write(f"# Fit window: {args.tau_min:.0f}-{args.tau_max:.0f} ps\n")
            f.write(f"# Bootstrap: {args.boot} replicates\n\n")
            
            f.write("Species\tD(m2/s)\tD_std(m2/s)\tCI_low(m2/s)\tCI_high(m2/s)\tR2\tn_fit\n")
            
            if results_urea:
                f.write(f"UREA\t{results_urea['D']:.6e}\t{results_urea['D_std']:.6e}\t"
                       f"{results_urea['ci'][0]:.6e}\t{results_urea['ci'][1]:.6e}\t"
                       f"{results_urea['R2']:.4f}\t{results_urea['n_fit']}\n")
            
            if results_water:
                f.write(f"WATER\t{results_water['D']:.6e}\t{results_water['D_std']:.6e}\t"
                       f"{results_water['ci'][0]:.6e}\t{results_water['ci'][1]:.6e}\t"
                       f"{results_water['R2']:.4f}\t{results_water['n_fit']}\n")
        
        print(f"\nResults written to: {args.output}")
    
    print("\nNotes:")
    print("  • Uses Einstein relation: D = <r²>/6t")
    print("  • Assumes trajectory is already unwrapped")
    print("  • Standard errors from bootstrap over molecules")    

    return 0

if __name__ == "__main__":
    sys.exit(main())
