#!/usr/bin/env python3
"""
Near-protein urea enrichment from PSF+DCD (NAMD/CHARMM)

Definition:
  - Near: urea COM within r_near of any protein heavy atom
  - Bulk: urea COM farther than r_bulk from all protein heavy atoms
  - Enrichment E = (N_near / V_near) / (N_bulk / V_bulk)
Volumes V_near, V_bulk are estimated via Monte-Carlo sampling inside the PBC box.

python urea_enrichment.py \
  --psf input.psf \
  --dcd step5.dcd \
  --protein_sel "protein and not name H*" \
  --urea_sel "resname URE" \
  --r-near 0.5 \
  --r-bulk 0.6 \
  --stride 1 \
  --nrand 20000 \
  --seed 7

"""

import argparse, sys, math
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance

A_TO_NM = 0.1

def parse_args():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--psf", required=True, help="PSF/Topology file")
    ap.add_argument("--dcd", required=True, help="DCD trajectory (wrapped; must include unit cell)")
    ap.add_argument("--protein_sel", default="protein and not name H*", help="Protein heavy-atom selection")
    ap.add_argument("--urea_sel", default="resname URE", help="Urea molecule selection (entire residue)")
    ap.add_argument("--r-near", type=float, default=0.4, help="Near cutoff (nm)")
    ap.add_argument("--r-bulk", type=float, default=1.2, help="Bulk cutoff (nm)")
    ap.add_argument("--stride", type=int, default=10, help="Analyze every Nth frame")
    ap.add_argument("--start", type=int, default=None, help="First frame index (0-based)")
    ap.add_argument("--stop", type=int, default=None, help="Stop frame index (exclusive)")
    ap.add_argument("--nrand", type=int, default=5000, help="Random points per analyzed frame for volume")
    ap.add_argument("--seed", type=int, default=13, help="RNG seed for reproducibility")
    return ap.parse_args()

def wrap_box(coords, L):
    """Wrap coords to [0, L) for orthorhombic boxes; coords, L in Å."""
    return coords - np.floor(coords / L) * L

def capped_pairs(A, B, cutoff, box):
    """
    Version-safe wrapper for MDAnalysis.lib.distances.capped_distance.
    Returns (i_idx, j_idx) regardless of whether the underlying function
    returns (i,j) or (i,j,d) or an array of pairs.
    """
    out = capped_distance(A, B, max_cutoff=cutoff, box=box, return_distances=False)
    if isinstance(out, tuple):
        # Most versions: (i, j) ; some: (i, j, d)
        i_idx = np.asarray(out[0], dtype=int)
        j_idx = np.asarray(out[1], dtype=int)
        return i_idx, j_idx
    # Fallback: array of pairs
    pairs = np.asarray(out)
    if pairs.size == 0:
        return np.array([], dtype=int), np.array([], dtype=int)
    return pairs[:, 0].astype(int), pairs[:, 1].astype(int)

def main():
    args = parse_args()
    u = mda.Universe(args.psf, args.dcd)
    prot = u.select_atoms(args.protein_sel)
    if prot.n_atoms == 0:
        sys.exit(f"Protein selection matched 0 atoms: {args.protein_sel}")
    urea = u.select_atoms(args.urea_sel)
    if urea.n_residues == 0:
        sys.exit(f"Urea selection matched 0 residues: {args.urea_sel}")

    rng = np.random.default_rng(args.seed)
    r_near_A = args.r_near * 10.0   # nm -> Å
    r_bulk_A = args.r_bulk * 10.0

    # Accumulators
    counts_near, counts_bulk = [], []
    vols_near, vols_bulk = [], []

    # Per-frame densities for SEM
    rho_near_frames, rho_bulk_frames = [], []

    # Frame iterator
    traj = u.trajectory[args.start:args.stop:args.stride]

    for ts in traj:
        # Box & sanity
        Lx, Ly, Lz = ts.dimensions[:3]
        if Lx <= 0 or Ly <= 0 or Lz <= 0:
            sys.exit("Non-positive box lengths; ensure your DCD contains unit cell info.")
        L = np.array([Lx, Ly, Lz], dtype=float)

        # Protein heavy atoms (Å), wrapped (idempotent if already wrapped)
        prot_pos = wrap_box(prot.positions.copy(), L)

        # Urea COMs (Å), PBC-aware & wrapped
        urea_COM = []
        for res in urea.residues:
            try:
                com = res.atoms.center_of_mass(pbc=True)  # prefer PBC-aware COM
            except TypeError:
                # Older MDAnalysis: fall back (trajectory is wrapped so this is usually fine)
                com = res.atoms.center_of_mass()
            urea_COM.append(com)
        urea_COM = wrap_box(np.asarray(urea_COM, dtype=float), L)

        # ---- classify urea molecules ----
        # Near: any protein atom within r_near_A
        i_u_near, _ = capped_pairs(urea_COM, prot_pos, r_near_A, ts.dimensions)
        is_near = np.zeros(len(urea_COM), dtype=bool)
        if i_u_near.size:
            is_near[np.unique(i_u_near)] = True

        # Bulk: NO protein atom within r_bulk_A
        i_u_within_bulk, _ = capped_pairs(urea_COM, prot_pos, r_bulk_A, ts.dimensions)
        is_within_bulk_exclusion = np.zeros(len(urea_COM), dtype=bool)
        if i_u_within_bulk.size:
            is_within_bulk_exclusion[np.unique(i_u_within_bulk)] = True
        is_bulk = ~is_within_bulk_exclusion

        n_near = int(is_near.sum())
        n_bulk = int(is_bulk.sum())

        # ---- Monte-Carlo volume estimation ----
        # Sample random points uniformly in [0, L) Å
        rand = rng.random((args.nrand, 3)) * L[None, :]

        # Points within r_near_A of any protein atom
        i_rand_near, _ = capped_pairs(rand, prot_pos, r_near_A, ts.dimensions)
        near_point_mask = np.zeros(args.nrand, dtype=bool)
        if i_rand_near.size:
            near_point_mask[np.unique(i_rand_near)] = True

        # Points within r_bulk_A of any protein atom (EXCLUDED from bulk)
        i_rand_bulk, _ = capped_pairs(rand, prot_pos, r_bulk_A, ts.dimensions)
        within_bulk_exclusion = np.zeros(args.nrand, dtype=bool)
        if i_rand_bulk.size:
            within_bulk_exclusion[np.unique(i_rand_bulk)] = True

        # Fractions -> volumes
        V_box_A3 = (Lx * Ly * Lz)  # Å^3
        frac_near = near_point_mask.mean()
        frac_bulk = (~within_bulk_exclusion).mean()
        # Guard against degenerate volumes
        if frac_near <= 0 or frac_bulk <= 0:
            # If either is zero, skip this frame
            continue
        V_near = frac_near * V_box_A3
        V_bulk = frac_bulk * V_box_A3

        # Densities (molecules per Å^3)
        rho_near = n_near / V_near
        rho_bulk = n_bulk / V_bulk

        # Accumulate
        counts_near.append(n_near); counts_bulk.append(n_bulk)
        vols_near.append(V_near);   vols_bulk.append(V_bulk)
        rho_near_frames.append(rho_near); rho_bulk_frames.append(rho_bulk)

    # ---- aggregate over frames ----
    if len(rho_near_frames) == 0:
        sys.exit("No valid frames accumulated (check r_near/r_bulk and that the trajectory has unit cell info).")

    counts_near = np.array(counts_near, int)
    counts_bulk = np.array(counts_bulk, int)
    vols_near = np.array(vols_near, float)
    vols_bulk = np.array(vols_bulk, float)
    rho_near_frames = np.array(rho_near_frames, float)
    rho_bulk_frames = np.array(rho_bulk_frames, float)

    # Frame-averaged densities and SEM
    rho_near_mean = np.nanmean(rho_near_frames)
    rho_bulk_mean = np.nanmean(rho_bulk_frames)

    def sem(x):
        x = x[np.isfinite(x)]
        return np.std(x, ddof=1) / np.sqrt(len(x)) if len(x) > 1 else np.nan

    rho_near_sem = sem(rho_near_frames)
    rho_bulk_sem = sem(rho_bulk_frames)

    # Enrichment ratio
    E = np.nan
    if np.isfinite(rho_near_mean) and np.isfinite(rho_bulk_mean) and rho_bulk_mean > 0:
        E = rho_near_mean / rho_bulk_mean

    # Pretty print (convert volumes to nm^3 for readability)
    A3_TO_NM3 = 1e-3
    print("\n=== Near-protein urea enrichment ===")
    print(f"Frames analyzed: {len(rho_near_frames)}   stride: {args.stride}")
    print(f"Definitions: near ≤ {args.r_near:.2f} nm ; bulk ≥ {args.r_bulk:.2f} nm from protein heavy atoms")
    print(f"Urea molecules (residues): {urea.n_residues}")
    print(f"Per-frame means:")
    print(f"  <N_near> = {counts_near.mean():.1f}    <V_near> = {np.mean(vols_near)*A3_TO_NM3:.3f} nm^3")
    print(f"  <N_bulk> = {counts_bulk.mean():.1f}    <V_bulk> = {np.mean(vols_bulk)*A3_TO_NM3:.3f} nm^3")
    print(f"  ρ_near = {rho_near_mean:.4e} ± {rho_near_sem:.4e}  molecules/Å^3")
    print(f"  ρ_bulk = {rho_bulk_mean:.4e} ± {rho_bulk_sem:.4e}  molecules/Å^3")
    print(f"\n  Enrichment E = ρ_near / ρ_bulk = {E:.3f}")

if __name__ == "__main__":
    sys.exit(main())
