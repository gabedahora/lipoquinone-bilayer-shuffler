
#!/usr/bin/env python3
import argparse
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

# -------------------------
# Utilities
# -------------------------
def infer_unit_scale_from_box(Lz, expected_Lz_nm=15.83588):
    """Return nm->internal scale and unit label (Å or nm) using box length heuristic."""
    if abs(Lz - expected_Lz_nm) < 3.0:
        return 1.0, "nm"
    if abs(Lz - expected_Lz_nm * 10.0) < 30.0:
        return 10.0, "Å"
    if Lz > 80.0:
        return 10.0, "Å"
    return 1.0, "nm"

def rodrigues_align(a, b):
    """
    Rotation matrix that aligns vector a to vector b (both 3D).
    Handles near-parallel and near-antiparallel cases.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    na = np.linalg.norm(a)
    nb = np.linalg.norm(b)
    if na < 1e-12 or nb < 1e-12:
        return np.eye(3)
    a /= na
    b /= nb

    v = np.cross(a, b)
    c = np.dot(a, b)

    if np.linalg.norm(v) < 1e-12:
        # parallel or antiparallel
        if c > 0:
            return np.eye(3)
        # 180° rotation around any axis orthogonal to a
        # pick an arbitrary orthogonal axis
        axis = np.array([1.0, 0.0, 0.0])
        if abs(a[0]) > 0.9:
            axis = np.array([0.0, 1.0, 0.0])
        v = np.cross(a, axis)
        v /= np.linalg.norm(v)
        # Rodrigues for 180°: R = -I + 2 vv^T
        return -np.eye(3) + 2.0 * np.outer(v, v)

    s = np.linalg.norm(v)
    vx = np.array([[0,    -v[2],  v[1]],
                   [v[2],  0,    -v[0]],
                   [-v[1], v[0],  0]])
    R = np.eye(3) + vx + (vx @ vx) * ((1 - c) / (s**2))
    return R

def rotz(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s, 0.0],
                     [s,  c, 0.0],
                     [0.0, 0.0, 1.0]])

def wrap_residue_as_whole(coords, ref_point, box):
    """
    Wrap a residue as a whole using a reference point (e.g., head center).
    Shifts all atoms by integer box vectors so ref_point ends up in [0, L).
    """
    Lx, Ly, Lz = box[:3]
    shift = np.zeros(3)
    shift[0] = -np.floor(ref_point[0] / Lx) * Lx
    shift[1] = -np.floor(ref_point[1] / Ly) * Ly
    shift[2] = -np.floor(ref_point[2] / Lz) * Lz
    return coords + shift

# -------------------------
# Argument parsing
# -------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Randomize UQ8 positions in bilayer without breaking molecules or placing in water.")
    p.add_argument("-i", "--input", default="step5_input.gro")
    p.add_argument("-o", "--output", default="shuffled.gro")
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--resname", default="UQ8")

    # Clash distances (nm)
    p.add_argument("--min-clash-prot-nm", type=float, default=0.30,
                   help="Min heavy-atom distance to protein and placed UQ8 (nm).")
    p.add_argument("--min-clash-lip-nm", type=float, default=0.18,
                   help="Min heavy-atom distance to lipids (nm). Lower helps packing.")
    # Protein groove exclusion (nm) using head/ring center to protein heavy atoms
    p.add_argument("--ring-protein-excl-nm", type=float, default=0.80,
                   help="Hard exclusion: head/ring center must be >= this from protein heavy atoms (nm).")

    # Z definitions (nm)
    p.add_argument("--z-margin-nm", type=float, default=0.50,
                   help="Margin inside phosphate planes for safe bounds (nm).")
    p.add_argument("--z-keepout-nm", type=float, default=0.20,
                   help="Extra keepout from phosphate planes for all heavy atoms (nm). Prevents water poking.")
    p.add_argument("--z-head-band-nm", type=float, default=0.40,
                   help="Gaussian sigma for head z around leaflet interface target (nm).")

    p.add_argument("--max-attempts", type=int, default=50000)

    # UQ8 head/tail selections: edit if needed after checking UQ8.itp
    p.add_argument("--head-names", nargs="*", default=["O1","O2","O3","O4","O5","O6"],
                   help="Atom names used to define UQ8 head/ring (defaults are guesses).")
    p.add_argument("--tail-last-n-heavy", type=int, default=6,
                   help="Fallback: use last N heavy atoms in residue as 'tail' if selection fails.")
    return p.parse_args()

# -------------------------
# Main
# -------------------------
def main():
    args = parse_args()
    if args.seed is not None:
        np.random.seed(args.seed)

    u = mda.Universe(args.input)
    box = u.dimensions.copy()
    Lx, Ly, Lz = box[:3]

    nm2int, unit = infer_unit_scale_from_box(Lz)
    min_prot = args.min_clash_prot_nm * nm2int
    min_lip  = args.min_clash_lip_nm  * nm2int
    ring_excl = args.ring_protein_excl_nm * nm2int
    z_margin = args.z_margin_nm * nm2int
    z_keepout = args.z_keepout_nm * nm2int
    z_band = args.z_head_band_nm * nm2int

    print(f"[INFO] Box: Lx={Lx:.3f}, Ly={Ly:.3f}, Lz={Lz:.3f} ({unit} internal)")
    print(f"[INFO] Thresholds ({unit}): min_prot={min_prot:.3f}, min_lip={min_lip:.3f}, ring_excl={ring_excl:.3f}")
    print(f"[INFO] Z params ({unit}): z_margin={z_margin:.3f}, z_keepout={z_keepout:.3f}, z_band={z_band:.3f}")

    ubq = u.select_atoms(f"resname {args.resname}")
    if ubq.n_atoms == 0:
        raise RuntimeError(f"No atoms found with resname {args.resname}")

    protein_heavy = u.select_atoms("protein and not name H*")

    lipid_resnames = ["POPE","POPG","DPPE","DPPG","TOCL2","DPPC","POPC"]
    lipids_heavy = u.select_atoms(f"resname {' '.join(lipid_resnames)} and not name H*")

    P = u.select_atoms("name P")
    if P.n_atoms == 0:
        raise RuntimeError("No phosphate atoms found with selection 'name P'. Check lipid naming.")

    # Center P mean z to box center to reduce PBC artifacts
    shift_z = (Lz/2.0) - P.positions[:,2].mean()
    u.atoms.positions[:,2] += shift_z

    # Leaflet detection by median split
    zP = P.positions[:,2]
    z_mid = np.median(zP)
    upper = zP[zP >= z_mid]
    lower = zP[zP <  z_mid]
    if len(upper) < 10 or len(lower) < 10:
        raise RuntimeError("Leaflet split failed; membrane may be mis-imaged.")

    upper_z = np.percentile(upper, 90)
    lower_z = np.percentile(lower, 10)
    z_center = 0.5*(upper_z + lower_z)

    # Safe bounds inside phosphate planes
    z_min_safe = lower_z + z_margin
    z_max_safe = upper_z - z_margin
    if z_min_safe >= z_max_safe:
        raise RuntimeError("Bad safe Z window. Reduce --z-margin-nm.")

    # Interface targets (where head should start): between P plane and midplane
    z_if_upper = upper_z - (upper_z - z_center)/3.0
    z_if_lower = lower_z + (z_center - lower_z)/3.0

    print(f"[INFO] Phosphate planes ({unit}): lower~{lower_z:.3f}, upper~{upper_z:.3f}, mid~{z_center:.3f}")
    print(f"[INFO] Safe Z ({unit}): z_min_safe={z_min_safe:.3f}, z_max_safe={z_max_safe:.3f}")
    print(f"[INFO] Head targets ({unit}): z_if_lower={z_if_lower:.3f}, z_if_upper={z_if_upper:.3f}")

    placed_uq8_heavy_indices = []

    for res in ubq.residues:
        res_atoms = res.atoms
        res_heavy = res_atoms.select_atoms("not name H*")
        if res_heavy.n_atoms == 0:
            continue

        # Define head selection by names; fallback to first ~6 heavy atoms if not found
        head = res_atoms.select_atoms("name " + " ".join(args.head_names))
        if head.n_atoms == 0:
            head = res_heavy[:min(6, res_heavy.n_atoms)]
            head_mode = "FALLBACK(first heavy)"
        else:
            head_mode = "head_names"

        # Define tail selection: last N heavy atoms (robust fallback)
        Ntail = min(args.tail_last_n_heavy, res_heavy.n_atoms)
        tail = res_heavy[-Ntail:]

        # Use head center as anchor
        head_center = head.positions.mean(axis=0)
        tail_center = tail.positions.mean(axis=0)

        # Reference coords relative to head_center (NOT COG)
        ref = res_atoms.positions - head_center
        # Head/tail vectors in reference frame
        v_ht = tail_center - head_center  # head->tail direction in original pose

        placed = False
        for attempt in range(1, args.max_attempts + 1):
            # Random XY
            x = np.random.rand() * Lx
            y = np.random.rand() * Ly

            # Choose leaflet for head placement
            if np.random.rand() < 0.5:
                z_head0 = z_if_lower
                target = np.array([0.0, 0.0, +1.0])  # tail should point +z toward midplane
                # Actually for LOWER leaflet, midplane is above => tail points +z
            else:
                z_head0 = z_if_upper
                target = np.array([0.0, 0.0, -1.0])  # tail should point -z toward midplane

            # Sample head z near interface and clamp to safe bounds
            z_head = np.random.normal(loc=z_head0, scale=z_band)
            z_head = min(max(z_head, z_min_safe), z_max_safe)

            # Align head->tail vector to target (+z or -z)
            R_align = rodrigues_align(v_ht, target)

            # Random azimuth about z to diversify
            theta = np.random.rand() * 2*np.pi
            R_azi = rotz(theta)

            # Apply rotations (azimuth after alignment)
            R = R_azi @ R_align

            coords_trial = ref @ R.T + np.array([x, y, z_head])
            # DO NOT atom-wise wrap here

            # Whole-residue wrap using head reference point
            head_trial_center = np.array([x, y, z_head])
            coords_trial = wrap_residue_as_whole(coords_trial, head_trial_center, box)
            head_trial_center = wrap_residue_as_whole(head_trial_center[None,:], head_trial_center, box)[0]

            # --- Prevent water placement: all heavy atoms must stay between phosphate planes with keepout
            heavy_mask = np.isin(res_atoms.indices, res_heavy.indices)
            coords_heavy_trial = coords_trial[heavy_mask]
            zmin = coords_heavy_trial[:,2].min()
            zmax = coords_heavy_trial[:,2].max()
            if (zmin < (lower_z + z_keepout)) or (zmax > (upper_z - z_keepout)):
                continue

            # --- Protein groove exclusion: head center distance to protein heavy
            if protein_heavy.n_atoms > 0:
                d_head_prot = distance_array(head_trial_center[None,:], protein_heavy.positions, box=box).min()
                if d_head_prot < ring_excl:
                    continue

            # --- Clash checks (separate)
            # 1) strict vs protein + already placed UQ8
            if protein_heavy.n_atoms > 0:
                if distance_array(coords_heavy_trial, protein_heavy.positions, box=box).min() < min_prot:
                    continue

            if placed_uq8_heavy_indices:
                placed_uq8_heavy = u.atoms[placed_uq8_heavy_indices]
                if distance_array(coords_heavy_trial, placed_uq8_heavy.positions, box=box).min() < min_prot:
                    continue

            # 2) softer vs lipids
            if lipids_heavy.n_atoms > 0:
                if distance_array(coords_heavy_trial, lipids_heavy.positions, box=box).min() < min_lip:
                    continue

            # Accept
            res_atoms.positions = coords_trial
            # Wrap as a residue (should be redundant now, but harmless)
            res_atoms.wrap(compound="residues")

            placed_uq8_heavy_indices.extend(res_heavy.indices.tolist())
            placed = True
            if attempt > 1:
                print(f"[INFO] Placed {args.resname} resid {res.resid} in {attempt} attempts (head mode: {head_mode})")
            break

        if not placed:
            print(f"[WARN] Failed to place {args.resname} resid {res.resid} after {args.max_attempts} attempts.")
            print("       Suggestions: lower --min-clash-lip-nm (0.18→0.15), "
                  "reduce --z-keepout-nm (0.20→0.15), increase --max-attempts, "
                  "or check head atom naming via --head-names.")

    with mda.coordinates.GRO.GROWriter(args.output, n_atoms=u.atoms.n_atoms) as w:
        w.write(u.atoms)

    print(f"[DONE] Wrote randomized configuration to {args.output}")

if __name__ == "__main__":
    main()
