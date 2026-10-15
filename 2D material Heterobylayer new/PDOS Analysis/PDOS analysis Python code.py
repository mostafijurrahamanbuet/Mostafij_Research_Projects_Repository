
# ─── Time formatter (global) ──────────────────────────────────────────────────
def fmt_time(seconds):
    seconds = int(seconds)
    if seconds < 60:
        return f"{seconds}s"
    elif seconds < 3600:
        return f"{seconds//60}m {seconds%60}s"
    else:
        return f"{seconds//3600}h {(seconds%3600)//60}m {seconds%60}s"


# ─── PDOS Function ────────────────────────────────────────────────────────────
def finding_PDOS_plain(path, initial_frame, Nc, n_atoms, dt,
                        current_block=1, total_block=1,
                        memmap_dir=r"C:\temp_memmap"):
    """
    Drop-in replacement for read_data() + finding_PDOS().
    Uses plain text (decompressed) dump file + numpy fast parsing.

    ADAPTED: this version expects "vx vy vz" only (3 columns per atom line),
    matching your LAMMPS dump command:
        dump vel_hBN hBN custom 50 <file> vx vy vz
    (no id/type columns -- your original friend's version expected 5 columns:
    "id type vx vy vz")

    Memory strategy  : memmap for power_* and vel_all (disk-backed)
    Peak RAM         : < 1 GB
    Speed vs .gz     : ~3-5x faster reading
    Speed vs Python  : ~3-5x faster parsing (numpy vs float())
    """

    os.makedirs(memmap_dir, exist_ok=True)

    N_fft        = 2**int(np.ceil(np.log2(2 * Nc)))
    n_freqs      = N_fft // 2 + 1
    window       = np.hanning(Nc)
    window_power = np.sum(window**2)

    # ── Create memory-mapped arrays on disk ───────────────────────────────
    print(f"\n{'='*60}")
    print(f"  Block {current_block}/{total_block} — Creating memmap files...")
    print(f"{'='*60}")

    vel_all = np.memmap(
        os.path.join(memmap_dir, f'vel_block{current_block}.dat'),
        dtype='float32', mode='w+',
        shape=(n_atoms, Nc, 3)
    )
    power_vx = np.memmap(
        os.path.join(memmap_dir, f'power_vx_block{current_block}.dat'),
        dtype='float32', mode='w+',
        shape=(n_atoms, n_freqs)
    )
    power_vy = np.memmap(
        os.path.join(memmap_dir, f'power_vy_block{current_block}.dat'),
        dtype='float32', mode='w+',
        shape=(n_atoms, n_freqs)
    )
    power_vz = np.memmap(
        os.path.join(memmap_dir, f'power_vz_block{current_block}.dat'),
        dtype='float32', mode='w+',
        shape=(n_atoms, n_freqs)
    )
    power_combined = np.memmap(
        os.path.join(memmap_dir, f'power_combined_block{current_block}.dat'),
        dtype='float32', mode='w+',
        shape=(n_atoms, n_freqs)
    )

    # ── Phase 1: Read frames ───────────────────────────────────────────────
    print(f"\n  [Phase 1/2] Reading {Nc} frames x {n_atoms} atoms...")

    lines_per_frame = n_atoms + 9
    skip_lines      = int(initial_frame) * lines_per_frame

    t_read_start = time.time()
    last_print_t = t_read_start
    frames_read  = 0
    first        = True
    last_frame   = 0

    with open(path, 'r') as f:

        # ── Skip to initial_frame efficiently ─────────────────────────────
        print(f"  Skipping to frame {int(initial_frame)}...")
        t_skip_start = time.time()
        for _ in range(skip_lines):
            f.readline()
        print(f"  Skip complete : {fmt_time(time.time() - t_skip_start)}")

        # ── Read Nc frames ─────────────────────────────────────────────────
        while frames_read < Nc:
            line = f.readline()
            if not line:
                print(f"  WARNING: reached end of file after only "
                      f"{frames_read}/{Nc} frames! Reduce Nc or "
                      f"initial_frame -- results for this block are INVALID.")
                break

            if line.strip() == "ITEM: TIMESTEP":
                last_frame = int(f.readline().split()[0])
                if first:
                    print(f"  Starting Frame : {last_frame}")
                    first = False

                # Skip 7 header lines
                for _ in range(7):
                    f.readline()

                # ── Read entire frame at once with numpy ───────────────────
                # Collect all atom lines for this frame
                raw_lines = ''.join(f.readline() for _ in range(n_atoms))

                # CHANGED: parse 3 columns (vx vy vz) instead of 5 (id type vx vy vz)
                # -- your dump file has no id/type columns, only velocities.
                frame_data = np.fromstring(raw_lines, sep=' ') \
                               .reshape(n_atoms, 3)

                # CHANGED: store all 3 columns directly (no [:, 2:] slice needed,
                # since there's no id/type to skip past anymore)
                vel_all[:, frames_read, :] = frame_data.astype('float32')

                frames_read += 1

                # ── Progress + ETA ─────────────────────────────────────────
                now        = time.time()
                time_since = now - last_print_t
                if time_since >= 30 or frames_read == Nc:
                    elapsed  = now - t_read_start
                    pct_done = frames_read / Nc
                    eta      = (elapsed / pct_done) * (1 - pct_done) if pct_done > 0 else 0
                    print(f"  Reading frames : {frames_read}/{Nc} "
                          f"({pct_done*100:.1f}%) | "
                          f"Elapsed: {fmt_time(elapsed)} | "
                          f"ETA: {fmt_time(eta)}")
                    last_print_t = now

    t_read_total = time.time() - t_read_start
    print(f"  Ending Frame   : {last_frame}")
    print(f"  Read complete  : {fmt_time(t_read_total)}")
    vel_all.flush()

    # ── Phase 2: FFT ──────────────────────────────────────────────────────
    print(f"\n  [Phase 2/2] Computing FFT for {n_atoms} atoms...")

    t_fft_start  = time.time()
    last_print_t = t_fft_start
    atom_times   = []

    for atom in range(n_atoms):
        t_atom_start = time.time()

        for arr, comp in [(power_vx, 0), (power_vy, 1), (power_vz, 2)]:
            # Pull single atom slice into RAM
            v  = np.array(vel_all[atom, :Nc, comp])
            C0 = np.mean(v**2)
            if C0 == 0:
                continue
            v_w   = v * window
            v_pad = np.pad(v_w, (0, N_fft - Nc), mode='constant')
            fft_v = np.fft.rfft(v_pad)
            arr[atom] = np.abs(fft_v)**2 / (window_power * C0)

        power_combined[atom] = (power_vx[atom] +
                                 power_vy[atom] +
                                 power_vz[atom])

        # Rolling average over last 200 atoms for stable ETA
        atom_times.append(time.time() - t_atom_start)
        if len(atom_times) > 200:
            atom_times.pop(0)

        now        = time.time()
        time_since = now - last_print_t
        if time_since >= 30 or atom == n_atoms - 1:
            elapsed      = now - t_fft_start
            atoms_done   = atom + 1
            atoms_left   = n_atoms - atoms_done
            avg_per_atom = sum(atom_times) / len(atom_times)
            eta          = avg_per_atom * atoms_left
            pct_done     = atoms_done / n_atoms
            print(f"  FFT atoms      : {atoms_done}/{n_atoms} "
                  f"({pct_done*100:.1f}%) | "
                  f"Elapsed: {fmt_time(elapsed)} | "
                  f"ETA: {fmt_time(eta)} | "
                  f"Speed: {1/avg_per_atom:.1f} atoms/s")
            last_print_t = now

    t_fft_total = time.time() - t_fft_start
    print(f"  FFT complete   : {fmt_time(t_fft_total)}")

    # Flush results to disk
    for arr in [power_vx, power_vy, power_vz, power_combined]:
        arr.flush()

    freqs         = np.fft.rfftfreq(N_fft, d=dt)
    pdos_vx       = np.mean(power_vx,       axis=0)
    pdos_vy       = np.mean(power_vy,       axis=0)
    pdos_vz       = np.mean(power_vz,       axis=0)
    pdos_combined = np.mean(power_combined, axis=0)

    # ── Convert to regular numpy arrays before closing memmap ────────────
    pdos_vx        = np.array(np.mean(power_vx,       axis=0))
    pdos_vy        = np.array(np.mean(power_vy,       axis=0))
    pdos_vz        = np.array(np.mean(power_vz,       axis=0))
    pdos_combined  = np.array(np.mean(power_combined, axis=0))
    out_power_vx       = np.array(power_vx)
    out_power_vy       = np.array(power_vy)
    out_power_vz       = np.array(power_vz)
    out_power_combined = np.array(power_combined)

    # ── Explicitly close + flush before deleting on Windows ──────────────
    vel_all._mmap.close()
    power_vx._mmap.close()
    power_vy._mmap.close()
    power_vz._mmap.close()
    power_combined._mmap.close()

    del vel_all, power_vx, power_vy, power_vz, power_combined

    # Force garbage collection — Windows needs this
    import gc
    gc.collect()

    # ── Now safe to delete files ──────────────────────────────────────────
    for fname in [f'vel_block{current_block}.dat',
                  f'power_vx_block{current_block}.dat',
                  f'power_vy_block{current_block}.dat',
                  f'power_vz_block{current_block}.dat',
                  f'power_combined_block{current_block}.dat']:
        fpath = os.path.join(memmap_dir, fname)
        if os.path.exists(fpath):
            try:
                os.remove(fpath)
            except PermissionError:
                # Wait and retry once — Windows sometimes needs a moment
                time.sleep(2)
                os.remove(fpath)

    # ── Block summary ─────────────────────────────────────────────────────
    t_block_total = t_read_total + t_fft_total
    print(f"\n  Block {current_block}/{total_block} complete")
    print(f"  Read : {fmt_time(t_read_total)} | "
          f"FFT  : {fmt_time(t_fft_total)} | "
          f"Total: {fmt_time(t_block_total)}")
    print(f"{'='*60}\n")

    return freqs, out_power_vx, out_power_vy, out_power_vz, out_power_combined, \
           pdos_vx, pdos_vy, pdos_vz, pdos_combined



# ─── Main Block ───────────────────────────────────────────────────────────────
# CHANGED: path now points to YOUR BP (phosphorene) dump file.
# Double-check this matches where the file actually sits on your machine --
# adjust the folder if your Python script isn't in the same folder as the .vc files.
path = r"D:\Another sum up of all simulation\DB3\3%\Y axis\dumpBP_hBN_NEMD.BP_DB3_y_axis_3_deg.vc"

# CHANGED: memmap_dir set to D:, which you confirmed has ~172 GB free -- plenty of room.
memmap_dir = r"D:\temp\mm"

# CHANGED: dt updated to match YOUR actual dump interval.
# Your LAMMPS dump writes every 50 steps, with timestep = 0.0001 ps,
# so time between frames = 50 * 0.0001 = 0.005 ps (NOT 0.00025, which was
# your friend's own simulation's setting -- it has nothing to do with yours).
dt = 0.005

# Read n_atoms from file header
with open(path, 'r') as f:
    for i in range(3):
        f.readline()
    n_atoms = int(f.readline().split()[0])
print(f"Number of Atoms : {n_atoms}")

# CHANGED: Nc and number_of_block set so the block/overlap scheme spans
# your FULL 400 ps trajectory (all 80,001 frames), per your request.
# WHY these numbers: Nc*(1+(blocks-1)/2) must not exceed your available
# frame count. With Nc=20000 and blocks=7: 20000*(1+6/2) = 80,000 frames
# needed, which fits within your 80,001-frame file almost exactly --
# covering 0 to 400 ps with 7 overlapping 100 ps blocks (50% overlap).
Nc              = 20000   # 20,000 frames x 0.005 ps = 100 ps per block
number_of_block = 7       # 7 overlapping blocks, 50% overlap, spans full 400 ps

# CHANGED: initial_frame starts at 0 -- the very beginning of the trajectory --
# since you want the full 400 ps run included, not just the steady-state tail.
# NOTE: this means the transient (non-steady-state) period at the start of
# your NEMD run is now mixed into the averaged spectrum, unlike the earlier
# steady-state-only version. That's a legitimate choice if you want a
# spectrum representative of the whole run -- just be aware the result
# reflects a blend of both regimes, not purely steady state.
initial_frame = 0

# Sanity check — print exact frames needed
total_frames = int(Nc * (1 + (number_of_block - 1) / 2))
print(f"Nc              : {Nc}")
print(f"Blocks          : {number_of_block}")
print(f"Frames needed   : {total_frames}  (starting at frame {initial_frame})")
print(f"Last frame used : {initial_frame + total_frames}  (file has 80001 total)")
print(f"Overlap         : 50%")

N_fft   = 2**int(np.ceil(np.log2(2 * Nc)))
n_freqs = N_fft // 2 + 1

pdos_vx_sum       = np.zeros(n_freqs)
pdos_vy_sum       = np.zeros(n_freqs)
pdos_vz_sum       = np.zeros(n_freqs)
pdos_combined_sum = np.zeros(n_freqs)

t_total_start = time.time()

for q in range(number_of_block):
    freqs, power_vx, power_vy, power_vz, power_combined, \
    pdos_vx, pdos_vy, pdos_vz, pdos_combined = finding_PDOS_plain(
        path, initial_frame, Nc, n_atoms, dt,
        current_block = q + 1,
        total_block   = number_of_block,
        memmap_dir    = memmap_dir
    )

    pdos_vx_sum       += pdos_vx
    pdos_vy_sum       += pdos_vy
    pdos_vz_sum       += pdos_vz
    pdos_combined_sum += pdos_combined

    initial_frame += Nc / 2

    # Overall progress
    blocks_done   = q + 1
    blocks_left   = number_of_block - blocks_done
    elapsed_total = time.time() - t_total_start
    avg_per_block = elapsed_total / blocks_done
    eta_total     = avg_per_block * blocks_left
    print(f">>> Overall : {blocks_done}/{number_of_block} blocks done | "
          f"Elapsed : {fmt_time(elapsed_total)} | "
          f"ETA : {fmt_time(eta_total)}\n")

pdos_vx       = pdos_vx_sum       / number_of_block
pdos_vy       = pdos_vy_sum       / number_of_block
pdos_vz       = pdos_vz_sum       / number_of_block
pdos_combined = pdos_combined_sum / number_of_block

# CHANGED: output filename renamed to reflect that this run is for BP (phosphorene)
# specifically (kept on D: since that's your drive with free space).
savepath = r"D:\PDOS\BP_DB3_y_axis_3_deg.npz"
os.makedirs(os.path.dirname(savepath), exist_ok=True)
np.savez(savepath,
         freqs    = freqs,
         vx       = pdos_vx,
         vy       = pdos_vy,
         vz       = pdos_vz,
         combined = pdos_combined)

total_time = time.time() - t_total_start
print(f"{'='*60}")
print(f"ALL BLOCKS COMPLETE")
print(f"Total time : {fmt_time(total_time)}")
print(f"Saved      : {savepath}")
print(f"{'='*60}")