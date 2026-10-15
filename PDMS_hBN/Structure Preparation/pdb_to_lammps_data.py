#!/usr/bin/env python3
import MDAnalysis as mda
import numpy as np
import os
import sys
from io import StringIO
import datetime

# =============================================================================
# CONFIGURATION & PARAMETERS
# =============================================================================
pdb_file_name = "50%Rectangular_Graphene_pdms_composite.pdb"
output_data_file = "50%Rectangular_Graphene_pdms_composite.data"

# According to your instructions:
# All atoms up to serial ID 3600 belong to the sheet (h-BN + Graphene).
# Atoms with serial ID > 3600 belong to the polymer matrix (PDMS).
DOPED_CARBON_ID_RANGE = (1, 3600)

# Atom properties for the full composite
# Logic:
# Type 1: PDMS Oxygen (O)
# Type 2: PDMS Silicon (Si)
# Type 3: PDMS Methyl Carbon (UA) -> Carbons with serial ID > 3600
# Type 4: Boron (hBN) -> Boron with serial ID <= 3600
# Type 5: Nitrogen (hBN) -> Nitrogen with serial ID <= 3600
# Type 6: Doped Carbon -> Carbon with serial ID <= 3600
atom_properties = {
    'O': {'type_id': 1, 'mass': 15.9994, 'charge': -0.3, 'desc': 'PDMS Oxygen'},
    'Si': {'type_id': 2, 'mass': 28.0855, 'charge': 0.3, 'desc': 'PDMS Silicon'},
    'C': {'type_id': 3, 'mass': 15.0344, 'charge': 0.0, 'desc': 'PDMS Methyl (UA)'},
    'B': {'type_id': 4, 'mass': 10.8110, 'charge': 0.0, 'desc': 'Boron (hBN)'},
    'N': {'type_id': 5, 'mass': 14.0070, 'charge': 0.0, 'desc': 'Nitrogen (hBN)'},
    'GC': {'type_id': 6, 'mass': 12.0110, 'charge': 0.0, 'desc': 'Doped/Sheet Carbon'}
}

print("=" * 80)
print(f"STARTING CONVERSION: {datetime.datetime.now()}")
print(f"Target PDB: {pdb_file_name}")
print(f"Output: {output_data_file}")
print("=" * 80)


# =============================================================================
# PDB PRE-PROCESSING (STRICT FORMAT RECONSTRUCTION)
# =============================================================================
def preprocess_pdb_content_robust(raw_pdb_lines):
    """
    Cleans and reformats HETATM/ATOM lines to strict PDB format for MDAnalysis.
    Uses serial ID <= 3600 to identify Graphene and h-BN sheet atoms, and
    serial ID > 3600 to identify polymer matrix (PDMS) atoms.
    """
    processed_lines = []
    skipped_lines = 0
    atom_count = 0

    print(f"Step 1: Re-processing PDB lines with sheet/matrix boundary (Cutoff: {DOPED_CARBON_ID_RANGE[1]})...")

    for line_num, line in enumerate(raw_pdb_lines, 1):
        if line.startswith("HETATM") or line.startswith("ATOM"):
            try:
                # Extract raw data using fixed column widths
                record_name = line[0:6].strip()
                atom_serial_str = line[6:11].strip()
                atom_serial = int(atom_serial_str) if atom_serial_str else 0
                atom_name_raw = line[12:16].strip().upper()
                res_name = line[17:20].strip().upper() or 'RES'

                # --- SHEET VS MATRIX CLASSIFICATION BY SERIAL ID ---
                if DOPED_CARBON_ID_RANGE[0] <= atom_serial <= DOPED_CARBON_ID_RANGE[1]:
                    # Sheet Atoms (h-BN + Graphene)
                    if 'C' in atom_name_raw:
                        atom_name = 'GC'
                    elif 'B' in atom_name_raw:
                        atom_name = 'B'
                    elif 'N' in atom_name_raw:
                        atom_name = 'N'
                    else:
                        atom_name = atom_name_raw
                else:
                    # Polymer Matrix Atoms (PDMS)
                    if 'SI' in atom_name_raw:
                        atom_name = 'Si'
                    elif 'O' in atom_name_raw:
                        atom_name = 'O'
                    elif 'C' in atom_name_raw:
                        atom_name = 'C'  # PDMS Methyl (UA)
                    else:
                        atom_name = atom_name_raw

                # Formatting cleanup
                if atom_name == 'SI': atom_name = 'Si'

                chain_id = line[21].strip() or 'A'
                res_seq_str = line[22:26].strip()
                res_seq = int(res_seq_str) if res_seq_str else 1

                # Coordinates
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())

                # Element for MDAnalysis internal parser
                element_symbol = atom_name if len(atom_name) <= 2 else atom_name[0]

                # Reconstruct strict PDB line
                new_line = (
                    f"{record_name:<6}{atom_serial:>5} {atom_name:<4} "
                    f"{res_name:<3} {chain_id:<1}{res_seq:>4}    "
                    f"{x:>8.3f}{y:>8.3f}{z:>8.3f}"
                    f"  1.00  0.00          {element_symbol:>2}"
                )
                processed_lines.append(new_line)
                atom_count += 1
            except (ValueError, IndexError):
                skipped_lines += 1
                continue
        elif any(line.startswith(tag) for tag in ["CONECT", "HEADER", "TITLE", "TER", "END"]):
            processed_lines.append(line.rstrip())

    print(f"   - Total Atoms Processed: {atom_count}")
    print(f"   - Malformed Lines Skipped: {skipped_lines}")
    return "\n".join(processed_lines)


# =============================================================================
# DATA EXTRACTION & ANALYSIS
# =============================================================================
if not os.path.exists(pdb_file_name):
    print(f"CRITICAL ERROR: File '{pdb_file_name}' not found.")
    sys.exit(1)

with open(pdb_file_name, 'r') as f:
    lines = f.readlines()

processed_pdb = preprocess_pdb_content_robust(lines)
u = mda.Universe(StringIO(processed_pdb), format='PDB')

# 1. Process Atoms
print("Step 2: Assigning Atom Types and Molecular IDs...")
atoms_data = []
molecule_map = {}
mol_counter = 1

for a in u.atoms:
    name = a.name.strip()
    res_key = (a.residue.resname, a.residue.resid)
    if res_key not in molecule_map:
        molecule_map[res_key] = mol_counter
        mol_counter += 1

    prop = atom_properties.get(name, {'type_id': 99, 'mass': 1.0, 'charge': 0.0})

    atoms_data.append({
        'id': a.id,
        'mol': molecule_map[res_key],
        'type': prop['type_id'],
        'q': prop['charge'],
        'x': a.position[0],
        'y': a.position[1],
        'z': a.position[2],
        'name': name
    })

atoms_data.sort(key=lambda x: x['id'])
id_to_name = {a['id']: a['name'] for a in atoms_data}

# 2. Process Bonds
print("Step 3: Generating Bond Topology from CONECT...")
bonds_list = []
bond_types = {}
bond_type_count = 0

for b in u.bonds:
    id1, id2 = b.atoms[0].id, b.atoms[1].id
    if id1 not in id_to_name or id2 not in id_to_name: continue

    n1, n2 = id_to_name[id1], id_to_name[id2]
    b_key = tuple(sorted((n1, n2)))

    if b_key not in bond_types:
        bond_type_count += 1
        bond_types[b_key] = bond_type_count

    bonds_list.append({
        'id': len(bonds_list) + 1,
        'type': bond_types[b_key],
        'a1': id1,
        'a2': id2,
        'key': b_key
    })

# 3. Process Angles
print("Step 4: Inferring Angles from Bonds...")
angles_list = []
angle_types = {}
angle_type_count = 0
adj = {a['id']: [] for a in atoms_data}
for b in bonds_list:
    adj[b['a1']].append(b['a2'])
    adj[b['a2']].append(b['a1'])

for mid in adj:
    neigh = adj[mid]
    if len(neigh) < 2: continue
    m_name = id_to_name[mid]
    for i in range(len(neigh)):
        for j in range(i + 1, len(neigh)):
            a1, a3 = neigh[i], neigh[j]
            n1, n3 = id_to_name[a1], id_to_name[a3]
            a_key = tuple(sorted((n1, n3))) + (m_name,)

            if a_key not in angle_types:
                angle_type_count += 1
                angle_types[a_key] = angle_type_count

            angles_list.append({
                'id': len(angles_list) + 1,
                'type': angle_types[a_key],
                'a1': a1, 'a2': mid, 'a3': a3,
                'key': a_key
            })

# 4. Process Dihedrals
print("Step 5: Inferring Dihedrals...")
dihedrals_list = []
dihedral_types = {}
dihedral_type_count = 0

for b in bonds_list:
    j, k = b['a1'], b['a2']
    nj, nk = id_to_name[j], id_to_name[k]
    for i in adj[j]:
        if i == k: continue
        ni = id_to_name[i]
        for l in adj[k]:
            if l == j or l == i: continue
            nl = id_to_name[l]
            d_key = tuple(sorted((ni, nl))) + tuple(sorted((nj, nk)))

            if d_key not in dihedral_types:
                dihedral_type_count += 1
                dihedral_types[d_key] = dihedral_type_count

            dihedrals_list.append({
                'id': len(dihedrals_list) + 1,
                'type': dihedral_types[d_key],
                'a1': i, 'a2': j, 'a3': k, 'a4': l,
                'key': d_key
            })

# =============================================================================
# BOX CALCULATION
# =============================================================================
coords = np.array([[a['x'], a['y'], a['z']] for a in atoms_data])
min_p = np.min(coords, axis=0) - 10.0
max_p = np.max(coords, axis=0) + 10.0
box = max_p - min_p

# =============================================================================
# WRITE LAMMPS DATA FILE
# =============================================================================
print(f"Step 6: Writing final data file: {output_data_file}")
with open(output_data_file, 'w') as f:
    f.write(f"LAMMPS data file - Doped hBN/PDMS Composite - Created {datetime.datetime.now()}\n\n")

    f.write(f"{len(atoms_data)} atoms\n")
    f.write(f"{len(bonds_list)} bonds\n")
    f.write(f"{len(angles_list)} angles\n")
    f.write(f"{len(dihedrals_list)} dihedrals\n\n")

    f.write(f"{len(atom_properties)} atom types\n")
    f.write(f"{len(bond_types)} bond types\n")
    f.write(f"{len(angle_types)} angle types\n")
    f.write(f"{dihedral_type_count} dihedral types\n\n")

    f.write(f"{-box[0] / 2:.6f} {box[0] / 2:.6f} xlo xhi\n")
    f.write(f"{-box[1] / 2:.6f} {box[1] / 2:.6f} ylo yhi\n")
    f.write(f"{-box[2] / 2:.6f} {box[2] / 2:.6f} zlo zhi\n\n")

    f.write("Masses\n\n")
    for name, p in sorted(atom_properties.items(), key=lambda x: x[1]['type_id']):
        f.write(f"{p['type_id']} {p['mass']:.4f}  # {name} ({p['desc']})\n")

    f.write("\nAtoms # full\n\n")
    for a in atoms_data:
        f.write(f"{a['id']} {a['mol']} {a['type']} {a['q']:.6f} {a['x']:.6f} {a['y']:.6f} {a['z']:.6f} # {a['name']}\n")

    f.write("\nBonds\n\n")
    for b in bonds_list:
        f.write(f"{b['id']} {b['type']} {b['a1']} {b['a2']} # {b['key'][0]}-{b['key'][1]}\n")

    f.write("\nAngles\n\n")
    for g in angles_list:
        f.write(f"{g['id']} {g['type']} {g['a1']} {g['a2']} {g['a3']} # {g['key'][0]}-{g['key'][2]}-{g['key'][1]}\n")

    f.write("\nDihedrals\n\n")
    for d in dihedrals_list:
        n1, n2, n3, n4 = id_to_name[d['a1']], id_to_name[d['a2']], id_to_name[d['a3']], id_to_name[d['a4']]
        f.write(f"{d['id']} {d['type']} {d['a1']} {d['a2']} {d['a3']} {d['a4']} # {n1}-{n2}-{n3}-{n4}\n")

# =============================================================================
# FINAL SUMMARY REPORT
# =============================================================================
print("\n" + "=" * 80)
print("FINAL REPORT")
print("=" * 80)
print(f"Total Atoms:     {len(atoms_data)}")
print(f"Total Bonds:     {len(bonds_list)}")
print(f"Total Angles:    {len(angles_list)}")
print(f"Total Dihedrals: {len(dihedrals_list)}")
print("-" * 30)
print("Atom Type Distribution:")
counts = {k: 0 for k in atom_properties.keys()}
for a in atoms_data:
    if a['name'] in counts:
        counts[a['name']] += 1
for k, v in counts.items():
    print(f"  - {k} (Type {atom_properties[k]['type_id']}): {v} atoms")

print("-" * 30)
print("Topology Info:")
print(f"  - Bond Types:     {len(bond_types)}")
print(f"  - Angle Types:    {len(angle_types)}")
print(f"  - Dihedral Types: {dihedral_type_count}")
print("=" * 80)
print(f"SUCCESS: {output_data_file} generated.")