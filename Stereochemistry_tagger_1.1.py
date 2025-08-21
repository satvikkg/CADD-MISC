# ==========================================================================================
# Script: SDF Preprocessor with Stereochemistry & Conformer Tagging (+ Optional SDF.GZ Export)
#
# DESCRIPTION:
# Processes an SDF from Schrödinger Maestro’s LigPrep (supports `.sdf`, `.sdf.gz`, `.sdfgz`)
# and produces deterministic, idempotent compound titles that include stereochemistry and
# conformer indices. Key behaviors:
#
#   1) Standardize names
#      - Trims whitespace and converts the SDF record name `_Name` to UPPERCASE for grouping.
#
#   2) Detect stereochemistry (deterministically)
#      - Explicitly calls RDKit AssignStereochemistry (works even with sanitize=False input).
#      - Builds chiral-center tags in stable atom-index order (e.g., `_4R_8S`).
#
#   3) Assign conformer identifiers
#      - Groups by base name + stereo signature and adds `_CONF1`, `_CONF2`, … as needed.
#
#   4) Export tagged metadata to CSV
#      - First column is the updated `Title` (mirrors `_Name`); subsequent columns are all SD tags.
#
#   5) Optional SDF.GZ output
#      - Writes a gzipped SDF (`.sdf.gz`) preserving original 3D coordinates and SD tags,
#        with the updated titles embedded (so Maestro shows the tagged names).
#
# KEY PROPERTIES:
#   - Idempotent: Skips re-tagging entries that already end with `_#R/_#S` (optionally `_CONF#`),
#                 regardless of case (e.g., `_7r_10s` is treated as tagged).
#   - Deterministic: Chiral-center suffix order is fixed by atom index.
#   - Non-destructive: Coordinates and SD tags are copied verbatim; only `_Name` and the `Title`
#                      SD tag are overwritten to keep Maestro and CSV in sync.
#
# INPUT:
#   - `.sdf`, `.sdf.gz`, or `.sdfgz` (the latter is copied to `.sdf.gz` for reading).
#
# OUTPUT:
#   - CSV: `<output_name>.csv` with `Title` + all SD tags found.
#   - (Optional) SDF.GZ: `<output_name>.sdf.gz` with updated titles (same 3D coords).
#
# DEPENDENCIES:
#   - RDKit, pandas
#
# USAGE EXAMPLE:
#   $ Stereochemistry_tagger.py
#       Enter the full path to your SDF file (.sdf, .sdf.gz, .sdfgz): /path/to/file.sdfgz
#       Enter output CSV filename (without extension): tagged_output
#       Also save updated SDF.GZ file with tagged titles? [y/n]: y
#
# NOTES:
#   - Reading uses `sanitize=False` to avoid altering structures; stereochemistry is then
#     assigned explicitly before detection.
#   - Overwriting `_Name` and `Title` ensures Maestro and CSV show identical tagged names.
#
# AUTHOR:
#   In Young Bae (with ChatGPT, OpenAI)
# ==========================================================================================

from rdkit import Chem
import pandas as pd
import os
import gzip
import shutil
import re
from collections import defaultdict

# ---------------------------------------------------------------------------
# User input & path normalization
# ---------------------------------------------------------------------------
# Accept .sdf / .sdf.gz / .sdfgz; normalize Windows backslashes and stray quotes.
sdf_path = input("Enter the full path to your SDF file (.sdf, .sdf.gz, .sdfgz): ").strip().replace("\\", "/").strip('"')

# If input is Maestro's single-extension gzip (.sdfgz), copy to .sdf.gz so RDKit can read it.
if sdf_path.endswith(".sdfgz"):
    new_path = sdf_path[:-6] + ".sdf.gz"
    shutil.copyfile(sdf_path, new_path)               # byte-for-byte copy; original remains untouched
    sdf_path = new_path
    print(f"🛠 Converted .sdfgz to .sdf.gz: {sdf_path}")

# ---------------------------------------------------------------------------
# Read molecules (non-sanitizing, pass-through of coords & tags)
# ---------------------------------------------------------------------------
# removeHs=False preserves explicit hydrogens; sanitize=False avoids graph modifications.
if sdf_path.endswith(".gz"):
    with gzip.open(sdf_path, 'rb') as gzfile:
        suppl = Chem.ForwardSDMolSupplier(gzfile, removeHs=False, sanitize=False)
        mols = [mol for mol in suppl if mol is not None]
else:
    suppl = Chem.SDMolSupplier(sdf_path, removeHs=False, sanitize=False)
    mols = [mol for mol in suppl if mol is not None]

# ---------------------------------------------------------------------------
# Normalize titles: trim + uppercase for consistent grouping
# ---------------------------------------------------------------------------
for mol in mols:
    if mol is not None and mol.HasProp('_Name'):
        # Trim leading/trailing whitespace then uppercase (e.g., " chembl123 " -> "CHEMBL123")
        mol.SetProp('_Name', mol.GetProp('_Name').strip().upper())

# ---------------------------------------------------------------------------
# Group molecules by current title
# ---------------------------------------------------------------------------
grouped = defaultdict(list)
for mol in mols:
    if mol is not None and mol.HasProp('_Name'):
        grouped[mol.GetProp('_Name')].append(mol)

# Idempotency: consider a record "already tagged" if it ends with _<num><R/S> repeated
# and an optional _CONF<num>. Case-insensitive so `_7r_10s` is also treated as tagged.
_ALREADY_TAGGED_RE = re.compile(r'(?:_\d+[RS])+(?:_CONF\d+)?$', re.IGNORECASE)

# ---------------------------------------------------------------------------
# Assign new titles (stereo + conformer tags)
# ---------------------------------------------------------------------------
new_title_map = {}   # maps id(mol) -> new title string

for title, group in grouped.items():
    # If this group's title already carries a stereo suffix (and maybe _CONF#), don't re-tag.
    if _ALREADY_TAGGED_RE.search(title):
        for mol in group:
            new_title_map[id(mol)] = title
        continue

    # Compute stereochemistry for each molecule in the group.
    # With sanitize=False, we must explicitly assign stereochemistry.
    stereo_tags = []
    for mol in group:
        Chem.AssignStereochemistry(mol, force=True, flagPossibleStereoCenters=True)
        # Find defined centers (skip unassigned '?'); return list of (atom_idx, 'R'/'S')
        centers = Chem.FindMolChiralCenters(
            mol, includeUnassigned=False, useLegacyImplementation=False
        )
        centers = sorted(centers, key=lambda x: x[0])      # stable order by atom index
        stereo_tags.append(tuple(centers))                 # hashable signature for grouping

    # Split the group by stereo signature so CONF indices apply within a stereo set
    stereo_to_mols = defaultdict(list)
    for i, tag in enumerate(stereo_tags):
        stereo_to_mols[tag].append(group[i])

    # Build final titles:
    #  - If stereo exists, append "_<idx><R/S>" blocks (uppercase R/S), e.g., "_4R_8S"
    #  - If multiple records share the same stereo, append "_CONF1", "_CONF2", ...
    #  - If no stereo exists but there are duplicates, append only "_CONF#"
    for stereo_tag, mol_list in stereo_to_mols.items():
        has_stereo = len(stereo_tag) > 0
        if has_stereo:
            # Force uppercase in case any source produced lowercase 'r'/'s'
            suffix = "_".join([f"{idx+1}{conf.upper()}" for idx, conf in stereo_tag])

        for i, mol in enumerate(mol_list, start=1):
            if has_stereo:
                new_title = f"{title}_{suffix}_CONF{i}" if len(mol_list) > 1 else f"{title}_{suffix}"
            else:
                new_title = f"{title}_CONF{i}" if len(mol_list) > 1 else title
            new_title_map[id(mol)] = new_title

# ---------------------------------------------------------------------------
# Flatten to a table and keep Maestro + CSV consistent
# ---------------------------------------------------------------------------
rows = []
for mol in mols:
    if mol is None:
        continue

    # Use the computed title when available; otherwise keep existing _Name or fallback
    new_title = new_title_map.get(id(mol), mol.GetProp('_Name') if mol.HasProp('_Name') else "UNKNOWN")

    # Overwrite both the record name and the Title SD tag so Maestro and CSV match
    mol.SetProp('_Name', new_title)
    mol.SetProp('Title', new_title)

    # Build a row of all SD tags for the CSV (Title first for readability)
    row = {'Title': new_title}
    for prop in mol.GetPropNames():
        # Skipping 'Title' here avoids redundant overwrite of the just-written value
        if prop == 'Title':
            continue
        row[prop] = mol.GetProp(prop)
    rows.append(row)

df = pd.DataFrame(rows)
# Ensure 'Title' stays as the first column
other_cols = [c for c in df.columns if c != 'Title']
df = df[['Title'] + other_cols]

# ---------------------------------------------------------------------------
# Write CSV next to the input file
# ---------------------------------------------------------------------------
output_name = input("Enter output CSV filename (without extension): ").strip()
if not output_name:
    output_name = "tagged_sdf_export"
output_dir = os.path.dirname(sdf_path)
csv_path = os.path.join(output_dir, f"{output_name}.csv")
df.to_csv(csv_path, index=False)
print(f"\n✅ Tagged metadata saved to CSV: {csv_path}")

# ---------------------------------------------------------------------------
# Optional: write updated SDF.GZ with the new titles
# ---------------------------------------------------------------------------
save_sdf = input("Also save updated SDF.GZ file with tagged titles? [y/n]: ").strip().lower().startswith("y")
if save_sdf:
    import tempfile

    sdf_gz_out_path = os.path.join(output_dir, f"{output_name}.sdf.gz")

    # RDKit's SDWriter doesn't write gzip directly. Write a temp .sdf, then gzip it.
    with tempfile.NamedTemporaryFile(delete=False, suffix=".sdf") as tmp:
        tmp_sdf = tmp.name

    writer = Chem.SDWriter(tmp_sdf)
    for mol in mols:
        writer.write(mol)
    writer.close()

    # Compress to .sdf.gz and remove the temp .sdf
    with open(tmp_sdf, "rb") as fin, gzip.open(sdf_gz_out_path, "wb") as fout:
        shutil.copyfileobj(fin, fout)
    os.remove(tmp_sdf)

    print(f"✅ Updated SDF.GZ file saved to: {sdf_gz_out_path}")
