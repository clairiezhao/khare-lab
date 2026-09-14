#!/usr/bin/env python3
"""
filter_metadata_repr.py  REPRIDS.txt  METADATA.csv  REPRMETADATA.csv

Filters a metadata CSV to only rows whose uniprot_id appears in the
representative-sequence ID list produced by MMseqs2.

Arguments (positional):
    1. REPRIDS.txt      - text file of representative UniProt IDs, one per line
    2. METADATA.csv     - input metadata CSV (columns: uniprot_id, afdb_id,
                          pdb_url, mean_plddt, length — or any superset)
    3. REPRMETADATA.csv - output CSV filtered to representative sequences
"""

import sys
import pandas as pd


def main():
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)

    repr_ids_path, metadata_path, output_path = sys.argv[1], sys.argv[2], sys.argv[3]

    # Load representative IDs
    with open(repr_ids_path) as fh:
        repr_ids = {line.strip() for line in fh if line.strip()}

    # Load metadata
    df = pd.read_csv(metadata_path)
    n_input = len(df)

    # Filter
    df_repr = df[df["uniprot_id"].isin(repr_ids)].copy()
    n_repr = len(df_repr)
    n_removed = n_input - n_repr

    # IDs in repr list but absent from metadata (informational)
    missing = repr_ids - set(df["uniprot_id"])

    df_repr.to_csv(output_path, index=False)

    print(f"Input metadata rows:        {n_input}")
    print(f"Representative IDs loaded:  {len(repr_ids)}")
    print(f"Rows kept:                  {n_repr}")
    print(f"Rows removed:               {n_removed}")
    if missing:
        print(f"IDs in repr list but not in metadata: {len(missing)}")
        for m in sorted(missing):
            print(f"  {m}")
    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()
