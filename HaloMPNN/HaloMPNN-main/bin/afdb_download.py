import requests
import sys
import os
import json
import pandas as pd

# First argument: input csv with columns pdb_url, afdb_id
# Second argument: output directory
# Example usage:
# afdb_download.py extreme_repr_afdb_20K_sample.csv afdb/

def download_pdb(pdb_url, output_file):
    """Download the PDB file from the provided URL.

    Args:
        pdb_url: The URL to download the PDB file from.
        output_file: The path to the output file.
    """

    response = requests.get(pdb_url)
    if response.status_code == 200:
        with open(output_file, "wb") as f:
            f.write(response.content)
    else:
        print(f"Failed to download PDB file from {pdb_url}.")


def main():
    input_csv=sys.argv[1]
    output_path=sys.argv[2]
    df = pd.read_csv(input_csv)
    for index, row in df.iterrows():
        if not os.path.exists(output_path + row['afdb_id'] + '.pdb'):
            download_pdb(row['pdb_url'], output_path + row['afdb_id'] + '.pdb')

if __name__ == "__main__":
    main()
