import os
import requests
import gemmi

SEARCH_API_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
DOWNLOAD_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"

RAW_DIR = "data/raw_pdb"
LATTICE_DIR = "data/lattice_pdb"

os.makedirs(RAW_DIR, exist_ok=True)
os.makedirs(LATTICE_DIR, exist_ok=True)


def get_already_processed_ids() -> set[str]:
    """Scans the output directory to find IDs that already have a completed lattice file."""
    processed = set()
    for fname in os.listdir(LATTICE_DIR):
        if fname.endswith("_lattice.pdb"):
            pdb_id = fname.split("_lattice.pdb")[0].upper()
            processed.add(pdb_id)
    return processed

def query_rcsb_page(start: int, rows: int = 1000) -> list[str]:
    """Fetches a single page of crystal protein PDB IDs from the Search API."""
    query_payload = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "exptl.method",
                        "operator": "exact_match",
                        "value": "X-RAY DIFFRACTION",
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "entity_poly.rcsb_entity_polymer_type",
                        "operator": "exact_match",
                        "value": "Protein",
                    },
                },
            ],
        },
        "request_options": {
            "paginate": {"start": start, "rows": rows}
        },
        "return_type": "entry",
    }

    response = requests.post(SEARCH_API_URL, json=query_payload)
    if response.status_code == 204:
        return []
    
    # If there's still an error, print the server's exact explanation
    if not response.ok:
        print("RCSB Error Response:", response.text)
        response.raise_for_status()

    results = response.json().get("result_set", [])
    return [entry["identifier"].upper() for entry in results]

def download_pdb(pdb_id: str) -> str:
    """Downloads the standard ASU PDB file if not already present."""
    filepath = os.path.join(RAW_DIR, f"{pdb_id}.pdb")
    if not os.path.exists(filepath):
        res = requests.get(DOWNLOAD_URL.format(pdb_id=pdb_id))
        res.raise_for_status()
        with open(filepath, "w", encoding="utf-8") as f:
            f.write(res.text)
    return filepath

# generate units in contact with center unit (15A)
# get general statistics:
# which residues are at contact points
# 1. residue frequency along entire dataset vs residue frequency along just contact points
# 2. how many chain contact points are there per ASU
# 3. what are the different types of symmetries (quantify different types of symmetries: c1, c2, d1, etc)
# 4. how much surface area each contact is
# solvent assessability surface area for each patch of contact
# training dataset: CATH families for proteins: recreate with just crystal proteins
# filter for redundancy in dataset
def generate_unit_cell_lattice(pdb_id: str, input_path: str) -> str:
    """Builds the full unit-cell assembly using space group symmetry operations."""
    structure = gemmi.read_structure(input_path)
    output_path = os.path.join(LATTICE_DIR, f"{pdb_id}_lattice.pdb")

    cell = structure.cell
    sg = structure.find_spacegroup()
    if not sg:
        raise ValueError(f"No valid space group found for {pdb_id}")

    lattice_struct = gemmi.Structure()
    lattice_struct.cell = cell
    lattice_struct.spacegroup_hm = structure.spacegroup_hm

    lattice_model = gemmi.Model("1")
    chain_idx = 0
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"

    for op in sg.operations():
        for model in structure:
            for chain in model:
                new_chain = chain.clone()
                new_chain.name = (
                    alphabet[chain_idx % len(alphabet)]
                    + str(chain_idx // len(alphabet) or "")
                )
                chain_idx += 1

                for residue in new_chain:
                    for atom in residue:
                        frac = cell.fractionalize(atom.pos)
                        sym_frac = op.apply_to_xyz([frac.x, frac.y, frac.z])
                        atom.pos = cell.orthogonalize(gemmi.Fractional(*sym_frac))

                lattice_model.add_chain(new_chain)

    lattice_struct.add_model(lattice_model)
    lattice_struct.write_minimal_pdb(output_path)
    return output_path


def process_unprocessed_structures(target_new_count: int = 5, page_size: int = 500):
    """
    Paginates through RCSB, skips already processed entries, and stops once
    target_new_count new structures have been processed.
    Set target_new_count=None to process all remaining structures.
    """
    processed_ids = get_already_processed_ids()
    print(f"Found {len(processed_ids)} already processed structures.")

    start_offset = 0
    newly_processed = 0

    while True:
        ids_in_page = query_rcsb_page(start=start_offset, rows=page_size)
        if not ids_in_page:
            print("No more entries returned by RCSB Search API.")
            break

        for pdb_id in ids_in_page:
            if pdb_id in processed_ids:
                continue

            try:
                raw_file = download_pdb(pdb_id)
                lattice_file = generate_unit_cell_lattice(pdb_id, raw_file)
                print(f"[{newly_processed + 1}] Successfully generated: {lattice_file}")
                processed_ids.add(pdb_id)
                newly_processed += 1

                if target_new_count is not None and newly_processed >= target_new_count:
                    print(f"Reached target batch size of {target_new_count} new entries.")
                    return

            except Exception as err:
                print(f"Skipping {pdb_id} due to error: {err}")

        start_offset += page_size


if __name__ == "__main__":
    # Test with 5 new structures per run; set to None when you want it to run continuously
    process_unprocessed_structures(target_new_count=5)