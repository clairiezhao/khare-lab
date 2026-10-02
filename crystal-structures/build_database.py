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

# create a new folder for each structure
# store each generated ASU as its own file
# 3. maybe write to file in folder
# 2. find all residues within 15A

import os
import requests
import gemmi

SEARCH_API_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
DOWNLOAD_URL = "https://files.rcsb.org/download/{pdb_id}.cif"

RAW_DIR = "data/raw_cif"
SHELL_DIR = "data/contact_shell_cif"

os.makedirs(RAW_DIR, exist_ok=True)
os.makedirs(SHELL_DIR, exist_ok=True)


def get_already_processed_ids() -> set[str]:
    """Scans the output directory to find IDs that already have a completed shell file."""
    processed = set()
    for fname in os.listdir(SHELL_DIR):
        if fname.endswith("_shell.cif"):
            pdb_id = fname.split("_shell.cif")[0].upper()
            processed.add(pdb_id)
    return processed


def query_rcsb_page(start: int, rows: int = 1000) -> list[str]:
    """Fetches a single page of crystal protein PDB IDs from the Search API v2."""
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
    
    if not response.ok:
        print("RCSB Error Response:", response.text)
        response.raise_for_status()

    results = response.json().get("result_set", [])
    return [entry["identifier"].upper() for entry in results]


def download_cif(pdb_id: str) -> str:
    """Downloads the standard ASU mmCIF file if not already present."""
    filepath = os.path.join(RAW_DIR, f"{pdb_id}.cif")
    if not os.path.exists(filepath):
        res = requests.get(DOWNLOAD_URL.format(pdb_id=pdb_id))
        res.raise_for_status()
        with open(filepath, "w", encoding="utf-8") as f:
            f.write(res.text)
    return filepath


def generate_first_contact_shell(
    input_path: str,
    output_path: str,
    cutoff_distance: float = 4.5,
) -> str:

    st = gemmi.read_structure(input_path)
    st.remove_waters()
    st.remove_hydrogens()

    model = st[0]
    cell = st.cell

    ns = gemmi.NeighborSearch(
        model,
        cell,
        cutoff_distance,
    ).populate(include_h=False)

    # Store unique symmetry images as:
    # (symmetry operation, PBC translation)
    images = set()

    for chain in model:
        for residue in chain:
            for atom in residue:

                for mark in ns.find_neighbors(
                    atom,
                    min_dist=0.1,
                    max_dist=cutoff_distance,
                ):

                    # Original atom corresponding to this Mark
                    original_atom = mark.to_cra(model).atom

                    # Find the actual PBC image of the symmetry-related
                    # atom that is near our central atom.
                    nearest = cell.find_nearest_pbc_image(
                        atom.pos,
                        original_atom.pos,
                        mark.image_idx,
                    )

                    # PBC translation
                    shift = (
                        nearest.pbc_shift[0],
                        nearest.pbc_shift[1],
                        nearest.pbc_shift[2],
                    )

                    key = (mark.image_idx, shift)

                    # Don't include the central ASU itself.
                    if mark.image_idx == 0 and shift == (0, 0, 0):
                        continue

                    images.add(key)

    # Create output structure.
    out = gemmi.Structure()
    out.cell = cell
    out.spacegroup_hm = st.spacegroup_hm

    out.add_model(gemmi.Model("1"))
    out_model = out[0]


    # Add central ASU.
    for chain in model:
        out_model.add_chain(chain.clone())

    # Add each contacting symmetry-related ASU.
    for image_idx, shift in images:

        transform = ns.get_image_transformation(image_idx)

        for chain in model:

            new_chain = out_model.add_chain(
                chain,
                unique_name=True,
            )

            for residue in new_chain:
                for atom in residue:

                    frac = transform.apply(
                        cell.fractionalize(atom.pos)
                    )

                    frac = gemmi.Fractional(
                        frac.x + shift[0],
                        frac.y + shift[1],
                        frac.z + shift[2],
                    )

                    atom.pos = cell.orthogonalize(frac)

    out.assign_serial_numbers(numbered_ter=False)
    out.make_mmcif_document().write_file(output_path)

    print(f"Cutoff: {cutoff_distance}")
    print(f"Number of contacting ASUs: {len(images)}")

    for image_idx, shift in sorted(images):
        print(image_idx, shift)
    
    return output_path


def process_unprocessed_structures(target_new_count: int = 5, page_size: int = 500):
    """Fetches, filters, and generates contact shells with automatic resumption."""
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
                raw_file = download_cif(pdb_id)
                output_file = os.path.join(SHELL_DIR, f"{pdb_id}_shell.cif")
                generate_first_contact_shell(raw_file, output_file, cutoff_distance=15)

                print(f"[{newly_processed + 1}] Successfully generated shell for {pdb_id}")
                processed_ids.add(pdb_id)
                newly_processed += 1

                if target_new_count is not None and newly_processed >= target_new_count:
                    print(f"Reached target batch size of {target_new_count} new entries.")
                    return

            except Exception as err:
                print(f"Skipping {pdb_id} due to error: {err}")

        start_offset += page_size


if __name__ == "__main__":
    # Processes 5 new entries per run; set target_new_count=None to run without an upper bound
    process_unprocessed_structures(target_new_count=10)