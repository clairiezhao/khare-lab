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
    """Computes and saves the central ASU surrounded by its first crystal contact shell."""
    structure = gemmi.read_structure(input_path)
    cell = structure.cell
    sg = structure.find_spacegroup()
    if not sg:
        raise ValueError("Structure lacks space group or unit cell parameters.")

    # Strip waters and hydrogens so contacts represent direct protein-protein packing
    structure.remove_waters()
    structure.remove_hydrogens()

    # Index ASU atoms into a periodic spatial grid
    ns = gemmi.NeighborSearch(structure[0], cell, max_radius=cutoff_distance)
    ns.populate(include_h=False)

    # Search through space group operations and direct neighboring unit cells (-1, 0, 1)
    identity_op = gemmi.Op("x,y,z")
    contact_transforms = []

    for sym_op in sg.operations():
        for ta in (-1, 0, 1):
            for tb in (-1, 0, 1):
                for tc in (-1, 0, 1):
                    # Skip identity at reference unit cell [0, 0, 0]
                    if sym_op == identity_op and ta == 0 and tb == 0 and tc == 0:
                        continue

                    # Check if any atom in this symmetry mate touches the reference ASU
                    has_contact = False
                    for chain in structure[0]:
                        if has_contact:
                            break
                        for res in chain:
                            if has_contact:
                                break
                            for atom in res:
                                # Transform Cartesian -> Fractional
                                frac = cell.fractionalize(atom.pos)
                                # Apply rotational/screw symmetry operation
                                sym_frac = sym_op.apply_to_xyz([frac.x, frac.y, frac.z])
                                # Add unit cell translation directly
                                final_frac = gemmi.Fractional(
                                    sym_frac[0] + ta,
                                    sym_frac[1] + tb,
                                    sym_frac[2] + tc,
                                )
                                sym_pos = cell.orthogonalize(final_frac)

                                # Keyword-only argument: radius=cutoff_distance
                                if len(ns.find_atoms(sym_pos, radius=cutoff_distance)) > 0:
                                    contact_transforms.append((sym_op, (ta, tb, tc)))
                                    has_contact = True
                                    break

    # Build the assembled contact shell
    shell_struct = gemmi.Structure()
    shell_struct.cell = cell
    shell_struct.spacegroup_hm = structure.spacegroup_hm
    shell_model = gemmi.Model("1")

    # Add reference ASU chains
    for chain in structure[0]:
        ref_chain = chain.clone()
        ref_chain.name = f"{chain.name}_ref"
        shell_model.add_chain(ref_chain)

    # Add interacting symmetry mates with unique chain identifiers
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
    chain_counter = 0

    for sym_op, (ta, tb, tc) in contact_transforms:
        for chain in structure[0]:
            sym_chain = chain.clone()
            suffix = chain_counter // len(alphabet)
            sym_chain.name = f"{alphabet[chain_counter % len(alphabet)]}{suffix if suffix > 0 else ''}"
            chain_counter += 1

            for res in sym_chain:
                for atom in res:
                    frac = cell.fractionalize(atom.pos)
                    sym_frac = sym_op.apply_to_xyz([frac.x, frac.y, frac.z])
                    final_frac = gemmi.Fractional(
                        sym_frac[0] + ta,
                        sym_frac[1] + tb,
                        sym_frac[2] + tc,
                    )
                    atom.pos = cell.orthogonalize(final_frac)

            shell_model.add_chain(sym_chain)

    shell_struct.add_model(shell_model)
    shell_struct.make_mmcif_document().write_file(output_path)
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
                generate_first_contact_shell(raw_file, output_file, cutoff_distance=4.5)

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
    process_unprocessed_structures(target_new_count=5)