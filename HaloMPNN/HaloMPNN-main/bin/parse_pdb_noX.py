from Bio.PDB import PDBParser
import numpy as np
import torch
import sys
import os


RES_NAMES = [
    'ALA','ARG','ASN','ASP','CYS',
    'GLN','GLU','GLY','HIS','ILE',
    'LEU','LYS','MET','PHE','PRO',
    'SER','THR','TRP','TYR','VAL'
]

RES_NAMES_1 = 'ARNDCQEGHILKMFPSTWYV'

to1letter = {aaa:a for a,aaa in zip(RES_NAMES_1,RES_NAMES)}
to3letter = {a:aaa for a,aaa in zip(RES_NAMES_1,RES_NAMES)}

ATOM_NAMES = [
    ("N", "CA", "C", "O", "CB"), # ala
    ("N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"), # arg
    ("N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"), # asn
    ("N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"), # asp
    ("N", "CA", "C", "O", "CB", "SG"), # cys
    ("N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"), # gln
    ("N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"), # glu
    ("N", "CA", "C", "O"), # gly
    ("N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"), # his
    ("N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"), # ile
    ("N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"), # leu
    ("N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"), # lys
    ("N", "CA", "C", "O", "CB", "CG", "SD", "CE"), # met
    ("N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"), # phe
    ("N", "CA", "C", "O", "CB", "CG", "CD"), # pro
    ("N", "CA", "C", "O", "CB", "OG"), # ser
    ("N", "CA", "C", "O", "CB", "OG1", "CG2"), # thr
    ("N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE2", "CE3", "NE1", "CZ2", "CZ3", "CH2"), # trp
    ("N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"), # tyr
    ("N", "CA", "C", "O", "CB", "CG1", "CG2") # val
]
        
idx2ra = {(RES_NAMES_1[i],j):(RES_NAMES[i],a) for i in range(20) for j,a in enumerate(ATOM_NAMES[i])}

aa2idx = {(r,a):i for r,atoms in zip(RES_NAMES,ATOM_NAMES) 
          for i,a in enumerate(atoms)}
aa2idx.update({(r,'OXT'):3 for r in RES_NAMES})

def parse_pdb(file:str, id:str):
    print(file)
    parser = PDBParser(QUIET=True)
    model = parser.get_structure(file=file, id='random')
    complete_structure = {}
    meta_data = {'seq': []}
    chains = []

    for chain in model.get_chains():
        chains.append(chain.get_id())
        
        full_seq, gaps_seq = [], []
        full_atoms_xyz = []
        bfac = []
        occ = []
        mask = []
        for residue in chain.get_residues():
            # adding the sequences for every chain
            # skip water HOH
            aa_3_name = residue.get_resname()
            if aa_3_name == 'HOH':
                continue
            
            try:
                aa_1_name = to1letter[aa_3_name]
            except KeyError:
                continue
            full_seq.append(aa_1_name)
            gaps_seq.append(aa_1_name)
            
            # range(14) because 14 different atom types, maximum is in
            # Tryptophan (see list)
            atoms_xyz = [[np.nan for _ in range(3)] for _ in range(14)]
            atoms_bfac = [np.nan for _ in range(14)]
            atoms_occ = [0 for _ in range(14)]
            atoms_mask = [False for _ in range(14)]
            for atom in residue.get_atoms():
                # index for final list, same order for every residue
                # from the parse_cif_noX.py script
                try:
                    atom_idx = aa2idx[aa_3_name, atom.get_id()]
                except KeyError:
                    continue
                x, y, z = atom.get_coord()
                
                atoms_xyz[atom_idx] = [np.round(n, 3) for n in [x,y,z]]
                atoms_bfac[atom_idx] = atom.bfactor
                atoms_occ[atom_idx] = atom.occupancy
                atoms_mask[atom_idx] = True
            
            full_atoms_xyz.append(atoms_xyz)
            bfac.append(atoms_bfac)
            occ.append(atoms_occ)
            mask.append(atoms_mask)
        
        chain_name = chain.get_id()
        tmp_d = {
                'seq': ''.join(gaps_seq),
                'xyz': np.array(full_atoms_xyz),
                'occ': np.array(occ),
                'bfac': np.array(bfac),
                'mask': np.array(mask)
        }
        complete_structure.update({chain_name: tmp_d})
        meta_data['seq'].append([''.join(full_seq), ''.join(gaps_seq)])
        
    # addint pseudo values for the other fields
    meta_data.update({'method' : 'biopython'})
    meta_data.update({'date': '2023-01-01'})
    meta_data.update({'resolution': 0.0})
    meta_data.update({'chains': chains})
    meta_data.update({'id': str(id)})
    meta_data.update({'asmb_chains': [','.join(chains)]})
    meta_data.update({'asmb_details': ['own_parsing']})
    meta_data.update({'asmb_method': ['none']})
    meta_data.update({'asmb_ids': ['1']})
    xform = [
        [[1. ,0. ,0. ,0.],
        [0. ,1. ,0. ,0.],
        [0. ,0. ,1. ,0.],
        [0. ,0. ,0. ,1.]]
    ]
    meta_data.update({'asmb_xform0': np.array(xform)})
    
    return complete_structure, meta_data        

IN = sys.argv[1]
OUT = sys.argv[2]
ID = str(sys.argv[3])

save_path = os.path.join(OUT, f'{ID}.pt')

if not os.path.exists(save_path):

    chains,metadata = parse_pdb(IN, ID)

    meta_pt = {}
    for k,v in metadata.items():
        if "asmb_xform" in k or k=="tm":
            v = torch.Tensor(v)
        meta_pt.update({k:v})

    torch.save(meta_pt, save_path)
