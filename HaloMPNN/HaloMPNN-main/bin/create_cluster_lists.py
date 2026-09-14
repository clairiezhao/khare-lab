#!/usr/bin/env python3
import torch
import pandas as pd
from tqdm import tqdm
from sklearn.model_selection import train_test_split
import os, sys
import logging

logging.basicConfig(format="[{asctime}|{levelname:7s}] {message}", style='{', 
        filemode='w', datefmt='%Y-%m-%d %H:%M:%S', level=logging.INFO)

def get_files_from_path(path:str, filesuffix:str):
    '''Return all files (with path) in a given path containing the filesuffix.'''
    return [os.path.join(path, file) for file in os.listdir(path) if file.endswith(filesuffix)]


def _get_pt_dict(input_path:str, debug:bool=False):
    logging.info('Find pt files for models.')
    pt_files = get_files_from_path(input_path, filesuffix='.pt')
    pt_files = [file for file in pt_files]

    if debug:
        pt_files = pt_files[:100]

    pt_dict = {}
    logging.info(f'Read pt files (n = {len(pt_files)})')
    for pt_file in pt_files:
        basename = pt_file.split('/')[-1].replace('.pt', '')
        torch_file = torch.load(pt_file)

        seqs = [torch_file['seq'][i][0] for i in range(len(torch_file['seq']))]
        seq_hash = str(torch_file['id']) + '_id'
        
        # chain_list = [f'{basename}_{letter}' for letter in torch_file['chains']]
        chain_list = [f'{basename}' for letter in torch_file['chains']]
        fasta_list = [f'{basename}' for letter in torch_file['chains']]
        
        # other properties
        date = torch_file['date']
        resolution = torch_file['resolution']
        
        for file, chain, seq in zip(fasta_list, chain_list, seqs):
            pt_dict.update({file : {'model' : chain, 'hash': seq_hash, 'seq': seq, 
                            'date': date, 'resolution': resolution}})

    return pt_dict


def _create_clusters(cluster_file:str, pt_dict:dict, outdir):
    cluster_df = pd.read_csv(filepath_or_buffer=cluster_file, sep='\t', names=['represent', 'member'])
    # display(cluster_df)

    cluster_names = cluster_df['represent'].unique()
    n_clusters = len(cluster_names)
    logging.info(f'N clusters: {n_clusters}')

    # iterate over all representatives and get the cluster members
    counter_cluster = 1
    member_dict = {'CHAINID': [], 'DEPOSITION': [], 'RESOLUTION': [], 'HASH': [], 'CLUSTER': [], 'SEQUENCE': []}
    logging.info('Create sequence list.')
    missing = 0
    for id_cluster in cluster_names:
        # store current representative and get all cluster members
        represent = id_cluster
        members = cluster_df.loc[cluster_df['represent'] == id_cluster, : ]['member']
        # print(pt_dict)
        # for every member create new entry in list
        for member in members:
            try:
                member_dict['CHAINID'].append(pt_dict[member]['model'])
                member_dict['DEPOSITION'].append(pt_dict[member]['date'])
                member_dict['RESOLUTION'].append(pt_dict[member]['resolution'])
                member_dict['HASH'].append(pt_dict[member]['hash'])
                member_dict['CLUSTER'].append(counter_cluster)
                member_dict['SEQUENCE'].append(pt_dict[member]['seq'])
            except KeyError:
                logging.error(f'Key {member} not in pt_dict... Why??')
                missing += 1

        # add counter + 1 for new cluster
        counter_cluster += 1
        
    if missing > 0:
        logging.warning(f'Missing key from pt_dict: {missing}')
    

    df_member = pd.DataFrame(data=member_dict)
    df_member.to_csv(path_or_buf=os.path.join(outdir, 'list_cluster.csv'), index=False)
    return df_member



if __name__ == '__main__':
    inputpt = sys.argv[1]
    clusterfile = sys.argv[2]
    outdir = sys.argv[3]
    
    pt_dict = _get_pt_dict(input_path=inputpt)
    df_member = _create_clusters(cluster_file=clusterfile, pt_dict=pt_dict, outdir=outdir)
    
    X_train, X_test, = train_test_split(df_member['CLUSTER'].unique(), test_size=0.4)
    X_test, X_val, = train_test_split(X_test, test_size=0.5)

    logging.info(f'Train: {X_train.shape[0]}, Valid: {X_val.shape[0]}, Test: {X_test.shape[0]}')
    pd.DataFrame(X_train).to_csv(os.path.join(outdir, 'train_cluster.txt'), index=False, header=False)
    pd.DataFrame(X_val).to_csv(os.path.join(outdir, 'valid_cluster.txt'), index=False, header=False)
    pd.DataFrame(X_test).to_csv(os.path.join(outdir, 'test_cluster.txt'), index=False, header=False)
    
    logging.info('Done.')
        
