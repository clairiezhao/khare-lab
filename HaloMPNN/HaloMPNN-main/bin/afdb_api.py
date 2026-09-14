import requests
import sys
import os
import json

def get_model(sequence):
    url=f"https://alphafold.ebi.ac.uk/api/prediction/{sequence}"
    headers={
    "accept":"application/json",
    "key":"AIzaSyCeurAJz7ZGjPQUtEaerUkBZ3TaBkXrY94",#ThereisnoofficialAPIaccess.ThisisaworkingfakeAPIkey.
    }
    response=requests.get(url,headers=headers)
    if response.status_code==200:
        data=response.json()
        if isinstance(data,list) and len(data)>0:
            latest_model=data[-1]
            return latest_model
    return None

def extract_metadata(uniprot_id,model):
    if model is None:
        return[uniprot_id,"NA","NA","NA","NA"]
    afdb_id=model.get("modelEntityId")
    pdb_url=model.get("pdbUrl")
    mean_plddt=str(model.get("globalMetricValue"))
    seqstart=model.get("sequenceStart")
    seqend=model.get("sequenceEnd")
    seqlength=str(int(seqend)-int(seqstart)+1)
    return[uniprot_id,afdb_id,pdb_url,mean_plddt,seqlength]

def main():
    uniprot_ids_file=sys.argv[1]
    output_file=sys.argv[2]
    with open(output_file,"w") as f:
        f.write("uniprot_id,afdb_id,pdb_url,mean_plddt,length"+"\n")
    with open(uniprot_ids_file,"r") as f:
        for uniprot_id in f:
            id = uniprot_id.rstrip("\n")
            model=get_model(id)
            metadata=extract_metadata(id,model)
            with open(output_file,"a") as f2:
                f2.write(",".join(metadata)+"\n")


if __name__ == "__main__":
    main()
