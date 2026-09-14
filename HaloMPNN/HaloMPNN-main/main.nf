#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ─────────────────────────────────────────────
//  Validate inputs
// ─────────────────────────────────────────────
if (!params.input_fasta) {
    error "Please provide an input FASTA file with --input_fasta"
}

// ─────────────────────────────────────────────
//  Processes
// ─────────────────────────────────────────────

process REMOVE_DUPLICATES {
    tag "seqkit rmdup"
    conda params.mmseqs_conda
    publishDir "${params.outdir}/dedup", mode: 'copy'

    input:
    path fasta

    output:
    path "input_nodups.fasta",  emit: nodups_fasta
    path "dupseqs.fasta",       emit: dup_seqs
    path "numdup.txt",          emit: num_dups

    script:
    """
    seqkit rmdup ${fasta} \\
        --id-regexp "\\|([A-Z0-9]+)\\|" \\
        -d dupseqs.fasta \\
        -D numdup.txt \\
        > input_nodups.fasta
    """
}

process EXTRACT_UNIPROT_IDS {
    tag "seqkit seq -> IDs"
    conda params.mmseqs_conda
    publishDir "${params.outdir}/dedup", mode: 'copy'

    input:
    path nodups_fasta

    output:
    path "uniprot_ids.txt", emit: uniprot_ids

    script:
    """
    seqkit seq ${nodups_fasta} \\
        -n -i \\
        --id-regexp "\\|([A-Z0-9]+)\\|" \\
        > uniprot_ids.txt
    """
}

process QUERY_AFDB_API {
    tag "AlphaFoldDB API"
    conda params.mmseqs_conda
    publishDir "${params.outdir}/metadata", mode: 'copy'

    input:
    path uniprot_ids

    output:
    path "metadata.csv", emit: metadata

    script:
    """
    python ${params.scripts_dir}/afdb_api.py \\
        ${uniprot_ids} \\
        metadata.csv
    """
}

process FILTER_LENGTH_PLDDT {
    tag "filter length=${params.min_length} pLDDT=${params.min_plddt}"
    conda params.mmseqs_conda
    publishDir "${params.outdir}/metadata", mode: 'copy'

    input:
    path metadata

    output:
    path "metadata_filtered.csv", emit: filtered_metadata
    path "filter_stats.txt",      emit: filter_stats

    script:
    """
    python ${params.scripts_dir}/filter_length_plddt.py \\
        ${params.min_length} \\
        ${params.min_plddt} \\
        --input  ${metadata} \\
        --output metadata_filtered.csv \\
        --stats  filter_stats.txt
    """
}

process MMSEQS_CLUSTER {
    tag "MMseqs2 cluster seqid=${params.seqid}"
    conda params.mmseqs_conda
    publishDir "${params.outdir}/cluster", mode: 'copy'

    // CPUs come from SLURM if running under Slurm; otherwise default to 4
    cpus { System.getenv("SLURM_CPUS_PER_TASK") ? System.getenv("SLURM_CPUS_PER_TASK") as int : 4 }

    input:
    path nodups_fasta

    output:
    path "cluster_assignments.tsv", emit: cluster_tsv
    path "repr_seqs.fasta",         emit: repr_fasta
    path "repr_ids.txt",            emit: repr_ids

    script:
    def cpus = System.getenv("SLURM_CPUS_PER_TASK") ?: task.cpus
    """
    mkdir -p tmp

    mmseqs createdb ${nodups_fasta} inputdb

    mmseqs cluster inputdb outputdb tmp \\
        --cov-mode 0 \\
        --min-seq-id ${params.seqid} \\
        --split-memory-limit ${params.split_memory_limit} \\
        --threads ${cpus}

    mmseqs createtsv inputdb inputdb outputdb cluster_assignments.tsv

    mmseqs createsubdb outputdb inputdb reprdb
    
    mmseqs convert2fasta reprdb repr_seqs.fasta

    seqkit seq repr_seqs.fasta \\
        -n -i \\
        --id-regexp "\\|([A-Z0-9]+)\\|" \\
        > repr_ids.txt
    """
}

process FILTER_METADATA_REPR {
    tag "filter metadata -> representatives"
    conda params.mmseqs_conda
    publishDir "${params.outdir}/metadata", mode: 'copy'

    input:
    path repr_ids
    path metadata

    output:
    path "repr_metadata.csv", emit: repr_metadata

    script:
    """
    python ${params.scripts_dir}/filter_metadata_repr.py \\
        ${repr_ids} \\
        ${metadata} \\
        repr_metadata.csv
    """
}

process DOWNLOAD_PDB_STRUCTURES {
    tag "AFDB download"
    conda params.mmseqs_conda
    publishDir "${params.datadir}", mode: 'copy'

    input:
    path repr_metadata

    output:
    path "pdb_files/", emit: pdb_dir

    script:
    """
    mkdir -p pdb_files
    python ${params.scripts_dir}/afdb_download.py \\
        ${repr_metadata} \\
        pdb_files/
    """
}


process PARSE_PDB_TO_PT {
    tag "parse_pdb_noX -> ${pdb_file.baseName}"
    conda params.aifold_conda
    publishDir "${params.datadir}", mode: 'copy'

    input:
    path pdb_file

    output:
    path "pt_files/${pdb_file.baseName}.pt", emit: pt_file   // one file, not the whole dir

    script:
    def id = pdb_file.baseName
    """
    mkdir -p pt_files
    python ${params.scripts_dir}/parse_pdb_noX.py \\
        ${pdb_file} \\
        pt_files/ \\
        ${id}
    """
}

process CREATE_CLUSTER_LIST {
    tag "create cluster TSV"
    conda params.mmseqs_conda
    publishDir "${params.outdir}/cluster", mode: 'copy'

    input:
    path pt_files

    output:
    path "cluster_list.tsv", emit: cluster_list

    script:
    """
    find . -name "*.pt" | xargs -I{} basename {} .pt | awk '{print \$1"\t"\$1}' > cluster_list.tsv
    """
}

process TRAIN_TEST_SPLIT {
    tag "create_cluster_lists.py"
    conda params.aifold_conda
    publishDir "${params.outdir}/cluster", mode: 'copy'

    input:
    path cluster_list

    output:
    path "data_files/", emit: data_files

    script:
    """
    mkdir -p data_files
    python ${params.scripts_dir}/create_cluster_lists.py \\
        "${params.datadir}/pt_files" \\
        ${cluster_list} \\
        'data_files'
    """
}

process MPNN_RETRAIN {
    tag "HyperMPNN retrain"
    conda params.aifold_conda

    publishDir "${params.outdir}/mpnn_training", mode: 'copy'

    input:
    path data_files

    output:
    path "train_output/", emit: train_output

    script:
    """
    mkdir -p train_output

    python3 ${params.scripts_dir}/training.py \\
        --path_for_training_data "${params.datadir}/pt_files" \\
        --path_for_data_files   ${data_files} \\
        --path_for_outputs      train_output \\
        --num_neighbors         ${params.num_neighbors} \\
        --backbone_noise        ${params.backbone_noise} \\
        --num_epochs            ${params.num_epochs} \\
        --batch_size            ${params.batch_size} \\
        --max_protein_length    ${params.max_protein_length} \\
        --num_examples_per_epoch ${params.num_examples_per_epoch} \\
        --rerun
    """
}

// ─────────────────────────────────────────────
//  Workflow
// ─────────────────────────────────────────────
workflow {

    // 1. Input channel
    input_fasta_ch = Channel.fromPath(params.input_fasta, checkIfExists: true)

    // 2. De-duplicate
    REMOVE_DUPLICATES(input_fasta_ch)

    // 3. Extract UniProt IDs
    EXTRACT_UNIPROT_IDS(REMOVE_DUPLICATES.out.nodups_fasta)

    // 4. Query AlphaFold DB API
    QUERY_AFDB_API(EXTRACT_UNIPROT_IDS.out.uniprot_ids)

    // 5. Filter by length & pLDDT
    FILTER_LENGTH_PLDDT(QUERY_AFDB_API.out.metadata)

    // 6. MMseqs2 clustering on de-duplicated sequences
    MMSEQS_CLUSTER(REMOVE_DUPLICATES.out.nodups_fasta)

    // 7. Filter metadata to representative sequences
    FILTER_METADATA_REPR(
        MMSEQS_CLUSTER.out.repr_ids,
        FILTER_LENGTH_PLDDT.out.filtered_metadata
    )

    // 8. Download PDB structures for representatives
    DOWNLOAD_PDB_STRUCTURES(FILTER_METADATA_REPR.out.repr_metadata)

    // 9. Parse each PDB -> PyTorch tensors (one process per file)
    pdb_files_ch = DOWNLOAD_PDB_STRUCTURES.out.pdb_dir
        .flatMap { dir -> file(dir).listFiles().findAll { it.name.endsWith('.pdb') } }

    PARSE_PDB_TO_PT(pdb_files_ch)

    pt_files_ch = PARSE_PDB_TO_PT.out.pt_file.collect()

    // 10. Create self-cluster TSV
    // CREATE_CLUSTER_LIST(PARSE_PDB_TO_PT.out.pt_dir)
    CREATE_CLUSTER_LIST(pt_files_ch)

    // 11. Train/valid/test split (aifold env)
    TRAIN_TEST_SPLIT(
        CREATE_CLUSTER_LIST.out.cluster_list
    )

    // 12. HyperMPNN retraining
    MPNN_RETRAIN(
        TRAIN_TEST_SPLIT.out.data_files
    )
}
