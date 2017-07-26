# emr_visualization
## Author: Edward Huang

## Preprocessing work for the LONI Parkinson's disease dataset.

1.  Call on the STITCH API to fetch batches of files containing chemical-drug
    interactions, as well as outputting a dictionary mapping WHO drug IDs to
    STITCH identifiers. The LONI Parkinson's dataset uses WHO drug IDs.
    
    ```bash
    python call_stitch_api.py
    ```

2.  Find interesting SNPs that are relevant to Parkinson's disease. Creates two
    files: ./data/ppmi/snp_healthy_freq.tsv, which shows the top SNPs ranked
    by their rareness in healthy populations. ./datappmi/snp_fisher_test.tsv,
    which shows the top SNPs as ranked by enrichments in gene mutations between
    PD patients and healthy patients.

    ```bash
    python snp_fisher.py [-h] -t {wgs,wes} -c {hard,ppmi} [-b {ignore,use}]
    ```
    Necessary runs:
    python snp_fisher.py -t wes -c hard -b ignore
    python snp_fisher.py -t wes -c ppmi

## Running ProSNet.

1.  Generate the low-dimensional vectors for each node in the network.

    ```bash
    python run_prosnet.py num_dim
    ```

2.  Build the patient feature matrices for baseline as well as normal run.

    ```bash
    python build_patient_feature_matrix.py [-h] -n {max,l1,l2} [-d NUM_DIM]
                                       [-s SIM_THRESH] [-w {before,after,both}]
                                       [-l {updrs,status,tau}] [-e {biospecimen}]
    ```
    Only -n and -l are required.

## Reducing dimensionality and visualizing the EMRs.

1.  
    ```bash
    python visualize_emrs.py [-h] [-n NORM_TYPE] [-d NUM_DIM] [-s SIM_THRESH]
                                [-w WHERE_NORM] [-p N_PCA_COMP] [-t TSNE_INIT]
                                [-r L_RATE] [-k KNN]
    ```
    -n, -p, -t required.

## Extra scripts for looking at the data.

1.  Writes out all symptom features to file, one on each new line.

    ```bash
    python write_symptoms_to_file.py
    ```