# emr_visualization
## Author: Edward Huang

## Preprocessing scripts for reading data.

### Preprocessing the MIMIC III data.

```bash
python process_mimic.py
```

### Preprocessing the LONI Parkinson's disease data.

```bash
python process_loni_parkinsons.py
```

## Obtaining external chemical-chemical network.

1.  Call on the STITCH API to fetch batches of files containing chemical-drug
    interactions, as well as outputting a dictionary mapping WHO drug IDs to
    STITCH identifiers. The LONI Parkinson's dataset uses WHO drug IDs.
    
    ```bash
    python call_stitch_api.py
    ```

2.  Download PPI network from HumanNet.


## Running ProSNet.

1.  Generate the low-dimensional vectors for each node in the network.

    ```bash
    python run_prosnet.py num_dim
    ```

2.  Build the patient feature matrices for baseline as well as normal run.

    ```bash
    python build_patient_feature_matrix.py [-h] [-n NORM_TYPE] [-d NUM_DIM]
                                                [-s SIM_THRESH] [-w WHERE_NORM]
    ```
    Only -n is required.

## Reducing dimensionality and visualizing the EMRs.

1.  
    ```bash
    python visualize_emrs.py [-h] [-n NORM_TYPE] [-d NUM_DIM] [-s SIM_THRESH]
                                [-w WHERE_NORM] [-p N_PCA_COMP] [-t TSNE_INIT]
    ```
    -n, -p, -t required.