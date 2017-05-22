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