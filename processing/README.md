## Processing

Run the following pipeline to bias correct CMIP6 data variables using ERA5, then calculate CFFDRS fire weather indices from the results. 

This README assumes you have downloaded copies of the CMIP6 and ERA5 source data (more explanation TBD).

### Setup

- Create a conda environment using the `environment.yml` (TBD)
```
conda env create -f environment.yml
```

- Set environment variables
```
export DATA_DIR=/path/to/source/data
export TMP_DIR=/path/for/temp/files
export OUT_DIR=/path/for/output/files
```

### Prep data

### Bias Correction

### Calculate CFFDRS indices