## Processing

Run the following pipeline to bias correct CMIP6 data variables using ERA5, then calculate CFFDRS fire weather indices from the results. 

This README assumes you have downloaded copies of the CMIP6 and ERA5 source data. Contact SNAP for information about obtaining source data for this processing pipeline.

### Setup

- Create a conda environment using the `environment.yml` (TBD)
```bash
conda env create -f environment.yml
```

- Set environment variables
```bash
export DATA_DIR=/path/to/source/data
export TMP_DIR=/path/for/temp/files
export OUT_DIR=/path/for/output/files
```

### Subsetting & standardization
- `00_unzip.py`
  - Unzip CMIP6 and ERA5 data.

- `01_process_era5.py`
  - Get daily summaries of desired variables, modify geographic extent, and standardize calendar. 

- `02_process_cmip6.py`
  - Modify geographic extent, standardize calendar, convert units (e.g., K -> degC).


### Bias correction

### CFFDRS calculation

### References