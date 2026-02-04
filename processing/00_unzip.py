from pathlib import Path
from config import ZIP_DIR, DATA_DIR, source_zips
import zipfile

def unzip_source_data():
    """Unzip source data files into the temporary directory."""

    if not ZIP_DIR or not DATA_DIR:
        raise EnvironmentError("ZIP_DIR and DATA_DIR environment variables must be set.")
    
    zip_path = Path(ZIP_DIR)
    data_path = Path(DATA_DIR)
    data_path.mkdir(parents=True, exist_ok=True)

    for zip_filename in source_zips:
        try:
            zip_file_path = zip_path / zip_filename
            print(f"Unzipping {zip_file_path} to {data_path}...")
            with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
                zip_ref.extractall(data_path)
        except Exception as e:
            print(f"Error unzipping {zip_filename}: {e}")
        
    print("Unzipping completed.")

if __name__ == "__main__":
    unzip_source_data()