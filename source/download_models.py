from pathlib import Path;import shutil,sys,os;
import sys;from pathlib import Path;sys.path.append(str(Path(__file__).parent))
from paths import dir_libraries
import os
import subprocess
import zipfile

def main():
    dirs=dir_libraries(mkdir=False)
    # https://drive.google.com/file/d/1vlk_jX2TfOQG8z3EAKCbqUwyaxJ9HHBG/view?usp=sharing
    # https://drive.google.com/file/d/10fHIyJ1aAhureFMS6p6L3L6Iuqfqbdtk/view?usp=sharing
    # https://drive.google.com/file/d/13cyrF25el7fpXdJpfgex3-uxi6DRaWoW/view?usp=sharing
    # https://drive.google.com/file/d/1f9EBJ3CAeisIOU7gN7W1MA4G2W9rUXoI/view?usp=sharing
    # Map local paths -> Google Drive file IDs
    FILES = {dirs.Root/'_DataArchive.zip':'1vlk_jX2TfOQG8z3EAKCbqUwyaxJ9HHBG',
    str(dirs.Data/'SNR_Models'/'BulkHold.pkl.zip'): "10fHIyJ1aAhureFMS6p6L3L6Iuqfqbdtk",
    str(dirs.Data/'SNR_Models'/'SNR_acausul.filter_V04_5s_bandwidth_100_bands.pkl.zip'): "13cyrF25el7fpXdJpfgex3-uxi6DRaWoW",
    str(dirs.Analysis/'BulkLoad.SR.Coherences_092625.pkl.zip'): "1f9EBJ3CAeisIOU7gN7W1MA4G2W9rUXoI",
    }
    # link=lambda FILE_ID:f'https://drive.google.com/file/d/{FILE_ID}/view?usp=sharing'
    # ensure gdown is installed
    try:
        subprocess.run(["gdown", "--version"], check=True, capture_output=True)
    except Exception:
        raise SystemExit(
            "gdown is required. Install with:\n\n    pip install gdown\n"
        )
    for local_path, file_id in FILES.items():
        folder = os.path.dirname(local_path)
        if folder:os.makedirs(folder, exist_ok=True)
        print(f"Downloading {local_path} ...")
        subprocess.run(
            ["gdown", "--id", file_id, "-O", local_path],
            check=True,
        )

        if Path(local_path).name=='_DataArchive.zip':
            with zipfile.ZipFile(local_path, "r") as zf:zf.extractall(Path(local_path).parent)

    print("All models downloaded.")
if __name__ == "__main__":
    main()
