import os
import sys
import subprocess
from glob import glob

# Add parent directory to sys.path and import config
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(sys.modules[__name__].__file__), "..")))
from oxkat import config as cfg

def safe_remove(patterns):
    """Remove files matching each pattern individually, reporting once per pattern."""
    for pattern in patterns:
        matched = glob(pattern)
        if not matched:
            print(f"No matches for pattern: {pattern}")
            continue

        for path in matched:
            try:
                subprocess.run(["rm", "-rf", path], check=True)
            except subprocess.CalledProcessError as e:
                print(f"Failed to remove {path}: {e}")

        print(f"Removed {len(matched)} file(s) matching pattern: {pattern}")

# Combined patterns
patterns = [
    f"{cfg.CWD}/*.log",
    f"{cfg.CWD}/*.txt",
    f"{cfg.CWD}/*.last",
    f"{cfg.CWD}/*scan*.ms",
    # f"{cfg.CWD}/*1024ch_*.ms",  # Optional
    f"{cfg.CWD}/*.parmdb",
    f"{cfg.CWD}/*.skel",
    f"{cfg.CWD}/*.html",
    f"{cfg.IMAGES}/*datablind-0*",
    f"{cfg.IMAGES}/*datablind-0*",
    f"{cfg.IMAGES}/*datablind-MFS-*image.fits",
    f"{cfg.IMAGES}/*datablind-MFS-*residual.fits",
    f"{cfg.IMAGES}/*datablind-MFS-*model.fits",
    f"{cfg.IMAGES}/*datablind-MFS-*dirty.fits",
    f"{cfg.IMAGES}/*datablind-MFS-*psf.fits",
    f"{cfg.IMAGES}/*datamask-0*",
    f"{cfg.IMAGES}/*datamask-MFS*",
    f"{cfg.IMAGES}/*uniform-0*",
    f"{cfg.IMAGES}/*uniform-MFS-*residual*",
    f"{cfg.IMAGES}/*uniform-MFS-*model*",
    f"{cfg.IMAGES}/*uniform-MFS-*dirty*",
    f"{cfg.IMAGES}/*uniform-MFS-*psf*",
    f"{cfg.IMAGES}/*intermask*",
    f"{cfg.IMAGES}/*pcalmask-0*",
    f"{cfg.IMAGES}/*pcalmask-MFS-*residual*",
    f"{cfg.IMAGES}/*pcalmask-MFS-*model*",
    f"{cfg.IMAGES}/*pcalmask-MFS-*dirty*",
    f"{cfg.IMAGES}/*pcalmask-MFS-*psf*",
    f"{cfg.IMAGES}/*snapmask*",
    f"{cfg.IMAGES}/*_mask*",
    f"{cfg.IMAGES}/*diagnostic-MFS*",
    f"{cfg.IMAGES}/*diagnostic-0*",
    f"{cfg.IMAGES}/*snapblind*",
    f"{cfg.IMAGES}/*psf*",
    f"{cfg.IMAGES}/*dirty*",
    f"{cfg.IMAGES}/*.pb.*",
    f"{cfg.IMAGES}/*.wt.*",
    f"{cfg.IMAGES}/*kernel*",
    f"{cfg.INTERVALS}/*psf*",
    f"{cfg.INTERVALS}/*modelsub*.fits",
    f"{cfg.INTERVALS}/*-psf*.fits",
    f"{cfg.RESULTS}/meerkat_gps_data",
]

# Run deletions
print("Cleaning directories...")
safe_remove(patterns)

