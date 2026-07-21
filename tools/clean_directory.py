import os
import sys
import subprocess
from glob import glob

# Add parent directory to sys.path and import config
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(sys.modules[__name__].__file__), "..")))
from oxkat import config as cfg

def safe_remove(patterns, exclude_patterns=None):
    """Remove files matching each pattern individually, reporting once per pattern.
    
    Args:
        patterns: List of glob patterns to match
        exclude_patterns: List of glob patterns to exclude from deletion
    """
    exclude_patterns = exclude_patterns or []
    
    for pattern in patterns:
        matched = glob(pattern)
        if not matched:
            print(f"No matches for pattern: {pattern}")
            continue

        # Filter out excluded files
        excluded_files = set()
        for exclude_pattern in exclude_patterns:
            excluded_files.update(glob(exclude_pattern))
        
        files_to_remove = [f for f in matched if f not in excluded_files]
        
        if not files_to_remove:
            print(f"No files to remove after exclusions for pattern: {pattern}")
            continue
        
        for path in files_to_remove:
            try:
                subprocess.run(["rm", "-rf", path], check=True)
            except subprocess.CalledProcessError as e:
                print(f"Failed to remove {path}: {e}")

        print(f"Removed {len(files_to_remove)} file(s) matching pattern: {pattern}")

# ================== COMMAND LINE ARGUMENTS ==================

if len(sys.argv) < 2:
    print("Usage: python3 clean_directory.py <source_name>")
    print("  source_name: Name of the source to match in filenames")
    print(f"\nReceived arguments: {sys.argv}")
    print(f"Number of arguments: {len(sys.argv)}")
    sys.exit(1)

try:
    source_name = sys.argv[1]
    print(f"Source name: {source_name}")
except IndexError as e:
    print(f"ERROR: Failed to get source_name from sys.argv[1]")
    print(f"sys.argv: {sys.argv}")
    print(f"sys.argv length: {len(sys.argv)}")
    print(f"Exception: {e}")
    sys.exit(1)

# Combined patterns
patterns = [
    f"{cfg.CWD}/casa*.log",
    f"{cfg.CWD}/*.txt",
    f"{cfg.CWD}/*.last",
    f"{cfg.CWD}/*.sh.sh",
    #f"{cfg.CWD}/*{source_name}*scan*.ms",
    #f"{cfg.CWD}/*1024ch*{source_name}*.ms*",  # Optional
    f"{cfg.CWD}/*{source_name}*.parmdb",
    f"{cfg.CWD}/*{source_name}*.skel",
    f"{cfg.CWD}/*.html",
    #f"{cfg.IMAGES}/{source_name}/*datablind*",
    f"{cfg.IMAGES}/{source_name}/*datamask*0*",
    f"{cfg.IMAGES}/{source_name}/*uniform-0*",
    f"{cfg.IMAGES}/{source_name}/*uniform-MFS-*residual*",
    f"{cfg.IMAGES}/{source_name}/*uniform-MFS-*model*",
    f"{cfg.IMAGES}/{source_name}/*uniform-MFS-*dirty*",
    f"{cfg.IMAGES}/{source_name}/*uniform-MFS-*psf*",
    f"{cfg.IMAGES}/{source_name}/*notaper-0*",
    f"{cfg.IMAGES}/{source_name}/*notaper-MFS-*residual*",
    f"{cfg.IMAGES}/{source_name}/*notaper-MFS-*model*",
    f"{cfg.IMAGES}/{source_name}/*notaper-MFS-*dirty*",
    f"{cfg.IMAGES}/{source_name}/*notaper-MFS-*psf*",
    f"{cfg.IMAGES}/{source_name}/*intermask*0*",
    f"{cfg.IMAGES}/{source_name}/*pcalmask-0*",
    f"{cfg.IMAGES}/{source_name}/*pcalmask-MFS-*residual*",
    f"{cfg.IMAGES}/{source_name}/*pcalmask-MFS-*model*",
    f"{cfg.IMAGES}/{source_name}/*pcalmask-MFS-*dirty*",
    f"{cfg.IMAGES}/{source_name}/*pcalmask-MFS-*psf*",
    f"{cfg.IMAGES}/{source_name}/*snapmask*",
    f"{cfg.IMAGES}/{source_name}/*_mask*",
    f"{cfg.IMAGES}/{source_name}/*diagnostic-MFS*",
    f"{cfg.IMAGES}/{source_name}/*diagnostic-0*",
    f"{cfg.IMAGES}/{source_name}/*snapblind*",
    f"{cfg.IMAGES}/{source_name}/*psf*",
    f"{cfg.IMAGES}/{source_name}/*dirty*",
    f"{cfg.IMAGES}/{source_name}/*.pb.*",
    f"{cfg.IMAGES}/{source_name}/*.wt.*",
    f"{cfg.IMAGES}/{source_name}/*kernel*",
    f"{cfg.INTERVALS}/*{source_name}*psf*",
    f"{cfg.INTERVALS}/*{source_name}*modelsub*.fits",
    f"{cfg.INTERVALS}/*{source_name}*-psf*.fits",
    # f"{cfg.RESULTS}/meerkat_gps_data",
]

# Patterns to exclude from deletion (e.g., datablind masks)
exclude_patterns = [
    f"{cfg.IMAGES}/{source_name}/*datablind*mask.fits",
]

# Run deletions
print("Cleaning directories...")
safe_remove(patterns, exclude_patterns)

