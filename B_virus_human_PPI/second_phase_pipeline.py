import subprocess
import sys

###################################

# READ ME !!!

# The script, <src/03_get_baits_from_clusters/03_03_extract_baits_from_significant_hits.py>
# Needs to be updated within the script itself to match the signficiant clusters found from phase 1 pipeline
# Internet is needed to run, <"src/04_stringDB/04_01_stringDB.py>

####################################

# List of Python scripts to run
scripts = [
    "src/03_get_baits_from_clusters/03_03_extract_baits_from_significant_hits.py",
    "src/04_stringDB/04_01_stringDB.py"
]

print("Starting the pipeline...")

for script in scripts:
    print(f"Running {script}...")
    try:
        result = subprocess.run(
            [sys.executable, script],
            check=True
        )
    except subprocess.CalledProcessError as e:
        print(f"Pipeline failed at {script} with return code {e.returncode}")
        sys.exit(1)

print("Pipeline completed successfully!")