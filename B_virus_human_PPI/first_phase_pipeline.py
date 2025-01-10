import subprocess
import sys

# List of Python scripts to run (relative paths)
scripts = [
    "src/01_create_graphs/01_01_preprocess_raw_data.py",
    "src/01_create_graphs/01_02_create_interaction_matrices.py",
    "src/01_create_graphs/01_03_create_empirical_graph.py",
    "src/01_create_graphs/01_04_create_random_graphs.py",
    "src/02_cluster_virus_families/02_01_run_all_algorithms.py",
    "src/02_cluster_virus_families/02_02_regex_to_compile_results.py",
    "src/02_cluster_virus_families/02_03_get_overlapping_signficant_result.py",
    "src/03_get_baits_from_clusters/03_01_get_baits_from_cluster_algorithms.py",
    "src/03_get_baits_from_clusters/03_02_collapse_on_human_accesssion.py",
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
