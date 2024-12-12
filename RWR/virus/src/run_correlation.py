import subprocess

# Path to the cophenetic_correlation.py script
script_path = "./cophenetic_correlation.py"

# Run the script 100 times with incremental arguments
for i in range(100):
    try:
        print(f"Running iteration {i}...")
        subprocess.run(["python", script_path, str(i)], check=True)
        print(f"Iteration {i} completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error in iteration {i}: {e}")
