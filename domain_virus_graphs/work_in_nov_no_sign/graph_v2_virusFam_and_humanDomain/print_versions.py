import importlib

# List of packages to check
packages = [
    "numpy",
    "igraph",
    "scipy",
    "nimfa",
    "networkx",
    "seaborn",
    "matplotlib",
    "pandas"
]

print("Installed Package Versions:\n")
for package in packages:
    try:
        # Dynamically import the package
        pkg = importlib.import_module(package)
        # Print the package name and version
        print(f"{package}: {pkg.__version__}")
    except ImportError:
        print(f"{package}: Not Installed")
    except AttributeError:
        print(f"{package}: Version information not available")
