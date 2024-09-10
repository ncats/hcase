from setuptools import setup, find_packages

setup(
    name="hcase",
    version="0.1.0",
    description="TO BE PROVIDED BY GERGELY.",
    packages=find_packages(include=["hcase", "hcase.*"]),  # Include only the hcase folder
    package_dir={"": "."},  # The hcase folder is in the root directory
    install_requires=[
        "pandas",
        "numpy",
        "hilbertcurve",
        "matplotlib",
        "networkx",
        "notebook",
        "scipy",
        "seaborn",
        "tqdm",
        "scikit-learn",
        "jupyter-contrib-nbextensions"
    ],
    extras_require={
        'rdkit': [],
    },
    long_description="""
    RDKit is not installable via pip. If you need RDKit functionality, 
    please install it separately via conda:
    
    `conda install -c conda-forge rdkit`
    """
)
