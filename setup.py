from setuptools import setup, find_packages

setup(
    name="hcase",
    version="0.1.0",
    description="The hcase library implements the novel Hilbert-Curve Assisted Structure Embedding (HCASE) method that is able to generate 2D maps where molecules are arranged in an intuitive manner following medicinal chemists' logic. Once a chemical space is established, compounds can embedded in it without affecting the positions of other, already embedded compounds, making this a unique method compared to existing space embedding and dimensional reduction techniques.",
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
