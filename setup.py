import setuptools

install_requires = ["h5py", "pandas", "tqdm", "numpy", "networkx", "time", "gurobipy", "graph-tool"]

setuptools.setup(
    name="PanTax",
    version="1.0.0",
    author="Weihai Zhang",
    url="https://github.com/LuoGroup2023/PanTax.git",
    packages=setuptools.find_packages(),   
    entry_points={
        "console_scripts": [
            "build_graph.py = pantax.build_graph:main",
            "extract_complete_genome.py = pantax.extract_complete_genome:main",
            "gaf_filter.py = pantax.gaf_filter:main",
            "genomes_cluster.py = pantax.genomes_cluster:main",
            "get_genomes_info.py = pantax.get_genomes_info:main",
            "otu_range.py = pantax.otu_range:main",
            "prepare_building.py = pantax.prepare_building:main",
            "read_classification.py = pantax.read_classification:main",
            "species_abundance_cal.py = pantax.species_abundance_cal:main",
            "strain_abundance_cal.py = pantax.strain_abundance_cal:main"
        ]
    },
    description="Strain-level taxonomic classification of metagenomic data using pangenome graphs",
    install_requires=install_requires
)