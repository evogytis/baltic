from setuptools import setup, find_packages

with open("README.md", "rb") as fh:
    long_description = fh.read().decode()

with open("requirements.txt") as fh:
    requirements = fh.read().splitlines()

setup(
    name="baltic",
    version="1.0.0",
    packages=find_packages(),
    url="https://github.com/evogytis/baltic",
    project_urls={
        "Documentation": "https://baltic.readthedocs.io/",
        "Gallery": "https://phylo-baltic.github.io/baltic-gallery/",
        "Source": "https://github.com/evogytis/baltic",
        "Releases": "https://github.com/evogytis/baltic/releases",
    },
    keywords = ['phylogeny', 'visualization'],
    license="gpl-3.0",
    author="Gytis Dudas",
    author_email="gytisdudas@gmail.com",
    description="Lightweight package for analyzing, manipulating and visualizing annotated phylogenetic trees",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=requirements,
    extras_require={
        "seq": ["biopython", "bcbio-gff"],
        "maps": ["cartopy"],
        "treemap": ["squarify"],
        "all": ["biopython", "bcbio-gff", "cartopy", "squarify"],
    },
    python_requires=">=3.10",
    include_package_data=False,
    zip_safe=False,
    # classifiers=[],
)
