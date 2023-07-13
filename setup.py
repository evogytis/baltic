from setuptools import setup, find_packages

with open("README.md", "rb") as fh:
    long_description = fh.read().decode()

with open("requirements.txt") as fh:
    requirements = fh.read().splitlines()

setup(
    name="baltic",
    version="0.2.2",
    packages=find_packages(),
    url="https://github.com/evogytis/baltic",
    download_url="https://github.com/evogytis/baltic/archive/v0.2.2.tar.gz",
    keywords = ['phylogeny', 'visualization'],
    license="gpl-3.0",
    author="Gytis Dudas",
    author_email="gytisdudas@gmail.com",
    description="Lightweight package for analyzing, manipulating and visualizing annotated phylogenetic trees",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=requirements,
    python_requires=">=3.5",
    include_package_data=False,
    zip_safe=False,
    # classifiers=[],
)
