from setuptools import setup, find_packages

with open("README.md", "rb") as fh:
    long_description = fh.read().decode()

with open("requirements.txt") as fh:
    requirements = fh.read().splitlines()

setup(
    name="baltic",
    version="0.1.0",
    packages=find_packages(),
    url="https://github.com/evogytis/baltic",
    # license="",
    author="Gytis Dudas",
    author_email="gytisdudas@gmail.com",
    description="Lightweight package for analyzing, manipulating and visualizing annotated phylogenetic trees",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=requirements,
    # python_requires=">=3.6.*",
    include_package_data=False,
    zip_safe=False,
    # classifiers=[],
)
