# setup.py  
# Author Roelof Rietbroek (roelof@geod.uni-bonn.de), 2019
import setuptools
from setuptools import find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

print(find_packages("pyfesom"))
setuptools.setup(
    name="pyfesom",
    author="FESOM team",
    author_email="??@awi.de",
    description="Collection of tools for basic handling of FESOM ocean model output",
    long_description=long_description,
    url="https://github.com/FESOM/pyfesom",
    packages=find_packages("."),
    package_dir={"":"."},
    classifiers=["Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT 1.0",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering",
        "Development Status :: 1 - Planning"]
    
)
