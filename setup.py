from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

with open('requirements.txt') as f:
        required = f.read().splitlines()

print(find_packages("pyfesom"))
setup(
    name="pyfesom",
    version="0.0.1",
    author="FESOM team",
    author_email="koldunovn@gmail.com",
    description="Collection of tools for basic handling of FESOM ocean model output",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=required,
    url="https://github.com/FESOM/pyfesom",
    packages=find_packages(include=['pyfesom']),
    package_dir={"":"."},
    classifiers=["Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT 1.0",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering",
        "Development Status :: 1 - Planning"]
    
)

