import setuptools

with open("README.md", "r") as fh:
        long_description = fh.read()

with open('requirements.txt') as f:
        required = f.read().splitlines()

setuptools.setup(
            name="pyfesom",
            version="0.0.1",
            author="Nikolai Kuldanov",
            author_email="pgierz@awi.de",
            description="Python Tools for working with FESOM",
            long_description=long_description,
            long_description_content_type="text/markdown",
            install_requires=required,
            url="https://github.com/pypa/sampleproject",
            packages=setuptools.find_packages(),
            classifiers=[
                    "Programming Language :: Python :: 3",
                    "Operating System :: OS Independent",
                    ],
            )
