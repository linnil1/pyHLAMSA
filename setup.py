import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyHLAMSA",
    version="0.2.0",
    author="linnil1",
    author_email="linnil1.886@gmail.com",
    description="A small and simple utility to read IMGT HLA data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/linnil1/pyHLAMSA",
    packages=setuptools.find_packages(),
    install_requires=[
        "biopython",
        "pysam"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
)
