import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyhlamsa",
    version="0.5.2",
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
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "License :: OSI Approved :: MIT License",
    ],
    python_requires=">=3.8",
)
