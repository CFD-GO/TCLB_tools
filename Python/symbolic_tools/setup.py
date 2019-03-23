import setuptools

with open("../../README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="test-tclb_tools-pkg",
    version="0.0.1",
    author="ggruszczynski",
    author_email="ggruszczynski@gmail.com",
    description="Various tools for the TCLB project.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/CFD-GO/TCLB_tools",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: GPL-3.0",
        "Operating System :: OS Independent",
    ],
)
