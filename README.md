# PAU4Chem

# Requirements

This code was written using Python 3.x, Anaconda 3, and operating system Ubuntu 18.04. The following Python libraries are required for running the code:

1. [bs4](https://anaconda.org/conda-forge/bs4)
2. [requests](https://anaconda.org/anaconda/requests)
3. [pandas](https://anaconda.org/anaconda/pandas)
4. [pyyaml](https://anaconda.org/anaconda/pyyaml/)

# How to use

## Web scraping module

PAU4Chem is a modular framework that uses web scraping for extracing the TRI information from the web and organizes it before the data engineering.
In order to run the web scraping module, navigate to the folder [extract](https://github.com/jodhernandezbe/PAU4Chem/tree/master/extract). Then, you execute the following command either on Windows CMD or Unix terminal:

```
python tri_web_scraper.py -Y TRI_Year -F TRI_File
```

The flag -Y represents the TRI reporting year the you want to get, while -F is the file from the TRI for retrieving the information (e.g., File 1a). Check [TRI Basic Plus Data Files Guides
](https://www.epa.gov/toxics-release-inventory-tri-program/tri-basic-plus-data-files-guides). PAU4Chem requires the files 1a, 2a, and 2b to run the data engineering.

## Data engineering module
