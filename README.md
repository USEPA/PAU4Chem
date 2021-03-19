# PAU4Chem

These Python codes build the PAU dataset for tracking and allocating the chemical flows in pollution abatement units (PAUs). The folder [datasets](https://github.com/jodhernandezbe/PAU4Chem/tree/master/transform/datasets) contains the outputs for each framework step. The [Chemicals_in_categories.csv](https://github.com/jodhernandezbe/PAU4Chem/blob/master/transform/chemicals/Chemicals_in_categories.csv) contains the chemicals for the TRI chemical categories.

<p align="center">
  <img src=https://github.com/jodhernandezbe/PAU4Chem/blob/master/PAU_dataset.svg width="80%">
</p>

## Requirements

This code was written using Python 3.x, Anaconda 3, and operating system Ubuntu 18.04.
Testing was completed in Windows 10 OS with Visual Studio Code Insiders Edition.
The following Python libraries are required for running the code:

1. [bs4](https://anaconda.org/conda-forge/bs4)
2. [requests](https://anaconda.org/anaconda/requests)
3. [pandas](https://anaconda.org/anaconda/pandas)
4. [pyyaml](https://anaconda.org/anaconda/pyyaml/)
5. [scipy](https://anaconda.org/conda-forge/scipy)
6. [numpy](https://anaconda.org/conda-forge/numpy)
7. [xlrd](https://anaconda.org/anaconda/xlrd)

## How to use

### Web scraping module

PAU4Chem is a modular framework that uses web scraping for extracing the TRI information from the web and organizes it before the data engineering.
In order to run the web scraping module, navigate to the folder [extract](https://github.com/jodhernandezbe/PAU4Chem/tree/master/extract). Then, you execute the following command either on Windows CMD or Unix terminal:

```bash
python tri_web_scraper.py -Y TRI_Year -F TRI_File
```

The flag -Y represents the TRI reporting year that you would like to get, while -F is the file from the TRI for retrieving the information (e.g., File 1a). Check [TRI Basic Plus Data Files Guides
](https://www.epa.gov/toxics-release-inventory-tri-program/tri-basic-plus-data-files-guides) for knowing more about the TRI files. PAU4Chem requires the files 1a, 1b, and 2b to run the data engineering.

The above command can work for retrieving many years and files at the same time, for example:

```bash
python tri_web_scraper.py -Y 2001 2002 2003 -F 1a 1b 2b
```

### Data engineering module

In order to run the steps in the data engineering, navigate to the folder [transform](https://github.com/jodhernandezbe/PAU4Chem/tree/master/transform). Execute the following command either on Windows CMD or Unix terminal:

```bash
python building_pau_db.py --help
```

You could see the following menu:

```bash
positional arguments:
  Option                What do you want to do:
                          [A]: Recover information from TRI.
                          [B]: File for statistics.
                          [C]: File for recycling.
                          [D]: Further cleaning of database.
                          [E]: Organizing file with flows (1987-2004).
                          [F]: Organizing file with substance prices (1987 - 2004).
                          [G]: Pollution abatement cost and expenditure (only 2004).
                          [H]: Pollution control unit positions (1987 - 2004).
                          [I]: Searching information for years after 2004

optional arguments:
  -h, --help            show this help message and exit
  -Y YEAR [YEAR ...], --Year YEAR [YEAR ...]
                        Records with up to how many PAUs you want to include?.
  -N_Bins N_BINS        Number of bins to split the middle waste flow values
```

The argument N_Bins is only used for running the option E that is for obtaining the waste input flows. Please, consider the information for the years provided for each option. For example, option E only works for years from 1987 to 2004. Execute each step starting with A and finishing with G. Additionally, start using the year 1987 and finish with 2018 (the most recent report). As example:

```bash
python building_pau_db.py B -Y 2004
```

However, you can run each option for many years at the same time:

```bash
python building_pau_db.py A -Y 1987 1988 1989
```

### Disclaimer

The views expressed in this article are those of the authors and do not necessarily represent the views or policies of
the U.S. Environmental Protection Agency. Any mention of trade names, products, or services does not imply an endorsement by the U.S.
Government or the U.S. Environmental Protection Agency. The U.S. Environmental Protection Agency does not endorse any commercial products, service, or enterprises.

### Acknowledgement

This research was supported in by an appointment for Jose D. Hernandez-Betancur to the Research Participation
Program at the Center for Environmental Solutions and Emergency Response, Office of Research and Development,
U.S. Environmental Protection Agency, administered by the Oak Ridge Institute for Science and Education through an Interagency Agreement No. DW-89-92433001 between the U.S. Department of Energy and the U.S. Environmental Protection Agency.
