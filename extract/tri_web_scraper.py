# tri_web_scraper.py (PAU4Chem)
# !/usr/bin/env python3
# coding=utf-8


"""TODO describe this file/module."""

# Importing libraries
from bs4 import BeautifulSoup
import requests
import zipfile
import time
import io
import os
import chardet
import codecs
import pandas as pd
import argparse
from common import config


class TRI_Scraper:
    """TODO describe this class."""

    def __init__(self, year, Files):
        """TODO describe this method.

        Args:
            year ([type]): [description]
            Files ([type]): [description]
        """
        self.year = year
        # Working Directory
        self._dir_path = os.path.dirname(os.path.realpath(__file__))
        # self._dir_path = os.getcwd() # if you are working on Jupyter Notebook
        self._config = config()['web_sites']['TRI']
        self._queries = self._config['queries']
        # Uniform Resource Locator (URL) of TRI Database
        self._url = self._config['url']
        self._TRI_File_Columns_Dictionary = {}  # TRI File Formats
        for File in Files:
            self._TRI_File_Columns_Dictionary[File] = []

    def _visit(self):
        """TODO describe this method."""
        html = requests.get(self._url).text
        time.sleep(15)
        self._soup = BeautifulSoup(html, 'html.parser')

    def _extracting_zip(self):
        """TODO describe this method."""
        self._visit()
        self._TRI_zip_options = {}
        for link in self._soup.find_all(self._queries['TRI_year_reported']):
            self._TRI_zip_options[link.text] = link.get(
                self._queries['TRI_zip'])

    def _Calling_TRI_columns(self):
        """TODO describe this method."""
        Path_columns = self._dir_path + '/../ancillary'
        for key in self._TRI_File_Columns_Dictionary.keys():
            inf_chardet = chardet.detect(open(
                Path_Columns + '/TRI_File_' + key + '_columns.txt',
                'rb').read())
            inf_encoding = inf_chardet['encoding']
            file_TRI_File_columns = codecs.open(
                Path_Columns + '/TRI_File_' + key + '_columns.txt', 'r',
                encoding=inf_encoding)
            TRI_File_columns = [
                TRI_File_columns.rstrip().replace('\n', '') for
                TRI_File_columns in file_TRI_File_columns.readlines()]
            self._TRI_File_Columns_Dictionary[key] = TRI_File_columns
            file_TRI_File_columns.close()

    # Method for Extracting  information according to TRI's year report
    def extacting_TRI_data_files(self):
        """TODO describe this method."""
        self._Calling_TRI_columns()
        self._extracting_zip()
        TRI_zip = self._TRI_zip_options[self.year]
        r_file = requests.get(TRI_zip)
        for key in self._TRI_File_Columns_Dictionary.keys():
            # Unzipping
            with zipfile.ZipFile(io.BytesIO(r_file.content)) as z:
                z.extract('US_' + key + '_' + self.year + '.txt',
                          self._dir_path + '/datasets')

            # Converting .txt to .csv
            df = pd.read_csv(
                self._dir_path + '/datasets/US_' + key + '_' +
                self.year + '.txt',
                header=None, encoding='ISO-8859-1', error_bad_lines=False,
                sep='\t', low_memory=False, skiprows=[0], lineterminator='\n',
                usecols=range(len(self._TRI_File_Columns_Dictionary[key])))
            df.columns = self._TRI_File_Columns_Dictionary[key]
            df.to_csv(self._dir_path + '/datasets/US_' + key + '_' +
                      self.year + '.csv',
                      sep=',', index=False)
            time.sleep(30)
            os.remove(self._dir_path + '/datasets/US_' + key + '_' +
                      self.year + '.txt')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS)

    parser.add_argument('-Y', '--Year', nargs='+',
                        help='What TRI year do you want to retrieve?.',
                        type=str)

    parser.add_argument(
        '-F', '--Files', nargs='+',
        help='What TRI Files do you want (e.g., 1a, 2a, etc).\n' +
        'Check:\nhttps://www.epa.gov/toxics-release-inventory-tri-program/' +
        'tri-basic-plus-data-files-guides',
        required=False)

    args = parser.parse_args()

    for Y in args.Year:
        Scraper = TRI_Scraper(Y, args.Files)
        Scraper.extacting_TRI_data_files()
