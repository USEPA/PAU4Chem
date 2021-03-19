# building_pau_db.py (PAU4Chem)
# !/usr/bin/env python3
# coding=utf-8


"""TODO describe this file/module."""

# Note:
# 1. Range of Influent Concentration was reported from 1987 through 2004
# 2. Treatment Efficiency Estimation was reported from 1987 through 2004

from scipy.stats import norm
from scipy.stats import lognorm
from itertools import combinations
from population import *

import os
import argparse
import numpy as np
import re
import time
import unicodedata
import yaml
import math

import warnings
import pandas as pd

warnings.simplefilter(action='ignore', category=FutureWarning)
pd.options.mode.chained_assignment = None


class PAU_DB:
    """TODO describe this class."""

    def __init__(self, Year):
        """TODO describe this method.

        Args:
            Year ([type]): [description]
        """
        # Working Directory
        self._dir_path = os.path.dirname(os.path.realpath(__file__))
        self.Year = Year
        # self._dir_path = os.getcwd() # if you are working on Jupyter Notebook

    def calling_tri_files(self):
        """TODO describe this method.

        Returns:
            [type]: [description]
        """
        TRI_Files = dict()
        for file in ['1a', '1b', '2b']:
            columns = pd.read_csv(
                self._dir_path + '/../ancillary/TRI_File_' + file +
                '_needed_columns.txt', header=None)
            columns = list(columns.iloc[:, 0])
            df = pd.read_csv(
                self._dir_path + '/../extract/datasets/US_' + file + '_' +
                str(self.Year) + '.csv', usecols=columns, low_memory=False)
            df = df.where(pd.notnull(df), None)
            TRI_Files.update({file: df})
        return TRI_Files

    def is_number(self, s):
        """TODO describe this method.

        Args:
            s ([type]): [description]

        Returns:
            [type]: [description]
        """
        try:
            float(s)
            return True
        except (TypeError, ValueError):
            pass
        try:
            unicodedata.numeric(s)
            return True
        except (TypeError, ValueError):
            pass
        return False

    def _efficiency_estimation_to_range(self, x):
        """TODO describe this method.

        Args:
            x ([type]): [description]

        Returns:
            [type]: [description]
        """
        if x != np.nan:
            x = np.abs(x)
            if (x >= 0.0) & (x <= 50.0):
                return 'E6'
            elif (x > 50.0) & (x <= 95.0):
                return 'E5'
            elif (x > 95.0) & (x <= 99.0):
                return 'E4'
            elif (x > 99.0) & (x <= 99.99):
                return 'E3'
            elif (x > 99.99) & (x <= 99.9999):
                return 'E2'
            elif (x > 99.9999):
                return 'E1'
        else:
            return None

    def _effic_est_empties_based_on_epa_regulation(self, classification,
                                                   HAP, RCRA):
        """TODO describe this method.

        Args:
            classification ([type]): [description]
            HAP ([type]): [description]
            RCRA ([type]): [description]

        Returns:
            [type]: [description]
        """
        if RCRA == 'YES':
            if classification == 'DIOXIN':
                result = np.random.uniform(99.9999, 100)
                if self.Year >= 2005:
                    result = self._efficiency_estimation_to_range(result)
            else:
                result = np.random.uniform(99.99, 100)
                if self.Year >= 2005:
                    result = self._efficiency_estimation_to_range(result)
            return result
        elif HAP == 'YES':
            result = np.random.uniform(95, 100)
            if self.Year >= 2005:
                result = self._efficiency_estimation_to_range(result)
            return result
        else:
            return None

    def _calling_SRS(self):
        """TODO describe this method.

        Returns:
            [type]: [description]
        """
        Acronyms = ['TRI', 'CAA', 'RCRA_F', 'RCRA_K', 'RCRA_P',
                    'RCRA_T', 'RCRA_U']
        Files = {Acronym: File for File in os.listdir(self._dir_path + '/srs')
                 for Acronym in Acronyms if Acronym in File}
        columns = ['ID', 'Internal Tracking Number']
        df = pd.read_csv(
            self._dir_path + '/srs/' + Files['TRI'], low_memory=False,
            usecols=['ID', 'Internal Tracking Number'],
            converters={
                'ID': lambda x: str(int(x)) if re.search('^\\d', x) else x},
            dtype={'Internal Tracking Number': 'int64'})
        df = df.assign(HAP=['NO']*df.shape[0], RCRA=['NO']*df.shape[0])
        del Files['TRI']
        for Acronym, File in Files.items():
            col = 'HAP'
            if Acronym in Acronyms[2:]:
                col = 'RCRA'
            ITN = pd.read_csv(self._dir_path + '/srs/' + File,
                              low_memory=False,
                              usecols=['Internal Tracking Number'],
                              dtype={'Internal Tracking Number': 'int64'})
            df.loc[df['Internal Tracking Number'].isin(
                ITN['Internal Tracking Number'].tolist()), col] = 'YES'
        df.drop(columns='Internal Tracking Number', inplace=True)
        df.rename(columns={'ID': 'CAS NUMBER'}, inplace=True)
        return df

    def _changin_management_code_for_2004_and_prior(self, x, m_n):
        """TODO describe this method.

        Args:
            x ([type]): [description]
            m_n ([type]): [description]

        Returns:
            [type]: [description]
        """
        Change = pd.read_csv(
            self._dir_path + '/../ancillary/Methods_TRI.csv',
            usecols=['Code 2004 and prior', 'Code 2005 and after'],
            low_memory=False)
        if list(x.values).count(None) != m_n:
            y = {v: 'T' if v in Change['Code 2004 and prior'].unique().tolist()
                 else 'F' for v in x.values}
            result = [Change.loc[Change['Code 2004 and prior'] == v,
                                 'Code 2005 and after'].tolist()[0]
                      if s == 'T' else None for v, s in y.items()]
            L = len(result)
            result = result + [None]*(m_n - L)
            return result
        else:
            return [None]*m_n

    def organizing(self):
        """TODO describe this method.

        Returns:
            [type]: [description]
        """
        dfs = self.calling_tri_files()
        df = dfs['2b'].where(pd.notnull(dfs['2b']), None)
        if self.Year >= 2005:
            df.drop(
                columns=df.iloc[:, list(range(18, 71, 13))].columns.tolist(),
                inplace=True)
        else:
            df.drop(
                columns=df.iloc[:, list(range(20, 73, 13))].columns.tolist(),
                inplace=True)
        df_PAUs = pd.DataFrame()
        Columns_0 = list(df.iloc[:, 0:8].columns)
        for i in range(5):
            Starting = 8 + 12*i
            Ending = Starting + 11
            Columns_1 = list(df.iloc[:, Starting:Ending + 1].columns)
            columns = Columns_0 + Columns_1
            df_aux = df[Columns]
            Columns_to_change = {col: re.sub(r'STREAM [1-5] - ', '', col)
                                 for col in Columns_1}
            df_aux.rename(columns=Columns_to_change, inplace=True)
            df_PAUs = pd.concat([df_PAUs, df_aux], ignore_index=True,
                                sort=True, axis=0)
            del Columns
        del df, df_aux
        cols = list(df_PAUs.iloc[:, 9:17].columns)
        df_PAUs.dropna(subset=cols, how='all', axis=0, inplace=True)
        if self.Year <= 2004:
            df_PAUs.dropna(
                subset=['WASTE STREAM CODE', 'RANGE INFLUENT CONCENTRATION',
                        'TREATMENT EFFICIENCY ESTIMATION'], how='any',
                axis=0, inplace=True)
            df_PAUs.reset_index(inplace=True, drop=True)
            df_PAUs['METHOD CODE - 2004 AND PRIOR'] = df_PAUs[cols].apply(
                lambda x: None if list(x).count(None) == len(cols)
                else ' + '.join(xx for xx in x if xx), axis=1)
            df_PAUs[cols] = df_PAUs.apply(lambda row: pd.Series(
                self._changin_management_code_for_2004_and_prior(
                    row[cols], len(cols))), axis=1)
            df_PAUs = df_PAUs.loc[pd.notnull(df_PAUs[cols]).any(axis=1)]
            df_PAUs['EFFICIENCY RANGE CODE'] = df_PAUs[
                'TREATMENT EFFICIENCY ESTIMATION'].apply(
                    lambda x: self._efficiency_estimation_to_range(float(x)))
            df_PAUs.rename(columns={
                'TREATMENT EFFICIENCY ESTIMATION': 'EFFICIENCY ESTIMATION'},
                           inplace=True)
            mask = pd.to_numeric(df_PAUs['RANGE INFLUENT CONCENTRATION'],
                                 errors='coerce').notnull()
            df_PAUs = df_PAUs[mask]
            df_PAUs['RANGE INFLUENT CONCENTRATION'] = df_PAUs[
                'RANGE INFLUENT CONCENTRATION'].apply(lambda x: abs(int(x)))
        else:
            df_PAUs.rename(columns={
                'TREATMENT EFFICIENCY RANGE CODE': 'EFFICIENCY RANGE CODE'},
                           inplace=True)
            df_PAUs.dropna(subset=['WASTE STREAM CODE',
                                   'EFFICIENCY RANGE CODE'],
                           how='any', axis=0, inplace=True)
        df_PAUs['METHOD CODE - 2005 AND AFTER'] = df_PAUs[cols].apply(
            lambda x: None if list(x).count(None) == len(cols) else ' + '.join(
                xx for xx in x if xx), axis=1)
        df_PAUs = df_PAUs.loc[pd.notnull(
            df_PAUs['METHOD CODE - 2005 AND AFTER'])]
        df_PAUs['TYPE OF MANAGEMENT'] = 'Treatment'
        df_PAUs.drop(columns=cols, inplace=True)
        df_PAUs.reset_index(inplace=True, drop=True)
        df_PAUs.loc[pd.isnull(df_PAUs['BASED ON OPERATING DATA?']),
                    'BASED ON OPERATING DATA?'] = 'NO'

        try:
            # On-site energy recovery
            df = dfs['1a'].iloc[:, list(range(12))]
            cols = [c for c in df.columns if 'METHOD' in c]
            df.dropna(subset=cols, how='all', axis=0, inplace=True)
            Columns_0 = list(df.iloc[:, 0:8].columns)
            Columns_1 = list(df.iloc[:, 8:].columns)
            dfs_energy = pd.DataFrame()
            for col in Columns_1:
                columns = Columns_0 + [col]
                df_aux = df[Columns]
                df_aux.rename(columns={col: re.sub(r' [1-4]', '', col)},
                              inplace=True)
                dfs_energy = pd.concat([dfs_energy, df_aux], ignore_index=True,
                                       sort=True, axis=0)
                del Columns

            del df, df_aux

            dfs_energy = dfs_energy.loc[pd.notnull(
                dfs_energy['ON-SITE ENERGY RECOVERY METHOD'])]
            dfs_energy['TYPE OF MANAGEMENT'] = 'Energy recovery'
            if self.Year <= 2004:
                dfs_energy['METHOD CODE - 2004 AND PRIOR'] = dfs_energy[
                    'ON-SITE ENERGY RECOVERY METHOD']
                dfs_energy['ON-SITE ENERGY RECOVERY METHOD'] = \
                    dfs_energy.apply(lambda row: pd.Series(
                        self._changin_management_code_for_2004_and_prior(
                            pd.Series(row['ON-SITE ENERGY RECOVERY METHOD']),
                            1)), axis=1)
                dfs_energy = dfs_energy.loc[pd.notnull(dfs_energy[
                    'ON-SITE ENERGY RECOVERY METHOD'])]
            dfs_energy.rename(columns={
                'ON-SITE ENERGY RECOVERY METHOD':
                    'METHOD CODE - 2005 AND AFTER'}, inplace=True)
            dfs_energy = dfs_energy.loc[pd.notnull(dfs_energy[
                'METHOD CODE - 2005 AND AFTER'])]
            df_PAUs = pd.concat([df_PAUs, dfs_energy], ignore_index=True,
                                sort=True, axis=0)
            del dfs_energy
        except ValueError as e:
            print(e)
            print('There is not information about energy recovery activities')

        try:
            # On-site recycling
            df = dfs['1a'].iloc[:, list(range(8)) + list(range(12, 19))]
            cols = [c for c in df.columns if 'METHOD' in c]
            df.dropna(subset=cols, how='all', axis=0, inplace=True)
            Columns_0 = list(df.iloc[:, 0:8].columns)
            Columns_1 = list(df.iloc[:, 8:].columns)
            dfs_recycling = pd.DataFrame()
            for col in Columns_1:
                columns = Columns_0 + [col]
                df_aux = df[Columns]
                df_aux.rename(columns={col: re.sub(r' [1-7]', '', col)},
                              inplace=True)
                dfs_recycling = pd.concat([dfs_recycling, df_aux],
                                          ignore_index=True,
                                          sort=True, axis=0)
                del Columns
            del df, df_aux
            dfs_recycling = dfs_recycling.loc[pd.notnull(dfs_recycling[
                'ON-SITE RECYCLING PROCESSES METHOD'])]
            dfs_recycling['TYPE OF MANAGEMENT'] = 'Recycling'
            if self.Year <= 2004:
                dfs_recycling['METHOD CODE - 2004 AND PRIOR'] = dfs_recycling[
                    'ON-SITE RECYCLING PROCESSES METHOD']
                dfs_recycling['ON-SITE RECYCLING PROCESSES METHOD'] = \
                    dfs_recycling.apply(lambda row: pd.Series(
                        self._changin_management_code_for_2004_and_prior(
                            pd.Series(row[
                                'ON-SITE RECYCLING PROCESSES METHOD']), 1)),
                                        axis=1)
                dfs_recycling = dfs_recycling.loc[pd.notnull(dfs_recycling[
                    'ON-SITE RECYCLING PROCESSES METHOD'])]
            dfs_recycling.rename(columns={
                'ON-SITE RECYCLING PROCESSES METHOD':
                    'METHOD CODE - 2005 AND AFTER'},
                                 inplace=True)
            dfs_recycling = dfs_recycling.loc[pd.notnull(
                dfs_recycling['METHOD CODE - 2005 AND AFTER'])]
            df_PAUs = pd.concat([df_PAUs, dfs_recycling], ignore_index=True,
                                sort=True, axis=0)
            del dfs_recycling
        except ValueError as e:
            print(e)
            print('There is not information about recycling activities')

        # Changing units
        df_PAUs = df_PAUs.loc[(df_PAUs.iloc[:0, :] != 'INV').all(axis=1)]
        df_PAUs.dropna(how='all', axis=0, inplace=True)
        if self.Year >= 2005:
            Change = pd.read_csv(
                self._dir_path + '/../ancillary/Methods_TRI.csv',
                usecols=['Code 2004 and prior', 'Code 2005 and after'],
                low_memory=False)
            Codes_2004 = Change.loc[(pd.notnull(
                Change['Code 2004 and prior'])) & (Change[
                    'Code 2005 and after'] != Change['Code 2004 and prior']),
                'Code 2004 and prior'].unique().tolist()
            idx = df_PAUs.loc[df_PAUs['METHOD CODE - 2005 AND AFTER'].isin(
                Codes_2004)].index.tolist()
            del Change, Codes_2004

            if len(idx) != 0:
                df_PAUs.loc[idx, 'METHOD CODE - 2005 AND AFTER'] = df_PAUs.loc[
                    idx
                ].apply(lambda row:
                        self._changin_management_code_for_2004_and_prior(
                            pd.Series(row['METHOD CODE - 2005 AND AFTER']), 1),
                        axis=1)

        # Adding methods name
        Methods = pd.read_csv(
            self._dir_path + '/../ancillary/Methods_TRI.csv',
            usecols=['Code 2004 and prior', 'Method 2004 and prior',
                     'Code 2005 and after', 'Method 2005 and after'])
        Methods.drop_duplicates(keep='first', inplace=True)

        # Adding chemical activities and uses
        df_PAUs['DOCUMENT CONTROL NUMBER'] = df_PAUs[
            'DOCUMENT CONTROL NUMBER'].apply(
                lambda x: str(int(float(x))) if self.is_number(x) else x)
        dfs['1b'].drop_duplicates(keep='first', inplace=True)
        dfs['1b']['DOCUMENT CONTROL NUMBER'] = dfs['1b'][
            'DOCUMENT CONTROL NUMBER'].apply(
                lambda x: str(int(float(x))) if self.is_number(x) else x)
        df_PAUs = pd.merge(
            df_PAUs, dfs['1b'],
            on=['TRIFID', 'DOCUMENT CONTROL NUMBER', 'CAS NUMBER'],
            how='inner')
        columns_DB_F = [
            'REPORTING YEAR', 'TRIFID', 'PRIMARY NAICS CODE', 'CAS NUMBER',
            'CHEMICAL NAME', 'METAL INDICATOR', 'CLASSIFICATION',
            'PRODUCE THE CHEMICAL', 'IMPORT THE CHEMICAL',
            'ON-SITE USE OF THE CHEMICAL',
            'SALE OR DISTRIBUTION OF THE CHEMICAL',
            'AS A BYPRODUCT', 'AS A MANUFACTURED IMPURITY',
            'USED AS A REACTANT',
            'ADDED AS A FORMULATION COMPONENT', 'USED AS AN ARTICLE COMPONENT',
            'REPACKAGING', 'AS A PROCESS IMPURITY', 'RECYCLING',
            'USED AS A CHEMICAL PROCESSING AID',
            'USED AS A MANUFACTURING AID',
            'ANCILLARY OR OTHER USE', 'WASTE STREAM CODE',
            'METHOD CODE - 2005 AND AFTER',
            'METHOD NAME - 2005 AND AFTER', 'TYPE OF MANAGEMENT',
            'EFFICIENCY RANGE CODE', 'BASED ON OPERATING DATA?']

        if self.Year <= 2004:
            Method = {row.iloc[0]: row.iloc[1] for index, row in Methods[[
                'Code 2004 and prior', 'Method 2004 and prior']].iterrows()}

            def _checking(x, M):
                if x:
                    return ' + '.join(M[xx] for xx in x.split(' + ') if xx
                                      and xx in M.keys())
                else:
                    return None

            df_PAUs = df_PAUs.loc[df_PAUs[
                'METHOD CODE - 2004 AND PRIOR'].str.contains(r'[A-Z]').where(
                    df_PAUs['METHOD CODE - 2004 AND PRIOR'].str.contains(
                        r'[A-Z]'), False)]
            df_PAUs['METHOD NAME - 2004 AND PRIOR'] = df_PAUs[
                'METHOD CODE - 2004 AND PRIOR'].apply(lambda x: _checking(
                    x, Method))
            df_PAUs = df_PAUs.loc[
                (df_PAUs['METHOD CODE - 2004 AND PRIOR'] != '') | (
                    pd.notnull(df_PAUs['METHOD CODE - 2004 AND PRIOR']))]
            columns_DB_F = [
                'REPORTING YEAR', 'TRIFID', 'PRIMARY NAICS CODE', 'CAS NUMBER',
                'CHEMICAL NAME', 'METAL INDICATOR', 'CLASSIFICATION',
                'PRODUCE THE CHEMICAL', 'IMPORT THE CHEMICAL',
                'ON-SITE USE OF THE CHEMICAL',
                'SALE OR DISTRIBUTION OF THE CHEMICAL', 'AS A BYPRODUCT',
                'AS A MANUFACTURED IMPURITY', 'USED AS A REACTANT',
                'ADDED AS A FORMULATION COMPONENT',
                'USED AS AN ARTICLE COMPONENT',
                'REPACKAGING', 'AS A PROCESS IMPURITY', 'RECYCLING',
                'USED AS A CHEMICAL PROCESSING AID',
                'USED AS A MANUFACTURING AID', 'ANCILLARY OR OTHER USE',
                'WASTE STREAM CODE', 'RANGE INFLUENT CONCENTRATION',
                'METHOD CODE - 2004 AND PRIOR', 'METHOD NAME - 2004 AND PRIOR',
                'METHOD CODE - 2005 AND AFTER', 'METHOD NAME - 2005 AND AFTER',
                'TYPE OF MANAGEMENT', 'EFFICIENCY RANGE CODE',
                'EFFICIENCY ESTIMATION', 'BASED ON OPERATING DATA?']
        Method = {row.iloc[0]: row.iloc[1] for index, row in Methods[[
            'Code 2005 and after', 'Method 2005 and after']].iterrows()}
        df_PAUs = df_PAUs.loc[df_PAUs[
            'METHOD CODE - 2005 AND AFTER'].str.contains(r'[A-Z]').where(
                df_PAUs['METHOD CODE - 2005 AND AFTER'].str.contains(
                    r'[A-Z]'), False)]
        df_PAUs['METHOD NAME - 2005 AND AFTER'] = df_PAUs[
            'METHOD CODE - 2005 AND AFTER'].apply(lambda x: ' + '.join(
                Method[xx] for xx in x.split(' + ') if xx and
                xx in Method.keys()))

        # Saving information
        df_PAUs['REPORTING YEAR'] = self.Year
        df_PAUs = df_PAUs[columns_DB_F]
        df_PAUs.to_csv(
            self._dir_path + '/datasets/intermediate_pau_datasets/PAUs_DB_' +
            str(self.Year) + '.csv', sep=',', index=False)

    def Building_database_for_statistics(self):
        """TODO describe this method."""
        columns = pd.read_csv(
            self._dir_path +
            '/../ancillary/TRI_File_2b_needed_columns_for_statistics.txt',
            header=None)
        columns = list(columns.iloc[:0, ])
        df = pd.read_csv(
            self._dir_path + '/../extract/datasets/US_2b_' +
            str(self.Year) + '.csv',
            usecols=columns, low_memory=False)
        df_statistics = pd.DataFrame()

        if self.Year >= 2005:
            df.drop(
                columns=df.iloc[:, list(range(12, 61, 12))].columns.tolist(),
                inplace=True)
            codes_incineration = ['A01', 'H040', 'H076', 'H122']
        else:
            df.drop(
                columns=df.iloc[:, list(range(13, 62, 12))].columns.tolist(),
                inplace=True)
            codes_incineration = ['A01', 'F01', 'F11', 'F19', 'F31',
                                  'F41', 'F42', 'F51', 'F61',
                                  'F71', 'F81', 'F82', 'F83',
                                  'F99']
        Columns_0 = list(df.iloc[:, 0:2].columns)
        for i in range(5):
            Columns_1 = list(
                df.iloc[:, [2 + 11*i, 11 + 11*i, 12 + 11*i]].columns)
            Treatmes = list(df.iloc[:, 3 + 11*i: 11 + 11*i].columns)
            columns = Columns_0 + Columns_1
            df_aux = df[Columns]
            df_aux['INCINERATION'] = 'NO'
            df_aux.loc[df[Treatmes].isin(codes_incineration).any(
                axis=1), 'INCINERATION'] = 'YES'
            df_aux['IDEAL'] = df[Treatmes].apply(
                lambda x: 'YES' if
                len(list(np.where(pd.notnull(x))[0])) == 1 else 'NO',
                axis=1)
            Columns_to_change = {col: re.sub(
                r'STREAM [1-5] - ', '', col) for col in Columns_1}
            df_aux.rename(columns=Columns_to_change, inplace=True)
            df_statistics = pd.concat(
                [df_statistics, df_aux], ignore_index=True, sort=True, axis=0)
            del Columns
        del df, df_aux
        if self.Year <= 2004:
            df_statistics.dropna(how='any', axis=0, inplace=True)
            mask = pd.to_numeric(df_statistics[
                'TREATMENT EFFICIENCY ESTIMATION'], errors='coerce').notnull()
            df_statistics = df_statistics[mask]
            df_statistics['TREATMENT EFFICIENCY ESTIMATION'] = df_statistics[
                'TREATMENT EFFICIENCY ESTIMATION'].astype(float)
            df_statistics['EFFICIENCY RANGE'] = df_statistics[
                'TREATMENT EFFICIENCY ESTIMATION'].apply(
                    lambda x: self._efficiency_estimation_to_range(float(x)))
            mask = pd.to_numeric(df_statistics['RANGE INFLUENT CONCENTRATION'],
                                 errors='coerce').notnull()
            df_statistics = df_statistics[mask]
            df_statistics['RANGE INFLUENT CONCENTRATION'] = df_statistics[
                'RANGE INFLUENT CONCENTRATION'].astype(int)
            df_statistics.rename(
                columns={
                    'TREATMENT EFFICIENCY ESTIMATION': 'EFFICIENCY ESTIMATION'
                }, inplace=True)
        else:
            df_statistics.rename(columns={
                'TREATMENT EFFICIENCY RANGE CODE': 'EFFICIENCY RANGE'},
                                 inplace=True)
            df_statistics.dropna(
                subset=['EFFICIENCY RANGE', 'WASTE STREAM CODE'], how='any',
                axis=0, inplace=True)
        df_statistics.rename(
            columns={'PRIMARY NAICS CODE': 'NAICS',
                     'CAS NUMBER': 'CAS',
                     'WASTE STREAM CODE': 'WASTE',
                     'RANGE INFLUENT CONCENTRATION': 'CONCENTRATION'},
            inplace=True)
        df_statistics.loc[
            df_statistics['INCINERATION'] == 'NO', 'IDEAL'] = None
        df_statistics.to_csv(
            self._dir_path + '/statistics/db_for_general/DB_for_Statistics_' +
            str(self.Year) + '.csv',
            sep=',', index=False)

    def Building_database_for_recycling_efficiency(self):
        """TODO describe this method."""
        def _division(row, elements_total):
            if row['ON-SITE - RECYCLED'] == 0.0:
                if row['CLASSIFICATION'] == 'TRI':
                    row['ON-SITE - RECYCLED'] = 0.5
                elif row['CLASSIFICATION'] == 'PBT':
                    row['ON-SITE - RECYCLED'] = 0.1
                else:
                    row['ON-SITE - RECYCLED'] = 0.0001
            values = [abs(v) for v in row[elements_total] if v != 0.0]
            cases = list()
            for n_elements_sum in range(1, len(values) + 1):
                comb = combinations(values, n_elements_sum)
                for comb_values in comb:
                    sumatory = sum(comb_values)
                    cases.append(row['ON-SITE - RECYCLED']/(row[
                        'ON-SITE - RECYCLED'] + sumatory)*100)
            try:
                if len(list(set(cases))) == 1 and cases[0] == 100:
                    return [100]*6 + [0] + [row['ON-SITE - RECYCLED']]
                else:
                    return [np.min(cases)] + np.quantile(
                        cases, [0.25, 0.5, 0.75]).tolist() + \
                           [np.max(cases), np.mean(cases),
                            np.std(cases)/np.mean(cases),
                            row['ON-SITE - RECYCLED']]
            except ValueError:
                a = np.empty((8))
                a[:] = np.nan
                return a.tolist()
        columns = pd.read_csv(
            self._dir_path +
            '/../ancillary/TRI_File_1a_needed_columns_for_statistics.txt',
            header=None)
        columns = list(columns.iloc[:0, ])
        df = pd.read_csv(
            self._dir_path + '/../extract/datasets/US_1a_' + str(self.Year) +
            '.csv', usecols=columns, low_memory=False)
        elements_total = list(set(df.iloc[:, 5:64].columns.tolist()) -
                              set(['ON-SITE - RECYCLED']))
        df.iloc[:, 5:64] = df.iloc[:, 5:64].where(
            pd.notnull(df.iloc[:, 5:64]), 0.0)
        df.iloc[:, 5:64] = df.iloc[:, 5:64].apply(
            pd.to_numeric, errors='coerce')
        cols = [c for c in df.columns if 'METHOD' in c]
        df['IDEAL'] = df[cols].apply(
            lambda x: 'YES' if len(list(np.where(pd.notnull(x))[0])) == 1
            else 'NO', axis=1)
        df = df.loc[df['IDEAL'] == 'YES']
        df['METHOD'] = df[cols].apply(lambda x: x.values[
            np.where(pd.notnull(x))[0]][0], axis=1)
        df.drop(columns=['IDEAL'] + cols, inplace=True)
        df = df.loc[df['METHOD'] != 'INV']
        df[['LOWER EFFICIENCY', 'Q1', 'Q2', 'Q3', 'UPPER EFFICIENCY',
            'MEAN OF EFFICIENCY', 'CV', 'ON-SITE - RECYCLED']] = df.apply(
                lambda x: pd.Series(_division(x, elements_total)), axis=1)
        df = df.loc[pd.notnull(df['UPPER EFFICIENCY'])]
        df['IQR'] = df.apply(lambda x: x['Q3'] - x['Q1'], axis=1)
        df['Q1 - 1.5xIQR'] = df.apply(
            lambda x: 0 if x['Q1'] - 1.5*x['IQR'] < 0
            else x['Q1'] - 1.5*x['IQR'], axis=1)
        df['Q3 + 1.5xIQR'] = df.apply(
            lambda x: 100 if x['Q3'] + 1.5*x['IQR'] > 100
            else x['Q3'] + 1.5*x['IQR'], axis=1)
        df['UPPER EFFICIENCY OUTLIER?'] = df.apply(
            lambda x: 'YES' if x['UPPER EFFICIENCY'] >
            x['Q3 + 1.5xIQR'] else 'NO', axis=1)
        df['LOWER EFFICIENCY OUTLIER?'] = df.apply(
            lambda x: 'YES' if x['LOWER EFFICIENCY'] <
            x['Q1 - 1.5xIQR'] else 'NO', axis=1)
        df['HIGH VARIANCE?'] = df.apply(
            lambda x: 'YES' if x['CV'] > 1 else 'NO', axis=1)
        df = df[[
            'TRIFID', 'PRIMARY NAICS CODE', 'CAS NUMBER', 'ON-SITE - RECYCLED',
            'UNIT OF MEASURE', 'LOWER EFFICIENCY', 'LOWER EFFICIENCY OUTLIER?',
            'Q1 - 1.5xIQR', 'Q1', 'Q2', 'Q3', 'Q3 + 1.5xIQR',
            'UPPER EFFICIENCY', 'UPPER EFFICIENCY OUTLIER?', 'IQR',
            'MEAN OF EFFICIENCY', 'CV', 'HIGH VARIANCE?', 'METHOD']]
        df.iloc[:, [5, 7, 8, 9, 10, 11, 12, 14, 15, 16]] = df.iloc[
            :, [5, 7, 8, 9, 10, 11, 12, 14, 15, 16]].round(4)
        df.to_csv(
            self._dir_path + '/statistics/db_for_solvents/DB_for_Solvents_' +
            str(self.Year) + '.csv', sep=',', index=False)

    def _searching_naics(self, x, naics):
        """TODO describe this method.

        Args:
            x ([type]): [description]
            naics ([type]): [description]

        Returns:
            [type]: [description]
        """
        # https://www.census.gov/programs-surveys/economic-census/guidance/understanding-naics.html
        values = {0: 'Nothing', 1: 'Nothing', 2: 'Sector', 3: 'Subsector',
                  4: 'Industry Group', 5: 'NAICS Industry',
                  6: 'National Industry'}
        naics = str(naics)
        x = str(x)
        equal = 0
        for idx, char in enumerate(naics):
            try:
                if char == x[idx]:
                    equal = equal + 1
                else:
                    break
            except IndexError:
                break
        return values[equal]

    def _phase_estimation_recycling(self, df_s, row):
        """TODO describe this method.

        Args:
            df_s ([type]): [description]
            row ([type]): [description]

        Returns:
            [type]: [description]
        """
        # Solvent recovery
        if row['METHOD CODE - 2005 AND AFTER'] == 'H20':
            phases = ['L']
        # Acid regeneration and other reactions
        elif row['METHOD CODE - 2005 AND AFTER'] == 'H39':
            phases = ['W']
        # Metal recovery
        elif row['METHOD CODE - 2005 AND AFTER'] == 'H10':
            phases = ['W', 'S']
            if self.Year <= 2004:
                Pyrometallurgy = ['R27', 'R28', 'R29']  # They work with scrap
                if row['METHOD CODE - 2004 AND PRIOR'] in Pyrometallurgy:
                    Phases = ['S']
                else:
                    Phases = ['W', 'S']
        naics_structure = [
            'National Industry', 'NAICS Industry', 'Industry Group',
            'Subsector', 'Sector', 'Nothing']
        df_cas = df_s.loc[df_s['CAS'] == row['CAS NUMBER'],
                          ['NAICS', 'WASTE', 'VALUE']]
        df_cas = df_cas.groupby(['NAICS', 'WASTE'], as_index=False).sum()
        df_cas.reset_index(inplace=True)
        if (not df_cas.empty):
            df_cas['NAICS STRUCTURE'] = df_cas.apply(
                lambda x: self._searching_naics(x['NAICS'],
                                                row['PRIMARY NAICS CODE']),
                axis=1)
            i = 0
            phase = None
            while i <= 5 and phase not in phases:
                structure = naics_structure[i]
                i = i + 1
                # for structure in naics_structure:
                df_naics = df_cas.loc[df_cas['NAICS STRUCTURE'] == structure]
                if (df_naics.empty):
                    phase = None
                    # continue
                else:
                    if (df_naics['WASTE'].isin(phases).any()):
                        df_phase = df_naics.loc[df_naics['WASTE'].isin(phases)]
                        row['NAICS STRUCTURE'] = structure
                        row['WASTE STREAM CODE'] = df_phase.loc[
                            df_phase['VALUE'].idxmax(), 'WASTE']
                    else:
                        row['NAICS STRUCTURE'] = structure
                        row['WASTE STREAM CODE'] = df_naics.loc[
                            df_naics['VALUE'].idxmax(), 'WASTE']
                    phase = row['WASTE STREAM CODE']
            return row
        else:
            row['NAICS STRUCTURE'] = None
            row['WASTE STREAM CODE'] = None
            return row

    def _concentration_estimation_recycling(self, df_s, cas, naics,
                                            phase, structure):
        """TODO describe this method.

        Args:
            df_s ([type]): [description]
            cas ([type]): [description]
            naics ([type]): [description]
            phase ([type]): [description]
            structure ([type]): [description]

        Returns:
            [type]: [description]
        """
        df_s = df_s[['NAICS', 'CAS', 'WASTE', 'CONCENTRATION', 'VALUE']]
        df_s = df_s.loc[(df_s['CAS'] == cas) &
                        (df_s['WASTE'] == phase)]
        df_s = df_s.groupby(['NAICS', 'CAS', 'WASTE', 'CONCENTRATION'],
                            as_index=False).sum()
        df_s['NAICS STRUCTURE'] = df_s.apply(
            lambda x: self._searching_naics(x['NAICS'],
                                            naics), axis=1)
        df = df_s.loc[(df_s['NAICS STRUCTURE'] == structure)]
        return df.loc[df['VALUE'].idxmax(), 'CONCENTRATION']

    def _recycling_efficiency(self, row, df_s):
        """TODO describe this method.

        Args:
            row ([type]): [description]
            df_s ([type]): [description]

        Returns:
            [type]: [description]
        """
        naics_structure = [
            'National Industry', 'NAICS Industry', 'Industry Group',
            'Subsector', 'Sector', 'Nothing']
        if self.Year <= 2004:
            code = row['METHOD CODE - 2004 AND PRIOR']
        else:
            code = row['METHOD CODE - 2005 AND AFTER']
        df_cas = df_s.loc[(df_s['CAS NUMBER'] == row['CAS NUMBER']) &
                          (df_s['METHOD'] == code)]
        if (not df_cas.empty):
            df_fid = df_cas.loc[df_cas['TRIFID'] == row['TRIFID']]
            if (not df_fid.empty):
                return df_fid['UPPER EFFICIENCY'].iloc[0]
            else:
                df_cas['NAICS STRUCTURE'] = df_cas.apply(
                    lambda x: self._searching_naics(x['PRIMARY NAICS CODE'],
                                                    row['PRIMARY NAICS CODE']),
                    axis=1)
                i = 0
                efficiency = None
                while (i <= 5) and (not efficiency):
                    structure = naics_structure[i]
                    i = i + 1
                    df_naics = df_cas.loc[
                        df_cas['NAICS STRUCTURE'] == structure]
                    if df_naics.empty:
                        efficiency = None
                    else:
                        efficiency = df_naics['UPPER EFFICIENCY'].median()
                return efficiency
        else:
            return None

    def _phase_estimation_energy(self, df_s, row):
        """TODO describe this method.

        Args:
            df_s ([type]): [description]
            row ([type]): [description]

        Returns:
            [type]: [description]
        """
        phases = ['S', 'L', 'A']
        if row['METHOD CODE - 2005 AND AFTER'] == 'U01':
            # Industrial Kilns (specially rotatory kilns) are used to
            # burn hazardous liquid and solid wastes
            phases = ['S', 'L']
        naics_structure = [
            'National Industry', 'NAICS Industry', 'Industry Group',
            'Subsector', 'Sector', 'Nothing']
        df_cas = df_s.loc[df_s['CAS'] == row['CAS NUMBER'],
                          ['NAICS', 'WASTE', 'VALUE', 'INCINERATION']]
        df_cas = df_cas.groupby(['NAICS', 'WASTE', 'INCINERATION'],
                                as_index=False).sum()
        df_cas.reset_index(inplace=True)
        if (not df_cas.empty):
            df_cas['NAICS STRUCTURE'] = df_cas.apply(
                lambda x: self._searching_naics(x['NAICS'],
                                                row['PRIMARY NAICS CODE']),
                axis=1)
            i = 0
            phase = None

            # for structure in naics_structure:
            while i <= 5 and phase not in phases:
                structure = naics_structure[i]
                i = i + 1
                df_naics = df_cas.loc[df_cas['NAICS STRUCTURE'] == structure]
                if df_naics.empty:
                    phase = None
                else:
                    df_incineration = df_naics.loc[
                        df_cas['INCINERATION'] == 'YES']
                    if df_incineration.empty:
                        if (df_naics['WASTE'].isin(phases).any()):
                            df_phase = df_naics.loc[
                                df_naics['WASTE'].isin(phases)]
                            row['NAICS STRUCTURE'] = structure
                            row['WASTE STREAM CODE'] = df_phase.loc[
                                df_phase['VALUE'].idxmax(), 'WASTE']
                        else:
                            row['NAICS STRUCTURE'] = structure
                            row['WASTE STREAM CODE'] = df_naics.loc[
                                df_naics['VALUE'].idxmax(), 'WASTE']
                        row['BY MEANS OF INCINERATION'] = 'NO'
                    else:
                        if (df_incineration['WASTE'].isin(phases).any()):
                            df_phase = df_incineration.loc[df_incineration[
                                'WASTE'].isin(phases)]
                            row['NAICS STRUCTURE'] = structure
                            row['WASTE STREAM CODE'] = df_phase.loc[
                                df_phase['VALUE'].idxmax(), 'WASTE']
                        else:
                            row['NAICS STRUCTURE'] = structure
                            row['WASTE STREAM CODE'] = df_incineration.loc[
                                df_incineration['VALUE'].idxmax(), 'WASTE']
                        row['BY MEANS OF INCINERATION'] = 'YES'
                    phase = row['WASTE STREAM CODE']
            return row
        else:
            row['NAICS STRUCTURE'] = None
            row['WASTE STREAM CODE'] = None
            row['BY MEANS OF INCINERATION'] = None
            return row

    def _concentration_estimation_energy(self, df_s, cas, naics, phase,
                                         structure, incineration):
        """TODO describe this method.

        Args:
            df_s ([type]): [description]
            cas ([type]): [description]
            naics ([type]): [description]
            phase ([type]): [description]
            structure ([type]): [description]
            incineration ([type]): [description]

        Returns:
            [type]: [description]
        """
        df_s = df_s[['NAICS', 'CAS', 'WASTE', 'CONCENTRATION',
                     'VALUE', 'INCINERATION']]
        df_s = df_s.loc[(df_s['CAS'] == cas) &
                        (df_s['WASTE'] == phase) &
                        (df_s['INCINERATION'] == incineration)]
        df_s = df_s.groupby(['NAICS', 'CAS', 'WASTE',
                             'CONCENTRATION', 'INCINERATION'],
                            as_index=False).sum()
        df_s['NAICS STRUCTURE'] = df_s.apply(
            lambda x: self._searching_naics(x['NAICS'],
                                            naics), axis=1)
        df = df_s.loc[(df_s['NAICS STRUCTURE'] == structure)]
        return df.loc[df['VALUE'].idxmax(), 'CONCENTRATION']

    def _energy_efficiency(self, df_s, row):
        """TODO describe this method.

        Args:
            df_s ([type]): [description]
            row ([type]): [description]

        Returns:
            [type]: [description]
        """
        if self.Year <= 2004:
            df_s = df_s[['NAICS', 'CAS', 'WASTE', 'INCINERATION', 'IDEAL',
                         'EFFICIENCY ESTIMATION']]
        else:
            df_s = df_s[['NAICS', 'CAS', 'WASTE', 'INCINERATION', 'IDEAL',
                         'EFFICIENCY RANGE', 'VALUE']]
            df_s['WASTE'] = df_s.groupby(['NAICS', 'CAS', 'WASTE', 'IDEAL',
                                          'INCINERATION', 'EFFICIENCY RANGE'],
                                         as_index=False).sum()
        df_s = df_s.loc[(df_s['CAS'] == row['CAS NUMBER']) &
                        (df_s['INCINERATION'] == 'YES') &
                        (df_s['IDEAL'] == 'YES')]
        if (not df_s.empty):
            df_s['NAICS STRUCTURE'] = df_s.apply(
                lambda x: self._searching_naics(x['NAICS'],
                                                row['PRIMARY NAICS CODE']),
                axis=1)
            df_structure = df_s.loc[
                df_s['NAICS STRUCTURE'] == row['NAICS STRUCTURE']]
            if not df_structure.empty:
                df_phase = df_structure.loc[
                    df_structure['WASTE'] == row['WASTE STREAM CODE']]
                if not df_phase.empty:
                    if self.Year <= 2004:
                        result = df_phase['EFFICIENCY ESTIMATION'].median()
                    else:
                        result = df_phase.loc[df_phase['VALUE'].idxmax(),
                                              'EFFICIENCY RANGE']
                else:
                    if self.Year <= 2004:
                        result = df_structure['EFFICIENCY ESTIMATION'].median()
                    else:
                        result = df_structure.loc[
                            df_structure['VALUE'].idxmax(), 'EFFICIENCY RANGE']
            else:
                df_phase = df_s.loc[df_s['WASTE'] == row['WASTE STREAM CODE']]
                if not df_phase.empty:
                    if self.Year <= 2004:
                        result = df_phase['EFFICIENCY ESTIMATION'].median()
                    else:
                        result = df_phase.loc[df_phase['VALUE'].idxmax(),
                                              'EFFICIENCY RANGE']
                else:
                    if self.Year <= 2004:
                        result = df_s['EFFICIENCY ESTIMATION'].median()
                    else:
                        result = df_s.loc[df_s['VALUE'].idxmax(),
                                          'EFFICIENCY RANGE']
        else:
            return None
        return result

    def cleaning_database(self):
        """TODO describe this method."""
        # Calling TRI restriction for metals
        Restrictions = pd.read_csv(
            self._dir_path +
            '/../ancillary/Metals_divided_into_4_groups_can_be_reported.csv',
            low_memory=False, usecols=['ID', "U01, U02, U03 (Energy recovery)",
                                       'H20 (Solvent recovey)'])
        Energy_recovery = Restrictions.loc[Restrictions[
            "U01, U02, U03 (Energy recovery)"] == 'NO', 'ID'].tolist()
        Solvent_recovery = Restrictions.loc[Restrictions[
            'H20 (Solvent recovey)'] == 'NO', 'ID'].tolist()

        # Calling PAU
        PAU = pd.read_csv(
            self._dir_path + '/datasets/intermediate_pau_datasets/PAUs_DB_' +
            str(self.Year) + '.csv',
            low_memory=False, converters={
                'CAS NUMBER': lambda x: x if re.search(r'^[A-Z]', x)
                else str(int(x))})
        columns_DB_F = PAU.columns.tolist()
        PAU['PRIMARY NAICS CODE'] = PAU['PRIMARY NAICS CODE'].astype('int')
        if self.Year <= 2004:
            grouping = ['TRIFID', 'METHOD CODE - 2004 AND PRIOR']
            PAU.sort_values(by=['PRIMARY NAICS CODE', 'TRIFID',
                            'METHOD CODE - 2004 AND PRIOR', 'CAS NUMBER'],
                            inplace=True)
        else:
            grouping = ['TRIFID', 'METHOD CODE - 2005 AND AFTER']
            PAU.sort_values(by=['PRIMARY NAICS CODE', 'TRIFID',
                            'METHOD CODE - 2005 AND AFTER', 'CAS NUMBER'],
                            inplace=True)
        # Calling database for statistics
        Statistics = pd.read_csv(
            self._dir_path + '/statistics/db_for_general/DB_for_Statistics_' +
            str(self.Year) + '.csv', low_memory=False,
            converters={'CAS': lambda x: x if re.search(r'^[A-Z]', x)
                        else str(int(x))})
        Statistics['NAICS'] = Statistics['NAICS'].astype('int')
        Statistics['VALUE'] = 1
        Statistics.sort_values(by=['NAICS', 'CAS'], inplace=True)
        # Treatment
        Efficiency_codes = ['E1', 'E2', 'E3', 'E4', 'E5', 'E6']
        df_N_PAU = PAU.loc[PAU['TYPE OF MANAGEMENT'] == 'Treatment']
        df_N_PAU = df_N_PAU.loc[df_N_PAU['EFFICIENCY RANGE CODE'].isin(
            Efficiency_codes)]
        # Recycling
        PAU_recycling = PAU.loc[PAU['TYPE OF MANAGEMENT'] == 'Recycling']
        if not PAU_recycling.empty:
            PAU_recycling = PAU_recycling.loc[~((
                PAU_recycling['METHOD CODE - 2005 AND AFTER'] == 'H20'
                ) & (PAU_recycling['CAS NUMBER'].isin(Solvent_recovery)))]
            PAU_recycling.reset_index(inplace=True, drop=True)
            PAU_recycling['BASED ON OPERATING DATA?'] = 'NO'
            # Calling database for recycling efficiency
            Recycling_statistics = pd.read_csv(
                self._dir_path +
                '/statistics/db_for_solvents/DB_for_Solvents_' +
                str(self.Year) + '.csv',
                low_memory=False, usecols=[
                    'TRIFID', 'PRIMARY NAICS CODE', 'CAS NUMBER',
                    'UPPER EFFICIENCY', 'UPPER EFFICIENCY OUTLIER?',
                    'METHOD', 'HIGH VARIANCE?'], converters={
                        'CAS NUMBER': lambda x: x if re.search(r'^[A-Z]', x)
                        else str(int(x))})
            Recycling_statistics['PRIMARY NAICS CODE'] = Recycling_statistics[
                'PRIMARY NAICS CODE'].astype('int')
            Recycling_statistics = Recycling_statistics.loc[
                (Recycling_statistics['UPPER EFFICIENCY OUTLIER?'] == 'NO') &
                (Recycling_statistics['HIGH VARIANCE?'] == 'NO')]
            Recycling_statistics.drop(
                columns=['UPPER EFFICIENCY OUTLIER?', 'HIGH VARIANCE?'],
                axis=1)
            efficiency_estimation = PAU_recycling.apply(
                lambda x: self._recycling_efficiency(x, Recycling_statistics),
                axis=1).round(4)
            PAU_recycling['EFFICIENCY RANGE CODE'] = \
                efficiency_estimation.apply(
                    lambda x: self._efficiency_estimation_to_range(x))
            PAU_recycling = PAU_recycling.loc[pd.notnull(
                PAU_recycling['EFFICIENCY RANGE CODE'])]
            PAU_recycling = PAU_recycling.apply(
                lambda x: self._phase_estimation_recycling(Statistics, x),
                axis=1)
            PAU_recycling = PAU_recycling.loc[pd.notnull(
                PAU_recycling['WASTE STREAM CODE'])]
            if self.Year <= 2004:
                PAU_recycling['EFFICIENCY ESTIMATION'] = efficiency_estimation
                PAU_recycling['RANGE INFLUENT CONCENTRATION'] = \
                    PAU_recycling.apply(
                        lambda x: self._concentration_estimation_recycling(
                            Statistics, x['CAS NUMBER'],
                            x['PRIMARY NAICS CODE'], x['WASTE STREAM CODE'],
                            x['NAICS STRUCTURE']), axis=1)
            PAU_recycling.drop(columns=['NAICS STRUCTURE'], inplace=True)
            df_N_PAU = pd.concat([df_N_PAU, PAU_recycling],
                                 ignore_index=True, sort=True, axis=0)
        else:
            pass

        # Energy recovery
        PAU_energy = PAU.loc[PAU['TYPE OF MANAGEMENT'] == 'Energy recovery']
        if not PAU_energy.empty:
            PAU_energy = PAU_energy.loc[~((PAU_energy[
                'METHOD CODE - 2005 AND AFTER'].isin(
                    ['U01', 'U02', 'U03'])) & (
                        PAU_energy['CAS NUMBER'].isin(Energy_recovery)))]
            PAU_energy.reset_index(inplace=True, drop=True)
            PAU_energy['BASED ON OPERATING DATA?'] = 'NO'
            PAU_energy = PAU_energy.apply(
                lambda x: self._phase_estimation_energy(Statistics, x), axis=1)
            PAU_energy = PAU_energy.loc[pd.notnull(
                PAU_energy['WASTE STREAM CODE'])]
            SRS = self._calling_SRS()
            if self.Year <= 2004:
                PAU_energy['RANGE INFLUENT CONCENTRATION'] = PAU_energy.apply(
                    lambda x: self._concentration_estimation_energy(
                        Statistics, x['CAS NUMBER'], x['PRIMARY NAICS CODE'],
                        x['WASTE STREAM CODE'], x['NAICS STRUCTURE'],
                        x['BY MEANS OF INCINERATION']), axis=1)
                PAU_energy.drop(columns=['BY MEANS OF INCINERATION'],
                                inplace=True)
                PAU_energy['EFFICIENCY ESTIMATION'] = PAU_energy.apply(
                    lambda x: self._energy_efficiency(Statistics, x),
                    axis=1).round(4)
                PAU_energy = pd.merge(PAU_energy, SRS,
                                      on='CAS NUMBER', how='left')
                PAU_energy['EFFICIENCY ESTIMATION'] = PAU_energy.apply(
                    lambda x: self._effic_est_empties_based_on_epa_regulation(
                        x['CLASSIFICATION'], x['HAP'], x['RCRA'])
                    if not x['EFFICIENCY ESTIMATION'] else
                    x['EFFICIENCY ESTIMATION'], axis=1)
                PAU_energy = PAU_energy.loc[pd.notnull(
                    PAU_energy['EFFICIENCY ESTIMATION'])]
                PAU_energy['EFFICIENCY RANGE CODE'] = PAU_energy[
                    'EFFICIENCY ESTIMATION'].apply(
                        lambda x: self._efficiency_estimation_to_range(
                            float(x)))

            else:
                PAU_energy.drop(columns=['BY MEANS OF INCINERATION'],
                                inplace=True)
                PAU_energy['EFFICIENCY RANGE CODE'] = PAU_energy.apply(
                    lambda x: self._energy_efficiency(Statistics, x), axis=1)
                PAU_energy = pd.merge(PAU_energy, SRS,
                                      on='CAS NUMBER', how='left')
                PAU_energy['EFFICIENCY RANGE CODE'] = PAU_energy.apply(
                    lambda x: self._effic_est_empties_based_on_epa_regulation(
                        x['CLASSIFICATION'], x['HAP'], x['RCRA'])
                    if not x['EFFICIENCY RANGE CODE'] else
                    x['EFFICIENCY RANGE CODE'], axis=1)
                PAU_energy = PAU_energy.loc[pd.notnull(
                    PAU_energy['EFFICIENCY RANGE CODE'])]
            PAU_energy.drop(columns=['NAICS STRUCTURE', 'HAP', 'RCRA'],
                            inplace=True)
            PAU_energy.loc[
                (PAU_energy['WASTE STREAM CODE'] == 'W') &
                (PAU_energy['TYPE OF MANAGEMENT'] == 'Energy recovery'),
                'WASTE STREAM CODE'] = 'L'
            df_N_PAU = pd.concat([df_N_PAU, PAU_energy],
                                 ignore_index=True, sort=True, axis=0)
        else:
            pass
        Chemicals_to_remove = ['MIXTURE', 'TRD SECRT']
        df_N_PAU = df_N_PAU.loc[~df_N_PAU['CAS NUMBER'].isin(
            Chemicals_to_remove)]
        df_N_PAU['CAS NUMBER'] = df_N_PAU['CAS NUMBER'].apply(
            lambda x: str(int(x)) if 'N' not in x else x)
        df_N_PAU = df_N_PAU[columns_DB_F]
        df_N_PAU.to_csv(
            self._dir_path + '/datasets/final_pau_datasets/PAUs_DB_filled_' +
            str(self.Year) + '.csv', sep=',', index=False)

        # Chemicals and groups
        Chemicals = df_N_PAU[['CAS NUMBER', 'CHEMICAL NAME']].drop_duplicates(
            keep='first')
        Chemicals['TYPE OF CHEMICAL'] = None
        Path_c = self._dir_path + '/chemicals/Chemicals.csv'
        if os.path.exists(Path_c):
            df_c = pd.read_csv(Path_c)
            for index, row in Chemicals.iterrows():
                if (df_c['CAS NUMBER'] != row['CAS NUMBER']).all():
                    df_c = df_c.append(
                        pd.Series(row, index=row.index.tolist()),
                        ignore_index=True)
            df_c.to_csv(Path_c, sep=',', index=False)
        else:
            Chemicals.to_csv(Path_c, sep=',', index=False)

    def _Calculating_possible_waste_feed_supply(self, Flow, Concentration,
                                                Efficiency):
        """TODO describe this method.

        Args:
            Flow ([type]): [description]
            Concentration ([type]): [description]
            Efficiency ([type]): [description]

        Returns:
            [type]: [description]
        """
        if Concentration == '1':
            percentanges_c = [1, 100]
        elif Concentration == '2':
            percentanges_c = [0.01, 1]
        elif Concentration == '3':
            percentanges_c = [0.0001, 0.01]
        elif Concentration == '4':
            percentanges_c = [0.0000001, 0.0001]
        elif Concentration == '5':
            percentanges_c = [0.000000001, 0.0000001]
        if Efficiency != 0.0:
            Chemical_feed_flow = 100*Flow/Efficiency
        else:
            Chemical_feed_flow = Flow
        Waste_flows = [100*Chemical_feed_flow/c for c in percentanges_c]
        Interval = tuple([min(Waste_flows), max(Waste_flows)])
        Middle = 0.5*(Waste_flows[0] + Waste_flows[-1])
        return Interval, Middle

    def Building_database_for_flows(self, nbins):
        """TODO describe this method.

        Args:
            nbins ([type]): [description]
        """

        def func(x):
            if x.first_valid_index() is None:
                return None
            else:
                return x[x.first_valid_index()]

        with open(self._dir_path + '/../ancillary/Flow_columns.yaml',
                  mode='r') as f:
            dictionary_of_columns = yaml.load(f, Loader=yaml.FullLoader)
        dictionary_of_columns = {
            key: [el.strip() for el in val['columns'].split(',')]
            for key, val in dictionary_of_columns['TRI_Files'].items()}
        dfs = dict()
        for file, columns in dictionary_of_columns.items():
            df = pd.read_csv(
                self._dir_path + '/../extract/datasets/US_{}_{}.csv'.format(
                    file, self.Year), usecols=columns, low_memory=False,
                dtype={'PRIMARY NAICS CODE': 'object'})
            dfs.update({file: df})
        # Energy recovery:
        cols_energy_methods = [col for col in dfs['1a'].columns if
                               'METHOD' in col and 'ENERGY' in col]
        df_energy = dfs['1a'][['TRIFID', 'CAS NUMBER', 'UNIT OF MEASURE',
                               'ON-SITE - ENERGY RECOVERY'] +
                              cols_energy_methods]
        df_energy = df_energy.loc[pd.notnull(
            df_energy[cols_energy_methods]).sum(axis=1) == 1]
        df_energy['METHOD CODE'] = df_energy[
            cols_energy_methods].apply(func, axis=1)
        df_energy.rename(
            columns={'ON-SITE - ENERGY RECOVERY': 'FLOW'}, inplace=True)
        df_energy.drop(columns=cols_energy_methods, inplace=True)
        del cols_energy_methods

        # Recycling:
        cols_recycling_methods = [col for col in dfs['1a'].columns if
                                  'METHOD' in col and 'RECYCLING' in col]
        df_recycling = dfs['1a'][['TRIFID', 'CAS NUMBER', 'UNIT OF MEASURE',
                                 'ON-SITE - RECYCLED'] +
                                 cols_recycling_methods]
        df_recycling = df_recycling.loc[pd.notnull(df_recycling[
            cols_recycling_methods]).sum(axis=1) == 1]
        df_recycling['METHOD CODE'] = df_recycling[
            cols_recycling_methods].apply(func, axis=1)
        df_recycling.rename(columns={'ON-SITE - RECYCLED': 'FLOW'},
                            inplace=True)
        df_recycling.drop(columns=cols_recycling_methods, inplace=True)
        del cols_recycling_methods

        # Treatment
        cols_treatment_methods = [col for col in dfs['2b'].columns if
                                  'METHOD' in col]
        dfs['2b'] = dfs['2b'].loc[pd.notnull(dfs['2b'][
            cols_treatment_methods]).sum(axis=1) == 1]
        dfs['2b']['METHOD CODE'] = dfs['2b'][
            cols_treatment_methods].apply(func, axis=1)
        dfs['2b'].drop(columns=cols_treatment_methods, inplace=True)
        cols_for_merging = ['TRIFID', 'DOCUMENT CONTROL NUMBER',
                            'CAS NUMBER', 'ON-SITE - TREATED',
                            'UNIT OF MEASURE']
        df_treatment = pd.merge(dfs['1a'][cols_for_merging], dfs['2b'],
                                how='inner',
                                on=['TRIFID',
                                    'DOCUMENT CONTROL NUMBER',
                                    'CAS NUMBER'])
        del dfs, cols_treatment_methods, cols_for_merging
        df_treatment.rename(columns={'ON-SITE - TREATED': 'FLOW'},
                            inplace=True)
        df_treatment.drop(columns=['DOCUMENT CONTROL NUMBER'], inplace=True)
        df_PAU_flows = pd.concat([df_treatment, df_recycling, df_energy],
                                 ignore_index=True,
                                 sort=True, axis=0)
        del df_treatment, df_recycling, df_energy
        Chemicals_to_remove = ['MIXTURE', 'TRD SECRT']
        df_PAU_flows = df_PAU_flows.loc[~df_PAU_flows['CAS NUMBER'].isin(
            Chemicals_to_remove)]
        df_PAU_flows['CAS NUMBER'] = df_PAU_flows['CAS NUMBER'].apply(
            lambda x: str(int(x)) if 'N' not in x else x)
        df_PAU_flows.loc[df_PAU_flows['UNIT OF MEASURE'] == 'Pounds',
                         'FLOW'] *= 0.453592
        df_PAU_flows.loc[df_PAU_flows['UNIT OF MEASURE'] == 'Grams',
                         'FLOW'] *= 10**-3
        df_PAU_flows['FLOW'] = df_PAU_flows['FLOW'].round(6)
        df_PAU_flows = df_PAU_flows.loc[df_PAU_flows['FLOW'] != 0.0]
        df_PAU_flows['UNIT OF MEASURE'] = 'kg'

        # Calling cleaned database
        columns_for_calling = [
            'TRIFID', 'CAS NUMBER', 'RANGE INFLUENT CONCENTRATION',
            'METHOD CODE - 2004 AND PRIOR', 'EFFICIENCY ESTIMATION',
            'PRIMARY NAICS CODE', 'WASTE STREAM CODE', 'TYPE OF MANAGEMENT']
        df_PAU_cleaned = pd.read_csv(
            self._dir_path +
            '/datasets/final_pau_datasets/PAUs_DB_filled_{}.csv'.format(
                self.Year),
            usecols=columns_for_calling,
            dtype={'PRIMARY NAICS CODE': 'object'})
        df_PAU_cleaned.rename(
            columns={'METHOD CODE - 2004 AND PRIOR': 'METHOD CODE'},
            inplace=True)

        # Merging
        df_PAU_flows = pd.merge(df_PAU_flows, df_PAU_cleaned, how='inner',
                                on=['TRIFID', 'CAS NUMBER', 'METHOD CODE'])
        df_PAU_flows[['RANGE INFLUENT CONCENTRATION', 'EFFICIENCY ESTIMATION',
                      'FLOW']] = df_PAU_flows[[
                          'RANGE INFLUENT CONCENTRATION',
                          'EFFICIENCY ESTIMATION', 'FLOW']].applymap(
                              lambda x: abs(x))
        df_PAU_flows['RANGE INFLUENT CONCENTRATION'] = df_PAU_flows[
            'RANGE INFLUENT CONCENTRATION'].apply(lambda x: str(int(x)))
        df_PAU_flows[['WASTE FLOW RANGE',
                      'MIDDLE WASTE FLOW']] = df_PAU_flows.apply(
                          lambda x: pd.Series(
                              self._Calculating_possible_waste_feed_supply(
                                  x['FLOW'], x['RANGE INFLUENT CONCENTRATION'],
                                  x['EFFICIENCY ESTIMATION'])), axis=1)
        Max_value = df_PAU_flows['MIDDLE WASTE FLOW'].max()
        Min_value = df_PAU_flows['MIDDLE WASTE FLOW'].min()
        Order_max = int(math.log10(Max_value)) - 1
        Order_min = math.ceil(math.log10(Min_value))
        Delta = (Order_max - Order_min)/(nbins - 2)
        Bin_values = [Min_value - 10**(math.log10(Min_value) - 1)]
        Bin_values = Bin_values + [10**(Order_min + Delta*n) for n in range(
            nbins - 1)]
        Bin_values = Bin_values + [Max_value]
        Bin_values.sort()
        Bin_labels = [str(val) for val in range(1, len(Bin_values))]
        df_PAU_flows['MIDDLE WASTE FLOW INTERVAL'] = pd.cut(
            df_PAU_flows['MIDDLE WASTE FLOW'], bins=Bin_values)
        df_PAU_flows['MIDDLE WASTE FLOW INTERVAL CODE'] = pd.cut(
            df_PAU_flows['MIDDLE WASTE FLOW'], bins=Bin_values,
            labels=Bin_labels, precision=0)
        df_PAU_flows.to_csv(
            self._dir_path +
            '/datasets/waste_flow/Waste_flow_to_PAUs_{}_{}.csv'.format(
                self.Year, nbins), sep=',', index=False)

    def Organizing_substance_prices(self):
        """TODO describe this method."""
        # Organizing information about prices
        df_prices = pd.read_csv(self._dir_path + '/prices/Chemical_Price.csv',
                                dtype={'CAS NUMBER': 'object'})
        df_prices = df_prices.loc[
            (pd.notnull(df_prices['PRICE'])) &
            (df_prices['PRICE'] != 'Not found')]

        File_exchange = [file for file in os.listdir(
            self._dir_path + '/prices') if 'Exchange' in file]
        df_exchange_rate = pd.read_csv(
            self._dir_path + '/prices/{}'.format(File_exchange[0]))
        Exchange_rate = {row['CURRENCY']: row['EXCHANGE RATE TO USD']
                         for idx, row in df_exchange_rate.iterrows()}
        del df_exchange_rate

        df_prices['PRICE'] = df_prices.apply(
            lambda x: Exchange_rate[x['CURRENCY']]*float(x['PRICE']), axis=1)
        df_prices['QUANTITY'] = df_prices['QUANTITY'].str.lower()
        df_prices['QUANTITY'] = df_prices['QUANTITY'].str.replace(' ', '')
        df_prices = df_prices[df_prices['QUANTITY'].str.contains(
            r'ton|[umk]{0,1}g')]
        idx = df_prices[
            ~df_prices['QUANTITY'].str.contains(r'x')].index.tolist()
        df_prices.loc[idx, 'QUANTITY'] = '1x' + df_prices.loc[idx, 'QUANTITY']
        df_prices[['TIMES', 'MASS', 'UNIT']] = df_prices[
            'QUANTITY'].str.extract('(\\d+)x(\\d+\\.?\\d*)(ton|[umk]{0,1}g)',
                                    expand=True)
        dictionary_mass = {'g': 1, 'mg': 0.001, 'kg': 1000,
                           'ug': 10**-6, 'ton': 907185}
        df_prices['QUANTITY'] = df_prices[['TIMES', 'MASS']].apply(
            lambda x: float(x.values[0])*float(x.values[1]),
            axis=1)
        df_prices['QUANTITY'] = df_prices[['QUANTITY', 'UNIT']].apply(
            lambda x: x.values[0]*dictionary_mass[x.values[1]], axis=1)
        df_prices['UNIT PRICE (USD/g)'] = df_prices[[
            'PRICE', 'QUANTITY']].apply(lambda x: x.values[0]/x.values[1],
                                        axis=1)
        df_prices.drop(columns=['COMPANY_NAME', 'COUNTRY', 'PURITY',
                                'CURRENCY', 'TIMES', 'MASS', 'UNIT',
                                'QUANTITY', 'PRICE'],
                       inplace=True)
        df_prices = df_prices.groupby('CAS NUMBER', as_index=False).min()

        # Calling PAU
        df_PAU = pd.read_csv(
            self._dir_path +
            f'/datasets/final_pau_datasets/PAUs_DB_filled_{self.Year}.csv',
            usecols=['TRIFID', 'CAS NUMBER', 'METHOD CODE - 2004 AND PRIOR'])
        df_PAU = df_PAU[
            ~df_PAU['METHOD CODE - 2004 AND PRIOR'].str.contains('+')]

        # Separating categories and chemicals
        categories = pd.read_csv(
            self._dir_path + '/chemicals/Chemicals_in_categories.csv',
            usecols=['CAS NUMBER', 'CATEGORY CODE'])
        categories['CAS NUMBER'] = categories[
            'CAS NUMBER'].str.replace('-', '')
        chemicals = pd.read_csv(
            self._dir_path + '/chemicals/Chemicals.csv',
            usecols=['CAS NUMBER'])
        chemicals = chemicals.loc[~chemicals['CAS NUMBER'].isin(
            list(categories['CATEGORY CODE'].unique())), 'CAS NUMBER'].tolist()
        df_PAU_chemicals = df_PAU.loc[df_PAU['CAS NUMBER'].isin(chemicals)]
        df_PAU_categories = df_PAU.loc[~df_PAU['CAS NUMBER'].isin(chemicals)]
        df_PAU_categories.rename(columns={'CAS NUMBER': 'CATEGORY CODE'},
                                 inplace=True)
        del chemicals, df_PAU

        # Merging prices with chemicals
        df_PAU_chemicals = pd.merge(df_PAU_chemicals, df_prices,
                                    how='inner',
                                    on='CAS NUMBER')

        # Calling CDR
        df_CDR = pd.read_csv(
            self._dir_path + '/cdr/Substances_by_facilities.csv',
            usecols=['STRIPPED_CHEMICAL_ID_NUMBER',
                     'PGM_SYS_ID'],
            dtype={'STRIPPED_CHEMICAL_ID_NUMBER': 'object'})
        df_CDR.rename(columns={'STRIPPED_CHEMICAL_ID_NUMBER': 'CAS NUMBER',
                               'PGM_SYS_ID': 'TRIFID'},
                      inplace=True)
        df_CDR = pd.merge(df_CDR, categories, on='CAS NUMBER', how='inner')
        df_PAU_categories = pd.merge(df_PAU_categories, df_CDR,
                                     on=['CATEGORY CODE', 'TRIFID'],
                                     how='inner')

        # Merging prices with categories
        df_PAU_categories = pd.merge(df_PAU_categories, df_prices,
                                     how='inner',
                                     on='CAS NUMBER')
        df_PAU_categories.drop(columns=['CATEGORY CODE'],
                               inplace=True)
        df_PAU = pd.concat([df_PAU_categories, df_PAU_chemicals],
                           ignore_index=True,
                           sort=True, axis=0)
        df_PAU.to_csv(
            self._dir_path +
            f'/datasets/chemical_price/Chemical_price_vs_PAU_{self.Year}.csv',
            sep=',', index=False)

    def pollution_abatement_cost_and_expenditure(self):
        """TODO describe this method."""
        # Calling PAU
        df_PAU = pd.read_csv(
            self._dir_path +
            '/datasets/waste_flow/Waste_flow_to_PAUs_2004_10.csv',
            low_memory=False,
            usecols=['PRIMARY NAICS CODE', 'WASTE STREAM CODE',
                     'TYPE OF MANAGEMENT', 'METHOD CODE', 'TRIFID',
                     'MIDDLE WASTE FLOW'],
            dtype={'PRIMARY NAICS CODE': 'object'})
        df_PAU['FLOW'] = df_PAU.groupby(
            ['TRIFID', 'METHOD CODE'])['MIDDLE WASTE FLOW'].transform('median')
        df_PAU.drop(columns=['TRIFID', 'METHOD CODE', 'MIDDLE WASTE FLOW'],
                    inplace=True)
        df_PAU.drop_duplicates(keep='first', inplace=True)
        df_PAU = df_PAU.groupby(['PRIMARY NAICS CODE', 'TYPE OF MANAGEMENT',
                                 'WASTE STREAM CODE'], as_index=False).agg(
                                     {'FLOW': ['mean', 'std']})
        df_PAU = pd.DataFrame.from_records(
            df_PAU.values, columns=['NAICS code', 'Activity', 'Media',
                                    'Mean flow', 'SD flow'])
        df_PAU.drop_duplicates(keep='first', inplace=True)
        df_PAU = df_PAU[df_PAU['NAICS code'].str.contains(r'^3[123]')]

        # Imputing coefficient of variation (<= 1)
        df_PAU['SD flow'] = df_PAU['SD flow'].fillna(df_PAU['Mean flow']*0.001)

        # Method of moments
        df_PAU['mu'] = df_PAU[['Mean flow', 'SD flow']].apply(
            lambda x: np.log(
                x.values[0]**2/(x.values[1]**2 + x.values[0]**2)**0.5),
            axis=1)
        df_PAU['theta_2'] = df_PAU[['Mean flow', 'SD flow']].apply(
            lambda x: np.log(x.values[1]**2/x.values[0]**2 + 1), axis=1)
        df_PAU_correlation = pd.read_csv(
            self._dir_path + '/../extract/datasets/US_1a_2004.csv',
            low_memory=False,
            usecols=['TRIFID', 'PRIMARY NAICS CODE', 'UNIT OF MEASURE',
                     'OFF-SITE - TOTAL TRANSFERRED FOR RECYCLING',
                     'OFF-SITE - TOTAL TRANSFERRED FOR ENERGY RECOVERY',
                     'OFF-SITE - TOTAL TRANSFERRED FOR TREATMENT',
                     'ON-SITE - RECYCLED', 'ON-SITE - ENERGY RECOVERY',
                     'ON-SITE - TREATED'])
        Flow_columns = ['OFF-SITE - TOTAL TRANSFERRED FOR RECYCLING',
                        'OFF-SITE - TOTAL TRANSFERRED FOR ENERGY RECOVERY',
                        'OFF-SITE - TOTAL TRANSFERRED FOR TREATMENT',
                        'ON-SITE - RECYCLED', 'ON-SITE - ENERGY RECOVERY',
                        'ON-SITE - TREATED']
        df_PAU_correlation.loc[
            df_PAU_correlation['UNIT OF MEASURE'] == 'Pounds',
            Flow_columns] *= 0.453592
        df_PAU_correlation.loc[df_PAU_correlation[
            'UNIT OF MEASURE'] == 'Grams', Flow_columns] *= 10**-3
        df_PAU_correlation['TREATMENT'] = df_PAU_correlation[
            ['OFF-SITE - TOTAL TRANSFERRED FOR TREATMENT',
             'ON-SITE - TREATED']].sum(axis=1)
        df_PAU_correlation['RECYCLING'] = df_PAU_correlation[
            ['OFF-SITE - TOTAL TRANSFERRED FOR RECYCLING',
             'ON-SITE - RECYCLED']].sum(axis=1)
        df_PAU_correlation['ENERGY RECOVERY'] = df_PAU_correlation[
            ['OFF-SITE - TOTAL TRANSFERRED FOR ENERGY RECOVERY',
             'ON-SITE - ENERGY RECOVERY']].sum(axis=1)
        df_PAU_correlation.drop(columns=[
            'UNIT OF MEASURE', 'OFF-SITE - TOTAL TRANSFERRED FOR RECYCLING',
            'OFF-SITE - TOTAL TRANSFERRED FOR ENERGY RECOVERY',
            'OFF-SITE - TOTAL TRANSFERRED FOR TREATMENT'], inplace=True)
        df_PAU_correlation = df_PAU_correlation.groupby(
            ['TRIFID', 'PRIMARY NAICS CODE'], as_index=False).sum()
        df_PAU_correlation.drop(columns=['TRIFID'], inplace=True)
        df_PAU_correlation = df_PAU_correlation.groupby(
            ['PRIMARY NAICS CODE'], as_index=False).sum()
        managements = {'TREATMENT': 'TREATED',
                       'RECYCLING': 'RECYCLED',
                       'ENERGY RECOVERY': 'ENERGY RECOVERY'}

        df_PAU_to_merge = pd.DataFrame()
        for management, on_site in managements.items():
            df_PAU_aux = pd.DataFrame()
            df_PAU_aux['NAICS code'] = df_PAU_correlation[
                'PRIMARY NAICS CODE'].apply(lambda x: str(int(float(x))))
            df_PAU_aux['% On-site flow'] = df_PAU_correlation[
                'ON-SITE - ' + on_site]*100/df_PAU_correlation[management]
            df_PAU_aux['Activity'] = management.capitalize()
            df_PAU_to_merge = pd.concat([df_PAU_to_merge, df_PAU_aux],
                                        ignore_index=True, sort=True, axis=0)
        df_PAU_correlation = df_PAU_to_merge.copy()
        df_PAU = pd.merge(df_PAU, df_PAU_correlation, how='left',
                          on=['Activity', 'NAICS code'])
        df_PAU.loc[df_PAU['% On-site flow'].isnull(), '% On-site flow'] = 0.0
        del df_PAU_aux, df_PAU_to_merge, df_PAU_correlation
        # U.S. Pollution Abatement Operating Costs - Survey 2005
        df_PAOC = pd.read_csv(
            self._dir_path +
            '/us_census_bureau/Pollution_Abatement_Operating_Costs_2005.csv',
            low_memory=False, header=None, skiprows=[0, 1],
            usecols=[0, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            names=['NAICS code', 'Total PAOC', 'Activity - treatment',
                   'Activity - prevention', 'Activity - recycling',
                   'Activity - disposal', 'Media - air', 'Media - water',
                   'Media - solid waste', 'RSE for total PAOC'])
        df_PAOC = df_PAOC.loc[pd.notnull(df_PAOC).all(axis=1)]
        # Substracting  prevention
        df_PAOC['Total PAOC'] = df_PAOC['Total PAOC'] - df_PAOC[
            'Activity - prevention'] - df_PAOC['Activity - disposal']
        df_PAOC.drop(columns=['Activity - prevention', 'Activity - disposal'],
                     inplace=True)

        # The media and the activity and supposed to be indepent events
        # Proportions activities
        col_activities = [col for col in df_PAOC.columns if 'Activity' in col]
        df_PAOC[col_activities] = df_PAOC[col_activities].div(
            df_PAOC[col_activities].sum(axis=1), axis=0)
        col_medias = [col for col in df_PAOC.columns if 'Media' in col]
        df_PAOC[col_medias] = df_PAOC[col_medias].div(df_PAOC[col_medias].sum(
            axis=1), axis=0)
        # U.S. Pollution Abatement Capital Expenditures - Survey 2005
        df_PACE = pd.read_csv(
            self._dir_path +
            '/us_census_bureau/' +
            'Pollution_Abatement_Capital_Expenditures_2005.csv',
            low_memory=False, header=None, skiprows=[0, 1],
            usecols=[0, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            names=['NAICS code', 'Total PACE', 'Activity - treatment',
                   'Activity - prevention', 'Activity - recycling',
                   'Activity - disposal', 'Media - air', 'Media - water',
                   'Media - solid waste', 'RSE for total PACE'])
        df_PACE = df_PACE.loc[pd.notnull(df_PACE).all(axis=1)]
        df_PACE['Total PACE'] = df_PACE['Total PACE'] - df_PACE[
            'Activity - prevention'] - df_PACE['Activity - disposal']
        df_PACE.drop(columns=['Activity - prevention', 'Activity - disposal'],
                     inplace=True)
        df_PACE[col_activities] = df_PACE[col_activities].div(
            df_PACE[col_activities].sum(axis=1), axis=0)
        df_PACE[col_medias] = df_PACE[col_medias].div(
            df_PACE[col_medias].sum(axis=1), axis=0)

        # Statistics of U.S. Businesses - Survey 2005
        # Sampled establishments in 2005
        df_SUSB = Organizing_sample(20378, self._dir_path)
        df_SUSB['Establishments (employees >= 20)'] = 1
        df_SUSB['Establishments (employees >= 20)'] = df_SUSB.groupby(
            'NAICS code')['Establishments (employees >= 20)'].transform('sum')
        df_SUSB['Total shipment'] = df_SUSB.groupby(
            'NAICS code')['Total shipment establishment'].transform('sum')
        df_SUSB['Info establishments'] = df_SUSB[[
            'Establishment', 'Total shipment establishment',
            'P-in-cluster']].apply(lambda x: x.tolist(), axis=1)
        df_SUSB.drop(
            columns=['P-cluster', 'Establishment', 'P-in-cluster',
                     'P-selected', 'Total shipment establishment', 'Unit'],
            inplace=True)
        df_SUSB = df_SUSB.groupby([
            'NAICS code', 'Total shipment',
            'Establishments (employees >= 20)'], as_index=False).agg(
                {'Info establishments': lambda x: {
                    val[0]: [val[1], val[2]] for idx, val in enumerate(x)}})

        # Joining sources from census
        df_PACE = pd.merge(df_PACE, df_SUSB, on='NAICS code', how='left')
        df_PAOC = pd.merge(df_PAOC, df_SUSB, on='NAICS code', how='left')

        # Searching higher naics levels (no in clusters but containing them)
        df_PACE = df_PACE.where(pd.notnull(df_PACE), None)
        df_PACE[['Establishments (employees >= 20)', 'Total shipment',
                 'Info establishments']] = df_PACE.apply(
                     lambda x: searching_establishments_by_hierarchy(
                         x['NAICS code'], df_SUSB)
                     if not x['Establishments (employees >= 20)']
                     else pd.Series(
                         [int(x['Establishments (employees >= 20)']),
                          x['Total shipment'],
                          x['Info establishments']]),
                     axis=1)
        df_PACE = df_PACE.loc[pd.notnull(df_PACE).all(axis=1)]
        df_PAOC = df_PAOC.where(pd.notnull(df_PAOC), None)
        df_PAOC[[
            'Establishments (employees >= 20)',
            'Total shipment', 'Info establishments']] = df_PAOC.apply(
                lambda x: searching_establishments_by_hierarchy(
                    x['NAICS code'], df_SUSB)
                if not x['Establishments (employees >= 20)']
                else pd.Series([int(
                    x['Establishments (employees >= 20)']),
                                x['Total shipment'],
                                x['Info establishments']]), axis=1)
        df_PAOC = df_PAOC.loc[pd.notnull(df_PAOC).all(axis=1)]
        del df_SUSB

        # Organizing by activity and media
        df_PACE_for_merging = pd.DataFrame()
        df_PAOC_for_merging = pd.DataFrame()
        Dictionary_relation = {
            'W': 'water', 'L': 'water', 'A': 'air', 'S': 'solid waste',
            'Treatment': 'Treatment', 'Energy recovery': 'Recycling',
            'Recycling': 'Recycling'}
        Medias = ['W', 'L', 'A', 'S']
        Activities = ['Treatment', 'Energy recovery', 'Recycling']
        # Inflation rate in the U.S. between 2005 and 2020 is 35.14%
        for Activity in Activities:
            for Media in Medias:
                df_PACE_aux = df_PACE[['NAICS code', 'Total PACE',
                                       'RSE for total PACE',
                                       'Establishments (employees >= 20)',
                                       'Info establishments',
                                       'Total shipment']]
                df_PACE_aux['Total PACE'] = \
                    df_PACE_aux['Total PACE'] * 1.3514 * 10**6
                df_PACE_aux['Media'] = Media
                df_PACE_aux['Activity'] = Activity
                df_PACE_aux['P-media'] = df_PACE[
                    'Media - ' + Dictionary_relation[Media].lower()]
                df_PACE_aux['P-activity'] = df_PACE[
                    'Activity - ' + Dictionary_relation[Activity].lower()]
                df_PACE_for_merging = pd.concat(
                    [df_PACE_for_merging, df_PACE_aux], ignore_index=True,
                    sort=True, axis=0)
                df_PAOC_aux = df_PAOC[['NAICS code', 'Total PAOC',
                                       'RSE for total PAOC',
                                       'Establishments (employees >= 20)',
                                       'Info establishments',
                                       'Total shipment']]
                df_PAOC_aux['Total PAOC'] = \
                    df_PAOC_aux['Total PAOC'] * 1.3514 * 10**6
                df_PAOC_aux['Media'] = Media
                df_PAOC_aux['Activity'] = Activity
                df_PAOC_aux['P-media'] = df_PAOC[
                    'Media - ' + Dictionary_relation[Media].lower()]
                df_PAOC_aux['P-activity'] = df_PAOC[
                    'Activity - ' + Dictionary_relation[Activity].lower()]
                df_PAOC_for_merging = pd.concat(
                    [df_PAOC_for_merging, df_PAOC_aux], ignore_index=True,
                    sort=True, axis=0)

        # Identifying probable establishment based on the
        # pobability of media and activity
        df_PAOC_for_merging = df_PAOC_for_merging.loc[pd.notnull(
            df_PAOC_for_merging[['P-media', 'P-activity']].all(axis=1))]
        df_PAOC_for_merging[[
            'Probable establishments by activity & media',
            'Info probable establishments']] = df_PAOC_for_merging[[
                'Info establishments', 'P-media', 'P-activity']].apply(
                    lambda x: selecting_establishment_by_activity_and_media(
                        x.values[0], x.values[1], x.values[2]), axis=1)
        df_PAOC_for_merging.drop(
            columns=['Info establishments',
                     'Establishments (employees >= 20)'],
            inplace=True)
        df_PACE_for_merging = df_PACE_for_merging.loc[
            pd.notnull(df_PACE_for_merging[['P-media', 'P-activity']].all(
                axis=1))]
        df_PACE_for_merging[[
            'Probable establishments by activity & media',
            'Info probable establishments']] = df_PACE_for_merging[[
                'Info establishments', 'P-media', 'P-activity']].apply(
                    lambda x: selecting_establishment_by_activity_and_media(
                        x.values[0], x.values[1], x.values[2]), axis=1)
        df_PACE_for_merging.drop(columns=['Info establishments',
                                          'Establishments (employees >= 20)'],
                                 inplace=True)

        # Joining census with TRI
        df_PACE = pd.merge(df_PACE_for_merging, df_PAU,
                           on=['NAICS code', 'Media', 'Activity'],
                           how='inner')
        df_PAOC = pd.merge(df_PAOC_for_merging, df_PAU,
                           on=['NAICS code', 'Media', 'Activity'],
                           how='inner')
        idx = df_PACE.loc[(df_PACE[['P-media', 'P-activity']].isnull()).all(
            axis=1)].index.tolist()

        df_PACE[['RSE for total PACE',
                 'Probable establishments by activity & media',
                 'P-media', 'P-activity', 'Total PACE', 'Total shipment',
                 'Info probable establishments']].iloc[idx] = df_PACE[[
                     'NAICS code', 'Media', 'Activity']].iloc[idx].apply(
                         lambda x: searching_census(
                             x.values[0], x.values[1], x.values[2],
                             df_PACE_for_merging), axis=1)

        df_PACE = df_PACE.loc[pd.notnull(df_PACE).all(axis=1)]
        idx = df_PAOC.loc[(df_PAOC[['P-media', 'P-activity']].isnull()).all(
            axis=1)].index.tolist()

        df_PAOC[['RSE for total PAOC',
                 'Probable establishments by activity & media', 'P-media',
                 'P-activity', 'Total PAOC', 'Total shipment',
                 'Info probable establishments']].iloc[idx] = df_PAOC[[
                     'NAICS code', 'Media', 'Activity']].iloc[idx].apply(
                         lambda x: searching_census(
                             x.values[0], x.values[1], x.values[2],
                             df_PAOC_for_merging), axis=1)
        df_PAOC = df_PAOC.loc[pd.notnull(df_PAOC).all(axis=1)]
        # Calculating total by media and activity
        df_PACE = df_PACE.groupby('NAICS code', as_index=False).apply(
            lambda x: normalizing_shipments(x))
        df_PAOC = df_PAOC.groupby('NAICS code', as_index=False).apply(
            lambda x: normalizing_shipments(x))
        # Calculating mass to activity and media
        # assuming lognormal distribution
        df_PAOC = df_PAOC[
            (df_PAOC[['P-media', 'P-activity']] != 0.0).all(axis=1)]
        df_PAOC[['Probable mass by activity & media',
                 'Total shipment by activity & media',
                 'Shipment/mass']] = df_PAOC[[
                     'mu', 'theta_2', 'Info probable establishments',
                     'P-media', 'P-activity']].apply(
                         lambda x: estimating_mass_by_activity_and_media(
                             x.values[0], x.values[1], x.values[2],
                             x.values[3], x.values[4]), axis=1)

        df_PACE = df_PACE[
            (df_PACE[['P-media', 'P-activity']] != 0.0).all(axis=1)]
        df_PACE[['Probable mass by activity & media',
                 'Total shipment by activity & media',
                 'Shipment/mass']] = df_PACE[[
                     'mu', 'theta_2', 'Info probable establishments',
                     'P-media', 'P-activity']].apply(
                         lambda x: estimating_mass_by_activity_and_media(
                             x.values[0], x.values[1], x.values[2],
                             x.values[3], x.values[4]), axis=1)

        # Assuming a normal distribution and a confidence level of 95%
        Z = norm.ppf(0.975)
        df_PAOC[['Total activity & media', 'Mean PAOC', 'SD PAOC',
                 'CI at 95% for Mean PAOC']] = df_PAOC.apply(
                     lambda x: mean_standard(
                         x['Probable establishments by activity & media'],
                         x['Shipment/mass'], x['RSE for total PAOC'],
                         x['Total PAOC'], x['Total shipment'], Z), axis=1)
        df_PAOC['Unit'] = 'USD/kg'
        df_PAOC = df_PAOC.loc[pd.notnull(df_PAOC).all(axis=1)]
        df_PAOC = df_PAOC.round(6)
        df_PAOC = df_PAOC.loc[df_PAOC['Mean PAOC'] != 0]
        df_PACE[['Total activity & media', 'Mean PACE',
                 'SD PACE', 'CI at 95% for Mean PACE']] = df_PACE.apply(
                     lambda x: mean_standard(
                         x['Probable establishments by activity & media'],
                         x['Shipment/mass'], x['RSE for total PACE'],
                         x['Total PACE'], x['Total shipment'], Z), axis=1)
        df_PACE['Unit'] = 'USD/kg'
        df_PACE = df_PACE.loc[pd.notnull(df_PACE).all(axis=1)]
        df_PACE = df_PACE.round(6)
        df_PACE = df_PACE.loc[df_PACE['Mean PACE'] != 0]
        #  Saving
        cols = ['NAICS code', 'Activity', 'Media',
                'Probable establishments by activity & media',
                'Probable mass by activity & media',
                'Total activity & media',
                'Total shipment by activity & media',
                'Shipment/mass', 'Mean PAOC',
                'SD PAOC', 'CI at 95% for Mean PAOC', 'Unit']
        df_PAOC = df_PAOC[cols]
        df_PAOC.to_csv(
            self._dir_path +
            '/datasets/pau_expenditure_and_cost/PAOC.csv', sep=',',
            index=False)
        cols = ['NAICS code', 'Activity', 'Media',
                'Probable establishments by activity & media',
                'Probable mass by activity & media',
                'Total activity & media',
                'Total shipment by activity & media',
                'Shipment/mass', 'Mean PACE',
                'SD PACE', 'CI at 95% for Mean PACE', 'Unit']
        df_PACE = df_PACE[cols]
        df_PACE.to_csv(
            self._dir_path + '/datasets/pau_expenditure_and_cost/PACE.csv',
            sep=',', index=False)

    def Pollution_control_unit_position(self):
        """TODO describe this method."""
        df_tri = pd.DataFrame()
        for Year in range(1987, 2005):
            df_tri_aux = pd.read_csv(
                self._dir_path +
                f'/datasets/intermediate_pau_datasets/PAUs_DB_{Year}.csv',
                usecols=['REPORTING YEAR', 'TRIFID',
                         'METHOD CODE - 2004 AND PRIOR'])
            df_tri_aux.drop_duplicates(keep='first', inplace=True)
            df_tri_aux.drop(columns=['REPORTING YEAR'], inplace=True)
            df_tri = pd.concat([df_tri, df_tri_aux], ignore_index=True,
                               sort=True, axis=0)
            del df_tri_aux
        df_tri.drop_duplicates(keep='first', inplace=True)
        df_tri.drop(columns=['TRIFID'], inplace=True)
        df_tri = df_tri[
            df_tri['METHOD CODE - 2004 AND PRIOR'].str.contains('+')]
        List_first = list()
        List_second = list()
        for idx, row in df_tri.iterrows():
            List_PAUs = row['METHOD CODE - 2004 AND PRIOR'].split(' + ')
            n = len(List_PAUs) - 1
            count = 0
            while count < n:
                List_first.append(List_PAUs[count])
                List_second.append(List_PAUs[count + 1])
                count = count + 1
        df_position = pd.DataFrame({'First': List_first,
                                    'Second': List_second})
        df_position.drop_duplicates(keep='first', inplace=True)
        df_position.to_csv(
            self._dir_path + '/datasets/pau_positions/Positions.csv',
            sep=',', index=False)

    def Searching_information_for_years_after_2004(self):
        """TODO describe this method."""
        df_tri_older = pd.DataFrame()
        for year in range(1987, 2005):
            df_tri_older_aux = pd.read_csv(
                f'{self._dir_path }/datasets/final_pau_datasets' +
                f'/PAUs_DB_filled_{year}.csv',
                usecols=['TRIFID', 'CAS NUMBER', 'WASTE STREAM CODE',
                         'METHOD CODE - 2004 AND PRIOR',
                         'METHOD NAME - 2004 AND PRIOR',
                         'METHOD CODE - 2005 AND AFTER',
                         'RANGE INFLUENT CONCENTRATION',
                         'EFFICIENCY ESTIMATION'], low_memory=False)
            df_tri_older = pd.concat(
                [df_tri_older, df_tri_older_aux], ignore_index=True,
                sort=True, axis=0)
            del df_tri_older_aux

        df_tri_older.drop_duplicates(keep='first', inplace=True)
        df_PAU = pd.read_csv(
            f'{self._dir_path }/datasets/final_pau_datasets' +
            f'/PAUs_DB_filled_{self.Year}.csv', low_memory=False)
        columns = df_PAU.columns.tolist()
        df_PAU.drop_duplicates(keep='first', inplace=True)
        df_PAU = pd.merge(df_PAU, df_tri_older, how='left',
                          on=['TRIFID',
                              'CAS NUMBER',
                              'WASTE STREAM CODE',
                              'METHOD CODE - 2005 AND AFTER'])
        del df_tri_older

        print(df_PAU.info())
        df_PAU.drop_duplicates(keep='first', subset=columns, inplace=True)
        df_PAU['EQUAL RANGE'] = df_PAU[
            ['EFFICIENCY RANGE CODE', 'EFFICIENCY ESTIMATION']].apply(
                lambda x: x.values[0] == self._efficiency_estimation_to_range(
                    x.values[1]) if x.values[1] else True, axis=1)

        idx = df_PAU.loc[((not df_PAU['EQUAL RANGE']) & (
            pd.notnull(df_PAU['EFFICIENCY ESTIMATION'])))].index.tolist()
        df_PAU.drop(columns=['EQUAL RANGE'], inplace=True)
        df_PAU_aux = df_PAU.loc[idx]
        df_PAU_aux.drop(columns=['METHOD CODE - 2004 AND PRIOR',
                                 'METHOD NAME - 2004 AND PRIOR',
                                 'RANGE INFLUENT CONCENTRATION',
                                 'EFFICIENCY ESTIMATION'],
                        inplace=True)
        df_PAU.drop(idx, inplace=True)
        df_PAU = pd.concat([df_PAU, df_PAU_aux], ignore_index=True,
                           sort=True, axis=0)
        del df_PAU_aux
        df_PAU.to_csv(
            f'{self._dir_path }/datasets/final_pau_datasets/' +
            f'PAUs_DB_filled_{self.Year}.csv', sep=',', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS)
    parser.add_argument(
        'Option',
        help='What do you want to do:\
            [A]: Recover information from TRI.\
            [B]: File for statistics. \
            [C]: File for recycling. \
            [D]: Further cleaning of database. \
            [E]: Organizing file with flows (1987-2004). \
            [F]: Organizing file with substance prices (1987 - 2004). \
            [G]: Pollution abatement cost and expenditure (only 2004). \
            [H]: Pollution control unit positions (1987 - 2004)\
            [I]: Searching information for years after 2004',
        type=str)

    parser.add_argument(
        '-Y', '--Year', nargs='+',
        help='Records with up to how many PAUs you want to include?.',
        type=str, required=False, default=[2018])

    parser.add_argument(
        '-N_Bins', help='Number of bins to split the middle waste flow values',
        type=int, default=10, required=False)

    args = parser.parse_args()
    start_time = time.time()

    for Year in args.Year:
        Building = PAU_DB(int(Year))
        if args.Option == 'A':
            Building.organizing()
        elif args.Option == 'B':
            Building.Building_database_for_statistics()
        elif args.Option == 'C':
            Building.Building_database_for_recycling_efficiency()
        elif args.Option == 'D':
            Building.cleaning_database()
        elif args.Option == 'E':
            Building.Building_database_for_flows(args.N_Bins)
        elif args.Option == 'F':
            Building.Organizing_substance_prices()
        elif args.Option == 'G':
            Building.pollution_abatement_cost_and_expenditure()
        elif args.Option == 'H':
            Building.Pollution_control_unit_position()
        elif args.Option == 'I':
            Building.Searching_information_for_years_after_2004()

    print('Execution time: %s sec' % (time.time() - start_time))
