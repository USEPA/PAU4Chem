# -*- coding: utf-8 -*-
#!/usr/bin/env python

# Note:
# 1. Range of Influent Concentration was reported from 1987 through 2004
# 2. Treatment Efficiency Estimation was reported from 1987 through 2004

import warnings
warnings.simplefilter(action = 'ignore', category = FutureWarning)
from scipy.stats import norm
from scipy.stats import lognorm
import pandas as pd
pd.options.mode.chained_assignment = None
import os
import argparse
import numpy as np
import re
import time
import unicodedata
from itertools import combinations
import yaml
import math
from Population import *

class PCU_DB:

    def __init__(self, Year):
        self._dir_path = os.path.dirname(os.path.realpath(__file__)) # Working Directory
        self.Year = Year
        #self._dir_path = os.getcwd() # if you are working on Jupyter Notebook

    def calling_TRI_Files(self):
        TRI_Files = dict()
        for file in ['1a', '1b', '2b']:
            columns = pd.read_csv(self._dir_path + '/Ancillary/TRI_File_' + file + '_needed_columns.txt',
                                header = None)
            columns =  list(columns.iloc[:,0])
            df = pd.read_csv(self._dir_path + '/Ancillary/US_' + file + '_' + str(self.Year) + '.csv',
                            usecols = columns,
                            low_memory = False)
            df = df.where(pd.notnull(df), None)
            TRI_Files.update({file: df})
        return TRI_Files


    def is_number(self, s):
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


    def _efficiency_estimation_empties_based_on_EPA_regulation(self, classification, HAP, RCRA):
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
        Acronyms = ['TRI', 'CAA', 'RCRA_F', 'RCRA_K', 'RCRA_P', 'RCRA_T', 'RCRA_U']
        Files = {Acronym: File for File in os.listdir(self._dir_path + '/SRS') for Acronym in Acronyms if Acronym in File}
        columns = ['ID', 'Internal Tracking Number']
        df = pd.read_csv(self._dir_path + '/SRS/' + Files['TRI'], low_memory = False,
                        usecols = ['ID', 'Internal Tracking Number'],
                        converters = {'ID': lambda x: str(int(x)) if re.search('^\d', x) else x},
                        dtype = {'Internal Tracking Number': 'int64'})
        df = df.assign(HAP = ['NO']*df.shape[0], RCRA = ['NO']*df.shape[0])
        del Files['TRI']
        for Acronym, File in Files.items():
            col = 'HAP'
            if Acronym in Acronyms[2:]:
                col = 'RCRA'
            ITN = pd.read_csv(self._dir_path + '/SRS/' + File,
                            low_memory = False,
                            usecols = ['Internal Tracking Number'],
                            dtype = {'Internal Tracking Number': 'int64'})
            df.loc[df['Internal Tracking Number'].isin( ITN['Internal Tracking Number'].tolist()), col] = 'YES'
        df.drop(columns = 'Internal Tracking Number', inplace = True)
        df.rename(columns = {'ID': 'CAS NUMBER'}, inplace = True)
        return df


    def _changin_management_code_for_2004_and_prior(self, x, m_n):
        Change = pd.read_csv(self._dir_path + '/Ancillary/Methods_TRI.csv',
                        usecols = ['Code 2004 and prior', 'Code 2005 and after'],
                        low_memory = False)
        if list(x.values).count(None) != m_n:
            y = {v:'T' if v in Change['Code 2004 and prior'].unique().tolist() else 'F' for v in x.values}
            result = [Change.loc[Change['Code 2004 and prior'] == v, 'Code 2005 and after'].tolist()[0] \
                    if s == 'T' else None for v, s in y.items()]
            L = len(result)
            result = result + [None]*(m_n - L)
            return result
        else:
            return [None]*m_n


    def organizing(self):
        dfs = self.calling_TRI_Files()
        df = dfs['2b'].where(pd.notnull(dfs['2b']), None)
        if self.Year >= 2005:
            df.drop(columns = df.iloc[:, list(range(18, 71, 13))].columns.tolist(), inplace = True)
        else:
            df.drop(columns = df.iloc[:, list(range(20, 73, 13))].columns.tolist(), inplace = True)
        df_PCUs = pd.DataFrame()
        Columns_0 = list(df.iloc[:, 0:8].columns)
        for i in range(5):
            Starting = 8 + 12*i
            Ending = Starting + 11
            Columns_1 = list(df.iloc[:, Starting:Ending + 1].columns)
            Columns = Columns_0 + Columns_1
            df_aux = df[Columns]
            Columns_to_change = {col: re.sub(r'STREAM [1-5] - ', '', col) for col in Columns_1}
            df_aux.rename(columns = Columns_to_change, inplace =  True)
            df_PCUs = pd.concat([df_PCUs, df_aux], ignore_index = True,
                                       sort = True, axis = 0)
            del Columns
        del df, df_aux
        cols =  list(df_PCUs.iloc[:, 9:17].columns)
        df_PCUs.dropna(subset = cols, how = 'all', axis = 0, inplace = True)
        if self.Year <= 2004:
            df_PCUs.dropna(subset = ['WASTE STREAM CODE', 'RANGE INFLUENT CONCENTRATION', \
                        'TREATMENT EFFICIENCY ESTIMATION'], how = 'any', axis = 0, inplace = True)
            df_PCUs.reset_index(inplace = True, drop = True)
            df_PCUs['METHOD CODE - 2004 AND PRIOR'] = df_PCUs[cols].apply(lambda x: None if  list(x).count(None) == len(cols) else ' + '.join(xx for xx in x if xx), axis = 1)
            df_PCUs[cols] = df_PCUs.apply(lambda row: pd.Series(self._changin_management_code_for_2004_and_prior(row[cols], len(cols))),
                                    axis =  1)
            df_PCUs = df_PCUs.loc[pd.notnull(df_PCUs[cols]).any(axis = 1)]
            df_PCUs['EFFICIENCY RANGE CODE'] = df_PCUs['TREATMENT EFFICIENCY ESTIMATION']\
                                      .apply(lambda x: self._efficiency_estimation_to_range(float(x)))
            df_PCUs.rename(columns = {'TREATMENT EFFICIENCY ESTIMATION': 'EFFICIENCY ESTIMATION'}, inplace = True)
            mask = pd.to_numeric(df_PCUs['RANGE INFLUENT CONCENTRATION'], errors='coerce').notnull()
            df_PCUs = df_PCUs[mask]
            df_PCUs['RANGE INFLUENT CONCENTRATION'] = df_PCUs['RANGE INFLUENT CONCENTRATION'].apply(lambda x: abs(int(x)))
        else:
            df_PCUs.rename(columns = {'TREATMENT EFFICIENCY RANGE CODE': 'EFFICIENCY RANGE CODE'}, inplace = True)
            df_PCUs.dropna(subset = ['WASTE STREAM CODE', 'EFFICIENCY RANGE CODE'],
                            how = 'any', axis = 0, inplace = True)
        df_PCUs['METHOD CODE - 2005 AND AFTER'] = df_PCUs[cols].apply(lambda x: None if  list(x).count(None) == len(cols) else ' + '.join(xx for xx in x if xx), axis = 1)
        df_PCUs = df_PCUs.loc[pd.notnull(df_PCUs['METHOD CODE - 2005 AND AFTER'])]
        df_PCUs['TYPE OF MANAGEMENT'] = 'Treatment'
        df_PCUs.drop(columns = cols, inplace = True)
        df_PCUs.reset_index(inplace =  True, drop =  True)
        df_PCUs.loc[pd.isnull(df_PCUs['BASED ON OPERATING DATA?']), 'BASED ON OPERATING DATA?'] = 'NO'
        try:
            # On-site energy recovery
            df = dfs['1a'].iloc[:, list(range(12))]
            cols = [c for c in df.columns if 'METHOD' in c]
            df.dropna(subset = cols, how = 'all', axis = 0, inplace = True)
            Columns_0 = list(df.iloc[:, 0:8].columns)
            Columns_1 = list(df.iloc[:, 8:].columns)
            dfs_energy = pd.DataFrame()
            for col in Columns_1:
                Columns = Columns_0 + [col]
                df_aux = df[Columns]
                df_aux.rename(columns = {col: re.sub(r' [1-4]', '', col)},
                                inplace =  True)
                dfs_energy = pd.concat([dfs_energy, df_aux], ignore_index = True,
                                           sort = True, axis = 0)
                del Columns
            del df, df_aux
            dfs_energy = dfs_energy.loc[pd.notnull(dfs_energy['ON-SITE ENERGY RECOVERY METHOD'])]
            dfs_energy['TYPE OF MANAGEMENT'] = 'Energy recovery'
            if self.Year <= 2004:
                dfs_energy['METHOD CODE - 2004 AND PRIOR'] = dfs_energy['ON-SITE ENERGY RECOVERY METHOD']
                dfs_energy['ON-SITE ENERGY RECOVERY METHOD'] = dfs_energy.apply(lambda row: \
                            pd.Series(self._changin_management_code_for_2004_and_prior(pd.Series(row['ON-SITE ENERGY RECOVERY METHOD']), 1)),
                                            axis =  1)
                dfs_energy = dfs_energy.loc[pd.notnull(dfs_energy['ON-SITE ENERGY RECOVERY METHOD'])]
            dfs_energy.rename(columns = {'ON-SITE ENERGY RECOVERY METHOD': 'METHOD CODE - 2005 AND AFTER'},
                            inplace =  True)
            dfs_energy = dfs_energy.loc[pd.notnull(dfs_energy['METHOD CODE - 2005 AND AFTER'])]
            df_PCUs = pd.concat([df_PCUs, dfs_energy], ignore_index = True,
                                   sort = True, axis = 0)
            del dfs_energy
        except ValueError as e:
            print('{}:\nThere is not information about energy recovery activities'.format(e))
        try:
            # On-site recycling
            df = dfs['1a'].iloc[:, list(range(8)) + list(range(12,19))]
            cols = [c for c in df.columns if 'METHOD' in c]
            df.dropna(subset = cols, how = 'all', axis = 0, inplace = True)
            Columns_0 = list(df.iloc[:, 0:8].columns)
            Columns_1 = list(df.iloc[:, 8:].columns)
            dfs_recycling = pd.DataFrame()
            for col in Columns_1:
                Columns = Columns_0 + [col]
                df_aux = df[Columns]
                df_aux.rename(columns = {col: re.sub(r' [1-7]', '', col)},
                                inplace =  True)
                dfs_recycling = pd.concat([dfs_recycling, df_aux], ignore_index = True,
                                           sort = True, axis = 0)
                del Columns
            del df, df_aux
            dfs_recycling = dfs_recycling.loc[pd.notnull(dfs_recycling['ON-SITE RECYCLING PROCESSES METHOD'])]
            dfs_recycling['TYPE OF MANAGEMENT'] = 'Recycling'
            if self.Year <= 2004:
                dfs_recycling['METHOD CODE - 2004 AND PRIOR'] = dfs_recycling['ON-SITE RECYCLING PROCESSES METHOD']
                dfs_recycling['ON-SITE RECYCLING PROCESSES METHOD'] = dfs_recycling.apply(lambda row: \
                            pd.Series(self._changin_management_code_for_2004_and_prior(pd.Series(row['ON-SITE RECYCLING PROCESSES METHOD']), 1)),
                                            axis =  1)
                dfs_recycling = dfs_recycling.loc[pd.notnull(dfs_recycling['ON-SITE RECYCLING PROCESSES METHOD'])]
            dfs_recycling.rename(columns = {'ON-SITE RECYCLING PROCESSES METHOD': 'METHOD CODE - 2005 AND AFTER'},
                            inplace =  True)
            dfs_recycling = dfs_recycling.loc[pd.notnull(dfs_recycling['METHOD CODE - 2005 AND AFTER'])]
            df_PCUs = pd.concat([df_PCUs, dfs_recycling], ignore_index = True,
                                   sort = True, axis = 0)
            del dfs_recycling
        except ValueError as e:
            print('{}:\nThere is not information about recycling activities'.format(e))
        # Changing units
        df_PCUs = df_PCUs.loc[(df_PCUs.iloc[:,0:] != 'INV').all(axis = 1)]
        df_PCUs.dropna(how = 'all', axis = 0, inplace = True)
        if self.Year >= 2005:
            Change = pd.read_csv(self._dir_path + '/Ancillary/Methods_TRI.csv',
                            usecols = ['Code 2004 and prior', 'Code 2005 and after'],
                            low_memory = False)
            Codes_2004 = Change.loc[(pd.notnull(Change['Code 2004 and prior'])) \
                        & (Change['Code 2005 and after'] != Change['Code 2004 and prior']),\
                        'Code 2004 and prior'].unique().tolist()
            idx = df_PCUs.loc[df_PCUs['METHOD CODE - 2005 AND AFTER'].isin(Codes_2004)].index.tolist()
            del Change, Codes_2004
            if len(idx) != 0:
                df_PCUs.loc[idx, 'METHOD CODE - 2005 AND AFTER'] = \
                df_PCUs.loc[idx]\
                            .apply(lambda row: self._changin_management_code_for_2004_and_prior(\
                            pd.Series(row['METHOD CODE - 2005 AND AFTER']),\
                            1),
                            axis = 1)
        # Adding methods name
        Methods = pd.read_csv(self._dir_path + '/Ancillary/Methods_TRI.csv',
                            usecols = ['Code 2004 and prior',
                                    'Method 2004 and prior',
                                    'Code 2005 and after',
                                    'Method 2005 and after'])
        Methods.drop_duplicates(keep =  'first', inplace = True)
        # Adding chemical activities and uses
        df_PCUs['DOCUMENT CONTROL NUMBER'] = df_PCUs['DOCUMENT CONTROL NUMBER'].apply(lambda x: str(int(float(x))) if self.is_number(x) else x)
        dfs['1b'].drop_duplicates(keep = 'first', inplace = True)
        dfs['1b']['DOCUMENT CONTROL NUMBER'] = dfs['1b']['DOCUMENT CONTROL NUMBER'].apply(lambda x: str(int(float(x))) if self.is_number(x) else x)
        df_PCUs = pd.merge(df_PCUs, dfs['1b'], on = ['TRIFID', 'DOCUMENT CONTROL NUMBER', 'CAS NUMBER'],
                                               how = 'inner')
        columns_DB_F = ['REPORTING YEAR', 'TRIFID', 'PRIMARY NAICS CODE', 'CAS NUMBER',
                         'CHEMICAL NAME', 'METAL INDICATOR', 'CLASSIFICATION',
                         'PRODUCE THE CHEMICAL', 'IMPORT THE CHEMICAL',
                         'ON-SITE USE OF THE CHEMICAL','SALE OR DISTRIBUTION OF THE CHEMICAL',
                         'AS A BYPRODUCT', 'AS A MANUFACTURED IMPURITY', 'USED AS A REACTANT',
                         'ADDED AS A FORMULATION COMPONENT', 'USED AS AN ARTICLE COMPONENT',
                         'REPACKAGING', 'AS A PROCESS IMPURITY', 'RECYCLING',
                         'USED AS A CHEMICAL PROCESSING AID', 'USED AS A MANUFACTURING AID',
                         'ANCILLARY OR OTHER USE',
                         'WASTE STREAM CODE', 'METHOD CODE - 2005 AND AFTER',
                         'METHOD NAME - 2005 AND AFTER', 'TYPE OF MANAGEMENT',
                         'EFFICIENCY RANGE CODE', 'BASED ON OPERATING DATA?']
        if self.Year <= 2004:
            Method = {row.iloc[0]: row.iloc[1] for index, row in Methods[['Code 2004 and prior', 'Method 2004 and prior']].iterrows()}
            def _checking(x, M):
                if x:
                    return ' + '.join(M[xx] for xx in x.split(' + ') if xx and xx and xx in M.keys())
                else:
                    return None
            df_PCUs = df_PCUs.loc[df_PCUs['METHOD CODE - 2004 AND PRIOR'].str.contains(r'[A-Z]').where(df_PCUs['METHOD CODE - 2004 AND PRIOR'].str.contains(r'[A-Z]'), False)]
            df_PCUs['METHOD NAME - 2004 AND PRIOR'] = df_PCUs['METHOD CODE - 2004 AND PRIOR'].apply(lambda x: _checking(x, Method))
            df_PCUs = df_PCUs.loc[(df_PCUs['METHOD CODE - 2004 AND PRIOR'] != '') | (pd.notnull(df_PCUs['METHOD CODE - 2004 AND PRIOR']))]
            columns_DB_F = ['REPORTING YEAR', 'TRIFID', 'PRIMARY NAICS CODE', 'CAS NUMBER',
                             'CHEMICAL NAME', 'METAL INDICATOR', 'CLASSIFICATION',
                             'PRODUCE THE CHEMICAL', 'IMPORT THE CHEMICAL', 'ON-SITE USE OF THE CHEMICAL',
                             'SALE OR DISTRIBUTION OF THE CHEMICAL', 'AS A BYPRODUCT',
                             'AS A MANUFACTURED IMPURITY', 'USED AS A REACTANT',
                             'ADDED AS A FORMULATION COMPONENT', 'USED AS AN ARTICLE COMPONENT',
                             'REPACKAGING', 'AS A PROCESS IMPURITY', 'RECYCLING',
                             'USED AS A CHEMICAL PROCESSING AID', 'USED AS A MANUFACTURING AID',
                             'ANCILLARY OR OTHER USE',
                             'WASTE STREAM CODE', 'RANGE INFLUENT CONCENTRATION',
                             'METHOD CODE - 2004 AND PRIOR', 'METHOD NAME - 2004 AND PRIOR',
                             'METHOD CODE - 2005 AND AFTER', 'METHOD NAME - 2005 AND AFTER',
                             'TYPE OF MANAGEMENT', 'EFFICIENCY RANGE CODE', 'EFFICIENCY ESTIMATION',
                             'BASED ON OPERATING DATA?']
        Method = {row.iloc[0]: row.iloc[1] for index, row in Methods[['Code 2005 and after', 'Method 2005 and after']].iterrows()}
        df_PCUs = df_PCUs.loc[df_PCUs['METHOD CODE - 2005 AND AFTER'].str.contains(r'[A-Z]').where(df_PCUs['METHOD CODE - 2005 AND AFTER'].str.contains(r'[A-Z]'), False)]
        df_PCUs['METHOD NAME - 2005 AND AFTER'] = df_PCUs['METHOD CODE - 2005 AND AFTER'].apply(lambda x: ' + '.join(Method[xx] for xx in x.split(' + ') if xx and xx in Method.keys()))
        # Saving information
        df_PCUs['REPORTING YEAR'] = self.Year
        df_PCUs = df_PCUs[columns_DB_F]
        df_PCUs.to_csv(self._dir_path + '/Datasets/Intermediate_PCU_datasets/PCUs_DB_' + str(self.Year) + '.csv',
                     sep = ',', index = False)


    def Building_database_for_statistics(self):
        columns = pd.read_csv(self._dir_path + '/Ancillary/TRI_File_2b_needed_columns_for_statistics.txt',
                             header = None)
        columns =  list(columns.iloc[:,0])
        df = pd.read_csv(self._dir_path + '/Ancillary/US_2b_' + str(self.Year) + '.csv',
                        usecols = columns,
                        low_memory = False)
        df_statistics = pd.DataFrame()
        if self.Year >= 2005:
            df.drop(columns = df.iloc[:, list(range(12, 61, 12))].columns.tolist(), inplace = True)
            codes_incineration = ['A01', 'H040', 'H076', 'H122']
        else:
            df.drop(columns = df.iloc[:, list(range(13, 62, 12))].columns.tolist(), inplace = True)
            codes_incineration = ['A01', 'F01', 'F11', 'F19', 'F31',
                                'F41', 'F42', 'F51', 'F61',
                                'F71', 'F81', 'F82', 'F83',
                                'F99']
        Columns_0 = list(df.iloc[:, 0:2].columns)
        for i in range(5):
            Columns_1 = list(df.iloc[:, [2 + 11*i, 11 + 11*i, 12 + 11*i]].columns)
            Treatmes = list(df.iloc[:, 3 + 11*i: 11 + 11*i].columns)
            Columns = Columns_0 + Columns_1
            df_aux = df[Columns]
            df_aux['INCINERATION'] = 'NO'
            df_aux.loc[df[Treatmes].isin(codes_incineration).any(axis = 1), 'INCINERATION'] = 'YES'
            df_aux['IDEAL'] = df[Treatmes].apply(lambda x: 'YES' if  \
                                            len(list(np.where(pd.notnull(x))[0])) == 1  \
                                            else 'NO',
                                            axis = 1)
            Columns_to_change = {col: re.sub(r'STREAM [1-5] - ', '', col) for col in Columns_1}
            df_aux.rename(columns = Columns_to_change, inplace =  True)
            df_statistics = pd.concat([df_statistics, df_aux], ignore_index = True,
                                       sort = True, axis = 0)
            del Columns
        del df, df_aux
        if self.Year <= 2004:
            df_statistics.dropna(how = 'any', axis = 0, inplace = True)
            mask = pd.to_numeric(df_statistics['TREATMENT EFFICIENCY ESTIMATION'], errors='coerce').notnull()
            df_statistics = df_statistics[mask]
            df_statistics['TREATMENT EFFICIENCY ESTIMATION'] = df_statistics['TREATMENT EFFICIENCY ESTIMATION'].astype(float)
            df_statistics['EFFICIENCY RANGE'] = df_statistics['TREATMENT EFFICIENCY ESTIMATION']\
                                  .apply(lambda x: self._efficiency_estimation_to_range(float(x)))
            mask = pd.to_numeric(df_statistics['RANGE INFLUENT CONCENTRATION'], errors='coerce').notnull()
            df_statistics = df_statistics[mask]
            df_statistics['RANGE INFLUENT CONCENTRATION'] = df_statistics['RANGE INFLUENT CONCENTRATION'].astype(int)
            df_statistics.rename(columns = {'TREATMENT EFFICIENCY ESTIMATION': 'EFFICIENCY ESTIMATION'},
                                inplace = True)
        else:
            df_statistics.rename(columns = {'TREATMENT EFFICIENCY RANGE CODE': 'EFFICIENCY RANGE'},
                                inplace = True)
            df_statistics.dropna(subset = ['EFFICIENCY RANGE', 'WASTE STREAM CODE'], how = 'any', axis = 0, inplace = True)
        df_statistics.rename(columns = {'PRIMARY NAICS CODE': 'NAICS',
                                    'CAS NUMBER': 'CAS',
                                    'WASTE STREAM CODE': 'WASTE',
                                    'RANGE INFLUENT CONCENTRATION': 'CONCENTRATION'},
                            inplace = True)
        df_statistics.loc[df_statistics['INCINERATION'] == 'NO', 'IDEAL'] = None
        df_statistics.to_csv(self._dir_path + '/Statistics/DB_for_general/DB_for_Statistics_' + str(self.Year) + '.csv',
                     sep = ',', index = False)


    def Building_database_for_recycling_efficiency(self):
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
                    cases.append(row['ON-SITE - RECYCLED']/(row['ON-SITE - RECYCLED'] + sumatory)*100)
            try:
                if len(list(set(cases))) == 1 and cases[0] == 100:
                    return [100]*6 + [0] + [row['ON-SITE - RECYCLED']]
                else:
                    return [np.min(cases)] + np.quantile(cases, [0.25, 0.5, 0.75]).tolist() + \
                           [np.max(cases), np.mean(cases), np.std(cases)/np.mean(cases),\
                            row['ON-SITE - RECYCLED']]
            except ValueError:
                a = np.empty((8))
                a[:] = np.nan
                return a.tolist()
        columns = pd.read_csv(self._dir_path + '/Ancillary/TRI_File_1a_needed_columns_for_statistics.txt',
                             header = None)
        columns =  list(columns.iloc[:,0])
        df = pd.read_csv(self._dir_path + '/Ancillary/US_1a_' + str(self.Year) + '.csv',
                        usecols = columns,
                        low_memory = False)
        elements_total = list(set(df.iloc[:, 5:64].columns.tolist()) - set(['ON-SITE - RECYCLED']))
        df.iloc[:, 5:64] = df.iloc[:, 5:64].where(pd.notnull(df.iloc[:, 5:64]), 0.0)
        df.iloc[:, 5:64] = df.iloc[:, 5:64].apply(pd.to_numeric, errors='coerce')
        cols = [c for c in df.columns if 'METHOD' in c]
        df['IDEAL'] = df[cols].apply(lambda x: 'YES' if  \
                                         len(list(np.where(pd.notnull(x))[0])) == 1  \
                                         else 'NO',
                                         axis = 1)
        df = df.loc[df['IDEAL'] == 'YES']
        df['METHOD'] = df[cols].apply(lambda x: x.values[np.where(pd.notnull(x))[0]][0], axis = 1)
        df.drop(columns = ['IDEAL'] + cols, inplace = True)
        df = df.loc[df['METHOD'] != 'INV']
        df[['LOWER EFFICIENCY', 'Q1', 'Q2','Q3', 'UPPER EFFICIENCY',
            'MEAN OF EFFICIENCY', 'CV', 'ON-SITE - RECYCLED']]\
             = df.apply(lambda x: pd.Series(_division(x, elements_total)), axis = 1)
        df = df.loc[pd.notnull(df['UPPER EFFICIENCY'])]
        df['IQR'] = df.apply(lambda x: x['Q3'] - x['Q1'], axis = 1)
        df['Q1 - 1.5xIQR'] = df.apply(lambda x: 0 if x['Q1'] - 1.5*x['IQR'] < 0 \
                                    else x['Q1'] - 1.5*x['IQR'], axis = 1)
        df['Q3 + 1.5xIQR'] = df.apply(lambda x: 100 if x['Q3'] + 1.5*x['IQR'] > 100 \
                                    else x['Q3'] + 1.5*x['IQR'], axis = 1)
        df['UPPER EFFICIENCY OUTLIER?'] = df.apply(lambda x: 'YES' if x['UPPER EFFICIENCY'] > \
                                x['Q3 + 1.5xIQR'] else 'NO', axis = 1)
        df['LOWER EFFICIENCY OUTLIER?'] = df.apply(lambda x: 'YES' if x['LOWER EFFICIENCY'] < \
                                x['Q1 - 1.5xIQR'] else 'NO', axis = 1)
        df['HIGH VARIANCE?'] = df.apply(lambda x: 'YES' if x['CV'] > 1 else 'NO', axis = 1)
        df = df[['TRIFID', 'PRIMARY NAICS CODE', 'CAS NUMBER', 'ON-SITE - RECYCLED', 'UNIT OF MEASURE', \
                'LOWER EFFICIENCY', 'LOWER EFFICIENCY OUTLIER?', 'Q1 - 1.5xIQR', 'Q1', \
                'Q2', 'Q3', 'Q3 + 1.5xIQR', 'UPPER EFFICIENCY', 'UPPER EFFICIENCY OUTLIER?', \
                'IQR', 'MEAN OF EFFICIENCY', 'CV', 'HIGH VARIANCE?', 'METHOD']]
        df.iloc[:, [5, 7, 8, 9, 10, 11, 12, 14, 15, 16]] = \
                df.iloc[:, [5, 7, 8, 9, 10, 11, 12, 14, 15, 16]].round(4)
        df.to_csv(self._dir_path + '/Statistics/DB_for_solvents/DB_for_Solvents_' + str(self.Year) + '.csv',
                      sep = ',', index = False)



    def _searching_naics(self, x, naics):
        # https://www.census.gov/programs-surveys/economic-census/guidance/understanding-naics.html
        values = {0:'Nothing',
                 1:'Nothing',
                 2:'Sector',
                 3:'Subsector',
                 4:'Industry Group',
                 5:'NAICS Industry',
                 6:'National Industry'}
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
        if row['METHOD CODE - 2005 AND AFTER'] == 'H20': # Solvent recovery
            phases = ['L']
        elif row['METHOD CODE - 2005 AND AFTER'] == 'H39': # Acid regeneration and other reactions
            phases = ['W']
        elif row['METHOD CODE - 2005 AND AFTER'] == 'H10': # Metal recovery
            phases = ['W', 'S']
            if self.Year <= 2004:
                Pyrometallurgy = ['R27', 'R28', 'R29'] # They work with scrap
                if row['METHOD CODE - 2004 AND PRIOR'] in Pyrometallurgy:
                    Phases = ['S']
                else:
                    Phases = ['W', 'S']
        naics_structure = ['National Industry', 'NAICS Industry', 'Industry Group',
                    'Subsector', 'Sector', 'Nothing']
        df_cas = df_s.loc[df_s['CAS'] == row['CAS NUMBER'], ['NAICS', 'WASTE', 'VALUE']]
        df_cas = df_cas.groupby(['NAICS', 'WASTE'], as_index = False).sum()
        df_cas.reset_index(inplace = True)
        if (not df_cas.empty):
            df_cas['NAICS STRUCTURE'] = df_cas.apply(lambda x: \
                            self._searching_naics(x['NAICS'], \
                            row['PRIMARY NAICS CODE']), \
                            axis = 1)
            i = 0
            phase = None
            while i <= 5 and (not phase in phases):
                structure = naics_structure[i]
                i = i + 1
            #for structure in naics_structure:
                df_naics = df_cas.loc[df_cas['NAICS STRUCTURE'] == structure]
                if (df_naics.empty):
                    phase = None
                    #continue
                else:
                    if (df_naics['WASTE'].isin(phases).any()):
                        df_phase = df_naics.loc[df_naics['WASTE'].isin(phases)]
                        row['NAICS STRUCTURE'] = structure
                        row['WASTE STREAM CODE'] = df_phase.loc[df_phase['VALUE'].idxmax(), 'WASTE']
                    else:
                        row['NAICS STRUCTURE'] = structure
                        row['WASTE STREAM CODE'] = df_naics.loc[df_naics['VALUE'].idxmax(), 'WASTE']
                    phase =  row['WASTE STREAM CODE']
            return row
        else:
            row['NAICS STRUCTURE'] = None
            row['WASTE STREAM CODE'] = None
            return row


    def _concentration_estimation_recycling(self, df_s, cas, naics, phase, structure):
        df_s = df_s[['NAICS', 'CAS', 'WASTE', 'CONCENTRATION', 'VALUE']]
        df_s = df_s.loc[(df_s['CAS'] == cas) & \
                    (df_s['WASTE'] == phase)]
        df_s = df_s.groupby(['NAICS', 'CAS', 'WASTE', 'CONCENTRATION'], as_index = False).sum()
        df_s['NAICS STRUCTURE'] = df_s.apply(lambda x: \
                        self._searching_naics(x['NAICS'], \
                        naics), \
                        axis = 1)
        df = df_s.loc[(df_s['NAICS STRUCTURE'] == structure)]
        return df.loc[df['VALUE'].idxmax(), 'CONCENTRATION']


    def _recycling_efficiency(self, row, df_s):
        naics_structure = ['National Industry', 'NAICS Industry', 'Industry Group',
                    'Subsector', 'Sector', 'Nothing']
        if self.Year <= 2004:
            code = row['METHOD CODE - 2004 AND PRIOR']
        else:
            code = row['METHOD CODE - 2005 AND AFTER']
        df_cas = df_s.loc[(df_s['CAS NUMBER'] == row['CAS NUMBER']) & (df_s['METHOD'] == code)]
        if (not df_cas.empty):
            df_fid = df_cas.loc[df_cas['TRIFID'] == row['TRIFID']]
            if (not df_fid.empty):
                return df_fid['UPPER EFFICIENCY'].iloc[0]
            else:
                df_cas['NAICS STRUCTURE'] = df_cas.apply(lambda x: \
                                self._searching_naics(x['PRIMARY NAICS CODE'], \
                                                    row['PRIMARY NAICS CODE']), \
                                axis = 1)
                i = 0
                efficiency = None
                while (i <= 5) and (not efficiency):
                    structure = naics_structure[i]
                    i = i + 1
                    df_naics = df_cas.loc[df_cas['NAICS STRUCTURE'] == structure]
                    if df_naics.empty:
                        efficiency = None
                    else:
                        efficiency =  df_naics['UPPER EFFICIENCY'].median()
                return efficiency
        else:
            return None


    def _phase_estimation_energy(self, df_s, row):
        phases = ['S', 'L', 'A']
        if row['METHOD CODE - 2005 AND AFTER'] == 'U01':
            phases = ['S', 'L'] # Industrial Kilns (specially rotatory kilns) are used to burn hazardous liquid and solid wastes
        naics_structure = ['National Industry', 'NAICS Industry', 'Industry Group',
                    'Subsector', 'Sector', 'Nothing']
        df_cas = df_s.loc[df_s['CAS'] == row['CAS NUMBER'], ['NAICS', 'WASTE', 'VALUE', 'INCINERATION']]
        df_cas = df_cas.groupby(['NAICS', 'WASTE', 'INCINERATION'], as_index = False).sum()
        df_cas.reset_index(inplace = True)
        if (not df_cas.empty):
            df_cas['NAICS STRUCTURE'] = df_cas.apply(lambda x: \
                            self._searching_naics(x['NAICS'], \
                            row['PRIMARY NAICS CODE']), \
                            axis = 1)
            i = 0
            phase = None
            #for structure in naics_structure:
            while i <= 5 and (not phase in phases):
                structure = naics_structure[i]
                i = i + 1
                df_naics = df_cas.loc[df_cas['NAICS STRUCTURE'] == structure]
                if df_naics.empty:
                    phase = None
                else:
                    df_incineration = df_naics.loc[df_cas['INCINERATION'] == 'YES']
                    if df_incineration.empty:
                        if (df_naics['WASTE'].isin(phases).any()):
                            df_phase = df_naics.loc[df_naics['WASTE'].isin(phases)]
                            row['NAICS STRUCTURE'] = structure
                            row['WASTE STREAM CODE'] = df_phase.loc[df_phase['VALUE'].idxmax(), 'WASTE']
                        else:
                            row['NAICS STRUCTURE'] = structure
                            row['WASTE STREAM CODE'] = df_naics.loc[df_naics['VALUE'].idxmax(), 'WASTE']
                        row['BY MEANS OF INCINERATION'] = 'NO'
                    else:
                        if (df_incineration['WASTE'].isin(phases).any()):
                            df_phase = df_incineration.loc[df_incineration['WASTE'].isin(phases)]
                            row['NAICS STRUCTURE'] = structure
                            row['WASTE STREAM CODE'] = df_phase.loc[df_phase['VALUE'].idxmax(), 'WASTE']
                        else:
                            row['NAICS STRUCTURE'] = structure
                            row['WASTE STREAM CODE'] = df_incineration.loc[df_incineration['VALUE'].idxmax(), 'WASTE']
                        row['BY MEANS OF INCINERATION'] = 'YES'
                    phase =  row['WASTE STREAM CODE']
            return row
        else:
            row['NAICS STRUCTURE'] = None
            row['WASTE STREAM CODE'] = None
            row['BY MEANS OF INCINERATION'] = None
            return row


    def _concentration_estimation_energy(self, df_s, cas, naics, phase, structure, incineration):
        df_s = df_s[['NAICS', 'CAS', 'WASTE', 'CONCENTRATION', \
                    'VALUE', 'INCINERATION']]
        df_s = df_s.loc[(df_s['CAS'] == cas) & \
                    (df_s['WASTE'] == phase) & \
                    (df_s['INCINERATION'] == incineration)]
        df_s = df_s.groupby(['NAICS', 'CAS', 'WASTE', 'CONCENTRATION', 'INCINERATION'],
                    as_index = False).sum()
        df_s['NAICS STRUCTURE'] = df_s.apply(lambda x: \
                        self._searching_naics(x['NAICS'], \
                        naics), \
                        axis = 1)
        df = df_s.loc[(df_s['NAICS STRUCTURE'] == structure)]
        return df.loc[df['VALUE'].idxmax(), 'CONCENTRATION']


    def _energy_efficiency(self, df_s, row):
        if self.Year <= 2004:
            df_s = df_s[['NAICS', 'CAS', 'WASTE', 'INCINERATION', 'IDEAL', 'EFFICIENCY ESTIMATION']]
        else:
            df_s = df_s[['NAICS', 'CAS', 'WASTE', 'INCINERATION', 'IDEAL', 'EFFICIENCY RANGE', 'VALUE']]
            df_s['WASTE'] = df_s.groupby(['NAICS', 'CAS', 'WASTE', 'IDEAL', 'INCINERATION', 'EFFICIENCY RANGE'],
                    as_index = False).sum()
        df_s = df_s.loc[(df_s['CAS'] == row['CAS NUMBER']) & \
                        (df_s['INCINERATION'] == 'YES') & \
                        (df_s['IDEAL'] == 'YES')]
        if (not df_s.empty):
            df_s['NAICS STRUCTURE'] = df_s.apply(lambda x: \
                        self._searching_naics(x['NAICS'], \
                        row['PRIMARY NAICS CODE']), \
                        axis = 1)
            df_structure = df_s.loc[df_s['NAICS STRUCTURE'] == row['NAICS STRUCTURE']]
            if (not df_structure.empty):
                df_phase = df_structure.loc[df_structure['WASTE'] \
                        == row['WASTE STREAM CODE']]
                if (not df_phase.empty):
                    if self.Year <= 2004:
                        result =  df_phase['EFFICIENCY ESTIMATION'].median()
                    else:
                        result = df_phase.loc[df_phase['VALUE'].idxmax(), 'EFFICIENCY RANGE']
                else:
                    if self.Year <= 2004:
                        result =  df_structure['EFFICIENCY ESTIMATION'].median()
                    else:
                        result = df_structure.loc[df_structure['VALUE'].idxmax(), 'EFFICIENCY RANGE']
            else:
                df_phase = df_s.loc[df_s['WASTE'] \
                        == row['WASTE STREAM CODE']]
                if (not df_phase.empty):
                    if self.Year <= 2004:
                        result =  df_phase['EFFICIENCY ESTIMATION'].median()
                    else:
                        result = df_phase.loc[df_phase['VALUE'].idxmax(), 'EFFICIENCY RANGE']
                else:
                    if self.Year <= 2004:
                        result =  df_s['EFFICIENCY ESTIMATION'].median()
                    else:
                        result = df_s.loc[df_s['VALUE'].idxmax(), 'EFFICIENCY RANGE']
        else:
            return None
        return result


    def cleaning_database(self):
        # Calling TRI restriction for metals
        Restrictions = pd.read_csv(self._dir_path + '/Ancillary/Metals_divided_into_4_groups_can_be_reported.csv',
                                    low_memory = False,
                                    usecols = ['ID',
                                               "U01, U02, U03 (Energy recovery)",
                                               'H20 (Solvent recovey)'])
        Energy_recovery = Restrictions.loc[Restrictions["U01, U02, U03 (Energy recovery)"] == 'NO', 'ID'].tolist()
        Solvent_recovery = Restrictions.loc[Restrictions['H20 (Solvent recovey)'] == 'NO', 'ID'].tolist()
        # Calling PCU
        PCU = pd.read_csv(self._dir_path + '/Datasets/Intermediate_PCU_datasets/PCUs_DB_' + str(self.Year) + '.csv',
                            low_memory = False,
                            converters = {'CAS NUMBER': lambda x: x if re.search(r'^[A-Z]', x) else str(int(x))})
        columns_DB_F = PCU.columns.tolist()
        PCU['PRIMARY NAICS CODE'] = PCU['PRIMARY NAICS CODE'].astype('int')
        if self.Year <= 2004:
            grouping = ['TRIFID', 'METHOD CODE - 2004 AND PRIOR']
            PCU.sort_values(by = ['PRIMARY NAICS CODE', 'TRIFID',
                            'METHOD CODE - 2004 AND PRIOR', 'CAS NUMBER'],
                            inplace = True)
        else:
            grouping = ['TRIFID', 'METHOD CODE - 2005 AND AFTER']
            PCU.sort_values(by = ['PRIMARY NAICS CODE', 'TRIFID',
                            'METHOD CODE - 2005 AND AFTER', 'CAS NUMBER'],
                            inplace = True)
        # Calling database for statistics
        Statistics = pd.read_csv(self._dir_path + '/Statistics/DB_for_general/DB_for_Statistics_' + str(self.Year) + '.csv',
                                low_memory = False,
                                converters = {'CAS': lambda x: x if re.search(r'^[A-Z]', x) else str(int(x))})
        Statistics['NAICS'] = Statistics['NAICS'].astype('int')
        Statistics['VALUE'] = 1
        Statistics.sort_values(by = ['NAICS', 'CAS'], inplace = True)
        # Treatment
        Efficiency_codes = ['E1', 'E2', 'E3', 'E4', 'E5', 'E6']
        df_N_PCU = PCU.loc[PCU['TYPE OF MANAGEMENT'] == 'Treatment']
        df_N_PCU = df_N_PCU.loc[df_N_PCU['EFFICIENCY RANGE CODE'].isin(Efficiency_codes)]
        # Recycling
        PCU_recycling = PCU.loc[PCU['TYPE OF MANAGEMENT'] == 'Recycling']
        if not PCU_recycling.empty:
            PCU_recycling =  PCU_recycling.loc[~ ((PCU_recycling['METHOD CODE - 2005 AND AFTER'] == 'H20') & (PCU_recycling['CAS NUMBER'].isin(Solvent_recovery)))]
            PCU_recycling.reset_index(inplace = True, drop = True)
            PCU_recycling['BASED ON OPERATING DATA?'] = 'NO'
            # Calling database for recycling efficiency
            Recycling_statistics = pd.read_csv(self._dir_path + '/Statistics/DB_for_solvents/DB_for_Solvents_' + str(self.Year) +  '.csv',
                                    low_memory = False,
                                    usecols = ['TRIFID', 'PRIMARY NAICS CODE', 'CAS NUMBER', \
                                               'UPPER EFFICIENCY', 'UPPER EFFICIENCY OUTLIER?',
                                               'METHOD', 'HIGH VARIANCE?'],
                                    converters = {'CAS NUMBER': lambda x: x if re.search(r'^[A-Z]', x) else str(int(x))})
            Recycling_statistics['PRIMARY NAICS CODE'] = Recycling_statistics['PRIMARY NAICS CODE'].astype('int')
            Recycling_statistics = Recycling_statistics\
                                    .loc[(Recycling_statistics['UPPER EFFICIENCY OUTLIER?'] == 'NO') &
                                         (Recycling_statistics['HIGH VARIANCE?'] == 'NO')]
            Recycling_statistics.drop(columns = ['UPPER EFFICIENCY OUTLIER?', 'HIGH VARIANCE?'], axis = 1)
            efficiency_estimation = \
                    PCU_recycling.apply(lambda x: self._recycling_efficiency(x, Recycling_statistics), axis = 1).round(4)
            PCU_recycling['EFFICIENCY RANGE CODE'] = \
                            efficiency_estimation.apply(lambda x: self._efficiency_estimation_to_range(x))
            PCU_recycling = PCU_recycling.loc[pd.notnull(PCU_recycling['EFFICIENCY RANGE CODE'])]
            PCU_recycling = \
                     PCU_recycling.apply(lambda x: \
                     self._phase_estimation_recycling(Statistics, x), axis = 1)
            PCU_recycling = PCU_recycling.loc[pd.notnull(PCU_recycling['WASTE STREAM CODE'])]
            if self.Year <= 2004:
                PCU_recycling['EFFICIENCY ESTIMATION'] = efficiency_estimation
                PCU_recycling['RANGE INFLUENT CONCENTRATION'] = \
                          PCU_recycling.apply(lambda x: \
                         self._concentration_estimation_recycling(Statistics, \
                                         x['CAS NUMBER'], \
                                         x['PRIMARY NAICS CODE'],\
                                         x['WASTE STREAM CODE'], \
                                         x['NAICS STRUCTURE']), \
                                         axis = 1)
            PCU_recycling.drop(columns = ['NAICS STRUCTURE'], inplace = True)
            df_N_PCU = pd.concat([df_N_PCU, PCU_recycling],
                                     ignore_index = True,
                                     sort = True, axis = 0)
        else:
            pass
        # Energy recovery
        PCU_energy = PCU.loc[PCU['TYPE OF MANAGEMENT'] == 'Energy recovery']
        if not PCU_energy.empty:
            PCU_energy =  PCU_energy.loc[~ ((PCU_energy['METHOD CODE - 2005 AND AFTER'].isin(['U01', 'U02', 'U03'])) & (PCU_energy['CAS NUMBER'].isin(Energy_recovery)))]
            PCU_energy.reset_index(inplace = True, drop = True)
            PCU_energy['BASED ON OPERATING DATA?'] = 'NO'
            PCU_energy = \
                     PCU_energy.apply(lambda x: \
                     self._phase_estimation_energy(Statistics, x), axis = 1)
            PCU_energy = PCU_energy.loc[pd.notnull(PCU_energy['WASTE STREAM CODE'])]
            SRS = self._calling_SRS()
            if self.Year <= 2004:
                PCU_energy['RANGE INFLUENT CONCENTRATION'] = \
                         PCU_energy.apply(lambda x: \
                         self._concentration_estimation_energy(Statistics, \
                                         x['CAS NUMBER'], \
                                         x['PRIMARY NAICS CODE'],\
                                         x['WASTE STREAM CODE'], \
                                         x['NAICS STRUCTURE'], \
                                         x['BY MEANS OF INCINERATION']), \
                         axis = 1)
                PCU_energy.drop(columns = ['BY MEANS OF INCINERATION'], inplace = True)
                PCU_energy['EFFICIENCY ESTIMATION'] = \
                        PCU_energy.apply(lambda x: \
                        self._energy_efficiency(Statistics, x), axis = 1).round(4)
                PCU_energy = pd.merge(PCU_energy, SRS, on = 'CAS NUMBER', how = 'left')
                PCU_energy['EFFICIENCY ESTIMATION'] = PCU_energy.apply(lambda x: \
                                    self._efficiency_estimation_empties_based_on_EPA_regulation(\
                                    x['CLASSIFICATION'], x['HAP'], x['RCRA']) \
                                    if not x['EFFICIENCY ESTIMATION'] else
                                    x['EFFICIENCY ESTIMATION'],
                                    axis =  1)
                PCU_energy = PCU_energy.loc[pd.notnull(PCU_energy['EFFICIENCY ESTIMATION'])]
                PCU_energy['EFFICIENCY RANGE CODE'] = PCU_energy['EFFICIENCY ESTIMATION']\
                                          .apply(lambda x: self._efficiency_estimation_to_range(float(x)))
            else:
                PCU_energy.drop(columns = ['BY MEANS OF INCINERATION'], inplace = True)
                PCU_energy['EFFICIENCY RANGE CODE'] = \
                        PCU_energy.apply(lambda x: \
                        self._energy_efficiency(Statistics, x), axis = 1)
                PCU_energy = pd.merge(PCU_energy, SRS, on = 'CAS NUMBER', how = 'left')
                PCU_energy['EFFICIENCY RANGE CODE'] = PCU_energy.apply(lambda x: \
                                    self._efficiency_estimation_empties_based_on_EPA_regulation(\
                                    x['CLASSIFICATION'], x['HAP'], x['RCRA']) \
                                    if not x['EFFICIENCY RANGE CODE'] else
                                    x['EFFICIENCY RANGE CODE'],
                                    axis =  1)
                PCU_energy = PCU_energy.loc[pd.notnull(PCU_energy['EFFICIENCY RANGE CODE'])]
            PCU_energy.drop(columns = ['NAICS STRUCTURE', 'HAP', 'RCRA'], inplace = True)
            PCU_energy.loc[(PCU_energy['WASTE STREAM CODE'] == 'W') & \
                           (PCU_energy['TYPE OF MANAGEMENT'] == 'Energy recovery'),\
                            'WASTE STREAM CODE'] = 'L'
            df_N_PCU = pd.concat([df_N_PCU, PCU_energy],
                                     ignore_index = True,
                                     sort = True, axis = 0)
        else:
            pass
        Chemicals_to_remove = ['MIXTURE', 'TRD SECRT']
        df_N_PCU = df_N_PCU.loc[~df_N_PCU['CAS NUMBER'].isin(Chemicals_to_remove)]
        df_N_PCU['CAS NUMBER'] = df_N_PCU['CAS NUMBER'].apply(lambda x: str(int(x)) if not 'N' in x else x)
        df_N_PCU = df_N_PCU[columns_DB_F]
        df_N_PCU.to_csv(self._dir_path + '/Datasets/Final_PCU_datasets/PCUs_DB_filled_' + str(self.Year) + '.csv',
                     sep = ',', index = False)
        # Chemicals and groups
        Chemicals = df_N_PCU[['CAS NUMBER', 'CHEMICAL NAME']].drop_duplicates(keep = 'first')
        Chemicals['TYPE OF CHEMICAL'] = None
        Path_c = self._dir_path + '/Chemicals/Chemicals.csv'
        if os.path.exists(Path_c):
            df_c = pd.read_csv(Path_c)
            for index, row in Chemicals.iterrows():
                if (df_c['CAS NUMBER'] != row['CAS NUMBER']).all():
                    df_c = df_c.append(pd.Series(row, index = row.index.tolist()), \
                                        ignore_index = True)
            df_c.to_csv(Path_c, sep = ',', index  = False)
        else:
            Chemicals.to_csv(Path_c, sep = ',', index = False)


    def _Calculating_possible_waste_feed_supply(self, Flow, Concentration, Efficiency):
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
            Chemical_feed_flow = 100*Flow/10**-4
        Waste_flows = [100*Chemical_feed_flow/c for c in percentanges_c]
        Interval = tuple([min(Waste_flows), max(Waste_flows)])
        Middle = 0.5*(Waste_flows[0] + Waste_flows[-1])
        return Interval, Middle


    def Building_database_for_flows(self, nbins):
        def func(x):
            if x.first_valid_index() is None:
                return None
            else:
                return x[x.first_valid_index()]
        with open(self._dir_path + '/Ancillary/Flow_columns.yaml', mode = 'r') as f:
            dictionary_of_columns = yaml.load(f, Loader=yaml.FullLoader)
        dictionary_of_columns = {key: [el.strip() for el in val['columns'].split(',')] for key, val in dictionary_of_columns['TRI_Files'].items()}
        dfs = dict()
        for file, columns in dictionary_of_columns.items():
            df = pd.read_csv(self._dir_path + '/Ancillary/US_{}_{}.csv'.format(file, self.Year),
                            usecols = columns,
                            low_memory = False,
                            dtype = {'PRIMARY NAICS CODE': 'object'})
            dfs.update({file: df})
        # Energy recovery:
        cols_energy_methods = [col for col in dfs['1a'].columns if 'METHOD' in col and 'ENERGY' in col]
        df_energy = dfs['1a'][['TRIFID', 'CAS NUMBER', 'UNIT OF MEASURE',
                               'ON-SITE - ENERGY RECOVERY'] + cols_energy_methods]
        df_energy = df_energy.loc[pd.notnull(df_energy[cols_energy_methods]).sum(axis = 1) == 1]
        df_energy['METHOD CODE'] = df_energy[cols_energy_methods].apply(func, axis = 1)
        df_energy.rename(columns = {'ON-SITE - ENERGY RECOVERY': 'FLOW'}, inplace = True)
        df_energy.drop(columns = cols_energy_methods, inplace = True)
        del cols_energy_methods
        # Recycling:
        cols_recycling_methods = [col for col in dfs['1a'].columns if 'METHOD' in col and 'RECYCLING' in col]
        df_recycling = dfs['1a'][['TRIFID', 'CAS NUMBER', 'UNIT OF MEASURE',
                                 'ON-SITE - RECYCLED'] + cols_recycling_methods]
        df_recycling = df_recycling.loc[pd.notnull(df_recycling[cols_recycling_methods]).sum(axis = 1) == 1]
        df_recycling['METHOD CODE'] = df_recycling[cols_recycling_methods].apply(func, axis = 1)
        df_recycling.rename(columns = {'ON-SITE - RECYCLED': 'FLOW'}, inplace = True)
        df_recycling.drop(columns = cols_recycling_methods, inplace = True)
        del cols_recycling_methods
        # Treatment
        cols_treatment_methods = [col for col in dfs['2b'].columns if 'METHOD' in col]
        dfs['2b'] = dfs['2b'].loc[pd.notnull(dfs['2b'][cols_treatment_methods]).sum(axis = 1) == 1]
        dfs['2b']['METHOD CODE'] = dfs['2b'][cols_treatment_methods].apply(func, axis = 1)
        dfs['2b'].drop(columns = cols_treatment_methods, inplace = True)
        cols_for_merging = ['TRIFID','DOCUMENT CONTROL NUMBER',
                            'CAS NUMBER', 'ON-SITE - TREATED',
                            'UNIT OF MEASURE']
        df_treatment = pd.merge(dfs['1a'][cols_for_merging], dfs['2b'],
                                how = 'inner',
                                on = ['TRIFID',
                                    'DOCUMENT CONTROL NUMBER',
                                    'CAS NUMBER'])
        del dfs, cols_treatment_methods, cols_for_merging
        df_treatment.rename(columns = {'ON-SITE - TREATED': 'FLOW'}, inplace = True)
        df_treatment.drop(columns = ['DOCUMENT CONTROL NUMBER'], inplace = True)
        df_PCU_flows = pd.concat([df_treatment, df_recycling, df_energy],
                                ignore_index = True,
                                sort = True, axis = 0)
        del df_treatment, df_recycling, df_energy
        Chemicals_to_remove = ['MIXTURE', 'TRD SECRT']
        df_PCU_flows = df_PCU_flows.loc[~df_PCU_flows['CAS NUMBER'].isin(Chemicals_to_remove)]
        df_PCU_flows['CAS NUMBER'] = df_PCU_flows['CAS NUMBER'].apply(lambda x: str(int(x)) if not 'N' in x else x)
        df_PCU_flows.loc[df_PCU_flows['UNIT OF MEASURE'] == 'Pounds', 'FLOW'] *= 0.453592
        df_PCU_flows.loc[df_PCU_flows['UNIT OF MEASURE'] == 'Grams', 'FLOW'] *= 10**-3
        df_PCU_flows['FLOW'] = df_PCU_flows['FLOW'].round(6)
        df_PCU_flows = df_PCU_flows.loc[~(df_PCU_flows['FLOW'] == 0)]
        df_PCU_flows['UNIT OF MEASURE'] = 'kg'
        # Calling cleaned database
        columns_for_calling = ['TRIFID', 'CAS NUMBER', 'RANGE INFLUENT CONCENTRATION',
                            'METHOD CODE - 2004 AND PRIOR', 'EFFICIENCY ESTIMATION',
                            'PRIMARY NAICS CODE', 'WASTE STREAM CODE', 'TYPE OF MANAGEMENT']
        df_PCU_cleaned = pd.read_csv(self._dir_path + '/Datasets/Final_PCU_datasets/PCUs_DB_filled_{}.csv'.format(self.Year),
                            usecols = columns_for_calling, dtype = {'PRIMARY NAICS CODE': 'object'})
        df_PCU_cleaned.rename(columns = {'METHOD CODE - 2004 AND PRIOR': 'METHOD CODE'}, inplace = True)
        # Merging
        df_PCU_flows = pd.merge(df_PCU_flows, df_PCU_cleaned, how = 'inner',
                                on = ['TRIFID', 'CAS NUMBER', 'METHOD CODE'])
        df_PCU_flows[['RANGE INFLUENT CONCENTRATION', 'EFFICIENCY ESTIMATION', 'FLOW']] = \
                    df_PCU_flows[['RANGE INFLUENT CONCENTRATION', 'EFFICIENCY ESTIMATION', 'FLOW']]\
                    .applymap(lambda x: abs(x))
        df_PCU_flows['RANGE INFLUENT CONCENTRATION'] = df_PCU_flows['RANGE INFLUENT CONCENTRATION'].apply(lambda x: str(int(x)))
        df_PCU_flows[['WASTE FLOW RANGE', 'MIDDLE WASTE FLOW']] = df_PCU_flows.apply(lambda x:
                            pd.Series(self._Calculating_possible_waste_feed_supply(
                                                        x['FLOW'],
                                                        x['RANGE INFLUENT CONCENTRATION'],
                                                        x['EFFICIENCY ESTIMATION'])),
                            axis = 1)
        Max_value = df_PCU_flows['MIDDLE WASTE FLOW'].max()
        Min_value = df_PCU_flows['MIDDLE WASTE FLOW'].min()
        Order_max = int(math.log10(Max_value)) - 1
        Order_min = math.ceil(math.log10(Min_value))
        Delta = (Order_max - Order_min)/(nbins - 2)
        Bin_values = [Min_value - 10**(math.log10(Min_value) - 1)]
        Bin_values = Bin_values + [10**(Order_min + Delta*n) for n in range(nbins - 1)]
        Bin_values = Bin_values + [Max_value]
        Bin_values.sort()
        Bin_labels = [str(val) for val in range(1, len(Bin_values))]
        df_PCU_flows['MIDDLE WASTE FLOW INTERVAL'] = pd.cut(df_PCU_flows['MIDDLE WASTE FLOW'],
                                                         bins = Bin_values)
        df_PCU_flows['MIDDLE WASTE FLOW INTERVAL CODE'] = pd.cut(df_PCU_flows['MIDDLE WASTE FLOW'],
                                                         bins = Bin_values,
                                                         labels = Bin_labels,
                                                         precision = 0)
        df_PCU_flows.to_csv(self._dir_path + '/Datasets/Waste_flow/Waste_flow_to_PCUs_{}_{}.csv'.format(self.Year, nbins), sep = ',', index = False)


    def Organizing_substance_prices(self):
        # Organizing information about prices
        df_scifinder = pd.read_csv(self._dir_path + '/SciFinder/Chemical_Price.csv',
                                dtype = {'CAS NUMBER': 'object'})
        df_scifinder = df_scifinder.loc[\
                            (pd.notnull(df_scifinder['PRICE'])) & \
                            (df_scifinder['PRICE'] != 'Not found')]
        df_scifinder
        File_exchange = [file for file in os.listdir(self._dir_path + '/SciFinder') if 'Exchange' in file]
        df_exchange_rate = pd.read_csv(self._dir_path + '/SciFinder/{}'.format(File_exchange[0]))
        Exchange_rate = {row['CURRENCY']:row['EXCHANGE RATE TO USD'] for idx, row in df_exchange_rate.iterrows()}
        del df_exchange_rate
        df_scifinder['PRICE'] = df_scifinder.apply(lambda x: \
                                    Exchange_rate[x['CURRENCY']] \
                                    *float(x['PRICE']),
                                    axis = 1)
        df_scifinder['QUANTITY'] = df_scifinder['QUANTITY'].str.lower()
        df_scifinder['QUANTITY'] = df_scifinder['QUANTITY'].str.replace(' ', '')
        df_scifinder = df_scifinder[df_scifinder['QUANTITY'].str.contains(r'ton|[umk]{0,1}g')]
        idx = df_scifinder[~df_scifinder['QUANTITY'].str.contains(r'x')].index.tolist()
        df_scifinder.loc[idx, 'QUANTITY'] = '1x' + df_scifinder.loc[idx, 'QUANTITY']
        df_scifinder[['TIMES', 'MASS', 'UNIT']] = \
                    df_scifinder['QUANTITY'].str.extract('(\d+)x(\d+\.?\d*)(ton|[umk]{0,1}g)',
                    expand = True)
        dictionary_mass = {'g':1, 'mg':0.001, 'kg':1000, 'ug':10**-6, 'ton':907185}
        df_scifinder['QUANTITY'] = df_scifinder[['TIMES', 'MASS']]\
                                .apply(lambda x: float(x.values[0])*float(x.values[1]),
                                axis = 1)
        df_scifinder['QUANTITY'] = df_scifinder[['QUANTITY', 'UNIT']]\
                            .apply(lambda x: x.values[0]*dictionary_mass[x.values[1]], axis = 1)
        df_scifinder['UNIT PRICE (USD/g)'] = df_scifinder[['PRICE', 'QUANTITY']]\
                                        .apply(lambda x: x.values[0]/x.values[1],
                                            axis = 1)
        df_scifinder.drop(columns = ['COMPANY_NAME', 'COUNTRY', 'PURITY',
                                    'CURRENCY', 'TIMES', 'MASS', 'UNIT',
                                    'QUANTITY', 'PRICE'],
                        inplace = True)
        df_scifinder = df_scifinder.groupby('CAS NUMBER', as_index = False).median()
        # Calling PCU
        df_PCU = pd.read_csv(self._dir_path + '/Datasets/Final_PCU_datasets/PCUs_DB_filled_{}.csv'.format(self.Year),
                            usecols = ['TRIFID', 'CAS NUMBER', 'METHOD CODE - 2004 AND PRIOR'])
        df_PCU = df_PCU[~df_PCU['METHOD CODE - 2004 AND PRIOR'].str.contains('\+')]
        # Separating categories and chemicals
        categories = pd.read_csv(self._dir_path + '/Chemicals/Chemicals_in_categories.csv',
                            usecols = ['CAS NUMBER', 'CATEGORY CODE'])
        categories['CAS NUMBER'] = categories['CAS NUMBER'].str.replace('-', '')
        chemicals = pd.read_csv(self._dir_path + '/Chemicals/Chemicals.csv',
                                usecols = ['CAS NUMBER'])
        chemicals = chemicals.loc[~chemicals['CAS NUMBER']\
                        .isin(list(categories['CATEGORY CODE'].unique())),\
                        'CAS NUMBER'].tolist()
        df_PCU_chemicals = df_PCU.loc[df_PCU['CAS NUMBER'].isin(chemicals)]
        df_PCU_categories = df_PCU.loc[~df_PCU['CAS NUMBER'].isin(chemicals)]
        df_PCU_categories.rename(columns = {'CAS NUMBER': 'CATEGORY CODE'},
                                inplace = True)
        del chemicals, df_PCU
        # Merging prices with chemicals
        df_PCU_chemicals = pd.merge(df_PCU_chemicals, df_scifinder,
                                    how = 'inner',
                                    on = 'CAS NUMBER')
        # Calling CDR
        df_CDR = pd.read_csv(self._dir_path +  '/CDR/Substances_by_facilities.csv',
                            usecols = ['STRIPPED_CHEMICAL_ID_NUMBER',
                                        'PGM_SYS_ID'],
                            dtype = {'STRIPPED_CHEMICAL_ID_NUMBER': 'object'}            )
        df_CDR.rename(columns = {'STRIPPED_CHEMICAL_ID_NUMBER': 'CAS NUMBER',
                                'PGM_SYS_ID': 'TRIFID'},
                     inplace = True)
        df_CDR = pd.merge(df_CDR, categories, on = 'CAS NUMBER',
                        how = 'inner')
        df_PCU_categories = pd.merge(df_PCU_categories, df_CDR,
                        on = ['CATEGORY CODE', 'TRIFID'],
                        how = 'inner')
        # Merging prices with categories
        df_PCU_categories = pd.merge(df_PCU_categories, df_scifinder,
                                    how = 'inner',
                                    on = 'CAS NUMBER')
        df_PCU_categories.drop(columns = ['CATEGORY CODE'],
                            inplace = True)
        df_PCU = pd.concat([df_PCU_categories, df_PCU_chemicals],
                            ignore_index = True,
                            sort = True, axis = 0)
        df_PCU.to_csv(self._dir_path + '/Datasets/Chemical_price/Chemical_price_vs_PCU_{}.csv'.format(self.Year),
                    sep = ',', index = False)


    def pollution_abatement_cost_and_expenditure(self):
        # Calling PCU
        df_PCU = pd.read_csv(self._dir_path + '/Datasets/Waste_flow/Waste_flow_to_PCUs_{}_10.csv'.format(2004),
                        low_memory = False,
                        usecols = ['PRIMARY NAICS CODE', 'WASTE STREAM CODE',
                                'TYPE OF MANAGEMENT', 'METHOD CODE', 'TRIFID',
                                'MIDDLE WASTE FLOW'],
                        dtype = {'PRIMARY NAICS CODE': 'object'})
        df_PCU['FLOW'] = df_PCU.groupby(['TRIFID', 'METHOD CODE'])['MIDDLE WASTE FLOW'].transform('median')
        df_PCU.drop(columns = ['TRIFID', 'METHOD CODE', 'MIDDLE WASTE FLOW'], inplace = True)
        df_PCU.drop_duplicates(keep = 'first', inplace = True)
        df_PCU = df_PCU.groupby(['PRIMARY NAICS CODE', 'TYPE OF MANAGEMENT',
                                 'WASTE STREAM CODE'], as_index = False)\
                                 .agg({'FLOW': ['mean', 'std']})
        df_PCU = pd.DataFrame.from_records(df_PCU.values, columns = ['NAICS code','Activity', 'Media', 'Mean flow','SD flow'])
        df_PCU.drop_duplicates(keep = 'first', inplace = True)
        df_PCU = df_PCU[df_PCU['NAICS code'].str.contains(r'^3[123]')]
        df_PCU['SD flow'] = df_PCU['SD flow'].fillna(df_PCU['Mean flow']*0.1) # Imputing coefficient of variation (<= 1)
        # Method of moments
        df_PCU['mu'] = df_PCU[['Mean flow', 'SD flow']].apply(lambda x: np.log(x.values[0]**2/(x.values[1]**2 + x.values[0]**2)**0.5) ,
                                                        axis = 1)
        df_PCU['theta_2'] = df_PCU[['Mean flow', 'SD flow']].apply(lambda x: np.log(x.values[1]**2/x.values[0]**2 + 1) ,
                                                        axis = 1)
        # U.S. Pollution Abatement Operating Costs - Survey 2005
        df_PAOC = pd.read_csv(self._dir_path + '/US_Census_Bureau/Pollution_Abatement_Operating_Costs_2005.csv',
                        low_memory = False, header = None, skiprows = [0,1],
                        usecols = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                        names = ['NAICS code', 'Total PAOC', 'Activity - treatment',
                                'Activity - prevention', 'Activity - recycling',
                                'Activity - disposal', 'Media - air', 'Media - water',
                                'Media - solid waste', 'RSE for total PAOC'])
        df_PAOC = df_PAOC.loc[pd.notnull(df_PAOC).all(axis = 1)]
        # Substracting  prevention
        df_PAOC['Total PAOC'] = df_PAOC['Total PAOC'] - df_PAOC['Activity - prevention'] - df_PAOC['Activity - disposal']
        df_PAOC.drop(columns = ['Activity - prevention', 'Activity - disposal'],
                    inplace = True)
        # The media and the activity and supposed to be indepent events
        # Proportions activities
        col_activities = [col for col in df_PAOC.columns if 'Activity' in col]
        df_PAOC[col_activities] = df_PAOC[col_activities].div(df_PAOC[col_activities].sum(axis = 1), axis = 0)
        col_medias = [col for col in df_PAOC.columns if 'Media' in col]
        df_PAOC[col_medias] = df_PAOC[col_medias].div(df_PAOC[col_medias].sum(axis = 1), axis = 0)
        # U.S. Pollution Abatement Capital Expenditures - Survey 2005
        df_PACE = pd.read_csv(self._dir_path + '/US_Census_Bureau/Pollution_Abatement_Capital_Expenditures_2005.csv',
                        low_memory = False, header = None, skiprows = [0,1],
                        usecols = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                        names = ['NAICS code', 'Total PACE', 'Activity - treatment',
                                'Activity - prevention', 'Activity - recycling',
                                'Activity - disposal', 'Media - air', 'Media - water',
                                'Media - solid waste', 'RSE for total PACE'])
        df_PACE = df_PACE.loc[pd.notnull(df_PACE).all(axis = 1)]
        df_PACE['Total PACE'] = df_PACE['Total PACE'] - df_PACE['Activity - prevention'] - df_PACE['Activity - disposal']
        df_PACE.drop(columns = ['Activity - prevention', 'Activity - disposal'],
                    inplace = True)
        df_PACE[col_activities] = df_PACE[col_activities].div(df_PACE[col_activities].sum(axis = 1), axis = 0)
        df_PACE[col_medias] = df_PACE[col_medias].div(df_PACE[col_medias].sum(axis = 1), axis = 0)
        for col_a in col_activities:
            for col_m in col_medias:
                df_PAOC['{} for {}'.format(col_a.replace('Activity - ', '').capitalize(),
                                           col_m.replace('Media - ', ''))] = \
                            df_PAOC[col_a]* df_PAOC[col_m] # P(media and activy) = P(media)*P(activity). Thet are independent events
                df_PACE['{} for {}'.format(col_a.replace('Activity - ', '').capitalize(),
                                           col_m.replace('Media - ', ''))] = \
                            df_PACE[col_a]* df_PACE[col_m] # P(media and activy) = P(media)*P(activity). Thet are independent events
        df_PAOC.drop(columns = col_medias + col_activities, inplace = True)
        df_PACE.drop(columns = col_medias + col_activities, inplace = True)
        # Statistics of U.S. Businesses - Survey 2005
        df_SUSB = Organizing_sample(20378, self._dir_path) # Sampled establishments in 2005
        df_SUSB['Establishments (employees >= 20)'] = 1
        df_SUSB['Establishments (employees >= 20)'] = \
                        df_SUSB.groupby('NAICS code')\
                        ['Establishments (employees >= 20)'].transform('sum')
        df_SUSB['Total shipment'] = \
                        df_SUSB.groupby('NAICS code')\
                        ['Total shipment establishment'].transform('sum')
        df_SUSB['Info establishments'] = df_SUSB[['Establishment', \
                                                'Total shipment establishment',  \
                                                'P-in-cluster']]\
                                                .apply(lambda x: x.tolist(), axis = 1)
        df_SUSB.drop(columns = ['P-cluster', 'Establishment', 'P-in-cluster',
                                'P-selected', 'Total shipment establishment', 'Unit'],
                    inplace = True)
        df_SUSB = df_SUSB.groupby(['NAICS code', 'Total shipment', \
                                    'Establishments (employees >= 20)'],
                                                          as_index = False)\
                                                .agg({'Info establishments': lambda x: {val[0]: [val[1], val[2]] for idx, val in enumerate(x)}})
        # Joining sources from census
        df_PACE = pd.merge(df_PACE, df_SUSB, on = 'NAICS code', how = 'left')
        df_PAOC = pd.merge(df_PAOC, df_SUSB, on = 'NAICS code', how = 'left')
        # Searching higher naics levels (no in clusters but containing them)
        df_PACE = df_PACE.where(pd.notnull(df_PACE), None)
        df_PACE[['Establishments (employees >= 20)', \
                'Total shipment',\
                'Info establishments']] = \
                        df_PACE.apply(lambda x: \
                        searching_establishments_by_hierarchy(x['NAICS code'], df_SUSB)
                        if not x['Establishments (employees >= 20)']
                        else pd.Series([int(x['Establishments (employees >= 20)']),\
                                            x['Total shipment'],\
                                            x['Info establishments']]),\
                        axis = 1)
        df_PACE = df_PACE.loc[pd.notnull(df_PACE).all(axis = 1)]
        df_PAOC = df_PAOC.where(pd.notnull(df_PAOC), None)
        df_PAOC[['Establishments (employees >= 20)', \
                'Total shipment',\
                'Info establishments']] = \
                        df_PAOC.apply(lambda x: \
                        searching_establishments_by_hierarchy(x['NAICS code'], df_SUSB)
                        if not x['Establishments (employees >= 20)']
                        else pd.Series([int(x['Establishments (employees >= 20)']),\
                                           x['Total shipment'],\
                                           x['Info establishments']]),\
                        axis = 1)
        df_PAOC = df_PAOC.loc[pd.notnull(df_PAOC).all(axis = 1)]
        del df_SUSB
        # Organizing by activity and media
        df_PACE_for_merging = pd.DataFrame()
        df_PAOC_for_merging = pd.DataFrame()
        Dictionary_relation = {'W': 'water', 'L': 'water', 'A': 'air', 'S': 'solid waste',
                            'Treatment': 'Treatment', 'Energy recovery': 'Recycling',
                            'Recycling': 'Recycling'}
        Medias = ['W', 'L', 'A', 'S']
        Activities = ['Treatment', 'Energy recovery', 'Recycling']
        # Inflation rate in the U.S. between 2005 and 2020 is 35.14%
        for Activity in Activities:
            for Media in Medias:
                Factor_col =  '{} for {}'.format(Dictionary_relation[Activity], Dictionary_relation[Media])
                df_PACE_aux = df_PACE[['NAICS code', 'Total PACE', \
                                       'RSE for total PACE', \
                                       'Establishments (employees >= 20)',\
                                       'Info establishments',\
                                       'Total shipment']]
                df_PACE_aux['Total PACE'] = df_PACE_aux['Total PACE']*1.3514*10**6
                df_PACE_aux['Media'] = Media
                df_PACE_aux['Activity'] = Activity
                df_PACE_aux['P-media_&_activiy'] = df_PACE[Factor_col]
                df_PACE_for_merging = pd.concat([df_PACE_for_merging, df_PACE_aux], ignore_index = True,
                                           sort = True, axis = 0)
                df_PAOC_aux = df_PAOC[['NAICS code', 'Total PAOC', \
                                       'RSE for total PAOC', \
                                       'Establishments (employees >= 20)',\
                                       'Info establishments',\
                                       'Total shipment']]
                df_PAOC_aux['Total PAOC'] = df_PAOC_aux['Total PAOC']*1.3514*10**6
                df_PAOC_aux['Media'] = Media
                df_PAOC_aux['Activity'] = Activity
                df_PAOC_aux['P-media_&_activiy'] = df_PAOC[Factor_col]
                df_PAOC_for_merging = pd.concat([df_PAOC_for_merging, df_PAOC_aux], ignore_index = True,
                                           sort = True, axis = 0)
        # Identifying probable establishment based on the pobability of media and activity
        df_PAOC_for_merging = df_PAOC_for_merging.loc[pd.notnull(df_PAOC_for_merging['P-media_&_activiy'])]
        df_PAOC_for_merging[['Probable establishments by activity & media',\
                            'Info probable establishments']] = \
                            df_PAOC_for_merging[['Info establishments', 'P-media_&_activiy']]\
                            .apply(lambda x: selecting_establishment_by_activity_and_media(\
                                                                        x.values[0],
                                                                        x.values[1]),
                            axis =  1)
        df_PAOC_for_merging.drop(columns = ['Info establishments', \
                                            'Establishments (employees >= 20)'],\
                                inplace = True)
        df_PACE_for_merging = df_PACE_for_merging.loc[pd.notnull(df_PACE_for_merging['P-media_&_activiy'])]
        df_PACE_for_merging[['Probable establishments by activity & media',\
                            'Info probable establishments']] = \
                            df_PACE_for_merging[['Info establishments', 'P-media_&_activiy']]\
                            .apply(lambda x: selecting_establishment_by_activity_and_media(\
                                                                        x.values[0],
                                                                        x.values[1]),
                            axis =  1)
        df_PACE_for_merging.drop(columns = ['Info establishments', \
                                            'Establishments (employees >= 20)'],\
                                inplace = True)
        # Calculating total by media and activity
        df_PACE_for_merging = df_PACE_for_merging.groupby('NAICS code',
                                                          as_index = False)\
                                                .apply(lambda x: normalizing_shipments(x))
        df_PAOC_for_merging = df_PAOC_for_merging.groupby('NAICS code',
                                                          as_index = False)\
                                                .apply(lambda x: normalizing_shipments(x))
        # Joining census with TRI
        df_PACE = pd.merge(df_PACE_for_merging, df_PCU,
                            on = ['NAICS code', 'Media', 'Activity'],
                            how = 'right')
        df_PAOC = pd.merge(df_PAOC_for_merging, df_PCU,
                            on = ['NAICS code', 'Media', 'Activity'],
                            how = 'right')
        idx = df_PACE.loc[df_PACE['P-media_&_activiy'].isnull()].index.tolist()
        df_PACE[['RSE for total PACE', \
                'Probable establishments by activity & media', \
                'P-media_&_activiy', 'Total PACE', 'Total shipment',\
                'Info probable establishments']].iloc[idx] = \
                df_PACE[['NAICS code', 'Media',
                        'Activity']].iloc[idx].apply(lambda x: searching_census(x.values[0],
                                                                            x.values[1],
                                                                            x.values[2],
                                                                            df_PACE_for_merging),
                                            axis = 1)
        df_PACE = df_PACE.loc[pd.notnull(df_PACE).all(axis = 1)]
        idx = df_PAOC.loc[df_PAOC['P-media_&_activiy'].isnull()].index.tolist()
        df_PAOC[['RSE for total PAOC', \
                'Probable establishments by activity & media', \
                'P-media_&_activiy', 'Total PAOC', 'Total shipment',\
                'Info probable establishments']].iloc[idx] = \
                df_PAOC[['NAICS code', 'Media',
                        'Activity']].iloc[idx].apply(lambda x: searching_census(x.values[0],
                                                                            x.values[1],
                                                                            x.values[2],
                                                                            df_PAOC_for_merging),
                                            axis = 1)
        df_PAOC = df_PAOC.loc[pd.notnull(df_PAOC).all(axis = 1)]
        # Calculating mass to activity and media assuming lognormal distribution
        df_PAOC = df_PAOC[df_PAOC['P-media_&_activiy'] != 0.0]
        df_PAOC['Probable mass by activity & media'] = \
                            df_PAOC[['mu', 'theta_2','Info probable establishments', 'P-media_&_activiy']]\
                            .apply(lambda x: estimating_mass_by_activity_and_media(x.values[0],
                                                                                   x.values[1],
                                                                                   x.values[2],
                                                                                   x.values[3]), axis =  1)
        df_PACE = df_PACE[df_PACE['P-media_&_activiy'] != 0.0]
        df_PACE['Probable mass by activity & media'] = \
                            df_PACE[['mu', 'theta_2','Info probable establishments', 'P-media_&_activiy']]\
                            .apply(lambda x: estimating_mass_by_activity_and_media(x.values[0],
                                                                                   x.values[1],
                                                                                   x.values[2],
                                                                                   x.values[3]), axis =  1)
        # Assuming a normal distribution and a confidence level of 95%
        Z = norm.ppf(0.975)
        df_PAOC[['Mean PAOC', 'SD PAOC', 'CI at 95% for Mean PAOC']] \
                             = df_PAOC.apply(lambda x: mean_standard(\
                                                 x, Z),
                                             axis = 1)
        df_PAOC['Unit'] = 'USD/kg'
        df_PAOC = df_PAOC.loc[pd.notnull(df_PAOC).all(axis = 1)]
        df_PAOC = df_PAOC.round(6)
        df_PAOC = df_PAOC.loc[df_PAOC['Mean PAOC'] != 0]
        df_PACE[['Mean PACE', 'SD PACE', 'CI at 95% for Mean PACE']] \
                              = df_PACE.apply(lambda x: mean_standard(\
                                                  x, Z),
                                              axis = 1)
        df_PACE['Unit'] = 'USD/kg'
        df_PACE = df_PACE.loc[pd.notnull(df_PACE).all(axis = 1)]
        df_PACE = df_PACE.round(6)
        df_PACE = df_PACE.loc[df_PACE['Mean PACE'] != 0]
        #  Saving
        cols = ['NAICS code', 'Activity', 'Media', \
                'Probable establishments by activity & media', \
                'Probable mass by activity & media', 'Mean PAOC', \
                'SD PAOC', 'CI at 95% for Mean PAOC','Unit']
        df_PAOC = df_PAOC[cols]
        df_PAOC.to_csv(self._dir_path + '/Datasets/PCU_expenditure_and_cost/PAOC.csv', sep = ',', index = False)
        cols = ['NAICS code', 'Activity', 'Media', \
                'Probable establishments by activity & media', \
                'Probable mass by activity & media', 'Mean PACE', \
                'SD PACE', 'CI at 95% for Mean PACE','Unit']
        df_PACE = df_PACE[cols]
        df_PACE.to_csv(self._dir_path + '/Datasets/PCU_expenditure_and_cost/PACE.csv', sep = ',', index = False)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(argument_default = argparse.SUPPRESS)

    parser.add_argument('Option',
                        help = 'What do you want to do:\
                        [A]: Recover information from TRI.\
                        [B]: File for statistics. \
                        [C]: File for recycling. \
                        [D]: Further cleaning of database. \
                        [E]: Organizing file with flows (1987-2004). \
                        [F]: Organizing file with substance prices (1987 - 2004). \
                        [G]: Pollution control cost and expenditure (only 2005).', \
                        type = str)

    parser.add_argument('-Y', '--Year', nargs = '+',
                        help = 'Records with up to how many PCUs you want to include?.',
                        type = str,
                        required = False,
                        default = [2018])

    parser.add_argument('-N_Bins',
                        help = 'Number of bins to split the middle waste flow values',
                        type = int,
                        default =  10,
                        required = False)


    args = parser.parse_args()
    start_time =  time.time()

    for Year in args.Year:
        Building = PCU_DB(int(Year))
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

    print('Execution time: %s sec' % (time.time() - start_time))
