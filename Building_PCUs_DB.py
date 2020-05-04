# -*- coding: utf-8 -*-
#!/usr/bin/env python

# Note:
# 1. Range of Influent Concentration was reported from 1987 through 2004
# 2. Treatment Efficiency Estimation was reported from 1987 through 2004

import warnings
warnings.simplefilter(action = 'ignore', category = FutureWarning)
import pandas as pd
pd.options.mode.chained_assignment = None
import os
import argparse
import numpy as np
import re
import time
import unicodedata
from itertools import combinations

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
        df_PCUs.to_csv(self._dir_path + '/Datasets/PCUs_DB_' + str(self.Year) + '.csv',
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
        df_statistics.to_csv(self._dir_path + '/Statistics/DB_for_Statistics_' + str(self.Year) + '.csv',
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
            values = [v for v in row[elements_total] if v != 0.0]
            cases = list()
            for n_elements_sum in range(1, len(values) + 1):
                comb = combinations(values, n_elements_sum)
                for comb_values in comb:
                    sumatory = sum(comb_values)
                    cases.append(row['ON-SITE - RECYCLED']/(row['ON-SITE - RECYCLED'] + sumatory)*100)
            try:
                if len(list(set(cases))) == 1 and cases[0] == 100:
                    return [100]*6 + [row['ON-SITE - RECYCLED']]
                else:
                    return [np.min(cases)] + np.quantile(cases, [0.25, 0.5, 0.75]).tolist() + [np.max(cases), np.mean(cases), row['ON-SITE - RECYCLED']]
            except ValueError:
                a = np.empty((7))
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
        df[['LOWER EFFICIENCY', '1ST QUARTILE OF EFFICIENCY', '2ND QUARTILE OF EFFICIENCY', \
            '3RD QUARTILE OF EFFICIENCY', 'UPPER EFFICIENCY', 'MEAN OF EFFICIENCY', 'ON-SITE - RECYCLED']]\
             = df.apply(lambda x: pd.Series(_division(x, elements_total)), axis = 1)
        df = df.loc[pd.notnull(df['UPPER EFFICIENCY'])]
        df['INTERQUARTILE RANGE'] = df.apply(lambda x: x['3RD QUARTILE OF EFFICIENCY'] \
                                 - x['1ST QUARTILE OF EFFICIENCY'], axis = 1)
        df['UPPER EFFICIENCY OUTLIER?'] = df.apply(lambda x: 'YES' if x['UPPER EFFICIENCY'] > \
                                1.5*x['INTERQUARTILE RANGE'] + x['3RD QUARTILE OF EFFICIENCY'] else 'NO',
                                axis = 1)
        df['LOWER EFFICIENCY OUTLIER?'] = df.apply(lambda x: 'YES' if x['LOWER EFFICIENCY'] < \
                                x['1ST QUARTILE OF EFFICIENCY'] - 1.5*x['INTERQUARTILE RANGE'] else 'NO',
                                axis = 1)
        df = df[['TRIFID', 'PRIMARY NAICS CODE', 'CAS NUMBER', 'ON-SITE - RECYCLED', 'UNIT OF MEASURE', \
                'LOWER EFFICIENCY', 'LOWER EFFICIENCY OUTLIER?', '1ST QUARTILE OF EFFICIENCY', \
                '2ND QUARTILE OF EFFICIENCY', '3RD QUARTILE OF EFFICIENCY', 'UPPER EFFICIENCY', \
                'UPPER EFFICIENCY OUTLIER?', 'MEAN OF EFFICIENCY', 'METHOD']]
        df.to_csv(self._dir_path + '/Statistics/DB_for_Solvents_' + str(self.Year) + '.csv',
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
        # Calling PCU
        PCU = pd.read_csv(self._dir_path + '/Datasets/PCUs_DB_' + str(self.Year) + '.csv',
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
        Statistics = pd.read_csv(self._dir_path + '/Statistics/DB_for_Statistics_' + str(self.Year) + '.csv',
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
            Recycling_statistics = pd.read_csv(self._dir_path + '/Statistics/DB_for_Solvents_' + str(self.Year) +  '.csv',
                                    low_memory = False,
                                    usecols = ['TRIFID', 'PRIMARY NAICS CODE', 'CAS NUMBER', \
                                               'UPPER EFFICIENCY', 'UPPER EFFICIENCY OUTLIER?', 'METHOD'],
                                    converters = {'CAS NUMBER': lambda x: x if re.search(r'^[A-Z]', x) else str(int(x))})
            Recycling_statistics['PRIMARY NAICS CODE'] = Recycling_statistics['PRIMARY NAICS CODE'].astype('int')
            Recycling_statistics = Recycling_statistics.loc[Recycling_statistics['UPPER EFFICIENCY OUTLIER?'] == 'NO']
            Recycling_statistics.drop(columns = ['UPPER EFFICIENCY OUTLIER?'], axis = 1)
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
        print(df_N_PCU.info())
        df_N_PCU.to_csv(self._dir_path + '/Datasets/PCUs_DB_filled_' + str(self.Year) + '.csv',
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



if __name__ == '__main__':

    parser = argparse.ArgumentParser(argument_default = argparse.SUPPRESS)

    parser.add_argument('Option',
                        help = 'What do you want to do:\
                        [A]: Recover information from TRI.\
                        [B]: File for statistics. \
                        [C]: File for recycling. \
                        [D]: Further cleaning of database.', \
                        type = str)

    parser.add_argument('-Y', '--Year', nargs = '+',
                        help = 'Records with up to how many PCUs you want to include?.',
                        type = str)


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

    print('Execution time: %s sec' % (time.time() - start_time))
