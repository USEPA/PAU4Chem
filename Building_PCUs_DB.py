# -*- coding: utf-8 -*-
#!/usr/bin/env python

# Note:
# 1. Range of Influent Concentration was reported from 1987 through 2004
# 2. Treatment Efficiency Estimation was reported from 1987 through 2004

import pandas as pd
pd.options.mode.chained_assignment = None
import os
import argparse
import numpy as np
import re
import time

class PCU_DB:

    def __init__(self, Year):
        self._dir_path = os.path.dirname(os.path.realpath(__file__)) # Working Directory
        self.Year = Year
        #self._dir_path = os.getcwd() # if you are working on Jupyter Notebook

    def calling_TRI_Files(self):
        TRI_Files = dict()
        for file in ['1a', '2b']:
            columns = pd.read_csv(self._dir_path + '/Ancillary/TRI_File_' + file + '_needed_columns.txt',
                                header = None)
            columns =  list(columns.iloc[:,0])
            df = pd.read_csv(self._dir_path + '/Ancillary/US_' + file + '_' + str(self.Year) + '.csv',
                            usecols = columns,
                            low_memory = False)
            df = df.where(pd.notnull(df), None)
            TRI_Files.update({file: df})
        return TRI_Files


    def _efficiency_estimation_to_range(self, x):
        if x != np.nan:
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
            if row['CLASSIFICATION'] == 'DIOXIN':
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


    # def _fuzzy_inference(self, df):
    #     N = dict()
    #     rows = df.shape[0]
    #     if rows == 1:
    #         return df
    #     else:
    #         naics = df['PRIMARY NAICS CODE'].unique().tolist()[0]
    #         cas = df['CAS NUMBER'].unique().tolist()[0]
    #         efficiency = df['EFFICIENCY RANGE CODE'].tolist()
    #         concentration = df['RANGE INFLUENT CONCENTRATION'].tolist()
    #         waste_code = df['WASTE STREAM CODE'].tolist()
    #         n_actual = len(waste_code)
    #         if self.Year <= 2004:
    #             efficiency = df['EFFICIENCY ESTIMATION'].tolist()
    #             #indices = [i if x == 0.0 else None for i, x in enumerate(df['EFFICIENCY ESTIMATION'].tolist())]
    #             indices = [i if x == 0.0 else None for i, x in enumerate(efficiency)]
    #             indices_f = {v: i for i, v in enumerate(set(range(len(waste_code))) - set(indices))}
    #             efficiency = [x for i, x in enumerate(efficiency) if not i in indices]
    #             concentration = [x for i, x in enumerate(concentration) if not i in indices]
    #             waste_code = [x for i, x in enumerate(waste_code) if not i in indices]
    #         n = len(waste_code)
    #         if n >= 2:
    #             # Based on position
    #             m_criteria = 1
    #             N_aux = np.empty((n, n))
    #             N_aux[:] = np.nan
    #             for i in range(n):
    #                 for j in range(i, n):
    #                     dist = j - i
    #                     if dist == 0:
    #                         N_aux[i][j] = 0
    #                     elif dist == 4:
    #                         N_aux[i][j] = 4
    #                     elif dist == 3:
    #                         N_aux[i][j] = 3
    #                     elif dist == 2:
    #                         N_aux[i][j] = 2
    #                     elif dist == 1:
    #                         N_aux[i][j] = 1
    #             N =  N_aux
    #             # Based on efficiency
    #             m_criteria = m_criteria + 1
    #             #dict_eff = {'E1':0, 'E2':0, 'E3':1, 'E4':2, 'E5':3, 'E6':4}
    #             N_aux = np.empty((n, n))
    #             N_aux[:] = np.nan
    #             for i in range(n):
    #                 for j in range(i, n):
    #                     #dist = dict_eff[efficiency[j]] - dict_eff[efficiency[i]]
    #                     if self.Year >= 2005:
    #                         dist =  (self._range_to_efficiency_estimation(efficiency[i]) - \
    #                                 self._range_to_efficiency_estimation(efficiency[j]))/100
    #                     else:
    #                         dist = (efficiency[i] - efficiency[j])/100
    #                     if dist == 0:
    #                         N_aux[i][j] = 0
    #                     elif (0.5 < dist) & (dist <= 1):#dist == 4:
    #                         N_aux[i][j] = 4
    #                     elif (0.245 < dist) & (dist <= 0.5):#dist == 3:
    #                         N_aux[i][j] = 3
    #                     elif (0.003 < dist) & (dist <= 0.245):#dist == 2:
    #                         N_aux[i][j] = 2
    #                     elif (0 < dist) & (dist <= 0.003):#dist == 1:
    #                         N_aux[i][j] = 1
    #                     elif (-1 <= dist) & (dist < -0.5):#dist == -4:
    #                         N_aux[i][j] = -4
    #                     elif (-0.5 <= dist) & (dist < -0.245):#dist == -3:
    #                         N_aux[i][j] = -3
    #                     elif (-0.245 <= dist) & (dist < -0.003):#dist == -2:
    #                         N_aux[i][j] = -2
    #                     elif (-0.003 <= dist) & (dist < 0):#dist == -1:
    #                         N_aux[i][j] = -1
    #             N =  np.concatenate((N, N_aux), axis = 0)
    #             # Based on concentration
    #             if self.Year <= 2004:
    #                 m_criteria = m_criteria + 1
    #                 #dict_con = {1:0, 2:1, 3:2, 4:3, 5:4}
    #                 N_aux = np.empty((n, n))
    #                 N_aux[:] = np.nan
    #                 for i in range(n):
    #                     for j in range(i, n):
    #                         dist =  (self._range_to_concentration_estimation(concentration[i]) - \
    #                                 self._range_to_concentration_estimation(concentration[j]))/100
    #                         #dist = dict_con[concentration[j]] - dict_con[concentration[i]]
    #                         if dist == 0:
    #                             N_aux[i][j] = 0
    #                         elif (0.25 < dist) & (dist <= 0.5):#dist == 4:
    #                             N_aux[i][j] = 4
    #                         elif (0.125 < dist) & (dist <= 0.25):#dist == 3:
    #                             N_aux[i][j] = 3
    #                         elif (0.004 < dist) & (dist <= 0.125):#dist == 2:
    #                             N_aux[i][j] = 2
    #                         elif (0 < dist) & (dist <= 0.004):#dist == 1:
    #                             N_aux[i][j] = 1
    #                         elif (-0.5 <= dist) & (dist < -0.25):#dist == -4:
    #                             N_aux[i][j] = -4
    #                         elif (-0.25 <= dist) & (dist < -0.125):#dist == -3:
    #                             N_aux[i][j] = -3
    #                         elif (-0.125 <= dist) & (dist < -0.004):#dist == -2:
    #                             N_aux[i][j] = -2
    #                         elif (-0.004 <= dist) & (dist < 0):#dist == -1:
    #                             N_aux[i][j] = -1
    #                 N =  np.concatenate((N, N_aux), axis = 0)
    #             # Based on statistics
    #             m_criteria = m_criteria + 1
    #             N_aux =  np.empty((n, n))
    #             N_aux[:] = np.nan
    #             df_statistics = pd.read_csv(self._dir_path + '/Statistics/DB_for_Statistics_' + str(self.Year) + '.csv',
    #                             usecols = ['CAS', 'NAICS', 'WASTE'],
    #                             low_memory = False,
    #                             converters = {'NAICS':  lambda x: str(int(float(x)))})
    #             df_statistics = df_statistics.loc[df_statistics['CAS'] == cas]
    #             if df_statistics.shape[0] != n:
    #                 df_statistics['NAICS STRUCTURE'] = df_statistics.apply(lambda x: self._searching_naics(\
    #                                     x['NAICS'], naics), axis = 1)
    #                 naics_structure = ['National Industry', 'NAICS Industry', 'Industry Group',
    #                                  'Subsector', 'Sector', 'Nothing']
    #                 count = {}
    #                 j = 0
    #                 while not count:
    #                     structure = naics_structure[j]
    #                     df_naics = df_statistics.loc[df_statistics['NAICS STRUCTURE'] ==  structure]
    #                     if (structure == 'National Industry') and (df_naics.shape[0] == n):
    #                         count = {i: waste_code.count(i)/len(waste_code) for i in list(set(waste_code))}
    #                     else:
    #                         values = df_naics['WASTE'].value_counts(normalize = True).keys().tolist()
    #                         counts = df_naics['WASTE'].value_counts(normalize = True).tolist()
    #                         count = {values[i]: counts[i] for i in range(len(values))}
    #                     j = j + 1
    #             else:
    #                 count = {i: waste_code.count(i)/len(waste_code) for i in list(set(waste_code))}
    #             for i in range(n):
    #                 for j in range(i, n):
    #                     dist =  count[waste_code[i]] - count[waste_code[j]]
    #                     if dist == 0.0:
    #                         N_aux[i][j] = 0
    #                     elif (dist > 0.0) and (dist <= 0.25):
    #                         N_aux[i][j] = 1
    #                     elif (dist > 0.25) and (dist <= 0.5):
    #                         N_aux[i][j] = 2
    #                     elif (dist > 0.5) and (dist <= 0.75):
    #                         N_aux[i][j] = 3
    #                     elif (dist > 0.75) and (dist <= 1.0):
    #                         N_aux[i][j] = 4
    #                     elif (dist < 0.0) and (dist >= -0.25):
    #                         N_aux[i][j] = -1
    #                     elif (dist < -0.25) and (dist >= -0.5):
    #                         N_aux[i][j] = -2
    #                     elif (dist < -0.5) and (dist >= -0.75):
    #                         N_aux[i][j] = -3
    #                     elif (dist < -0.75) and (dist >= -1.0):
    #                         N_aux[i][j] = -4
    #             N =  np.concatenate((N, N_aux), axis = 0)
    #             # Definition of variables
    #             W = np.zeros((1, n)) # initial value of weights vector (desired output)
    #             w = np.zeros((m_criteria, n)) # initial value of weithts matrix
    #             phi = np.zeros((m_criteria, 1)) # diversification degree
    #             eps = np.zeros((m_criteria, 1)) # entropy
    #             theta = np.zeros((1, m_criteria)) # criteria' uncertainty degrees
    #             sumphi = 0 # Initial value of diversification degree
    #             # Triangular fuzzy numbers (TFN)
    #             # TFN is a vector with n segments and each one has 3 different numbers
    #             # The segments are the linguistic scales
    #             # The 3 differente numbers in the segments are a triangular fuzzy number (l,m,u)
    #             # the first segment: equal importance; the second one: moderate importance of one over another;
    #             # the third one: strong importance of one over another; the fourth one: very strong importance of one over another
    #             # the fifth one: Absolute importance of one over another
    #             TFN = np.array([1, 1, 1 ,2/3, 1, 3/2, 3/2, 2, 5/2, 5/2, 3, 7/2, 7/2, 4, 9/2])
    #             for k in range(1, m_criteria + 1):
    #                 a = np.zeros((1, n*n*3))  # Comparison matrix (In this case is a vector because of computational memory
    #                 for i in range(k*n-(n-1), k*n + 1):
    #                     for j in range(i-n*(k-1), n + 1):
    #                         # This is the position of the third element of the segment for
    # 			            # a*(i,j) (upper triangular part)
    #                         jj = 3*(n*((i-n*(k-1))-1)+j)
    #                         # This is the position of the thrid element of the segment for
    # 			            # a*(j,i) (lower triangular part)
    #                         jjj = 3*(n*(j-1) + i-n*(k-1))
    #                         if N[i - 1][j - 1] == -4:
    #                             a[0][jjj-3:jjj] =  TFN[12:15]
    #                             a[0][jj-3:jj] = np.array([TFN[14]**-1, TFN[13]**-1, TFN[12]**-1])
    #                         elif N[i - 1][j - 1] == -3:
    #                             a[0][jjj-3:jjj] =  TFN[9:12]
    #                             a[0][jj-3:jj] = np.array([TFN[11]**-1, TFN[10]**-1, TFN[9]**-1])
    #                         elif N[i - 1][j - 1] == -2:
    #                             a[0][jjj-3:jjj] = TFN[6:9]
    #                             a[0][jj-3:jj] = np.array([TFN[8]**-1, TFN[7]**-1, TFN[6]**-1])
    #                         elif N[i - 1][j - 1] == -1:
    #                             a[0][jjj-3:jjj] =  TFN[3:6]
    #                             a[0][jj-3:jj] = np.array([TFN[5]**-1, TFN[4]**-1, TFN[3]**-1])
    #                         elif N[i - 1][j - 1] == 0:
    #                             a[0][jj-3:jj] =  TFN[0:3]
    #                             a[0][jjj-3:jjj] = TFN[0:3]
    #                         elif N[i - 1][j - 1] == 1:
    #                             a[0][jj-3:jj] =  TFN[3:6]
    #                             a[0][jjj-3:jjj] = np.array([TFN[5]**-1, TFN[4]**-1, TFN[3]**-1])
    #                         elif N[i - 1][j - 1] == 2:
    #                             a[0][jj-3:jj] = TFN[6:9]
    #                             a[0][jjj-3:jjj] = np.array([TFN[8]**-1, TFN[7]**-1, TFN[6]**-1])
    #                         elif N[i - 1][j - 1] == 3:
    #                             a[0][jj-3:jj] =  TFN[9:12]
    #                             a[0][jjj-3:jjj] = np.array([TFN[11]**-1, TFN[10]**-1, TFN[9]**-1])
    #                         elif N[i - 1][j - 1] == 4:
    #                             a[0][jj-3:jj] =  TFN[12:15]
    #                             a[0][jjj-3:jjj] = np.array([TFN[14]**-1, TFN[13]**-1, TFN[12]**-1])
    #                 # (2) fuzzy synthetic extension
    #                 A = np.zeros((n,3))
    #                 B = np.zeros((1,3))
    #                 for i in range(1, n + 1):
    #                     for j in range(1, n + 1):
    #                         jj = 3*(n*(i-1)+j)
    #                         A[i - 1][:] = A[i - 1][:] + a[0][jj-3:jj]
    #                     B = B + A[i - 1][:]
    #                 BB = np.array([B[0][2]**-1, B[0][1]**-1, B[0][0]**-1])
    #                 S = A*BB
    #                 # (3) Degree of possibility
    #                 for i in range(n):
    #                     V = np.zeros(n)
    #                     for j in range(n):
    #                         if S[i][1] >= S[j][1]:
    #                             V[j] = 1
    #                         elif S[j][0] >= S[i][2]:
    #                             V[j] = 0
    #                         else:
    #                             V[j] = (S[j][0] - S[i][2])/((S[i][1] - S[i][2]) - (S[j][1] - S[j][0]))
    #                     w[k - 1][i] = np.min(V)
    #                 # (4) Weight of each flow for a criterium
    #                 w[k - 1][:] = (np.sum(w[k - 1][:])**-1)*w[k - 1][:]
    #                 # (5) Criteria' uncertainty degrees
    #                 for i in range(n):
    #                     if w[k - 1][i] != 0:
    #                         eps[k - 1][0] = eps[k - 1][0] - w[k - 1][i]*np.log(w[k - 1][i])
    #                 eps[k - 1][0] = eps[k - 1][0]/np.log(n)
    #                 phi[k - 1][0] = 1 + eps[k - 1]
    #                 sumphi = sumphi + phi[k - 1][0]
    #             # (6) Final weight of all flows
    #             for i in range(n):
    #                 for k in range(m_criteria):
    #                     theta[0][k] = phi[k][0]/sumphi
    #                     W[0][i] = W[0][i] + w[k][i]*theta[0][k]
    #             W = (np.sum(W)**-1)*W # Final weights
    #             if self.Year <= 2004:
    #                 W = np.array([0.0 if i in indices else W[0][indices_f[i]] for i in range(n_actual)])
    #         else:
    #             W = np.array([0.0 if i in indices else 1.0 for i in range(n_actual)])
    #         df = df.assign(WEIGHT = W)
    #         df['MANAGED FLOW'] = df.apply(lambda x: x['MANAGED FLOW']*x['WEIGHT'], axis = 1)
    #         df.drop(columns = 'WEIGHT', inplace = True)
    #         return df


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
            df.drop(columns = df.iloc[:, list(range(17, 70, 13))].columns.tolist(), inplace = True)
        else:
            df.drop(columns = df.iloc[:, list(range(19, 72, 13))].columns.tolist(), inplace = True)
        df_PCUs = pd.DataFrame()
        Columns_0 = list(df.iloc[:, 0:7].columns)
        for i in range(5):
            Starting = 7 + 12*i
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
        cols =  list(df_PCUs.iloc[:, 8:16].columns)
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
        # On-site energy recovery
        df = dfs['1a'].iloc[:, list(range(11))]
        cols = [c for c in df.columns if 'METHOD' in c]
        df.dropna(subset = cols, how = 'all', axis = 0, inplace = True)
        Columns_0 = list(df.iloc[:, 0:7].columns)
        Columns_1 = list(df.iloc[:, 7:].columns)
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
        # On-site recycling
        df = dfs['1a'].iloc[:, list(range(7)) + list(range(11,18))]
        cols = [c for c in df.columns if 'METHOD' in c]
        df.dropna(subset = cols, how = 'all', axis = 0, inplace = True)
        Columns_0 = list(df.iloc[:, 0:7].columns)
        Columns_1 = list(df.iloc[:, 7:].columns)
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
        # Changing units
        df_PCUs = df_PCUs.loc[(df_PCUs.iloc[:,0:] != 'INV').all(axis = 1)]
        df_PCUs.dropna(how = 'all', axis = 0, inplace = True)
        # Adding methods name
        Methods = pd.read_csv(self._dir_path + '/Ancillary/Methods_TRI.csv',
                            usecols = ['Code 2004 and prior',
                                    'Method 2004 and prior',
                                    'Code 2005 and after',
                                    'Method 2005 and after'])
        Methods.drop_duplicates(keep =  'first', inplace = True)
        columns_DB_F = ['REPORTING YEAR', 'TRIFID', 'PRIMARY NAICS CODE', 'CAS NUMBER',
                         'CHEMICAL NAME', 'METAL INDICATOR', 'CLASSIFICATION',
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
            codes_incineration = ['H040', 'H076', 'H122']
        else:
            df.drop(columns = df.iloc[:, list(range(13, 62, 12))].columns.tolist(), inplace = True)
            codes_incineration = ['F01', 'F11', 'F19', 'F31',
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
        phases = ['W', 'L']
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
                    #continue
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
        PCU = pd.read_csv(self._dir_path + '/Datasets/PCUs_DB_' + str(self.Year) + '.csv',
                            low_memory = False,
                            converters = {'CAS NUMBER': lambda x: x if re.search(r'^[A-Z]', x) else str(int(x))})
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
        PCU_energy =  PCU_recycling.loc[~ ((PCU_recycling['METHOD CODE - 2005 AND AFTER'] == 'H20') & (PCU_recycling['CAS NUMBER'].isin(Solvent_recovery)))]
        PCU_recycling.reset_index(inplace = True, drop = True)
        PCU_recycling['BASED ON OPERATING DATA?'] = 'NO'
        efficiency_estimation = [p[0] for p in np.random.uniform(75, 99.5, size=(PCU_recycling.shape[0],1))]
        efficiency_range = [self._efficiency_estimation_to_range(x) for x in efficiency_estimation]
        PCU_recycling['EFFICIENCY RANGE CODE'] = pd.Series(efficiency_range)
        PCU_recycling = \
                 PCU_recycling.apply(lambda x: \
                 self._phase_estimation_recycling(Statistics, x), axis = 1)
        PCU_recycling = PCU_recycling.loc[pd.notnull(PCU_recycling['WASTE STREAM CODE'])]
        if self.Year <= 2004:
            PCU_recycling['EFFICIENCY ESTIMATION'] = pd.Series(efficiency_estimation)
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
        # Energy recovery
        PCU_energy = PCU.loc[PCU['TYPE OF MANAGEMENT'] == 'Energy recovery']
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
                    self._energy_efficiency(Statistics, x), axis = 1)
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
        print(PCU_energy.info())
        df_N_PCU = pd.concat([df_N_PCU, PCU_energy],
                                 ignore_index = True,
                                 sort = True, axis = 0)
        Chemicals_to_remove = ['MIXTURE', 'TRD SECRT']
        df_N_PCU = df_N_PCU.loc[~df_N_PCU['CAS NUMBER'].isin(Chemicals_to_remove)]
        df_N_PCU['CAS NUMBER'] = df_N_PCU['CAS NUMBER'].apply(lambda x: str(int(x)) if not 'N' in x else x)
        if self.Year <= 2004:
            columns_DB_F = ['REPORTING YEAR', 'TRIFID', 'PRIMARY NAICS CODE', 'CAS NUMBER',
                             'CHEMICAL NAME', 'METAL INDICATOR', 'CLASSIFICATION',
                             'WASTE STREAM CODE', 'RANGE INFLUENT CONCENTRATION',
                             'METHOD CODE - 2004 AND PRIOR', 'METHOD NAME - 2004 AND PRIOR',
                             'METHOD CODE - 2005 AND AFTER', 'METHOD NAME - 2005 AND AFTER',
                             'TYPE OF MANAGEMENT', 'EFFICIENCY RANGE CODE',
                             'EFFICIENCY ESTIMATION', 'BASED ON OPERATING DATA?']
        else:
            columns_DB_F = ['REPORTING YEAR', 'TRIFID', 'PRIMARY NAICS CODE', 'CAS NUMBER',
                            'CHEMICAL NAME', 'METAL INDICATOR', 'CLASSIFICATION',
                            'WASTE STREAM CODE', 'METHOD CODE - 2005 AND AFTER',
                            'METHOD NAME - 2005 AND AFTER', 'TYPE OF MANAGEMENT',
                            'EFFICIENCY RANGE CODE', 'BASED ON OPERATING DATA?']
        df_N_PCU = df_N_PCU[columns_DB_F]
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
                        [C]: Further cleaning of database.', \
                        type = str)

    parser.add_argument('-Y', '--Year',
                        help = 'Records with up to how many PCUs you want to include?.',
                        type = str,
                        required = False,
                        default = 2004)


    args = parser.parse_args()
    Year = args.Year

    start_time =  time.time()
    Building = PCU_DB(int(Year))

    if args.Option == 'A':
        Building.organizing()
    elif args.Option == 'B':
        Building.Building_database_for_statistics()
    elif args.Option == 'C':
        Building.cleaning_database()

    print('Execution time: %s sec' % (time.time() - start_time))
