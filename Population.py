import pandas as pd
import numpy as np
from scipy.stats import norm
import os, bisect
from Building_PCUs_DB import *

def Calling_US_census(dir_path):
    # 2008 Annual Survey of Manufactures (ASM):
    # Link: https://www.census.gov/data/tables/2008/econ/asm/2008-asm.html
    path_ASM_2008 = dir_path + '/US_Census_Bureau/ASM_2008.xlsx'
    df_ASM_2008 = pd.read_excel(path_ASM_2008, sheet_name = 'ASM_2008',
                        usecols = ['NAICS code',
                                'Year', 'Total value of shipments ($1,000)',
                                'Relative standard error of total value of shipments (%)'],
                        dtype = {'NAICS code': 'object'})
    df_ASM_2008 = df_ASM_2008[df_ASM_2008['Year'] == 2008]
    df_ASM_2008.drop(columns = ['Year'], inplace = True)
    df_ASM_2008['NAICS code'] = df_ASM_2008['NAICS code'].apply(lambda x: str(x).strip())
    # 2008 SUSB Annual Datasets by Establishment Industry
    # Link: https://www.census.gov/data/datasets/2008/econ/susb/2008-susb.html
    path_SUSB_2008 =  dir_path + '/US_Census_Bureau/SUSB_2008.csv'
    df_SUSB_2008 = pd.read_csv(path_SUSB_2008, usecols = ['NAICS code', 'ESTB', 'ENTRSIZEDSCR'],
                       dtype = {'NAICS code': 'object'})
    df_SUSB_2008['NAICS code'] = df_SUSB_2008['NAICS code'].apply(lambda x: str(x).strip())
    df_SUSB_2008 = df_SUSB_2008[df_SUSB_2008['NAICS code'].str.contains(r'^3[123]', na = False)]
    df_SUSB_2008 = df_SUSB_2008[df_SUSB_2008['ENTRSIZEDSCR'] == 'Total']
    df_SUSB_2008.drop(columns = ['ENTRSIZEDSCR'], inplace = True)
    # Merging
    Merged = pd.merge(df_SUSB_2008, df_ASM_2008, on = 'NAICS code', how = 'inner')
    Merged[['Total value of shipments ($1,000)',
            'Relative standard error of total value of shipments (%)']] = \
                Merged[['Total value of shipments ($1,000)',
                        'Relative standard error of total value of shipments (%)']]\
                .applymap(float)
    Merged[['Mean value of shipments ($1,000)', 'SD value of shipments ($1,000)']] = \
                        Merged.apply(lambda x:
                            pd.Series([x.values[2]/x.values[1],
                                       x.values[3]*x.values[2]/(100*x.values[1]**0.5)]),
                            axis = 1)
    Merged = Merged[['NAICS code',
                     'Mean value of shipments ($1,000)',
                     'SD value of shipments ($1,000)']]
    return Merged


def Probability_establishments_within_cluster(naics, establishment, df):
    values = {'Sector': 2,
                'Subsector': 3,
                'Industry Group': 4,
                'NAICS Industry': 5}
    df_interest = df.loc[df['NAICS code'] == naics]
    if df_interest.empty:
        PCU_class = PCU_DB(2008)
        df['NAICS structure'] = df['NAICS code'].apply(lambda x: PCU_class._searching_naics(x, naics))
        df['NAICS structure'] = df['NAICS structure'].map(values)
        Max =  df['NAICS structure'].max()
        df_interest = df[df['NAICS structure'] == Max]
    mean = df_interest['Mean value of shipments ($1,000)'].iloc[0]
    sd = df_interest['SD value of shipments ($1,000)'].iloc[0]
    # measure-of-size (MOS) (e.g., value of shipments, number of employees, etc.),
    # which was highly correlated with pollution abatement operating costs
    MOS = norm.rvs(size = establishment, scale = sd, loc = mean)
    MOS = [M if M > 0 else abs(M - mean) + mean for M in MOS] # Reflection respect of the mean to avoid negative values
    Best = max(MOS)
    Worst = min(MOS)
    MOS.sort()
    MOS_std = {str(idx + 1):(val - Worst)/(Best - Worst) for idx, val in enumerate(MOS)}
    return MOS_std


def Probability_cluster_being_sampled(naics, establishment, total_establishments, n_clusters, df_census, df_tri):
    values = {'Sector': 2,
                'Subsector': 3,
                'Industry Group': 4,
                'NAICS Industry': 5}
    df_interest = df_tri.loc[df_tri['NAICS code'] == naics]
    if df_interest.empty:
        PCU_class = PCU_DB(2008)
        df_tri['NAICS structure'] = df_tri['NAICS code'].apply(lambda x: PCU_class._searching_naics(x, naics))
        df_tri['NAICS structure'] = df_tri['NAICS structure'].map(values)
        Max =  df_tri['NAICS structure'].max()
        df_interest = df_tri[df_tri['NAICS structure'] == Max]
        P_cluster = df_interest['% establishments without PAA'].mean()
    else:
        P_cluster = df_interest['% establishments without PAA'].iloc[0]
    Pro_establishment = Probability_establishments_within_cluster(naics, establishment, df_census)
    Pro_establishment_accumulated = dict()
    sum = 0.0
    for key, val in Pro_establishment.items():
        sum = sum + val
        Pro_establishment_accumulated.update({key: [sum, val]})
    return pd.Series([P_cluster/100, Pro_establishment_accumulated])


def calling_TRI_for_prioritization_sectors(dir_path):
    # The survey prioritized the clusters based on PACE 1994
    Columns = ['ON-SITE - TOTAL WASTE MANAGEMENT',
               #'OFF-SITE - TOTAL TRANSFERRED FOR FURTHER WASTE MANAGEMENT',
               #'OFF-SITE - TOTAL POTW TRANSFER',
               'ON-SITE - TOTAL LAND RELEASES',
               'PRIMARY NAICS CODE',
               'TRIFID', 'UNIT OF MEASURE']
    df_TRI_1994 = pd.read_csv(dir_path + '/Ancillary/US_1a_1994.csv',
                                low_memory = False, usecols =  Columns,
                                dtype = {'PRIMARY NAICS CODE': 'object'})
    df_TRI_1994 = df_TRI_1994[df_TRI_1994['PRIMARY NAICS CODE'].str.contains(r'^3[123]', na = False)]
    Flow_columns = ['ON-SITE - TOTAL WASTE MANAGEMENT',
               #'OFF-SITE - TOTAL TRANSFERRED FOR FURTHER WASTE MANAGEMENT',
               #'OFF-SITE - TOTAL POTW TRANSFER',
               'ON-SITE - TOTAL LAND RELEASES']
    df_TRI_1994.loc[df_TRI_1994['UNIT OF MEASURE'] == 'Pounds', Flow_columns] *= 0.453592
    df_TRI_1994.loc[df_TRI_1994['UNIT OF MEASURE'] == 'Grams', Flow_columns] *= 10**-3
    df_TRI_1994 = df_TRI_1994.groupby(['TRIFID', 'PRIMARY NAICS CODE'],
                                        as_index = False).sum()
    df_TRI_1994['Any pollution abatement?'] = 'No'
    df_TRI_1994.loc[(df_TRI_1994[Flow_columns] != 0.0).any(axis = 1),\
                    'Any pollution abatement?'] = 'Yes'
    df_TRI_1994.drop(columns = ['TRIFID'] + Flow_columns,
                    inplace = True)
    df_TRI_1994.rename(columns = {'PRIMARY NAICS CODE': 'NAICS code'},
                       inplace = True)
    df_TRI_1994['Number of establishments'] = 1
    df_TRI_1994 = df_TRI_1994.groupby(['Any pollution abatement?', 'NAICS code'],
                        as_index = False).sum()
    df_TRI_1994['Number of establishments in cluster'] = \
                df_TRI_1994.groupby('NAICS code',
                                    as_index =  False)['Number of establishments']\
                            .transform('sum')
    df_TRI_1994 = df_TRI_1994[df_TRI_1994['Any pollution abatement?'] == 'No']
    df_TRI_1994.drop(columns = ['Any pollution abatement?'], inplace = True)
    df_TRI_1994['% establishments without PAA'] = df_TRI_1994[['Number of establishments', \
                                                            'Number of establishments in cluster']] \
                                                .apply(lambda x: 100*x.values[0]/x.values[1],
                                                        axis = 1)
    return df_TRI_1994


def Organizing_sample(n_sampled_establishments, dir_path):
    sampled_clusters = pd.read_csv(dir_path + '/US_Census_Bureau/Selected_clusters_2005.txt',
                                    header = None, index_col = False)
    sampled_clusters = [str(val) for val in sampled_clusters.ix[:,0]]
    n_sampled_clusters = len(sampled_clusters)
    # Statistics of U.S. Businesses - Survey 2005
    # Note: 1. The PAOC and PACE only have information for establishments with greather or equal to 20 employees
    #       2. The PAOC and PACE are on establishments
    #       3. The industry sectors surveyed were NAICS codes 31-33
    # Source: https://www.census.gov/prod/2008pubs/ma200-05.pdf
    df_SUSB_2005 = pd.read_csv(dir_path + '/US_Census_Bureau/Statistics_of_US_businesses_2004.csv',
                    low_memory = False, header = None,
                    usecols = [1, 4, 11],
                    names = ['NAICS code', 'Establishments (employees >= 20)', 'Employment size'])
    df_SUSB_2005 = df_SUSB_2005[df_SUSB_2005['NAICS code'].str.contains(r'^3[123]')]
    df_SUSB_2005['Establishments (employees >= 20)'] = pd.to_numeric( \
                            df_SUSB_2005['Establishments (employees >= 20)'],
                            errors='coerce')
    df_SUSB_2005 = df_SUSB_2005[pd.notnull(df_SUSB_2005['Establishments (employees >= 20)'])]
    df_SUSB_2005['Establishments (employees >= 20)'] = \
                df_SUSB_2005['Establishments (employees >= 20)'].astype('int')
    row_names = ['20-99 employees', '100-499 employees', '500 + employees']
    df_SUSB_2005 = df_SUSB_2005[df_SUSB_2005['Employment size'].isin(row_names)]
    df_SUSB_2005.drop(columns = ['Employment size'],
                inplace = True)
    df_SUSB_2005 = df_SUSB_2005.groupby('NAICS code', as_index = False).sum()
    df_SUSB_2005 = df_SUSB_2005.loc[df_SUSB_2005['NAICS code'].isin(sampled_clusters)]
    N_total_establishments = df_SUSB_2005['Establishments (employees >= 20)'].sum()
    # Calling information from 1994 TRI
    df_TRI_1994 = calling_TRI_for_prioritization_sectors(dir_path)
    # Calling information from 2008 census
    df_CENSUS_2008 = Calling_US_census(dir_path)
    df_SUSB_2005[['P-cluster', 'P-establishment']] = \
            df_SUSB_2005.apply(lambda x: Probability_cluster_being_sampled(
                                            x.values[0], x.values[1],
                                            N_total_establishments,
                                            n_sampled_clusters,
                                            df_CENSUS_2008,
                                            df_TRI_1994),
                               axis = 1)
    #df_SUSB_2005['P-cluster'] = 1 - df_SUSB_2005['P-cluster']
    df_SUSB_2005.sort_values(by = ['P-cluster'], inplace = True)
    df_SUSB_2005 = df_SUSB_2005.reset_index()
    df_SUSB_2005['P-cluster accumulated'] = df_SUSB_2005['P-cluster'].cumsum()
    Accumulated = df_SUSB_2005['P-cluster accumulated'].max()
    List_aux = df_SUSB_2005['P-cluster accumulated'].tolist()
    NAICS_list = list()
    P_cluster_list = list()
    Establishment_list = list()
    MOS_list = list()
    n_sampled = 0
    while n_sampled < n_sampled_establishments:
        rnd_cluster = np.random.uniform(0, Accumulated)
        idx = bisect.bisect_left(List_aux, rnd_cluster)
        naics = df_SUSB_2005['NAICS code'].iloc[idx]
        P_cluster =  df_SUSB_2005['P-cluster'].iloc[idx]
        List_estab_value = [val[0] for val in list(df_SUSB_2005['P-establishment'].iloc[idx].values())]
        List_estab_MOS = [val[1] for val in list(df_SUSB_2005['P-establishment'].iloc[idx].values())]
        List_estab_number = list(df_SUSB_2005['P-establishment'].iloc[idx].keys())
        rnd_establishment = np.random.uniform(0, max(List_estab_value))
        pos = bisect.bisect_left(List_estab_value, rnd_establishment)
        key = List_estab_number[pos]
        MOS = List_estab_MOS[pos]
        Pro_establishment_accumulated = dict()
        sum = 0.0
        for p, v in enumerate(List_estab_MOS):
            if p != pos:
                sum = sum + v
                Pro_establishment_accumulated.update({List_estab_number[p]: [sum, v]})
        df_SUSB_2005.iloc[idx]['P-establishment'] = Pro_establishment_accumulated
        NAICS_list.append(naics)
        P_cluster_list.append(P_cluster)
        Establishment_list.append(key)
        MOS_list.append(MOS)
        n_sampled = n_sampled + 1
    df_result = pd.DataFrame({'NAICS code': NAICS_list,
                            'P-cluster': P_cluster_list,
                            'Establishment': Establishment_list,
                            'MOS': MOS_list})
    return df_result


def searching_census(naics, media, activity, df):
    values = {'Sector': 2,
             'Subsector': 3,
             'Industry Group': 4,
             'NAICS Industry': 5}
    df = df.loc[df['NAICS code'].apply(lambda x: True if len(x) <= 5 else False)]
    PCU_class = PCU_DB(2008)
    df['NAICS structure'] = df['NAICS code'].apply(lambda x: PCU_class._searching_naics(x, naics))
    df['NAICS structure'] = df['NAICS structure'].map(values)
    Max =  df['NAICS structure'].max()
    df = df[df['NAICS structure'] == Max]
    df =  df.loc[(df['Activity'] == activity) & \
                 (df['Media'] == media)]
    dictionary = {col1: col2 for col1 in ['RSE', 'establishments', 'Factor', \
                                        'Total'] for col2 in df.columns if col1 in col2}
    s =  pd.Series([df[col].values for col in dictionary.values()])
    return s


def mean_standard(x, inflation, confidence):
    try:
        establishments = x.iloc[6]
        factor = x.iloc[1]
        rse = x.iloc[4]
        total = x.iloc[5]*inflation
        Mean = total*factor/(establishments)*10**6
        SD = (rse*total*factor/(100*(establishments)**0.5))*10**6
        CI = [Mean - confidence*SD/(establishments)**0.5,
              Mean + confidence*SD/(establishments)**0.5]
        return pd.Series([Mean, SD, CI])
    except KeyError:
        return pd.Series([None]*3)


def selecting_establishment_by_activity_and_media(establishments, probability):
    selected_establishments = 0
    for i in range(int(establishments)):
        rnd = np.random.rand()
        if probability > rnd:
            selected_establishments = selected_establishments + 1
    if probability != 0 and selected_establishments == 0:
        selected_establishments = 1
    return selected_establishments


if __name__ == '__main__':
    dir_path = os.path.dirname(os.path.realpath(__file__))
    df_result = Organizing_sample(20378, dir_path)
    df_result.to_csv(dir_path + '/US_Census_Bureau/Sampled_establishments.csv', index = False)
