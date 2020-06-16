import pandas as pd
import numpy as np
from scipy.stats import norm
import os, bisect
from Building_PCUs_DB import PCU_DB

def Calling_US_census():
    dir_path = os.path.dirname(os.path.realpath(__file__))
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
    MOS = [M if M > 0 else abs(M - mean) + mean for M in MOS]
    MOS.sort()
    MOS = {str(idx + 1):val for idx, val in enumerate(MOS)}
    return MOS


def Probability_cluster_being_sampled(naics, establishment, total_establishments, n_clusters, df):
    P_cluster = n_clusters*establishment/total_establishments
    Pro_establishment = Probability_establishments_within_cluster(naics, establishment, df)
    Pro_establishment_accumulated = dict()
    sum = 0.0
    for key, val in Pro_establishment.items():
        sum = sum + val
        Pro_establishment_accumulated.update({key: [sum, val]})
    return pd.Series([P_cluster, Pro_establishment_accumulated])


def Organizing_sample(n_sampled_establishments):
    dir_path = os.path.dirname(os.path.realpath(__file__))
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
    # Calling information from 2008 census
    df_CENSUS_2008 = Calling_US_census()
    df_SUSB_2005[['P-cluster', 'P-establishment']] = \
            df_SUSB_2005.apply(lambda x: Probability_cluster_being_sampled(
                                            x.values[0], x.values[1],
                                            N_total_establishments,
                                            n_sampled_clusters,
                                            df_CENSUS_2008),
                               axis = 1)
    df_SUSB_2005.sort_values(by = ['P-cluster'], inplace = True,
                            ascending = True)
    df_SUSB_2005['P-cluster accumulated'] = df_SUSB_2005['P-cluster'].cumsum()
    List_aux = df_SUSB_2005['P-cluster accumulated'].tolist()
    NAICS_list = list()
    Establishment_list = list()
    MOS_list = list()
    for i in range(1, n_sampled_establishments + 1):
        rnd_cluster = np.random.uniform(0, n_sampled_clusters)
        idx = bisect.bisect_left(List_aux, rnd_cluster)
        naics = df_SUSB_2005['NAICS code'].iloc[idx]
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
        Establishment_list.append(key)
        MOS_list.append(MOS)
    df_result = pd.DataFrame({'NAICS code': NAICS_list,
                  'Establishment': Establishment_list,
                  'MOS': MOS_list})
    print(len(list(df_result['NAICS code'].unique())))
    return df_result


if __name__ == '__main__':
    dir_path = os.path.dirname(os.path.realpath(__file__))
    df_result = Organizing_sample(20378)
    df_result.to_csv(dir_path + '/US_Census_Bureau/Sampled_establishments.csv', index = False)
