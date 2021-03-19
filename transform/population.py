# population.py (PAU4Chem)
# !/usr/bin/env python3
# coding=utf-8


"""TODO describe this file/module."""

import pandas as pd
import numpy as np
from scipy.stats import lognorm
import os
import bisect
from building_pau_db import *
import re


def Calling_US_census(dir_path):
    """TODO describe this method.

    Args:
        dir_path (str): [description]

    Returns:
        [type]: [description]
    """
    # 2008 Annual Survey of Manufactures (ASM):
    # Link: https://www.census.gov/data/tables/2008/econ/asm/2008-asm.html
    path_ASM_2008 = dir_path + '/us_census_bureau/ASM_2008.xlsx'
    df_ASM_2008 = pd.read_excel(
        path_ASM_2008, sheet_name='ASM_2008',
        usecols=['NAICS code',
                 'Year', 'Total value of shipments ($1,000)',
                 'Relative standard error of total value of shipments (%)'],
        dtype={'NAICS code': 'object'})
    df_ASM_2008 = df_ASM_2008[df_ASM_2008['Year'] == 2008]
    df_ASM_2008.drop(columns=['Year'], inplace=True)
    df_ASM_2008['NAICS code'] = df_ASM_2008['NAICS code'].apply(
        lambda x: str(x).strip())
    # 2008 SUSB Annual Datasets by Establishment Industry
    # Link: https://www.census.gov/data/datasets/2008/econ/susb/2008-susb.html
    path_SUSB_2008 = dir_path + '/us_census_bureau/SUSB_2008.csv'
    df_SUSB_2008 = pd.read_csv(
        path_SUSB_2008, usecols=['NAICS code', 'ESTB', 'ENTRSIZEDSCR'],
        dtype={'NAICS code': 'object'})
    df_SUSB_2008['NAICS code'] = df_SUSB_2008['NAICS code'].apply(
        lambda x: str(x).strip())
    df_SUSB_2008 = df_SUSB_2008[df_SUSB_2008['NAICS code'].str.contains(
        r'^3[123]', na=False)]
    df_SUSB_2008 = df_SUSB_2008[df_SUSB_2008['ENTRSIZEDSCR'] == 'Total']
    df_SUSB_2008.drop(columns=['ENTRSIZEDSCR'], inplace=True)
    # Merging
    Merged = pd.merge(df_SUSB_2008, df_ASM_2008, on='NAICS code', how='inner')
    Merged[[
        'Total value of shipments ($1,000)',
        'Relative standard error of total value of shipments (%)']] = Merged[[
            'Total value of shipments ($1,000)',
            'Relative standard error of total value of shipments (%)'
        ]].applymap(float)

    Merged[['Mean value of shipments ($1,000)',
            'SD value of shipments ($1,000)']] = Merged.apply(
                lambda x: pd.Series([
                    x.values[2]/x.values[1],
                    x.values[3]*x.values[2]/(100*x.values[1]**0.5)]),
                axis=1)
    Merged = Merged[['NAICS code',
                     'Mean value of shipments ($1,000)',
                     'SD value of shipments ($1,000)']]
    return Merged


def Probability_establishments_within_cluster(naics, establishment, df):
    """TODO describe this method.

    Args:
        naics ([type]): [description]
        establishment ([type]): [description]
        df ([type]): [description]

    Returns:
        [type]: [description]
    """
    values = {'Sector': 2,
              'Subsector': 3,
              'Industry Group': 4,
              'NAICS Industry': 5}
    df_interest = df.loc[df['NAICS code'] == naics]
    if df_interest.empty:
        PAU_class = PAU_DB(2008)
        df['NAICS structure'] = df['NAICS code'].apply(
            lambda x: PAU_class._searching_naics(x, naics))
        df['NAICS structure'] = df['NAICS structure'].map(values)
        Max = df['NAICS structure'].max()
        df_interest = df[df['NAICS structure'] == Max]
    mean = df_interest['Mean value of shipments ($1,000)'].iloc[0]
    sd = df_interest['SD value of shipments ($1,000)'].iloc[0]
    # measure-of-size (MOS) (e.g., value of shipments,
    #                        number of employees, etc.),
    # which was highly correlated with pollution abatement operating costs
    # Method of moments
    mu = np.log(mean**2/(sd**2 + mean**2)**0.5)
    theta_2 = np.log(sd**2/mean**2 + 1)
    MOS = lognorm.rvs(s=theta_2**0.5,
                      scale=np.exp(mu),
                      size=int(establishment))
    Best = max(MOS)
    # For avoiding 0 probability
    Worst = min(MOS) - 10**(np.log10(min(MOS)) - 2)
    MOS.sort()
    # High values of MOS represent a possible high value of PAA.
    # Establishments with high values of PAOC and PACE had a
    # probability of 1 of being selected
    MOS_std = {
        str(idx + 1): [(val - Worst)/(Best - Worst), val*10**3]
        for idx, val in enumerate(MOS)}
    return MOS_std


def Probability_cluster_being_sampled(naics, establishment,
                                      total_establishments, n_clusters,
                                      df_census, df_tri):
    """TODO describe this method.

    Args:
        naics ([type]): [description]
        establishment ([type]): [description]
        total_establishments ([type]): [description]
        n_clusters ([type]): [description]
        df_census ([type]): [description]
        df_tri ([type]): [description]

    Returns:
        [type]: [description]
    """
    values = {'Sector': 2,
              'Subsector': 3,
              'Industry Group': 4,
              'NAICS Industry': 5}
    df_interest = df_tri.loc[df_tri['NAICS code'] == naics]
    if df_interest.empty:
        PAU_class = PAU_DB(2008)
        df_tri['NAICS structure'] = df_tri['NAICS code'].apply(
            lambda x: PAU_class._searching_naics(x, naics))
        df_tri['NAICS structure'] = df_tri['NAICS structure'].map(values)
        Max = df_tri['NAICS structure'].max()
        df_interest = df_tri[df_tri['NAICS structure'] == Max]
        P_cluster = df_interest['% establishments without PAA'].mean()
    else:
        P_cluster = df_interest['% establishments without PAA'].iloc[0]
    Pro_establishment = Probability_establishments_within_cluster(
        naics, establishment, df_census)
    Pro_establishment_accumulated = dict()
    Shipment_value_establishment = dict()
    sum = 0.0
    for key, val in Pro_establishment.items():
        sum = sum + val[0]
        Pro_establishment_accumulated.update({key: [sum, val[0]]})
        Shipment_value_establishment.update({key: val[1]})
    return pd.Series([P_cluster/100, Pro_establishment_accumulated,
                      Shipment_value_establishment])


def calling_TRI_for_prioritization_sectors(dir_path):
    """TODO describe this method.

    Args:
        dir_path ([type]): [description]

    Returns:
        [type]: [description]
    """
    # The survey prioritized the clusters based on PACE 1994
    columns = ['ON-SITE - TOTAL WASTE MANAGEMENT',
               'ON-SITE - TOTAL LAND RELEASES',
               'PRIMARY NAICS CODE', 'TRIFID', 'UNIT OF MEASURE']
    df_TRI_1994 = pd.read_csv(dir_path + '/../extract/datasets/US_1a_1994.csv',
                              low_memory=False, usecols=Columns,
                              dtype={'PRIMARY NAICS CODE': 'object'})
    df_TRI_1994 = df_TRI_1994[df_TRI_1994[
        'PRIMARY NAICS CODE'].str.contains(r'^3[123]', na=False)]
    Flow_columns = ['ON-SITE - TOTAL WASTE MANAGEMENT',
                    'ON-SITE - TOTAL LAND RELEASES']
    df_TRI_1994.loc[df_TRI_1994[
        'UNIT OF MEASURE'] == 'Pounds', Flow_columns] *= 0.453592
    df_TRI_1994.loc[df_TRI_1994[
        'UNIT OF MEASURE'] == 'Grams', Flow_columns] *= 10**-3
    df_TRI_1994 = df_TRI_1994.groupby(
        ['TRIFID', 'PRIMARY NAICS CODE'], as_index=False).sum()
    df_TRI_1994['Any pollution abatement?'] = 'No'
    df_TRI_1994.loc[(df_TRI_1994[Flow_columns] != 0.0).any(axis=1),
                    'Any pollution abatement?'] = 'Yes'
    df_TRI_1994.drop(columns=['TRIFID'] + Flow_columns,
                     inplace=True)
    df_TRI_1994.rename(columns={'PRIMARY NAICS CODE': 'NAICS code'},
                       inplace=True)
    df_TRI_1994['Number of establishments'] = 1
    df_TRI_1994 = df_TRI_1994.groupby(
        ['Any pollution abatement?', 'NAICS code'], as_index=False).sum()
    df_TRI_1994['Number of establishments in cluster'] = df_TRI_1994.groupby(
        'NAICS code', as_index=False)[
            'Number of establishments'].transform('sum')
    df_TRI_1994 = df_TRI_1994[df_TRI_1994['Any pollution abatement?'] == 'No']
    df_TRI_1994.drop(columns=['Any pollution abatement?'], inplace=True)
    df_TRI_1994['% establishments without PAA'] = df_TRI_1994[
        ['Number of establishments', 'Number of establishments in cluster']
    ].apply(lambda x: 100*x.values[0]/x.values[1], axis=1)
    return df_TRI_1994


def Organizing_sample(n_sampled_establishments, dir_path):
    """TODO describe this method.

    Args:
        n_sampled_establishments ([type]): [description]
        dir_path ([type]): [description]

    Returns:
        [type]: [description]
    """
    sampled_clusters = pd.read_csv(
        dir_path + '/us_census_bureau/Selected_clusters_2005.txt',
        header=None, index_col=False)
    sampled_clusters = [str(val) for val in sampled_clusters.iloc[:, 0]]
    n_sampled_clusters = len(sampled_clusters)

    # Statistics of U.S. Businesses - Survey 2005
    # Note: 1. The PAOC and PACE only have information for establishments with
    #          greather or equal to 20 employees
    #       2. The PAOC and PACE are on establishments
    #       3. The industry sectors surveyed were NAICS codes 31-33
    # Source: https://www.census.gov/prod/2008pubs/ma200-05.pdf

    df_SUSB_2005 = pd.read_csv(
        dir_path + '/us_census_bureau/Statistics_of_US_businesses_2004.csv',
        low_memory=False, header=None, usecols=[1, 4, 11],
        names=['NAICS code', 'Establishments (employees >= 20)',
               'Employment size'])
    df_SUSB_2005 = df_SUSB_2005[df_SUSB_2005[
        'NAICS code'].str.contains(r'^3[123]')]
    df_SUSB_2005['Establishments (employees >= 20)'] = pd.to_numeric(
        df_SUSB_2005['Establishments (employees >= 20)'], errors='coerce')
    df_SUSB_2005 = df_SUSB_2005[pd.notnull(df_SUSB_2005[
        'Establishments (employees >= 20)'])]
    df_SUSB_2005['Establishments (employees >= 20)'] = \
        df_SUSB_2005['Establishments (employees >= 20)'].astype('int')
    row_names = ['20-99 employees', '100-499 employees', '500 + employees']
    df_SUSB_2005 = df_SUSB_2005[df_SUSB_2005[
        'Employment size'].isin(row_names)]
    df_SUSB_2005.drop(columns=['Employment size'],
                      inplace=True)
    df_SUSB_2005 = df_SUSB_2005.groupby('NAICS code', as_index=False).sum()
    df_SUSB_2005 = df_SUSB_2005.loc[
        df_SUSB_2005['NAICS code'].isin(sampled_clusters)]
    N_total_establishments = df_SUSB_2005[
        'Establishments (employees >= 20)'].sum()
    # Calling information from 1994 TRI
    df_TRI_1994 = calling_TRI_for_prioritization_sectors(dir_path)
    # Calling information from 2008 census
    df_CENSUS_2008 = Calling_US_census(dir_path)
    df_SUSB_2005[['P-cluster', 'P-establishment',
                  'Shipment-establishment']] = df_SUSB_2005.apply(
                      lambda x: Probability_cluster_being_sampled(
                          x.values[0], x.values[1], N_total_establishments,
                          n_sampled_clusters, df_CENSUS_2008, df_TRI_1994),
                      axis=1)
    df_SUSB_2005.sort_values(by=['P-cluster'], inplace=True)
    df_SUSB_2005 = df_SUSB_2005.reset_index()
    df_SUSB_2005['P-cluster accumulated'] = df_SUSB_2005['P-cluster'].cumsum()
    Accumulated = df_SUSB_2005['P-cluster accumulated'].max()
    List_aux = df_SUSB_2005['P-cluster accumulated'].tolist()
    NAICS_list = list()
    P_cluster_list = list()
    Establishment_list = list()
    MOS_list = list()
    P_selected_list = list()
    Shipment_list = list()
    n_sampled = 0
    while n_sampled < n_sampled_establishments:
        P_selected = 0
        # The PACE survey demanded a minimum probability of 0.05
        while P_selected <= 0.05:
            rnd_cluster = np.random.uniform(0, Accumulated)
            idx = bisect.bisect_left(List_aux, rnd_cluster)
            naics = df_SUSB_2005['NAICS code'].iloc[idx]
            P_cluster = df_SUSB_2005['P-cluster'].iloc[idx]
            List_estab_value = [val[0] for val in list(
                df_SUSB_2005['P-establishment'].iloc[idx].values())]
            List_estab_MOS = [val[1] for val in list(
                df_SUSB_2005['P-establishment'].iloc[idx].values())]
            List_estab_shipment = [val for val in list(
                df_SUSB_2005['Shipment-establishment'].iloc[idx].values())]
            List_estab_number = list(
                df_SUSB_2005['P-establishment'].iloc[idx].keys())
            rnd_establishment = np.random.uniform(0, max(List_estab_value))
            pos = bisect.bisect_left(List_estab_value, rnd_establishment)
            key = List_estab_number[pos]
            MOS = List_estab_MOS[pos]
            Shipment = List_estab_shipment[pos]
            # P(selected) = P(select in cluster)*P(cluster be selected)
            # they are independent events
            P_selected = MOS*P_cluster
        Pro_establishment_accumulated = dict()
        Shipment_value_establishment = dict()
        sum = 0.0
        for p, v in enumerate(List_estab_MOS):
            if p != pos:
                sum = sum + v
                Pro_establishment_accumulated.update(
                    {List_estab_number[p]: [sum, v]})
                Shipment_value_establishment.update(
                    {List_estab_number[p]: List_estab_shipment[p]})

        df_SUSB_2005.iloc[idx]['P-establishment'] = \
            Pro_establishment_accumulated
        df_SUSB_2005.iloc[idx]['Shipment-establishment'] = \
            Shipment_value_establishment

        NAICS_list.append(naics)
        P_cluster_list.append(P_cluster)
        P_selected_list.append(P_selected)
        Establishment_list.append(key)
        MOS_list.append(MOS)
        Shipment_list.append(Shipment)
        n_sampled = n_sampled + 1

    df_result = pd.DataFrame({'NAICS code': NAICS_list,
                              'P-cluster': P_cluster_list,
                              'Establishment': Establishment_list,
                              'P-in-cluster': MOS_list,
                              'P-selected': P_selected_list,
                              'Total shipment establishment': Shipment_list})

    # Inflation rate from 2008 to 2020 is 19.09%:
    df_result['Total shipment establishment'] = df_result[
        'Total shipment establishment']*1.1909
    df_result['Unit'] = 'USD'
    print('Total value of shipments for NAICS 31-33 ' +
          'in 2008 5,486,265,506,000 USD')
    print('Contibution of the selected establishments to the total ' +
          f'value of shipments [%]: ', round(100*df_result[
              'Total shipment establishment'].sum()/(5486265506000*1.1909), 4))
    return df_result


def searching_census(naics, media, activity, df):
    """TODO describe this method.

    Args:
        naics ([type]): [description]
        media ([type]): [description]
        activity ([type]): [description]
        df ([type]): [description]

    Returns:
        [type]: [description]
    """
    values = {'Sector': 2,
              'Subsector': 3,
              'Industry Group': 4,
              'NAICS Industry': 5}
    df = df.loc[df['NAICS code'].apply(
        lambda x: True if len(x) <= 5 else False)]
    PAU_class = PAU_DB(2008)
    df['NAICS structure'] = df['NAICS code'].apply(
        lambda x: PAU_class._searching_naics(x, naics))
    df['NAICS structure'] = df['NAICS structure'].map(values)
    Max = df['NAICS structure'].max()
    if Max == 2:
        df = df[df['NAICS code'] == '31–33']
    else:
        df = df[df['NAICS structure'] == Max]
        df['Length'] = df['NAICS code'].str.len()
        df = df.loc[df.groupby(['Activity', 'Media'])['Length'].idxmin()]
    df = df.loc[
        (df['Activity'] == activity) & (df['Media'] == media)]
    dictionary = {col1: col2 for col1 in [
        'RSE', ' activity & media', 'P-media', 'P-activity',
        'Total PA', 'Total s', 'Info'] for col2 in df.columns if col1 in col2}
    s = pd.Series([df[col].values for col in dictionary.values()])
    return s


def mean_standard(establishments, shipment_flow, rse,
                  total, shipment, confidence):
    """TODO describe this method.

    Args:
        establishments ([type]): [description]
        shipment_flow ([type]): [description]
        rse ([type]): [description]
        total ([type]): [description]
        shipment ([type]): [description]
        confidence ([type]): [description]

    Returns:
        [type]: [description]
    """
    if total >= shipment:
        ratio = 0.0044
    else:
        ratio = total/shipment
    try:
        total = ratio*shipment_flow
        Mean = total/establishments
        SD = (rse*total/(100*(establishments)**0.5))
        CI = [Mean - confidence*SD/(establishments)**0.5,
              Mean + confidence*SD/(establishments)**0.5]
        return pd.Series([total, Mean, SD, CI])
    except KeyError:
        return pd.Series([None]*4)


def selecting_establishment_by_activity_and_media(info_establishments,
                                                  probability_activity,
                                                  probability_media):
    """TODO describe this method.

    Args:
        info_establishments ([type]): [description]
        probability_activity ([type]): [description]
        probability_media ([type]): [description]

    Returns:
        [type]: [description]
    """
    # P(activity and media and establishment) =
    # P(activity and media) * P(establishment). They are independent events
    selected_establishments = 0
    Info_selected_establishments = dict()
    for estab, info in info_establishments.items():
        rnd_a = np.random.rand()
        if probability_activity > rnd_a:
            rnd_m = np.random.rand()
            if probability_media > rnd_m:
                probability_establishment = info[1]
                rnd_establishment = np.random.rand()
                if probability_establishment > rnd_establishment:
                    selected_establishments = selected_establishments + 1
                    Info_selected_establishments.update({estab: info})

    not_zero = probability_activity * probability_media != 0
    if not_zero and selected_establishments == 0:
        selected_establishments = 1
        maximum = max([val[0] for val in info_establishments.values()])
        Info_selected_establishments.update(
            {key: val for key, val in info_establishments.items()
             if val[0] == maximum})

    return pd.Series([selected_establishments, Info_selected_establishments])


def estimating_mass_by_activity_and_media(mu, theta_2, info_establishments,
                                          P_media, P_activity):
    """TODO describe this method.

    Args:
        mu ([type]): [description]
        theta_2 ([type]): [description]
        info_establishments ([type]): [description]
        P_media ([type]): [description]
        P_activity ([type]): [description]

    Returns:
        [type]: [description]
    """
    prob = [1 - vals[1]*P_activity*P_media
            for vals in info_establishments.values()]  # 1 - cdf
    shipments = [vals[0] for vals in info_establishments.values()]
    try:
        vals = lognorm.isf(prob,
                           s=theta_2**0.5,
                           scale=np.exp(mu))
    except ValueError:
        vals = lognorm.isf(prob,
                           s=10**-9,
                           scale=np.exp(mu))
    shipment_mass = sum([shipments[idx]/val for idx, val in enumerate(vals)])
    total_mass = sum(vals)
    total_shipment = sum(shipments)
    return pd.Series([total_mass, total_shipment, shipment_mass])


def searching_establishments_by_hierarchy(naics, df):
    """TODO describe this method.

    Args:
        naics ([type]): [description]
        df ([type]): [description]

    Returns:
        [type]: [description]
    """
    try:
        if (naics == '31–33') | (len(naics) == 3):
            naics = '3[123]'
            n = '2,4'
        elif len(naics) == 4:
            n = '1,2'
        elif len(naics) == 5:
            n = '1'
        regex = re.compile(r'^%s[0-9]{%s}' % (naics, n))
        df_result = df.loc[df['NAICS code'].apply(
            lambda x: True if re.match(regex, x) else False)]
        if not df_result.empty:
            estab = int(df_result['Establishments (employees >= 20)'].sum())
            shipment = df_result['Total shipment'].sum()
            Dictionary_establishments = dict()
            i = 0
            for idx, row in df_result.iterrows():
                for values in row['Info establishments'].values():
                    i = i + 1
                    Dictionary_establishments.update({str(i): values})
            return pd.Series([estab, shipment, Dictionary_establishments])
        else:
            return pd.Series([None]*3)
    except UnboundLocalError:
        return pd.Series([None]*3)


def normalizing_shipments(df):
    """TODO describe this method.

    Args:
        df ([type]): [description]
    """

    def function_aux(info_estab, P_m, P_a, P_on, estab, Total):
        vals = info_estab[estab]
        info_estab.update({estab: [vals[0]*P_m*P_a*P_on*0.01/Total, vals[1]]})
        return info_estab

    estabs = set(
        [key for idx, row in df.iterrows() for key in
         row['Info probable establishments'].keys()])
    for estab in estabs:
        idx = df.loc[df['Info probable establishments'].apply(
            lambda x: True if estab in x.keys() else False)].index.tolist()
        Total = (df.loc[idx, 'P-media']*df.loc[idx, 'P-activity']).sum()

        df.loc[idx]['Info probable establishments'] = df.loc[
            idx][['Info probable establishments', 'P-media', 'P-activity',
                  '% On-site flow']].apply(
                      lambda x: function_aux(x.values[0],
                                             x.values[1],
                                             x.values[2],
                                             x.values[3],
                                             estab,
                                             Total),
                      axis=1)
    return df


if __name__ == '__main__':
    dir_path = os.path.dirname(os.path.realpath(__file__))
    df_result = Organizing_sample(20378, dir_path)
    df_result.to_csv(
        dir_path + '/us_census_bureau/Sampled_establishments.csv',
        index=False)
