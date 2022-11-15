# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 18:08:36 2022

@author: José Euclides de Melo Ferraz
"""
import math
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from sklearn.cluster import DBSCAN
from scipy.stats import t

#------------------------------------------------------------------------------
# Function to calculate haversine distance
#------------------------------------------------------------------------------
def haversine_distance(r, lat1, long1, lat2, long2, degree = True):
    """
    Calculates the distance between two points along a maximum circle of
    a spherical surface.

    Parameters
    ----------
    r : float
        Radius of the spherical surface.
    lat1 : float
        Latitude of the first first point.
    long1 : float
        Longitude of the first second point.
    lat2 : float
        Latitude of the second point.
    long2 : TYPE
        Longitude of the second point.
    degree : Bool
        If true, assume coordinates in degrees.
        If false, assume coordinates in radians.

    Returns
    -------
    d : float

    """
    #Convert coordinates to radians, if necessary
    if degree:
        m = np.pi / 180 
        lat1 =  m * lat1
        long1 = m * long1
        lat2 = m * lat2
        long2 = m * long2

    x1 = math.sin((lat2-lat1)/2)**2
    x2 = math.cos(lat1) * math.cos(lat2) * math.sin((long2-long1)/2)**2
    d = 2 * r * math.asin(np.sqrt(x1 + x2))
    return d


#------------------------------------------------------------------------------
# Function to perform test
#-=----------------------------------------------------------------------------
def test_difference(df):
    """
    Tests the difference in proportion of votes to candidate 13 between older 
    models and 2020 model (UE2020), using the cluster of municipalities as the 
    experimental units.

    Parameters
    ----------
    df : dataframe
        Dataframe with the results from the 2020 Brazilian Election, as in 
        file `VOTOS_T1E2.xlsx`, but with the added column 'cluster'.

    Returns
    -------
    df_summary : dataframe
        DataFrame with the following columns:
        - `Statistic`: ["Sample size", "Mean", "Standard deviation", "t-statistic", "p-value"]
        - `First_ound` : statistics for the test of favoring to candidate 13
            in older urns x model 2020, in the first round.
        - `Second_ound` : statistics for the test of favoring to candidate 13
            in older urns x model 2020, in the second round.
        
    df_out : dataframe
        Dataframe with the following columns:
        - `REGION` :  region_list
        - `CLUSTER` : cluster_filter_list
        - `VOTE_13_NON_UE2020_T1` : vote_13_urn_other_T1_list
        - `VOTE_22_NON_UE2020_T1` : vote_22_urn_other_T1_list
        - `VOTE_TOTAL_NON_UE2020_T1` : vote_total_urn_other_T1_list
        - `VOTE_13_UE2020_T1` : vote_13_urn_2020_T1_list
        - `VOTE_22_UE2020_T1` : vote_22_urn_2020_T1_list
        - `VOTE_TOTAL_2020_T1` : vote_total_urn_2020_T1_list
        - `PROP_13_NON_UE2020_T1` :  proportion_13_urn_other_T1_list
        - `PROP_13_UE2020_T1` : proportion_13_urn_2020_T1_list
        - `DIF_13_T1`: dif_13_T1_list
        - `PROP_22_NON_UE2020_T1` : proportion_22_urn_other_T1_list
        - `PROP_22_UE2020_T1` : proportion_22_urn_2020_T1_list
        - `DIF_22_T1` : dif_22_T1_list
        - `VOTE_13_NON_UE2020_T2` : vote_13_urn_other_T2_list
        - `VOTE_22_NON_UE2020_T2`: vote_22_urn_other_T2_list
        - `VOTE_TOTAL_NON_UE2020_T2`: vote_total_urn_other_T2_list
        - `VOTE_13_UE2020_T2` : vote_13_urn_2020_T2_list
        - `VOTE_22_UE2020_T2` : vote_22_urn_2020_T2_list
        - `VOTE_TOTAL_2020_T2` : vote_total_urn_2020_T2_list
        - `PROP_13_NON_UE2020_T2` : proportion_13_urn_other_T2_list
        - `PROP_13_UE2020_T2` : proportion_13_urn_2020_T2_list
        - `DIF_13_T2` : dif_13_T2_list
        - `PROP_22_NON_UE2020_T2` : proportion_22_urn_other_T2_list
        - `PROP_22_UE2020_T2` : proportion_22_urn_2020_T2_list
        - `DIF_22_T2` : dif_22_T2_list            

    """
    #Extract urn models and clusters
    cluster_list_aux = df["CLUSTER"].to_list()
    urn_model_list_aux = df["LOG_MODELO"].to_list()
    
    cluster_list = list(set(cluster_list_aux))
    urn_model_list = list(set(urn_model_list_aux))
    
    cluster_list.sort()
    urn_model_list.sort()
    
    #Eliminate non-eletronic urn from the urn list
    urn_model_list_aux = urn_model_list
    urn_model_list = [urn for urn in urn_model_list_aux if urn != "-"]
    
    #Sum votes for 13 on each cluster, distriminating between UE2020 and other urn models
    cluster_filter_list = []
    region_list = []
    vote_13_urn_other_T1_list = []
    vote_13_urn_2020_T1_list = []
    vote_22_urn_other_T1_list = []
    vote_22_urn_2020_T1_list = []
    vote_total_urn_other_T1_list = []
    vote_total_urn_2020_T1_list = []
    proportion_13_urn_other_T1_list = []
    proportion_13_urn_2020_T1_list = []
    proportion_22_urn_other_T1_list = []
    proportion_22_urn_2020_T1_list = []
    dif_13_T1_list = []
    dif_22_T1_list = []
    vote_13_urn_other_T2_list = []
    vote_13_urn_2020_T2_list = []
    vote_22_urn_other_T2_list = []
    vote_22_urn_2020_T2_list = []
    vote_total_urn_other_T2_list = []
    vote_total_urn_2020_T2_list = []
    proportion_13_urn_other_T2_list = []
    proportion_13_urn_2020_T2_list = []
    proportion_22_urn_other_T2_list = []
    proportion_22_urn_2020_T2_list = []
    dif_13_T2_list = []
    dif_22_T2_list = []
    
    for cluster in cluster_list:
        #Identify cluster
        df_cluster = df.loc[df["CLUSTER"] == cluster]
        
        #Identify urns with model not-2020 in cluster
        df_other = df_cluster.loc[(df_cluster["LOG_MODELO"] != 'UE2020') & (df_cluster["LOG_MODELO"] != '-')]
        
        #Identify urns with model 2020 in cluster
        df_2020 = df_cluster.loc[df_cluster["LOG_MODELO"] == 'UE2020']
        
        #Select cluster with both types of urn model
        if len(df_other) > 0 and len(df_2020) > 0:
            region = df_2020.iloc[0]["REGIAO"]
            
            #First round
            sum_13_urn_other_T1 = df_other["T1QT13"].sum()
            sum_13_urn_2020_T1 = df_2020["T1QT13"].sum()
            sum_22_urn_other_T1 = df_other["T1QT22"].sum()
            sum_22_urn_2020_T1 = df_2020["T1QT22"].sum()
            sum_total_urn_other_T1 = df_other["T1QTVAL"].sum()
            sum_total_urn_2020_T1 = df_2020["T1QTVAL"].sum()
            proportion_13_urn_other_T1 = sum_13_urn_other_T1 / sum_total_urn_other_T1
            proportion_13_urn_2020_T1 = sum_13_urn_2020_T1 / sum_total_urn_2020_T1
            proportion_22_urn_other_T1 = sum_22_urn_other_T1 / sum_total_urn_other_T1
            proportion_22_urn_2020_T1 = sum_22_urn_2020_T1 / sum_total_urn_2020_T1
            dif_13_T1 = proportion_13_urn_other_T1 - proportion_13_urn_2020_T1
            dif_22_T1 = proportion_22_urn_other_T1 - proportion_22_urn_2020_T1
            
            #Second round
            sum_13_urn_other_T2 = df_other["T2QT13"].sum()
            sum_13_urn_2020_T2 = df_2020["T2QT13"].sum()
            sum_22_urn_other_T2 = df_other["T2QT22"].sum()
            sum_22_urn_2020_T2 = df_2020["T2QT22"].sum()        
            sum_total_urn_other_T2 = df_other["T2QTVAL"].sum()
            sum_total_urn_2020_T2 = df_2020["T2QTVAL"].sum()
            proportion_13_urn_other_T2 = sum_13_urn_other_T2 / sum_total_urn_other_T2
            proportion_13_urn_2020_T2 = sum_13_urn_2020_T2 / sum_total_urn_2020_T2
            proportion_22_urn_other_T2 = sum_22_urn_other_T2 / sum_total_urn_other_T2
            proportion_22_urn_2020_T2 = sum_22_urn_2020_T2 / sum_total_urn_2020_T2
            dif_13_T2 = proportion_13_urn_other_T2 - proportion_13_urn_2020_T2
            dif_22_T2 = proportion_22_urn_other_T2 - proportion_22_urn_2020_T2
            
            #Append variables
            cluster_filter_list.append(cluster)
            region_list.append(region)
            
            vote_13_urn_other_T1_list.append(sum_13_urn_other_T1)
            vote_22_urn_other_T1_list.append(sum_22_urn_other_T1)
            vote_total_urn_other_T1_list.append(sum_total_urn_other_T1)
            vote_13_urn_2020_T1_list.append(sum_13_urn_2020_T1)
            vote_22_urn_2020_T1_list.append(sum_22_urn_2020_T1)
            vote_total_urn_2020_T1_list.append(sum_total_urn_2020_T1)
            proportion_13_urn_other_T1_list.append(proportion_13_urn_other_T1)
            proportion_13_urn_2020_T1_list.append(proportion_13_urn_2020_T1)
            proportion_22_urn_other_T1_list.append(proportion_22_urn_other_T1)
            proportion_22_urn_2020_T1_list.append(proportion_22_urn_2020_T1)
            dif_13_T1_list.append(dif_13_T1)
            dif_22_T1_list.append(dif_22_T1)
            
            vote_13_urn_other_T2_list.append(sum_13_urn_other_T2)
            vote_22_urn_other_T2_list.append(sum_22_urn_other_T2)
            vote_total_urn_other_T2_list.append(sum_total_urn_other_T2)
            vote_13_urn_2020_T2_list.append(sum_13_urn_2020_T2)
            vote_22_urn_2020_T2_list.append(sum_22_urn_2020_T2)
            vote_total_urn_2020_T2_list.append(sum_total_urn_2020_T2)
            proportion_13_urn_other_T2_list.append(proportion_13_urn_other_T2)
            proportion_13_urn_2020_T2_list.append(proportion_13_urn_2020_T2)
            proportion_22_urn_other_T2_list.append(proportion_22_urn_other_T2)
            proportion_22_urn_2020_T2_list.append(proportion_22_urn_2020_T2)
            dif_13_T2_list.append(dif_13_T2)
            dif_22_T2_list.append(dif_22_T2)
            
    #Prepare output dataframe
    df_out = pd.DataFrame()
    df_out["region"] = region_list
    df_out["cluster"] = cluster_filter_list
    
    df_out["vote_13_non_UE2020_t1"] = vote_13_urn_other_T1_list
    df_out["vote_22_non_UE2020_t1"] = vote_22_urn_other_T1_list
    df_out["vote_total_non_UE2020_t1"] = vote_total_urn_other_T1_list
    df_out["vote_13_UE2020_t1"] = vote_13_urn_2020_T1_list
    df_out["vote_22_UE2020_t1"] = vote_22_urn_2020_T1_list
    df_out["vote_total_2020_t1"] = vote_total_urn_2020_T1_list
    df_out["prop_13_non_UE2020_t1"] = proportion_13_urn_other_T1_list
    df_out["prop_13_UE2020_t1"] = proportion_13_urn_2020_T1_list
    df_out["diff_13_t1"] = dif_13_T1_list
    df_out["prop_22_non_UE2020_t1"] = proportion_22_urn_other_T1_list
    df_out["prop_22_UE2020_t1"] = proportion_22_urn_2020_T1_list
    df_out["diff_22_t1"] = dif_22_T1_list
    
    df_out["vote_13_non_UE2020_t2"] = vote_13_urn_other_T2_list
    df_out["vote_22_non_UE2020_t2"] = vote_22_urn_other_T2_list
    df_out["vote_total_non_UE2020_t2"] = vote_total_urn_other_T2_list
    df_out["vote_13_UE2020_t2"] = vote_13_urn_2020_T2_list
    df_out["vote_22_UE2020_t2"] = vote_22_urn_2020_T2_list
    df_out["vote_total_2020_t2"] = vote_total_urn_2020_T2_list
    df_out["prop_13_non_UE2020_t2"] = proportion_13_urn_other_T2_list
    df_out["prop_13_UE2020_t2"] = proportion_13_urn_2020_T2_list
    df_out["diff_13_t2"] = dif_13_T2_list
    df_out["prop_22_non_UE2020_t2"] = proportion_22_urn_other_T2_list
    df_out["prop_22_UE2020_t2"] = proportion_22_urn_2020_T2_list
    df_out["diff_22_t2"] = dif_22_T2_list
    
    #Sort df_out
    df_out.sort_values(by = ["region", "cluster"], ascending = [True, True], inplace = True)
    
    #Calculate test statistic for first round - candidate 13
    # H0 : d1 = 0
    # H1: d1 > 0
    d_13_T1_H0 = 0.0
    d_13_T1 = df_out["diff_13_t1"].to_list()
    d_13_T1_mean = np.mean(d_13_T1)
    d_13_T1_std = np.std(d_13_T1)
    n = len(df_out)
    d_13_T1_t = np.sqrt(n)*(d_13_T1_mean - d_13_T1_H0)/d_13_T1_std
    d_13_T1_pvalue = 1 - t.cdf(d_13_T1_t, n-1)
    
    #Calculate test statistic for second round - candidate 13
    # H0 : d2 = 0
    # H1: d2 > 0
    d_13_T2_H0 = 0.0
    d_13_T2 = df_out["diff_13_t2"].to_list()
    d_13_T2_mean = np.mean(d_13_T2)
    d_13_T2_std = np.std(d_13_T2)
    n = len(df_out)
    d_13_T2_t = np.sqrt(n)*(d_13_T2_mean - d_13_T2_H0)/d_13_T2_std
    d_13_T2_pvalue = 1 - t.cdf(d_13_T2_t, n-1)

    #Create dataframe with summary statistics
    df_summary = pd.DataFrame()
    stat_list = ["Sample size", "Mean", "Standard deviation", "t-statistic", "p-value"]
    first_round_list = [n, d_13_T1_mean, d_13_T1_std, d_13_T1_t, d_13_T1_pvalue]
    second_round_list = [n, d_13_T2_mean, d_13_T2_std, d_13_T2_t, d_13_T2_pvalue]
    df_summary["statistic"] = stat_list
    df_summary["first_round"] = first_round_list
    df_summary["second_round"] = second_round_list
    
    return df_summary, df_out

#------------------------------------------------------------------------------
#Define oprational folder
folder = "D:\\Dados\\Users\\Euclides\\Diversos\\Brazil_was_stolen"

#Read data file
df = pd.read_excel(os.path.join(folder, "VOTOS_T1E2.xlsx"))
urn_city_name_list = df["NM_MUNICIPIO"].to_list()
urn_city_name_list_lower = [x.lower() for x in urn_city_name_list]

#Read data file with municipalities and respective geographical coordinates
df_cities = pd.read_csv(os.path.join(folder, "municipios.csv"))
city_name_list = df_cities["nome"].to_list()
city_name_list_lower = [x.lower() for x in city_name_list]

#Read data file with states
df_states = pd.read_csv(os.path.join(folder, "estados.csv"))

#Extract coordinates from municipalities
coords = np.matrix(df_cities[["latitude", "longitude"]])

#Define average radius of Earth
r = 6371.088 #km

#Define list of cuttoff distances, in km
cutoff_list = [400.0, 300.0, 200.0, 100.0, 50.0, 20.0, 10.0, 5.0, 0.1]

#Fill state name and region in df_cities
codigo_uf_list = df_states["codigo_uf"].to_list()
uf_in_df_cities_list = []
region_symbol_in_df_cities_list = []
for i, row in df_cities.iterrows():
    uf_ind = codigo_uf_list.index(row["codigo_uf"])
    uf = df_states.iloc[uf_ind]["uf"]
    uf_in_df_cities_list.append(uf)
df_cities["uf"] = uf_in_df_cities_list

#Find urn of the city in the city_name_list
city_index_list = []
missing = -1
for city in urn_city_name_list_lower:
    if city in city_name_list_lower:
        city_index_list.append(city_name_list_lower.index(city))
    else:
        city_index_list.append(missing)

#Add city index in dataframe
df["city_index"] = city_index_list

#Filter out records of unlocated cities
df = df.loc[df["city_index"] != missing]

#Perform analysis
columns = ["cutoff", "n_t1", "diff_mean_t1", "diff_std_t2", "diff_t_t1", "diff_pvalue_t1",\
           "n_t2", "diff_mean_t2", "diff_std_t2", "diff_t_t2", "diff_pvalue_t2"]
df_cutoff = pd.DataFrame(columns = columns)
for cutoff in cutoff_list:
    #Get clusters based on the haversine distance
    epsilon = cutoff / r
    db = DBSCAN(eps=epsilon, min_samples=1, algorithm='ball_tree', metric='haversine').fit(np.degrees(coords))
    cluster_label_list = db.labels_
    clusters = list(set(cluster_label_list))
    
    #Append cluster in the city dataframe
    df_cities_aux = df_cities.copy()
    df_cities_aux["cluster"] = cluster_label_list
    
    #Assign cluster to the urn dataset
    city_index_list = df["city_index"].to_list()
    urn_city_cluster_list = [cluster_label_list[x] for x in city_index_list]
    df["CLUSTER"] = urn_city_cluster_list

    #Separate cluster in NE from clusters in other regions
    df_NE = df.loc[df["REGIAO"] == "NE"]

    #Perform analysis
    df_summary_cluster_NE, df_out_cluster_NE = test_difference(df_NE)
    
    #For each cluster, find the components (cities)
    comp_list = []    
    for i, row_i in df_out_cluster_NE.iterrows():
        df_components = df_cities_aux.loc[df_cities_aux["cluster"] == row_i["cluster"]]
        comp_i_str = ""
        for j, row_j in df_components.iterrows():
            comp_i_str += f'{row_j["nome"]} / {row_j["uf"]};'
        comp_list.append(comp_i_str[0:-1])
    df_out_cluster_NE["cluster_components"] = comp_list
    
    #Append data in output dataframe
    line = [cutoff] + df_summary_cluster_NE["first_round"].to_list() + df_summary_cluster_NE["second_round"].to_list()
    df_cutoff.loc[len(df_cutoff)] = line
    
#Save last dataframes (at city level)
df_summary_cluster_NE.to_csv(os.path.join(folder, "df_summary_cluster_NE.csv"), sep=";", decimal=",", encoding = 'latin1')
df_out_cluster_NE.to_csv(os.path.join(folder, "df_out_cluster_NE.csv"), sep=";", decimal=",", encoding = 'latin1')
df_cutoff.to_csv(os.path.join(folder, "df_cutoff.csv"), sep=";", decimal=",")
    
#Plot histogram d_13_T1 discrimating NE x other regions - candidate 13
d_13_t1 = df_out_cluster_NE["diff_13_t1"]
plt.hist(d_13_t1)
plt.title("Hist. diff. urn non-2020 x urn 2020, for candidate 13, 1st round, NE cities")
plt.savefig(os.path.join(folder, "histogram_cluster_NE_13_T1.jpg"))
plt.show()
plt.close()

d_13_t2 = df_out_cluster_NE["diff_13_t2"]
plt.hist(d_13_t2)
plt.title("Hist. diff. urn non-2020 x urn 2020, for candidate 13, 2nd round, NE cities")
plt.savefig(os.path.join(folder, "histogram_cluster_NE_13_T2.jpg"))
plt.show()
plt.close()

#Plot p-valyues x cutoff
plt.plot(cutoff_list, df_cutoff["diff_pvalue_t1"].to_list(), cutoff_list, df_cutoff["diff_pvalue_t2"].to_list())
plt.legend(["1st round", "2nd round"])
plt.title("p-values from test H0: dmean_13 = 0 x H1: dmean_13 > 0")
plt.xlabel("cutoff (km)")
plt.ylabel("p-value")
plt.savefig(os.path.join(folder, "pvalue_cluster_NE.jpg"))
plt.show()
plt.close()


    