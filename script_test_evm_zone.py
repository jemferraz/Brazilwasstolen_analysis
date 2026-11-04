# -*- coding: utf-8 -*-
"""
Created on Sat Nov  5 12:05:08 2022

"""
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy.stats import t

#Define parameters
OUTPUT_FILE_ROOT = "summary_zone"
ALPHA = 0.05 #Confidence level for the intervals

#Read data file
df_all = pd.read_excel("data/VOTOS_T1E2.xlsx")

#Choose a filtering option:
# 0. All regions
# 1. N + NE
# 2. All regions but N + NE
region_option = 1
if region_option == 0:
    df = df_all
    suffix = ""
elif region_option == 1: 
    df = df_all.loc[(df_all["REGIAO"] == "N")|(df_all["REGIAO"] == "NE")]    
    suffix = "_NNE"
elif region_option == 2:
    df = df_all.loc[(df_all["REGIAO"] != "N")&(df_all["REGIAO"] != "NE")]
    suffix = "_other"
else:
    region_option = 0
    df = df_all
    suffix = ""    
    
OUTPUT_FILE = f'{OUTPUT_FILE_ROOT}{suffix}.txt'  

    
#Create unique identifier of the electoral zone by composing id_municipio + zone
id_city_list = df["ID_MUNICIPIO"].to_list()
nr_zone_list = df["NR_ZONA"].to_list()
n = len(df)
zone_list_aux = [f'{id_city_list[i]}_{nr_zone_list[i]}' for i in range(n)]
df["ID_MUNICIPIO_NR_ZONE"] = zone_list_aux

zone_list = list(set(zone_list_aux))
zone_list.sort()

#Extract evm models and ditricts
evm_model_list_aux = df["LOG_MODELO"].to_list()
evm_model_list = list(set(evm_model_list_aux))
evm_model_list.sort()

#Eliminate non-eletronic evm from the evm list
evm_model_list_aux = evm_model_list
evm_model_list = [evm for evm in evm_model_list_aux if evm != "-"]

#Sum votes for 13 and 22 on each zone in the NE, distriminating between UE2020 and other evm models
zone_filter_list = []
region_list = []
uf_list = []
city_list = []
vote_13_evm_other_T1_list = []
vote_13_evm_2020_T1_list = []
vote_total_evm_other_T1_list = []
vote_total_evm_2020_T1_list = []
proportion_13_evm_other_T1_list = []
proportion_13_evm_2020_T1_list = []
dif_13_T1_list = []
vote_13_evm_other_T2_list = []
vote_13_evm_2020_T2_list = []
vote_total_evm_other_T2_list = []
vote_total_evm_2020_T2_list = []
proportion_13_evm_other_T2_list = []
proportion_13_evm_2020_T2_list = []
dif_13_T2_list = []

for zone in zone_list:
    #Identify district
    df_zone = df.loc[(df["ID_MUNICIPIO_NR_ZONE"] == zone) & (df["REGIAO"] == 'NE')]
    
    #Identify evms with model not-2020 in district
    df_other = df_zone.loc[(df_zone["LOG_MODELO"] != 'UE2020') & (df_zone["LOG_MODELO"] != '-')]
    
    #Identify evms with model 2020 in district
    df_2020 = df_zone.loc[df_zone["LOG_MODELO"] == 'UE2020']
    
    #If district has both types of evm...
    if len(df_other) > 0 and len(df_2020) > 0:
        region = df_2020.iloc[0]["REGIAO"]
        uf = df_2020.iloc[0]["UF"]
        city = df_2020.iloc[0]["NM_MUNICIPIO"]
        zone_only = df_2020.iloc[0]["NR_ZONA"]
        
        #First round
        sum_13_evm_other_T1 = df_other["T1QT13"].sum()
        sum_13_evm_2020_T1 = df_2020["T1QT13"].sum()
        sum_total_evm_other_T1 = df_other["T1QTVAL"].sum()
        sum_total_evm_2020_T1 = df_2020["T1QTVAL"].sum()
        proportion_13_evm_other_T1 = sum_13_evm_other_T1 / sum_total_evm_other_T1
        proportion_13_evm_2020_T1 = sum_13_evm_2020_T1 / sum_total_evm_2020_T1
        dif_13_T1 = proportion_13_evm_other_T1 - proportion_13_evm_2020_T1
        
        #Second round
        sum_13_evm_other_T2 = df_other["T2QT13"].sum()
        sum_13_evm_2020_T2 = df_2020["T2QT13"].sum()    
        sum_total_evm_other_T2 = df_other["T2QTVAL"].sum()
        sum_total_evm_2020_T2 = df_2020["T2QTVAL"].sum()
        proportion_13_evm_other_T2 = sum_13_evm_other_T2 / sum_total_evm_other_T2
        proportion_13_evm_2020_T2 = sum_13_evm_2020_T2 / sum_total_evm_2020_T2
        dif_13_T2 = proportion_13_evm_other_T2 - proportion_13_evm_2020_T2
        
        #Append variables
        zone_filter_list.append(zone_only)
        city_list.append(city)
        uf_list.append(uf)
        region_list.append(region)
        
        vote_13_evm_other_T1_list.append(sum_13_evm_other_T1)
        vote_total_evm_other_T1_list.append(sum_total_evm_other_T1)
        vote_13_evm_2020_T1_list.append(sum_13_evm_2020_T1)
        vote_total_evm_2020_T1_list.append(sum_total_evm_2020_T1)
        proportion_13_evm_other_T1_list.append(proportion_13_evm_other_T1)
        proportion_13_evm_2020_T1_list.append(proportion_13_evm_2020_T1)
        dif_13_T1_list.append(dif_13_T1)
                
        vote_13_evm_other_T2_list.append(sum_13_evm_other_T2)
        vote_total_evm_other_T2_list.append(sum_total_evm_other_T2)
        vote_13_evm_2020_T2_list.append(sum_13_evm_2020_T2)
        vote_total_evm_2020_T2_list.append(sum_total_evm_2020_T2)
        proportion_13_evm_other_T2_list.append(proportion_13_evm_other_T2)
        proportion_13_evm_2020_T2_list.append(proportion_13_evm_2020_T2)
        dif_13_T2_list.append(dif_13_T2)
        
#Prepare output dataframe
df_out = pd.DataFrame()
df_out["REGIAO"] = region_list
df_out["UF"] = uf_list
df_out["NM_MUNICIPIO"] = city_list
df_out["NR_ZONA"] = zone_filter_list

df_out["VOTE_13_NON_UE2020_T1"] = vote_13_evm_other_T1_list
df_out["VOTE_TOTAL_NON_UE2020_T1"] = vote_total_evm_other_T1_list
df_out["VOTE_13_UE2020_T1"] = vote_13_evm_2020_T1_list
df_out["VOTE_TOTAL_2020_T1"] = vote_total_evm_2020_T1_list
df_out["PROP_13_NON_UE2020_T1"] = proportion_13_evm_other_T1_list
df_out["PROP_13_UE2020_T1"] = proportion_13_evm_2020_T1_list
df_out["DIF_13_T1"] = dif_13_T1_list

df_out["VOTE_13_NON_UE2020_T2"] = vote_13_evm_other_T2_list
df_out["VOTE_TOTAL_NON_UE2020_T2"] = vote_total_evm_other_T2_list
df_out["VOTE_13_UE2020_T2"] = vote_13_evm_2020_T2_list
df_out["VOTE_TOTAL_2020_T2"] = vote_total_evm_2020_T2_list
df_out["PROP_13_NON_UE2020_T2"] = proportion_13_evm_other_T2_list
df_out["PROP_13_UE2020_T2"] = proportion_13_evm_2020_T2_list
df_out["DIF_13_T2"] = dif_13_T2_list


#Sort df_out
df_out.sort_values(by = ["REGIAO", "UF", "NM_MUNICIPIO", "NR_ZONA"], ascending = [True, True, True, True], inplace = True)
df_out.to_csv(f"df_out_zone{suffix}.csv", sep=";", decimal=",", encoding='latin1')

#Calculate test statistic for first round
# H0 : d1 = 0
# H1: d1 > 0
d1_H0 = 0.0
d1 = df_out["DIF_13_T1"].to_list()
d1_mean = np.mean(d1)
d1_std = np.std(d1)
n = len(df_out)
t1 = np.sqrt(n)*(d1_mean - d1_H0)/d1_std
pvalue1 = 1 - t.cdf(t1, n-1)
#Interval for the estimates of the differences 
t_ALPHA = t.ppf(1 - ALPHA/2, n-1)
d1_lower = d1_mean - t_ALPHA*d1_std / np.sqrt(n)
d1_upper = d1_mean + t_ALPHA*d1_std / np.sqrt(n)
 

#Plot histogram d1
plt.hist(d1)
plt.title("Histogram difference EVM older than UE2020 x EVM UE2020,\n for candidate 13, 1st round, by UF")
plt.savefig(f"img/histogram_Zone_T1{suffix}.jpg")
plt.show()
plt.close()


#Calculate test statistic for second round
# H0 : d2 = 0
# H1: d2 > 0
d2_H0 = 0.0
d2 = df_out["DIF_13_T2"].to_list()
d2_mean = np.mean(d2)
d2_std = np.std(d2)
n = len(df_out)
t2 = np.sqrt(n)*(d2_mean - d2_H0)/d2_std
pvalue2 = 1 - t.cdf(t2, n-1)
#Interval for the estimates of the differences
d2_lower = d2_mean - t_ALPHA*d2_std / np.sqrt(n)
d2_upper = d2_mean + t_ALPHA*d2_std / np.sqrt(n)


#Plot histogram d2
plt.hist(d2)
plt.title("Histogram difference EVM older than UE2020 x EVM UE2020,\n for candidate 13, 2nd round, by UF")
plt.savefig(f"img/histogram_Zone_T2{suffix}.jpg")
plt.show()
plt.close()


#Create report
with open(OUTPUT_FILE, "w") as file:
    file.write( "============================================================================================\n")
    file.write( "                    Summary report for Electoral Zone analysis                              \n")
    file.write( "============================================================================================\n\n")
    file.write( "--------------------------------------------------------------------------------------------\n")
    file.write( "Differences between proportion of votes for candidate 13 in EVMs older than UE2020 x UE2020\n")
    file.write( "--------------------------------------------------------------------------------------------\n")
    file.write( "1.Estimates for the first round:\n")
    file.write(f" - Mean difference:                      {(d1_mean*100):>12,.1f}% \n")
    file.write(f" - Standard deviation difference:        {(d1_std*100):>12,.1f}% \n")
    file.write(f" - Lower bound for {(1-ALPHA)*100}% confidence level:{(d1_lower*100):>11,.1f}% \n")
    file.write(f" - Upper bound for {(1-ALPHA)*100}% confidence level:{(d1_upper*100):>11,.1f}% \n")
    file.write(f" - t_statistic:                          {(t1):>12,.2f} \n")
    file.write(f" - p-value:                              {(pvalue1):>12,.2e} \n")
    file.write(f" - Degrees of freedom:                   {(n-1):>12,.0f} \n")
    file.write( "2.Estimates for the second round:\n")
    file.write(f" - Mean difference:                      {(d2_mean*100):>12,.1f}% \n")
    file.write(f" - Standard deviation difference:        {(d2_std*100):>12,.1f}% \n")
    file.write(f" - Lower bound for {(1-ALPHA)*100}% confidence level:{(d2_lower*100):>11,.1f}% \n")
    file.write(f" - Upper bound for {(1-ALPHA)*100}% confidence level:{(d2_upper*100):>11,.1f}% \n")
    file.write(f" - t_statistic:                          {(t2):>12,.2f} \n")
    file.write(f" - p-value:                              {(pvalue2):>12,.2e} \n")
    file.write(f" - Degrees of freedom:                   {(n-1):>12,.0f} \n\n")