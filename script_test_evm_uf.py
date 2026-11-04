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
OUTPUT_FILE_ROOT = "summary_UF"
ALPHA = 0.05 #Confidence level for the intervals

#Read data file
df_all = pd.read_excel("data/VOTOS_T1E2.xlsx")

#Choose a filtering option:
# 0. All regions
# 1. N + NE
# 2. All regions but N + NE
region_option = 0
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

#Extract evm models and ditricts
uf_list_aux = df["UF"].to_list()
evm_model_list_aux = df["LOG_MODELO"].to_list()

uf_list = list(set(uf_list_aux))
evm_model_list = list(set(evm_model_list_aux))

uf_list.sort()
evm_model_list.sort()

#Eliminate non-eletronic evm from the evm list
evm_model_list_aux = evm_model_list
evm_model_list = [evm for evm in evm_model_list_aux if evm != "-"]

#Count votes for 13 on each district, distriminating between UE2020 and other evm models
uf_filter_list = []
region_list = []
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

for uf in uf_list:
    #Identify district
    df_uf = df.loc[df["UF"] == uf]
    
    #Identify evms with model not-2020 in district
    df_other = df_uf.loc[(df_uf["LOG_MODELO"] != 'UE2020') & (df_uf["LOG_MODELO"] != '-')]
    
    #Identify evms with model 2020 in district
    df_2020 = df_uf.loc[df_uf["LOG_MODELO"] == 'UE2020']
    
    #Select uf with both types of evm
    if len(df_other) > 0 and len(df_2020) > 0:
        region = df_2020.iloc[0]["REGIAO"]
        
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
        sum_22_evm_2020_T2 = df_2020["T2QT22"].sum()        
        sum_total_evm_other_T2 = df_other["T2QTVAL"].sum()
        sum_total_evm_2020_T2 = df_2020["T2QTVAL"].sum()
        proportion_13_evm_other_T2 = sum_13_evm_other_T2 / sum_total_evm_other_T2
        proportion_13_evm_2020_T2 = sum_13_evm_2020_T2 / sum_total_evm_2020_T2
        dif_13_T2 = proportion_13_evm_other_T2 - proportion_13_evm_2020_T2
        
        #Append variables
        uf_filter_list.append(uf)
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
df_out["UF"] = uf_filter_list

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
df_out.sort_values(by = ["REGIAO", "UF"], ascending = [True, True], inplace = True)
df_out.to_csv("df_out_uf.csv", sep=";", decimal=",")

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
plt.savefig(f"img/histogram_UF_T1{suffix}.jpg")
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
plt.savefig(f"img/histogram_UF_T2{suffix}.jpg")
plt.show()
plt.close()

if region_option == 0:
    #Seggregate data between N, NE x other regions
    d1_nne = []
    d1_other = []
    d2_nne = []
    d2_other = []
    for i, row in df_out.iterrows():
        if row["REGIAO"] in ["N", "NE"]:
            d1_nne.append(row["DIF_13_T1"])
            d2_nne.append(row["DIF_13_T2"])         
        else:
            d1_other.append(row["DIF_13_T1"])
            d2_other.append(row["DIF_13_T2"])
            
    #Plot histogram d1 discrimating NE x other regions
    plt.hist([d1_other, d1_nne])
    plt.legend(["Other regions", "N, NE"])
    plt.title("Histogram difference EVM older than UE2020 x EVM UE2020,\n for candidate 13, 1st round, by UF")
    plt.savefig("img/histogram_UF_T1_region.jpg")
    plt.show()
    plt.close()
            
    #Plot histogram d2 discrimating NE x other regions
    plt.hist([d2_other, d2_nne])
    plt.legend(["Other regions", "N, NE"])
    plt.title("Histogram difference EVM older than UE2020 x EVM UE2020,\n for candidate 13, 2nd round, by UF")
    plt.savefig("img/histogram_UF_T2_region.jpg")
    plt.show()
    plt.close()

#Estimate amount of votes transferred from the other candidates to 13
df_older_2020 = df.loc[(df["LOG_MODELO"] < 'UE2020') & (df["LOG_MODELO"] != '-')]
votes_T1 = df["T1QTVAL"].sum()
votes_T2 = df["T2QTVAL"].sum()
votes_older_UE2020_T1 = df_older_2020["T1QTVAL"].sum() 
votes_older_UE2020_T2 = df_older_2020["T2QTVAL"].sum()
#Estimates of transferred votes T1
transf_votes_T1_mean = votes_older_UE2020_T1 * d1_mean
transf_votes_T1_lower = votes_older_UE2020_T1 * d1_lower
transf_votes_T1_upper = votes_older_UE2020_T1 * d1_upper
#Estimates of transferred votes T2
transf_votes_T2_mean = votes_older_UE2020_T2 * d2_mean
transf_votes_T2_lower = votes_older_UE2020_T2 * d2_lower
transf_votes_T2_upper = votes_older_UE2020_T2 * d2_upper

#Estimate of actual votes for candidates
votes_13_T1 = df["T1QT13"].sum()
votes_22_T1 = df["T1QT22"].sum()
votes_13_T2 = df["T2QT13"].sum()
votes_22_T2 = df["T2QT22"].sum()

actual_votes_13_T1_mean = votes_13_T1 - transf_votes_T1_mean
actual_votes_13_T1_lower = votes_13_T1 - transf_votes_T1_upper
actual_votes_13_T1_upper = votes_13_T1 - transf_votes_T1_lower
actual_votes_22_T1_mean = votes_22_T1 + transf_votes_T1_mean
actual_votes_22_T1_lower = votes_22_T1 + transf_votes_T1_lower
actual_votes_22_T1_upper = votes_22_T1 + transf_votes_T1_upper

actual_votes_13_T2_mean = votes_13_T2 - transf_votes_T2_mean
actual_votes_13_T2_lower = votes_13_T2 - transf_votes_T2_upper
actual_votes_13_T2_upper = votes_13_T2 - transf_votes_T2_lower
actual_votes_22_T2_mean = votes_22_T2 + transf_votes_T2_mean
actual_votes_22_T2_lower = votes_22_T2 + transf_votes_T2_lower
actual_votes_22_T2_upper = votes_22_T2 + transf_votes_T2_upper

#Estimate (gross, due to Jensen's inequality) of the proportion of vote change to 13
p_change_T1_mean = transf_votes_T1_mean / actual_votes_22_T1_mean
p_change_T1_lower = transf_votes_T1_lower / actual_votes_22_T1_upper
p_change_T1_upper = transf_votes_T1_upper / actual_votes_22_T1_lower
p_change_T2_mean = transf_votes_T2_mean / actual_votes_22_T2_mean
p_change_T2_lower = transf_votes_T2_lower / actual_votes_22_T2_upper
p_change_T2_upper = transf_votes_T2_upper / actual_votes_22_T2_lower


#Create report
with open(OUTPUT_FILE, "w") as file:
    file.write( "============================================================================================\n")
    file.write( "                          Summary report for UF analysis                                    \n")
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
    file.write( "--------------------------------------------------------------------------------------------\n")
    file.write( "Transferred votes from the other candidates to 13 \n")
    file.write( "--------------------------------------------------------------------------------------------\n")
    file.write( "1.Estimates for first round:\n")
    file.write(f" - Total valid votes:                     {votes_T1:>12,.0f} {(votes_T1/votes_T1*100):>10,.1f}%\n")
    file.write(f" - Total votes in EVMs older than UE2020: {votes_older_UE2020_T1:>12,.0f} {(votes_older_UE2020_T1/votes_T1*100):>10,.1f}%\n")
    file.write(f" - Mean transferred votes:                {transf_votes_T1_mean:>12,.0f} {(transf_votes_T1_mean/votes_T1*100):>10,.1f}%\n")
    file.write(f" - Lower bound transferred votes:         {transf_votes_T1_lower:>12,.0f} {(transf_votes_T1_lower/votes_T1*100):>10,.1f}%\n")
    file.write(f" - Upper bound transferred votes:         {transf_votes_T1_upper:>12,.0f} {(transf_votes_T1_upper/votes_T1*100):>10,.1f}%\n\n")
    file.write( "2.Estimates for second round:\n")
    file.write(f" - Total valid votes:                     {votes_T2:>12,.0f} {(votes_T2/votes_T2*100):>10,.1f}%\n")
    file.write(f" - Total votes in EVMs older than UE2020: {votes_older_UE2020_T2:>12,.0f} {(votes_older_UE2020_T2/votes_T2*100):>10,.1f}%\n")
    file.write(f" - Mean transferred votes:                {transf_votes_T2_mean:>12,.0f} {(transf_votes_T2_mean/votes_T2*100):>10,.1f}%\n")
    file.write(f" - Lower bound transferred votes:         {transf_votes_T2_lower:>12,.0f} {(transf_votes_T2_lower/votes_T2*100):>10,.1f}%\n")
    file.write(f" - Upper bound transferred votes:         {transf_votes_T2_upper:>12,.0f} {(transf_votes_T2_upper/votes_T2*100):>10,.1f}%\n\n")
    file.write( "--------------------------------------------------------------------------------------------\n")
    file.write( "Estimates of actual votes for 13 \n")
    file.write( "--------------------------------------------------------------------------------------------\n")
    file.write( "1.First round:\n")
    file.write(f" - Total votes in the dataset             {votes_13_T1:>12,.0f} {(votes_13_T1/votes_T1*100):>10,.1f}%\n")
    file.write(f" - Mean:                                  {actual_votes_13_T1_mean:>12,.0f} {(actual_votes_13_T1_mean/votes_T1*100):>10,.1f}%\n")
    file.write(f" - Lower:                                 {actual_votes_13_T1_lower:>12,.0f} {(actual_votes_13_T1_lower/votes_T1*100):>10,.1f}%\n")
    file.write(f" - Upper bound:                           {actual_votes_13_T1_upper:>12,.0f} {(actual_votes_13_T1_upper/votes_T1*100):>10,.1f}%\n\n")
    file.write( "2.Second round:\n")
    file.write(f" - Total votes in the dataset             {votes_13_T2:>12,.0f} {(votes_13_T2/votes_T2*100):>10,.1f}%\n")
    file.write(f" - Mean:                                  {actual_votes_13_T2_mean:>12,.0f} {(actual_votes_13_T2_mean/votes_T2*100):>10,.1f}%\n")
    file.write(f" - Lower:                                 {actual_votes_13_T2_lower:>12,.0f} {(actual_votes_13_T2_lower/votes_T2*100):>10,.1f}%\n")
    file.write(f" - Upper bound:                           {actual_votes_13_T2_upper:>12,.0f} {(actual_votes_13_T2_upper/votes_T2*100):>10,.1f}%\n\n")
    file.write( "--------------------------------------------------------------------------------------------\n")
    file.write( "Estimates of actual votes for 22 \n")
    file.write( "--------------------------------------------------------------------------------------------\n")
    file.write( "1.First round:\n")
    file.write(f" - Total votes in the dataset             {votes_22_T1:>12,.0f} {(votes_22_T1/votes_T1*100):>10,.1f}%\n")
    file.write(f" - Mean:                                  {actual_votes_22_T1_mean:>12,.0f} {(actual_votes_22_T1_mean/votes_T1*100):>10,.1f}%\n")
    file.write(f" - Lower:                                 {actual_votes_22_T1_lower:>12,.0f} {(actual_votes_22_T1_lower/votes_T1*100):>10,.1f}%\n")
    file.write(f" - Upper bound:                           {actual_votes_22_T1_upper:>12,.0f} {(actual_votes_22_T1_upper/votes_T1*100):>10,.1f}%\n\n")
    file.write( "2.Second round:\n")
    file.write(f" - Total votes in the dataset             {votes_22_T2:>12,.0f} {(votes_22_T2/votes_T2*100):>10,.1f}%\n")
    file.write(f" - Mean:                                  {actual_votes_22_T2_mean:>12,.0f} {(actual_votes_22_T2_mean/votes_T2*100):>10,.1f}%\n")
    file.write(f" - Lower:                                 {actual_votes_22_T2_lower:>12,.0f} {(actual_votes_22_T2_lower/votes_T2*100):>10,.1f}%\n")
    file.write(f" - Upper bound:                           {actual_votes_22_T2_upper:>12,.0f} {(actual_votes_22_T2_upper/votes_T2*100):>10,.1f}%\n\n")
    file.write( "--------------------------------------------------------------------------------------------\n")
    file.write( "Estimate of proportion of vote change from 22 to 13 (gross, due to Jensen's inequality) \n")
    file.write( "--------------------------------------------------------------------------------------------\n")
    file.write( "1.First round:\n")
    file.write(f" - Mean:                                  {(p_change_T1_mean*100):>12,.1f}%\n")
    file.write(f" - Lower:                                 {(p_change_T1_lower*100):>12,.1f}%\n")
    file.write(f" - Upper bound:                           {(p_change_T1_upper*100):>12,.1f}%\n\n")
    file.write( "2.Second round:\n")
    file.write(f" - Mean:                                  {(p_change_T2_mean*100):>12,.1f}%\n")
    file.write(f" - Lower:                                 {(p_change_T2_lower*100):>12,.1f}%\n")
    file.write(f" - Upper bound:                           {(p_change_T2_upper*100):>12,.1f}%\n\n")
    
    