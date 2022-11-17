# -*- coding: utf-8 -*-
"""
Created on Sat Nov  5 12:05:08 2022

@author: José Euclides de Melo Ferraz
"""
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy.stats import t

#Define oprational folder
folder = "D:\\Dados\\Users\\Euclides\\Diversos\\Brazil_was_stolen"

#Read data file
df = pd.read_excel(os.path.join(folder, "VOTOS_T1E2.xlsx"))

#Extract urn models and ditricts
uf_list_aux = df["UF"].to_list()
urn_model_list_aux = df["LOG_MODELO"].to_list()

uf_list = list(set(uf_list_aux))
urn_model_list = list(set(urn_model_list_aux))

uf_list.sort()
urn_model_list.sort()

#Eliminate non-eletronic urn from the urn list
urn_model_list_aux = urn_model_list
urn_model_list = [urn for urn in urn_model_list_aux if urn != "-"]

#Count votes for 13 on each district, distriminating between UE2020 and other urn models
uf_filter_list = []
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

for uf in uf_list:
    #Identify district
    df_uf = df.loc[df["UF"] == uf]
    
    #Identify urns with model not-2020 in district
    df_other = df_uf.loc[(df_uf["LOG_MODELO"] != 'UE2020') & (df_uf["LOG_MODELO"] != '-')]
    
    #Identify urns with model 2020 in district
    df_2020 = df_uf.loc[df_uf["LOG_MODELO"] == 'UE2020']
    
    #Select uf with both types of urn
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
        uf_filter_list.append(uf)
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
df_out["REGIAO"] = region_list
df_out["UF"] = uf_filter_list

df_out["VOTE_13_NON_UE2020_T1"] = vote_13_urn_other_T1_list
df_out["VOTE_22_NON_UE2020_T1"] = vote_22_urn_other_T1_list
df_out["VOTE_TOTAL_NON_UE2020_T1"] = vote_total_urn_other_T1_list
df_out["VOTE_13_UE2020_T1"] = vote_13_urn_2020_T1_list
df_out["VOTE_22_UE2020_T1"] = vote_22_urn_2020_T1_list
df_out["VOTE_TOTAL_2020_T1"] = vote_total_urn_2020_T1_list
df_out["PROP_13_NON_UE2020_T1"] = proportion_13_urn_other_T1_list
df_out["PROP_13_UE2020_T1"] = proportion_13_urn_2020_T1_list
df_out["DIF_13_T1"] = dif_13_T1_list
df_out["PROP_22_NON_UE2020_T1"] = proportion_22_urn_other_T1_list
df_out["PROP_22_UE2020_T1"] = proportion_22_urn_2020_T1_list
df_out["DIF_22_T1"] = dif_22_T1_list

df_out["VOTE_13_NON_UE2020_T2"] = vote_13_urn_other_T2_list
df_out["VOTE_22_NON_UE2020_T2"] = vote_22_urn_other_T2_list
df_out["VOTE_TOTAL_NON_UE2020_T2"] = vote_total_urn_other_T2_list
df_out["VOTE_13_UE2020_T2"] = vote_13_urn_2020_T2_list
df_out["VOTE_22_UE2020_T2"] = vote_22_urn_2020_T2_list
df_out["VOTE_TOTAL_2020_T2"] = vote_total_urn_2020_T2_list
df_out["PROP_13_NON_UE2020_T2"] = proportion_13_urn_other_T2_list
df_out["PROP_13_UE2020_T2"] = proportion_13_urn_2020_T2_list
df_out["DIF_13_T2"] = dif_13_T2_list
df_out["PROP_22_NON_UE2020_T2"] = proportion_22_urn_other_T2_list
df_out["PROP_22_UE2020_T2"] = proportion_22_urn_2020_T2_list
df_out["DIF_22_T2"] = dif_22_T2_list


#Sort df_out
df_out.sort_values(by = ["REGIAO", "UF"], ascending = [True, True], inplace = True)
df_out.to_csv(os.path.join(folder, "df_out_uf_13_22.csv"), sep=";", decimal=",")

#Define probability of type I error
n = len(df_out)
alpha = 0.05
t_inv = t.isf(alpha/2, n-1)

#Calculate test statistic for first round - candidate 13
# H0 : d1 = 0
# H1: d1 > 0
d_13_T1_H0 = 0.0
d_13_T1 = df_out["DIF_13_T1"].to_list()
d_13_T1_mean = np.mean(d_13_T1)
d_13_T1_std = np.std(d_13_T1)
n = len(df_out)
t_13_T1 = np.sqrt(n)*(d_13_T1_mean - d_13_T1_H0)/d_13_T1_std
pvalue_13_T1 = 1 - t.cdf(t_13_T1, n-1)
d_13_T1_mean_lb = d_13_T1_mean - t_inv * d_13_T1_std / np.sqrt(n)
d_13_T1_mean_ub = d_13_T1_mean + t_inv * d_13_T1_std / np.sqrt(n)

#Plot histogram d_13_T1
plt.hist(d_13_T1)
plt.title("Hist. diff. urn non-2020 x urn 2020, for candidate 13, 1st round, by UF")
plt.savefig(os.path.join(folder, "histogram_UF_13_T1.jpg"))
plt.show()
plt.close()

#Calculate test statistic for second round - candidate 13
# H0 : d2 = 0
# H1: d2 > 0
d_13_T2_H0 = 0.0
d_13_T2 = df_out["DIF_13_T2"].to_list()
d_13_T2_mean = np.mean(d_13_T2)
d_13_T2_std = np.std(d_13_T2)
n = len(df_out)
t_13_T2 = np.sqrt(n)*(d_13_T2_mean - d_13_T2_H0)/d_13_T2_std
pvalue_13_T2 = 1 - t.cdf(t_13_T2, n-1)
pvalue_13_T2 = 1 - t.cdf(t_13_T2, n-1)
d_13_T2_mean_lb = d_13_T2_mean - t_inv * d_13_T2_std / np.sqrt(n)
d_13_T2_mean_ub = d_13_T2_mean + t_inv * d_13_T2_std / np.sqrt(n)

#Plot histogram d_13_T2
plt.hist(d_13_T2)
plt.title("Hist. diff. urn non-2020 x urn 2020, for candidate 13, 2nd round, by UF")
plt.savefig(os.path.join(folder, "histogram_UF_13_T2.jpg"))
plt.show()
plt.close()


#Calculate number of votes for differences in proportion - candidate 13
votes_total_urn_other_T1 = df_out["VOTE_TOTAL_NON_UE2020_T1"].sum()
votes_total_urn_other_T2 = df_out["VOTE_TOTAL_NON_UE2020_T2"].sum()
votes_d_13_T1_mean = d_13_T1_mean * votes_total_urn_other_T1
votes_d_13_T1_mean_lb = d_13_T1_mean_lb * votes_total_urn_other_T1
votes_d_13_T1_mean_ub = d_13_T1_mean_ub * votes_total_urn_other_T1
votes_d_13_T2_mean = d_13_T2_mean * votes_total_urn_other_T2
votes_d_13_T2_mean_lb = d_13_T2_mean_lb * votes_total_urn_other_T2
votes_d_13_T2_mean_ub = d_13_T2_mean_ub * votes_total_urn_other_T2


#Calculate test statistic for first round - candidate 22
# H0 : d1 = 0
# H1: d1 > 0
d_22_T1_H0 = 0.0
d_22_T1 = df_out["DIF_22_T1"].to_list()
d_22_T1_mean = np.mean(d_22_T1)
d_22_T1_std = np.std(d_22_T1)
t_22_T1 = np.sqrt(n)*(d_22_T1_mean - d_22_T1_H0)/d_22_T1_std
pvalue_22_T1 = 1 - t.cdf(t_22_T1, n-1)
d_22_T1_mean_lb = d_22_T1_mean - t_inv * d_22_T1_std / np.sqrt(n)
d_22_T1_mean_ub = d_22_T1_mean + t_inv * d_22_T1_std / np.sqrt(n)

#Plot histogram d_13_T1
plt.hist(d_22_T1)
plt.title("Hist. diff. urn non-2020 x urn 2020, for candidate 22, 1st round, by UF")
plt.savefig(os.path.join(folder, "histogram_UF_22_T1.jpg"))
plt.show()
plt.close()

#Calculate test statistic for second round - candidate 22
# H0 : d2 = 0
# H1: d2 > 0
d_22_T2_H0 = 0.0
d_22_T2 = df_out["DIF_22_T2"].to_list()
d_22_T2_mean = np.mean(d_22_T2)
d_22_T2_std = np.std(d_22_T2)
t_22_T2 = np.sqrt(n)*(d_22_T2_mean - d_22_T2_H0)/d_22_T2_std
pvalue_22_T2 = 1 - t.cdf(t_22_T2, n-1)
d_22_T2_mean_lb = d_22_T2_mean - t_inv * d_22_T2_std / np.sqrt(n)
d_22_T2_mean_ub = d_22_T2_mean + t_inv * d_22_T2_std / np.sqrt(n)

#Plot histogram d_13_T2
plt.hist(d_22_T2)
plt.title("Hist. diff. urn non-2020 x urn 2020, for candidate 22, 2nd round, by UF")
plt.savefig(os.path.join(folder, "histogram_UF_22_T2.jpg"))
plt.show()
plt.close()


#Calculate number of votes for differences in proportion - candidate 22
votes_total_urn_other_T1 = df_out["VOTE_TOTAL_NON_UE2020_T1"].sum()
votes_total_urn_other_T2 = df_out["VOTE_TOTAL_NON_UE2020_T2"].sum()
votes_d_22_T1_mean = d_22_T1_mean * votes_total_urn_other_T1
votes_d_22_T1_mean_lb = d_22_T1_mean_lb * votes_total_urn_other_T1
votes_d_22_T1_mean_ub = d_22_T1_mean_ub * votes_total_urn_other_T1
votes_d_22_T2_mean = d_22_T2_mean * votes_total_urn_other_T2
votes_d_22_T2_mean_lb = d_22_T2_mean_lb * votes_total_urn_other_T2
votes_d_22_T2_mean_ub = d_22_T2_mean_ub * votes_total_urn_other_T2


#------------------------------------------------------------------------------
#Segregate data between NE and other regions
d_13_T1_ne = []
d_13_T1_other = []
d_13_T2_ne = []
d_13_T2_other = []
d_22_T1_ne = []
d_22_T1_other = []
d_22_T2_ne = []
d_22_T2_other = []

for i, row in df_out.iterrows():
    if row["REGIAO"] == "NE":
        d_13_T1_ne.append(row["DIF_13_T1"])
        d_13_T2_ne.append(row["DIF_13_T2"])
        d_22_T1_ne.append(row["DIF_22_T1"])
        d_22_T2_ne.append(row["DIF_22_T2"])
    else:
        d_13_T1_other.append(row["DIF_13_T1"])
        d_13_T2_other.append(row["DIF_13_T2"])
        d_22_T1_other.append(row["DIF_22_T1"])
        d_22_T2_other.append(row["DIF_22_T2"])
        
        
#Plot histogram d_13_T1 discrimating NE x other regions - candidate 13
plt.hist([d_13_T1_other, d_13_T1_ne])
plt.legend(["Other regions", "NE"])
plt.title("Hist. diff. urn non-2020 x urn 2020, for candidate 13, 1st round, by UF")
plt.savefig(os.path.join(folder, "histogram_UF_13_T1_region.jpg"))
plt.show()
plt.close()
        
#Plot histogram d_13_T2 discrimating NE x other regions - candidate 13
plt.hist([d_13_T2_other, d_13_T2_ne])
plt.legend(["Other regions", "NE"])
plt.title("Hist. diff. urn non-2020 x urn 2020, for candidate 13, 2nd round, by UF")
plt.savefig(os.path.join(folder, "histogram_UF_13_T2_region.jpg"))
plt.show()
plt.close()


#Plot histogram d_22_T1 discrimating NE x other regions - candidate 22
plt.hist([d_22_T1_other, d_22_T1_ne])
plt.legend(["Other regions", "NE"])
plt.title("Hist. diff. urn non-2020 x urn 2020, for candidate 22, 1st round, by UF")
plt.savefig(os.path.join(folder, "histogram_UF_22_T1_region.jpg"))
plt.show()
plt.close()
        
#Plot histogram d_22_T2 discrimating NE x other regions - candidate 22
plt.hist([d_22_T2_other, d_22_T2_ne])
plt.legend(["Other regions", "NE"])
plt.title("Hist. diff. urn non-2020 x urn 2020, for candidate 22, 2nd round, by UF")
plt.savefig(os.path.join(folder, "histogram_UF_22_T2_region.jpg"))
plt.show()
plt.close()
