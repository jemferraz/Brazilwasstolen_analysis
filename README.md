# Analysis of the `BrazilWasStolen` main claim

It was published by an Argentine individual in the website [BraziWasStolen.com](https://brazilwasstolen.com/en) a 76 page document with an analysis of the 2022 Brazilian Election for Presidency. The main claim in this document is that there is an anomalous difference in the proportion of votes for candidate 13 in the electronic ballot models manufactured prior to 2020 (old models: UE2015, UE2013, etc.) than on those manufactured in the 2020 model (UE2020).

The objective of this work is to check this claim by means of a formal statistical test. In the next section we present the data that was used. Next, we explain the statistical procedure, and the corresponding results for different levels of aggregation of the electronic ballot machines. In the last section we present our conclusions.


## Data

The data was downloaded from [BraziWasStolen.com](https://brazilwasstolen.com/en). However, the site is currently inaccessible from Brazil, due to interference of the Electoral Supreme Court (TSE). In order to make the analysis possible to other researchers, we provided the Excel spreadsheet `VOTOS_T1E2.xlsx` in the folder `data`.

Some people were worried that this particular dataset might have not come from the [official TSE website](https://dadosabertos.tse.jus.br/). However a simple check on a sample of the number of votes for the candidates, at any aggregation level (state, municipality, electoral zone) indicates that the values in the spreadsheet are consistent with the ones presented at the official website. 

Nonetheless, we tried to download the original data from the [official TSE website](https://dadosabertos.tse.jus.br/), only to find that the data had been tampered: the data files inside the zip files have dates from 2022-11-07 22:00, one week after the election, which, by itself, is highly suspicious.


## Methods

The assertion that the proportion of votes to candidate 13 in the old models is higher than the proportion of votes to candidate 13 in the 2020 model can be tested by means of the so called paired comparison design (see, for example, Montgomery(1983)). The basic idea is to compare the treatments (proportions of votes in the old and new models) in the same experimental unit, which here represents an aggregation of the electronic ballot machines that contains both types of machines, old and 2020 models. 

Let's call *d* the difference between those proportions for a particular experimental unit. Under the null hypothesis of no difference in proportions, the mean *d* should be zero, against the alternative hypothesis that the mean *d* should be positive. That is, 

$$ H_{0}: \mu_{d} = 0 $$

$$ H_{1}: \mu_{d} > 0 $$


The test statistics for this hypothesis is

$$ t = \sqrt{n} \frac{\overline{d}}{s_{d}} $$,

which has a standard t-student distribution with *n*-1 degrees of freedom under the null. 


## Results

We performed our analysis for three aggregation levels:
- Units of the federation (UF): the states plus the federal district, where the capital of the country, Brasilia, is located.
- Clusters of municipalities (cities), based on their distance from each other.
- Electoral zones.

### UFs as experimental units

For this case we developed the Python script `script_test_urns_uf.py`.

Table 1 and figure 1 below illustrate the results for the first and second rounds of the election, with a total of 27 UFs. Since the histograms seem to indicate a bi- or tri-modal distribution, we also plotted the observed differences separating the UFs of the North-eastern region from the other UFs in the country.

![table_1](/img/table_UF.jpg)
*Table 1 – Proportions for candidate 13 in the old and new models for all UFs, and respective differences*

![figure_1](/img/histograms_UF.jpg)
*Figure 1 – Histograms for the differences in the first round (left) and in the second round (right), for all UFs (top) and discriminating UFs in the NE from other UFs (bottom)*

In the first round, the sample average difference was 7.17%, with a sample standard deviation of 8.01%, resulting in a t-statistic of 4.65, and corresponding p-value of 4.21e-5.  In the second round, the sample average difference was 6.27%, with a sample standard deviation of 8.46%, resulting in a t-statistic of 3.85, and corresponding p-value of 3.44e-4. Therefore, we reject H0, for both the first and the second rounds, at the usual significance levels. Or in other words, we found evidence beyond a reasonable doubt, that the proportion of votes for candidate 13 in the old models was higher than that in the 2020 model, as stated by the Argentine website "BrazilWasStolen".

Since the total number of votes in the old models in our sample was 64 660 376 and 64 832 379, for the first and second rounds, respectively, we estimate that the corresponding average difference in votes to be 4 637 092 and 4 065 286 in favour of candidate 13 in the old models.


### Clusters of cities in the North-eastern region as experimental units

After the analysis for the UFs in the previous section, some people raised the hypothesis that the observed difference might be due to a non-random allocation of the old machine models to places where candidate 13 would naturally have had more votes. That is possible, though also strange: why would the TSE not allocate the machine models randomly in the UFs?

In effect, this hypothesis corresponds to a supposed "heterogeneity" in the previous experimental units (UFs). Thus, in order to obtain more "homogeneous" groups we could aggregate the machine models based on their haversine distance from each other by means of Cluster Analysis, and sample those clusters with both old and 2020 models. The haversine distance is the distance along a maximum circle of a spherical surface between two points, given their geographical coordinates (latitude and longitude). The cluster are formed by specifying a maximum (cut-off) distance between members of the cluster.

By varying the cut-off distance, from very large (say, 400 km) to very small (say, 01. km), we should observe every more homogeneous groups, if the machine models had been non randomly distributed, weakening the test, and causing the failure of the rejection the original H0. However, if, on the contrary, we still reject the original H0 after decreasing the cut-off distance, we would fail to infer any heterogeneity in the allocation of the machine models.

Script `script_test_urns_cluster_NE.py` implements the analysis described above for clusters formed from machine models in municipalities from the North-eastern region of the country, since that region was the source of the difference observed in the analysis for the UFs. The results obtained by varying the cut-off distance from 400 km to 0.1 km are presented in the table 2 and figure 2 below. In table 3 we present the results for the most disaggregated level (at 0.1 km) which corresponds to only one member (municipality) per cluster.

![table_2](/img/table_cluster.jpg)

*Table 2 – Summary of the p-values for the test, as a function of the cut-off distance (km)*

![figure_2](/img/fig_pvalue_cluster_NE.jpg)

*Figure 2 – p-values for the test, as a function of the cut-off distance (km)*

![table_3](/img/table_cities_NE.jpg)
*Table 3 – Proportions for candidate 13 in the old and new models for the municipalities in the NE region with both types of models, and respective differences. These clusters are generated at the most disaggregated level (cut-off = 0.1 km).*



We can observe from the chart that, for all cut-off levels, the p-values of the H0 test remains well below the usual significance level (5%), implying that the observed difference in favour to candidate 13 in the old machine models is maintained at all aggregation levels of the machines, in the North-eastern region of the country. Therefore, we failed to find any indication of a non-random allocation of the machine models.


### Electoral zones in the North-eastern region as experimental units

Eventually, the last level we could aggregate machine models is by electoral zone. Script `script_test_urns_zone_NE.py` performs the analysis for the electoral zones in the North-eastern region. The results are presented in table 4 below.

![table_4](/img/table_zone_NE.jpg)
*Table 4 – Proportions for candidate 13 in the old and new models for the electoral zones in the NE region with both types of models, and respective differences.*

In the first round, the sample average difference was -4.79%, with a sample standard deviation of 6.81%, resulting in a t-statistic of -2.11, and corresponding p-value of 0.966.  In the second round, the sample average difference was -5.11%%, with a sample standard deviation of 6.34%, resulting in a t-statistic of -2.42, and corresponding p-value of 0.979.

Therefore, we failed to reject H0 at the usual significance levels, but that is to be expected due to the very weak power of the test, since the sample is small (only 8 experimental units), biased (all electoral zones are from municipalities in the state of Maranhão, MA), and irrelevant (the total number of votes in the old models for the sampled experimental units is only 8181 and 8158 votes, for the first and second rounds, respectively).

Therefore, there is no compelling evidence that the observed difference in favour of candidate 13 in the old electronic ballot models is due to a non-random allocation of the machine models to particular electoral zones.



## Conclusions

We confirm the assertion of "BrazilWasStolen" that the proportion of votes for candidate 13 in the old electronic ballot models is significantly higher than those in the 2020 model, beyond a reasonable doubt.

We found no compelling evidence that this was the result of a non-random allocation of the old models to places where candidate 13 would have more votes. On the contrary, we found compelling evidence that the observed effect is maintained at every aggregation level of the municipalities.

Therefore, the actual reason for the observed significant difference must be something else.

## References

Montgomery, Douglas C., **Design and Analysis of Experiments**, Second Edition, John Wiley & Sons, 1983.
"# Brazilwasstolen_analysis" 
