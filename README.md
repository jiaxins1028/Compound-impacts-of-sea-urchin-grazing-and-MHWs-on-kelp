# Compound-impacts-of-sea-urchin-grazing-and-MHWs-on-kelp
Data and Code for the research of the Compound impacts of sea urchin grazing and marine heatwaves (MHWs) on kelp collapse along east coast of Australia.
Due to their size, some raw datasets must be accessed directly from the providers. However, all datasets are freely available, and the instructions below explain exactly how to download them and where to place the files.  MHW calculations are performed using xmhw: https://github.com/coecms/xmhw.

| Raw Data | Source |
|-----:|---------------|
| Long-term Reef Monitoring Program kelp cover (*Ecklonia radiata*) and urchin density (*Centrostephanus rodgersii*) | Australia’s National Reef Monitoring Network. Raw data and instructions can be accessed from the Australian Temperate Reef Collaboration (ATRC) and Reef Life Survey (RLS). [https://catalogue-imos.aodn.org.au/geonetwork/srv/eng/catalog.search#/metadata/ec424e4f-0f55-41a5-a3f2-726bc4541947] |
| Sea surface temperature (SST) | High Resolution NOAA Optimum Interpolation 1/4 Degree Daily SST (OISST) Analysis, Version 2.1 (Huang et al. 2021; Reynolds et al. 2007) [https://www.ncei.noaa.gov/metadata/geoportal/rest/metadata/item/gov.noaa.ncdc:C00844/html] |
| Global climate simulation  outputs of SST | The Coupled Model Intercomparison Project, phase 6 (CMIP6; O'Neill et al. 2016)  available at [https://doi.org/10.25914/6009627c7af03] |


| Data File | Description |
|-----:|---------------|
| ercr_atrcrls_id.csv | Dataset of *Ecklonia* and *Centrostephanus* at survey_id level after quality control. The identical survey_id represents a survey conducted on the same transect.|


| Code File | Description |
|-----:|---------------|
| binomial_glm_projection.ipynb | Binomial GLM projection of kelp collapse probablity under different MHW intensity and urchin increases and CMIP6 SSP585 scenario (Davis et al. 2023) |
| boxplots_timeseries.ipynb | Time period series of *Ecklonia* cover and *Centrostephanus* density, and MHW maximum intensities. |
| kum_glm.pdf | Backward stepwise model selection of Binomial GLM regression of kelp collapse probablity with urchin density and MHW metrics. |
| mhw_metrics_each_survey.ipynb | Compiling *Ecklonia* cover and *Centrostephanus* density, and MHW metrics annually |
| nb_glm_kelp_urchin_mhwintensity.R | Relationship bwteern *Ecklonia* cover and *Centrostephanus* density, under different MHW maximum intensities levels. |
| urchin_distribution_collapse_propotion.ipynb | Latitudinal distribution of *Centrostephanus* density and derived kelp resilience propotion. |

# References

Hobday, A. J., Alexander, L. V., Perkins, S. E., Smale, D. A., Straub, S. C., Oliver, E. C. J., Benthuysen, J. A., Burrows, M. T., Donat, M. G., Feng, M., Holbrook, N. J., Moore, P. J., Scannell, H. A., Sen Gupta, A., & Wernberg, T. (2016). A hierarchical approach to defining marine heatwaves. Progress in Oceanography, 141, 227-238.

Huang, B., Liu, C., Banzon, V., Freeman, E., Graham, G., Hankins, B., Smith, T., & Zhang, H.-M. (2021). Improvements of the Daily Optimum Interpolation Sea Surface Temperature (DOISST) Version 2.1. Journal of Climate, 34(8), 2923-2939. 

Institute for Marine and Antarctic Studies (IMAS); Parks Victoria; Department of Primary Industries (DPI), New South Wales Government; Parks and Wildlife Tasmania; Department for Environment and Water (DEWNR), South Australia; Department of Biodiversity, Conservation and Attractions (DBCA), Western Australia; Integrated Marine Observing System (IMOS), 2024, IMOS - National Reef Monitoring Network Sub-Facility – Benthic cover data (in situ surveys), database provided, 10/10/2024.

O'Neill, B.C., Tebaldi, C., van Vuuren, D.P., Eyring, V., Friedlingstein, P., Hurtt, G. et al. (2016). The Scenario Model Intercomparison Project (ScenarioMIP) for CMIP6. Geosci. Model Dev., 9, 3461-3482.

Reynolds, R.W., Smith, T.M., Liu, C., Chelton, D.B., Casey, K.S. & Schlax, M.G. (2007). Daily High-Resolution-Blended Analyses for Sea Surface Temperature. Journal of Climate, 20, 5473-5496.
