# Covid 19 Physiology

In collaboration with: [Thomas Hillen](http://www.math.ualberta.ca/~thillen/) and [Jay Newby](https://newby-jay.github.io/)

---

## Physiology of Covid-19 research

This research focuses on within-host dynamics of the SARS-CoV-2. Why are people with underlying health conditions particularly at risk? Why some recovered patients experience long-term effects? How does the get to different parts of the body, and what it does to different organs? How does the within host dynamics change from variant to variant? This is a compiling and exploratory research of within host dynamics of SARS-CoV-2.

- Risk of COVID-19 using machine learning
- 

[List of acronyms and abbreviations](./references/acronyms.md)

[List of papers](./references/covid-19_papers.md)

[COVID-19 news](./references/news.md)

### Notebooks

- [Model for the covid-19 in body organs](./julia/08%20SARS-Cov2%20Physiology%20model.ipynb) [(html preview)](https://htmlpreview.github.io/?https://github.com/ccontrer/20-08Covid19Physiology/blob/master/julia/08%20SARS-Cov2%20Physiology%20model.html)

--- 

## Virus Load Function research

A virus load curve is a time series of the amount of virus (measured in infections particles per mL) within a host.  The virus infects cells and replicatess as the inmune system fights the infection. The dynamics between the virus, healthy cells, infected cells, and the immune system can be described by a system of differential equations ([review of viral kinetic models](https://link.springer.com/article/10.1007%2Fs10928-014-9363-3)). Many COVID-19 researchers focus on the within-host dynamics of the virus and require a viral load input for their models. Here we propose a personalizable and practical model for the SARS-CoV-2 viral load that can be used as input for models requiring a viral load.

- Hillen, T., Contreras, C., & Newby, J. M. (2021). Personalized Virus Load Curves of SARS-CoV-2 Infection. [Preprint](https://www.medrxiv.org/content/10.1101/2021.01.21.21250268v1).

[List of acronyms and abbreviations](./references/acronyms.md)

[List of references](./references/virus_load.md)

### Notebooks

- [Virus load function: data fitting](./julia/02%20Fitting%20the%20virus%20load%20function.ipynb)
- [Viral kinetic ODE model: numerical solution](./julia/01%20Solving%20the%20virus%20target%20model.ipynb)
- [Viral kinetic ODE model: data fitting](./julia/03%20Fitting%20the%20virus%20target%20model.ipynb)
- [Comparing the viral load function with the ODE model](./julia/04%20Comparing%20the%20virus%20load%20function%20with%20the%20virus-target%20model.ipynb)
- Fitting models to data:
  - [Influenza data](./julia/05%20Fitting%20models%20to%20influenza%20data.ipynb)
  - [Macaque monkey SARS-Cov2 data](./julia/06%20Fitting%20models%20to%20macaque%20monkey%20SARS-Cov2%20data.ipynb)
  - [Human SARS-Cov2 data](./julia/07%20Fitting%20models%20to%20human%20SARS-Cov2%20data.ipynb)

---

### Resources:

- [Covid-19 Physiology Group](https://sites.google.com/ualberta.ca/cov-pg/home) 
- [Treatment guidelines](https://www.idsociety.org/practice-guideline/covid-19-guideline-treatment-and-management/)
- [Data sets](references/datasets.md)
- [WHO Quick links](https://www.who.int/emergencies/diseases/novel-coronavirus-2019)
- [PhysiCell simulator](http://physicell.org/)

### Visualization tools

- [John Hopkins University dashboard](https://gisanddata.maps.arcgis.com/apps/opsdashboard/index.html#/bda7594740fd40299423467b48e9ecf6)
- [WHO dashboard](https://covid19.who.int/table)
- [nCoV2019.live (@AviSchiffmann)](https://ncov2019.live/)
- [Tracker Canada](https://covid19tracker.ca/)
- [Covid Trends (aatsihb)](https://aatishb.com/covidtrends/)
- [Next strain](https://nextstrain.org/groups/neherlab/ncov/S.N501?c=gt-S_501,69&p=grid&r=country)