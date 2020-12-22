# Math 499 Projects

## PhysiCell: extension and simulations 

This project involves the [PhysiCell simulator](https://nanohub.org/tools/pc4covid19) [1] udated by the [SARS-Cov-2 Tissue Simulation Coalition](http://physicell.org/covid19/). The goal of this open-source software is to simulate the dynamics of SARS-Cov-2 within host using the expertise of researcher from different fields including virology, inmunology, cell biology, mathematical biology and cloud computing [2]. The simulator currently has implemented chemotaxis, viral infection, cell death, inflamatory response and immune response. Dr. Hillen, Dr. Newby and Dr. Carlos are working on testing the possibility that SARS-Cov-2 affects differents organs and tissues of the body via the inmmune system. The goal of this project is to 1) implement a module to incorporate our models into the PhysiCell simulator and 2) use the simulator to test our assumptions. 
  
**Required skills and knowledge**: programming (C, Java, python), mathematical modelling (enzyme kinetic, cell biology).

[1] Ghaffarizadeh, A., Heiland, R., Friedman, S. H., Mumenthaler, S. M., & Macklin, P. (2018). PhysiCell: An open source physics-based cell simulator for 3-D multicellular systems. PLOS Computational Biology, 14(2), e1005991. https://doi.org/10.1371/journal.pcbi.1005991

[2] Getz, M., Wang, Y., An, G., Becker, A., Cockrell, C., & Collier, N. (2020). Rapid community-driven development of a SARS- CoV-2 tissue simulator. 1–55. https://doi.org/10.1101/2020.04.02.019075.

---

## Classification and ranking COVID-19 patients using machine learning

The goal of this project is to use machine learning to classify patients at risk. This is done by assigning a ranking or risk measument to patients based on thier age, sex, medical history, co-morbidities and other clinical features. The main idea is to reproduce the ideas and results in these two articles: [1] [2], possibly for the case of Canada and using different clasification methods. A preliminary task is obtain and mine data, posibly from [Kaggle](https://www.kaggle.com/allen-institute-for-ai/CORD-19-research-challenge), [Statistics Canada](https://www150.statcan.gc.ca/n1/en/type/data?MM=1), or [GenOMICC](https://genomicc.org/).

**Required skills and knowledge**: programming (python, julia, R), statistics, machine learning.

[1] Knight, S. R., Ho, A., Pius, R., Buchan, I., Carson, G., Drake, T. M., Dunning, J., Fairfield, C. J., Gamble, C., Green, C. A., Gupta, R., Halpin, S., Hardwick, H. E., Holden, K. A., Horby, P. W., Jackson, C., McLean, K. A., Merson, L., Nguyen-Van-Tam, J. S., … Harrison, E. M. (2020). Risk stratification of patients admitted to hospital with covid-19 using the ISARIC WHO Clinical Characterisation Protocol: Development and validation of the 4C Mortality Score. The BMJ, 370(September), 1–13. https://doi.org/10.1136/bmj.m3339

[2] Sudre, C., Lee, K., Ni Lochlainn, M., Varsavsky, T., Murray, B., Graham, M., Menni, C., Modat, M., Bowyer, R., Nguyen, L., Drew, D. A., Joshi, A., Ma, W., Guo, C. G., Lo, C. H., Ganesh, S., Buwe, A., Capdevila Pujol, J., Lavigne du Cadet, J., … Ourselin, S. (2020). Symptom clusters in Covid19: A potential clinical prediction tool from the COVID Symptom study app. 1–12. https://doi.org/10.1101/2020.06.12.20129056