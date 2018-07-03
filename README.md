
---
title: A comparison of models for polygon data, point data or both for malaria mapping
---



1. Use Madagascar + Indonesia as case studies.
2. do cross validation on both pr and polygon. no good way to combine so just have to report both in every case
3. can't see good way to compare stacked generalisation with linear polygon so have to do all linear.
4. if this paper falls out easily can expand to ml and polygon as covariates. or project for ra?
5. probably want random cv and spatial cv?
6. might need more countries but do first draft like this.
7. spatial only. but need to choose year/year range.
8. compare Pearson and spearman. maybe Mae as well.
9. how to handle spatial priors? try 3 values and pick best CV in each case.
10. gives country X model X CV X metric X prevalence vs incidence.


The script `run_analysis.R` is the main script that runs the full analysis.
The other scripts contain function definitions used by `run_analysis.R`.
Currently I am reading data directly from Z, mostly from the GBD2017 folder.

