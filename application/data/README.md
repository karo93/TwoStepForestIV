## Data Source

The empirical application builds on the ICU pseudo-data provided by [Keele et al. (2020)](https://academic.oup.com/jrsssa/article/183/4/1501/7056438), based on the SPOTlight study, a prospective cohort of 13,011 adult patients from 48 UK hospitals (2010–2011).

Since the original patient-level data cannot be publicly released, Keele et al. constructed a pseudo-dataset by resampling observations from the original study to approximate their distribution.The pseudo-data used in this application were obtained from the replication materials accompanying Keele et al. (2020). A related version of the ICU pseudo-data is available through the [R package ivmodel](https://cran.r-project.org/web/packages/ivmodel/index.html). Researchers should consult the original publication and its associated replication materials to obtain the required data.

## Data Preprocessing

Patients older than 93 years are excluded, resulting in a final sample of 12,903 observations.

The analysis focuses on:

* Outcome: 28-day mortality (dead28).
* Treatment: ICU admission (icu_bed).
* Instrument: ICU bed availability (open_bin), coded as 1 if fewer than four ICU beds were available.
* Covariates: Patient demographics, comorbidities, clinical severity scores (ICNARC, NEWS, SOFA), and indicators of observed and recommended levels of care.

The hospital identifier (site) is excluded from the main analysis but included in an additional robustness check.
