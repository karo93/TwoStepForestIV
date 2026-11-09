# Data Source

The analysis builds on the pseudo-dataset provided by Keele et al., based on the SPOTlight study, a prospective cohort of 13,011 adult patients from 48 UK hospitals (2010–2011).

Since the original patient-level data cannot be publicly released, Keele et al. constructed a pseudo-dataset by resampling the original observations to approximate their distribution.

The complete pseudo-dataset is available from the original data source provided by Keele et al. [add link].

The R package ivmodel provides related data and functionality for instrumental variable analysis.

# Data Preprocessing

Patients older than 93 years are excluded, resulting in a final sample of 12,903 patients.

The analysis focuses on:

* Outcome: 28-day mortality (dead28).
* Treatment: ICU admission (icu_bed).
* Instrument: ICU bed availability (open_bin), coded as 1 if fewer than four ICU beds were available.
* Covariates: Patient demographics, comorbidities, clinical severity scores (ICNARC, NEWS, SOFA), and indicators of observed and recommended levels of care.

The hospital identifier (site) is excluded from the main analysis but considered in an additional robustness check.
