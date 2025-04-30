# Dataset Description

**WP_50k_MS_20240807.csv**
* This file contains the Woods Point aftershock sequence (WPAS), which is the observed seismicity following the mainshock (21 September 2021) until 7 August 2024. The data cover a 50‑km radius from the mainshock's epicentre (SRC location).

**WP_50k_from2000_beforeMS.csv**
* This file comprises the pre‑mainshock background seismicity recorded from 2000 until the mainshock, within a 50‑km radius of the mainshock's epicentre (SRC location).

**other_SCRs_study.csv**
* This file is derived from Table 3 in Ebel (2009) (see also Appendix 4). Only events recorded as Mw or Ml are included since our study focuses exclusively on these two magnitude types. 

**dataset1_Mw.csv, dataset2_Mw.csv, dataset3_Mw.csv, dataset4_Mw.csv, dataset5_Mw.csv, dataset6_Mw.csv**
* These six datasets represent the WPAS over different spatio-temporal domains (see Appendix 3). The maximum spatial domains are found in dataset3_Mw.csv and dataset6_Mw.csv, which cover areas within a 20‑km radius of the mainshock's epicentre. These files contain converted parameters generated using the code from the Parameter_preparation module. The columns, from left to right, are: latitude, longitude, depth, Ml, Mw, Mo, RA, and RA_scaling.

**dataset6_August.csv**
* This file represents the WPAS within a 20‑km radius of the mainshock epicentre, covering events from the mainshock until 7 August 2024. dataset6_Mw.csv covers the same spatial domain as dataset6_August.csv, but its temporal window spans from the mainshock until 30 April 2024.

**WP_50k_AS.txt**
* This catalogue contains WPAS data (local magnitude) recorded from the occurrence of the mainshock until August 2024 within a 50‑km radius of the mainshock epicentre. Columns (from left to right) are: occurrence time (year/month/day), occurrence time (hour/minute/second), latitude, longitude, depth (km), and magnitude.

**WP_50k_AS_Mw.txt**
* This file has the same spatial and temporal domain as WP_50k_AS.txt, but the magnitudes are in moment magnitude (Mw) format.

# Environment Setup

**Matlab**

*Mc value estimation
  1. MAXC
  2. MAXC + Boostrapping
  3. MBS
  4. MBS + Boostrapping

*Gutenberg–Richter (G-R) Law
  1. Parameters Estimation (a and b values)
  2. Parameters Comparison with other SCRs (Ebel, 2009)
  3. Epistemic Uncertainty in the selection of Background Seismicity Models

*Seismicity Spatial Distribution
  1. Pre-mainshock Seismicity
  2. Post-mainshock Seismicity (WPAS)

*Seismic Moment Monthly Released Rate
  1. Pre-mainshock Seismicity
  2. Post-mainshock Seismicity (WPAS)

*Fault Planar Plane Fitting
  1. Fault Plane Orientation
  2. Rupture Area Identification
  3. Seismic Momenet Released from WPAS relate to the Fault Plane

**Python**
*Modified Omori's Law
  1. Utsu_Omori's Law Parameters (c, k and p values) Estimation
  2. Reasenberg and Jones Omori's Law Parameter (a value) Estimation

*Anticipated Apparent Duration
