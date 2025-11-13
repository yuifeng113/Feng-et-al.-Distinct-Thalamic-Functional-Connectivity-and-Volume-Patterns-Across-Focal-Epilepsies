* Encoding: UTF-8.

*Overall patients vs. controls.
GLM Ant_ipsi Ant_contra Lat_ipsi Lat_contra Med_ipsi Med_contra Post_ipsi Post_contra BY Controls
  /WSFACTOR=nuclei 8 Polynomial 
  /METHOD=SSTYPE(3)
  /PLOT=PROFILE(Controls*nuclei) TYPE=LINE ERRORBAR=SE(2) MEANREFERENCE=NO 
    YAXIS=AUTO
  /EMMEANS=TABLES(Controls*nuclei) COMPARE(Controls) ADJ(BONFERRONI)
  /EMMEANS=TABLES(Controls*nuclei) COMPARE(nuclei) ADJ(BONFERRONI)
  /PRINT=DESCRIPTIVE ETASQ 
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=nuclei
  /DESIGN=Controls.

*TLE vs. FLE vs. PQE.
USE ALL.
COMPUTE filter_$=(Controls = 0).
VARIABLE LABELS filter_$ 'Controls = 0 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.
GLM Ant_ipsi Ant_contra Lat_ipsi Lat_contra Med_ipsi Med_contra Post_ipsi Post_contra BY Focal
  /WSFACTOR=nuclei 4 Polynomial side 2 Polynomial 
  /METHOD=SSTYPE(3)
  /PLOT=PROFILE( Focal*side  Focal*side*nuclei  nuclei*Focal Focal*nuclei) TYPE=LINE ERRORBAR=SE(2) MEANREFERENCE=NO 
    YAXIS=AUTO
   /EMMEANS=TABLES( Focal) COMPARE ADJ(BONFERRONI)
   /EMMEANS=TABLES( Focal*side) COMPARE(Focal) ADJ(BONFERRONI)
  /EMMEANS=TABLES( Focal*side) COMPARE(side) ADJ(BONFERRONI)
  /EMMEANS=TABLES( Focal*nuclei) COMPARE(Focal) ADJ(BONFERRONI)
  /EMMEANS=TABLES( Focal*nuclei*side) COMPARE(Focal) ADJ(BONFERRONI)
  /EMMEANS=TABLES( Focal*nuclei*side) COMPARE(nuclei) ADJ(BONFERRONI)
  /EMMEANS=TABLES( Focal*nuclei*side) COMPARE(side) ADJ(BONFERRONI)
  /EMMEANS=TABLES(nuclei) COMPARE ADJ(BONFERRONI)
  /EMMEANS=TABLES( nuclei*side) COMPARE(side) ADJ(BONFERRONI)
  /EMMEANS=TABLES( nuclei*side) COMPARE(nuclei) ADJ(BONFERRONI)
  /PRINT=DESCRIPTIVE ETASQ 
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=nuclei side nuclei*side
  /DESIGN= Focal.

*One-sample t-tests were used to compare z-scores of each epilepsy group (overall focal epilepsy, TLE, FLE, and PQE) to zero (the control group mean) 
    for each thalamic nuclei group (anterior, lateral, medial, posterior) on both the ipsilateral and contralateral sides. 
  T-TEST
  /TESTVAL=0
  /MISSING=ANALYSIS
  /VARIABLES=Ant_ipsi Lat_ipsi Med_ipsi Post_ipsi Ant_contra 
    Lat_contra Med_contra Post_contra
  /ES DISPLAY(TRUE)
  /CRITERIA=CI(.95).

* Ipsi vs contra. Paired sample t tests.
T-TEST PAIRS=Ant_ipsi Lat_ipsi Med_ipsi Post_ipsi WITH Ant_contra Lat_contra Med_contra Post_contra 
    (PAIRED)
  /ES DISPLAY(TRUE) STANDARDIZER(SD)
  /CRITERIA=CI(.9500)
  /MISSING=ANALYSIS.

*Congenital vs. acquired etiologies.
USE ALL.
COMPUTE filter_$=(congenital_acquired_aetiology = 'A' |     congenital_acquired_aetiology='C').
VARIABLE LABELS filter_$ "congenital_acquired_aetiology = 'A' |     "+
    "congenital_acquired_aetiology='C' (FILTER)".
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.
GLM Ant_ipsi Ant_contra Lat_ipsi Lat_contra Med_ipsi Med_contra Post_ipsi Post_contra BY congenital_acquired_aetiology
  /WSFACTOR=nuclei 4 Polynomial side 2 Polynomial 
  /METHOD=SSTYPE(3)
  /PLOT=PROFILE(congenital_acquired_aetiology*side congenital_acquired_aetiology*side*nuclei nuclei*congenital_acquired_aetiology congenital_acquired_aetiology*nuclei) TYPE=LINE ERRORBAR=SE(2) MEANREFERENCE=NO 
    YAXIS=AUTO
  /EMMEANS=TABLES(congenital_acquired_aetiology*side) COMPARE(congenital_acquired_aetiology) ADJ(BONFERRONI)
  /EMMEANS=TABLES(congenital_acquired_aetiology*side) COMPARE(side) ADJ(BONFERRONI)
  /EMMEANS=TABLES(congenital_acquired_aetiology*nuclei) COMPARE(congenital_acquired_aetiology) ADJ(BONFERRONI)
  /EMMEANS=TABLES(congenital_acquired_aetiology*nuclei*side) COMPARE(congenital_acquired_aetiology) ADJ(BONFERRONI)
  /EMMEANS=TABLES(congenital_acquired_aetiology*nuclei*side) COMPARE(nuclei) ADJ(BONFERRONI)
  /EMMEANS=TABLES(congenital_acquired_aetiology*nuclei*side) COMPARE(side) ADJ(BONFERRONI)
  /PRINT=DESCRIPTIVE ETASQ 
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=nuclei side nuclei*side
  /DESIGN=congenital_acquired_aetiology.

*FBTCS.
USE ALL.
COMPUTE filter_$=(types_focal_to_bilateral_tonic_clonic = 'yes' | 
    types_focal_to_bilateral_tonic_clonic='no').
VARIABLE LABELS filter_$ "types_focal_to_bilateral_tonic_clonic = 'yes' | "+
    "types_focal_to_bilateral_tonic_clonic='no' (FILTER)".
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

GLM Ant_ipsi Ant_contra Lat_ipsi Lat_contra Med_ipsi Med_contra Post_ipsi Post_contra BY types_focal_to_bilateral_tonic_clonic
  /WSFACTOR=nuclei 4 Polynomial side 2 Polynomial 
  /METHOD=SSTYPE(3)
  /PLOT=PROFILE(types_focal_to_bilateral_tonic_clonic*nuclei types_focal_to_bilateral_tonic_clonic*side types_focal_to_bilateral_tonic_clonic*side*nuclei) TYPE=LINE ERRORBAR=SE(2) MEANREFERENCE=NO 
    YAXIS=AUTO
  /EMMEANS=TABLES(types_focal_to_bilateral_tonic_clonic*side) COMPARE(types_focal_to_bilateral_tonic_clonic) ADJ(BONFERRONI)
  /EMMEANS=TABLES(types_focal_to_bilateral_tonic_clonic*side) COMPARE(side) ADJ(BONFERRONI)
  /EMMEANS=TABLES(types_focal_to_bilateral_tonic_clonic*nuclei) COMPARE(types_focal_to_bilateral_tonic_clonic) ADJ(BONFERRONI)
  /EMMEANS=TABLES(types_focal_to_bilateral_tonic_clonic*nuclei) COMPARE(nuclei) ADJ(BONFERRONI)
  /EMMEANS=TABLES(types_focal_to_bilateral_tonic_clonic*nuclei*side) COMPARE(types_focal_to_bilateral_tonic_clonic) ADJ(BONFERRONI)
  /EMMEANS=TABLES(types_focal_to_bilateral_tonic_clonic*nuclei*side) COMPARE(nuclei) ADJ(BONFERRONI)
  /EMMEANS=TABLES(types_focal_to_bilateral_tonic_clonic*nuclei*side) COMPARE(side) ADJ(BONFERRONI)
  /PRINT=DESCRIPTIVE ETASQ 
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=nuclei side nuclei*side
  /DESIGN=types_focal_to_bilateral_tonic_clonic.

*TLE-HS vs. other TLE.
GLM Ant_ipsi Ant_contra Lat_ipsi Lat_contra Med_ipsi Med_contra Post_ipsi Post_contra BY mts
  /WSFACTOR=nuclei 4 Polynomial side 2 Polynomial 
  /METHOD=SSTYPE(3)
  /PLOT=PROFILE(mts*side mts*side*nuclei nuclei*mts mts*nuclei) TYPE=LINE ERRORBAR=SE(2) MEANREFERENCE=NO 
    YAXIS=AUTO
  /EMMEANS=TABLES(mts*side) COMPARE(mts) ADJ(BONFERRONI)
  /EMMEANS=TABLES(mts*side) COMPARE(side) ADJ(BONFERRONI)
  /EMMEANS=TABLES(mts*nuclei) COMPARE(mts) ADJ(BONFERRONI)
  /EMMEANS=TABLES(mts*nuclei*side) COMPARE(mts) ADJ(BONFERRONI)
  /EMMEANS=TABLES(mts*nuclei*side) COMPARE(nuclei) ADJ(BONFERRONI)
  /EMMEANS=TABLES(mts*nuclei*side) COMPARE(side) ADJ(BONFERRONI)
  /PRINT=DESCRIPTIVE ETASQ 
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=nuclei side nuclei*side
  /DESIGN=mts.

*Seizure free vs. Not seizure free.
USE ALL.
COMPUTE filter_$=(sf_latest = 1 | sf_latest =0).
VARIABLE LABELS filter_$ "sf_latest = 1 | sf_latest = 0 (FILTER)".
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

GLM Ant_ipsi Ant_contra Lat_ipsi Lat_contra Med_ipsi Med_contra Post_ipsi Post_contra BY sf_latest
  /WSFACTOR=nuclei 4 Polynomial side 2 Polynomial 
  /METHOD=SSTYPE(3)
  /PLOT=PROFILE(sf_latest*side sf_latest*side*nuclei) TYPE=LINE ERRORBAR=SE(2) MEANREFERENCE=NO 
    YAXIS=AUTO
  /EMMEANS=TABLES(sf_latest*side) COMPARE(sf_latest) ADJ(BONFERRONI)
  /EMMEANS=TABLES(sf_latest*side) COMPARE(side) ADJ(BONFERRONI)
  /EMMEANS=TABLES(sf_latest*nuclei*side) COMPARE(sf_latest) ADJ(BONFERRONI)
  /EMMEANS=TABLES(sf_latest*nuclei*side) COMPARE(nuclei) ADJ(BONFERRONI)
  /EMMEANS=TABLES(sf_latest*nuclei*side) COMPARE(side) ADJ(BONFERRONI)
  /PRINT=DESCRIPTIVE ETASQ 
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=nuclei side nuclei*side
  /DESIGN=sf_latest.

* Disease duration.
CORRELATIONS
  /VARIABLES=Ant_ipsi Lat_ipsi Med_ipsi Post_ipsi Ant_contra Lat_contra Med_contra Post_contra duration_epilepsy_fmri_years
  /PRINT=TWOTAIL NOSIG FULL
  /MISSING=PAIRWISE.

*NS - Volume correlation.
CORRELATIONS
  /VARIABLES=Ant_ipsi NS_Ant_ipsi_native_newharmo
  /PRINT=TWOTAIL NOSIG FULL
  /MISSING=PAIRWISE.
 CORRELATIONS
  /VARIABLES=Lat_ipsi NS_Lat_ipsi_native_newharmo
  /PRINT=TWOTAIL NOSIG FULL
  /MISSING=PAIRWISE.
CORRELATIONS
  /VARIABLES=Med_ipsi NS_Med_ipsi_native_newharmo
  /PRINT=TWOTAIL NOSIG FULL
  /MISSING=PAIRWISE.
CORRELATIONS
  /VARIABLES=Post_ipsi NS_Post_ipsi_native_newharmo
  /PRINT=TWOTAIL NOSIG FULL
  /MISSING=PAIRWISE.

CORRELATIONS
  /VARIABLES=Ant_contra NS_Ant_contra_native_newharmo
  /PRINT=TWOTAIL NOSIG FULL
  /MISSING=PAIRWISE.
CORRELATIONS
  /VARIABLES=Lat_contra NS_Lat_contra_native_newharmo
  /PRINT=TWOTAIL NOSIG FULL
  /MISSING=PAIRWISE.
CORRELATIONS
  /VARIABLES=Med_contra NS_Med_contra_native_newharmo
  /PRINT=TWOTAIL NOSIG FULL
  /MISSING=PAIRWISE.
CORRELATIONS
  /VARIABLES=Post_contra NS_Post_contra_native_newharmo
  /PRINT=TWOTAIL NOSIG FULL
  /MISSING=PAIRWISE.