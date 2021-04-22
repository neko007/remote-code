DSET /home/zzhzhao/data/CMORPH_sta/2001/SURF_CLI_CHN_MERGE_CMP_PRE_HOUR_GRID_0.10-%y4%m2%d2%h2.grd
*
UNDEF -999.0
*
OPTIONS   little_endian  template
*
TITLE  China Hourly Merged Precipitation Analysis
*
xdef  700 linear  70.0  0.10
*
ydef  440 linear  15.0  0.10 
*
ZDEF     1 LEVELS 1  
*
TDEF 8760 LINEAR 00:00Z01Jan2001 1hr 
*
VARS 2                           
crain      0 99  CH01   combined analysis (mm/Hour)
gsamp      0 99  CH02   gauge numbers
ENDVARS
