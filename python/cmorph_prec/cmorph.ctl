DSET /home/zzhzhao/data/CMORPH_sta/2008/SEVP_CLI_CHN_MERGE_CMP_PRE_HOUR_GRID_0.10-%y4%m2%d2%h2.grd
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
TDEF 8784 LINEAR 00Z01jan2008 1hr 
*
VARS 2                           
crain      1 99  CH01   combined analysis (mm/Hour)
gsamp      1 99  CH02   gauge numbers
ENDVARS
