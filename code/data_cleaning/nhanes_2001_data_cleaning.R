# NHANES data cleaning for pairwise KMR
# Nate Wiecha
# Using approximate data set from Gibson et al (2019), "An overview of methods..."

# setwd("H:/Pairwise_Additive_GP_Tests/data/nhanes_2001")
setwd("~/GitHub/Pairwise_Additive_GP_Tests/data/nhanes_2001")
rm(list=ls())

library(haven)
library(tidyverse)

# Telomere measurements
telo_dat <- read_xpt("TELO_B.xpt") %>%
  mutate(ln_ltl=log(TELOMEAN)) %>%
  select(-TELOSTD, -TELOMEAN) %>%
  drop_na()

# exposures
exp_dat <- read_xpt("L28POC_B.xpt") %>%
  select(SEQN, PCB74=LBX074LA,
         PCB99=LBX099LA,
         PCB138=LBX138LA,
         PCB153=LBX153LA,
         PCB170=LBX170LA,
         PCB180=LBX180LA,
         PCB187=LBX187LA,
         PCB194=LBX194LA,
         PCB126=LBXPCBLA,
         PCB169=LBXHXCLA,
         PCB118=LBX118LA,
         hxcdd123678=LBXD03LA,# 1,2,3,6,7,8-hxcdd
         hpcdd1234678=LBXD05LA, # 1,2,3,4,6,7,8-hpcdd
         ocdd12346789=LBXD07LA , # 1,2,3,4,6,7,8,9-ocdd
         pcndf23478=LBXF03LA, # 2,3,4,7,8-pncdf 
         hxcdf123478=LBXF04LA, # 1,2,3,4,7,8-hxcdf 
         hxcdf123678=LBXF05LA, # 1,2,3,6,7,8-hxcdf
         hcxdf1234678=LBXF08LA # 1,2,3,4,6,7,8-hxcdf
         )  %>%
  mutate(ln_PCB74=log(PCB74),
         ln_PCB99=log(PCB99),
         ln_PCB138=log(PCB138),
         ln_PCB153=log(PCB153),
         ln_PCB170=log(PCB170),
         ln_PCB180=log(PCB180),
         ln_PCB187=log(PCB187),
         ln_PCB194=log(PCB194),
         ln_PCB126=log(PCB126),
         ln_PCB169=log(PCB169),
         ln_PCB118=log(PCB118),
         ln_hxcdd123678=log(hxcdd123678),# 1,2,3,6,7,8-hxcdd
         ln_hpcdd1234678=log(hpcdd1234678), # 1,2,3,4,6,7,8-hpcdd
         ln_ocdd12346789=log(ocdd12346789) , # 1,2,3,4,6,7,8,9-ocdd
         ln_pcndf23478=log(pcndf23478), # 2,3,4,7,8-pncdf 
         ln_hxcdf123478=log(hxcdf123478), # 1,2,3,4,7,8-hxcdf 
         ln_hxcdf123678=log(hxcdf123678), # 1,2,3,6,7,8-hxcdf
         ln_hcxdf1234678=log(hcxdf1234678)) %>% # 1,2,3,4,6,7,8-hxcdf)
  select(-PCB74, -PCB99, -PCB138, -PCB153, -PCB170,
         -PCB180, -PCB187, -PCB194,  -PCB126, -PCB169,
         -PCB118, -hxcdd123678, -hpcdd1234678, -ocdd12346789,
         -pcndf23478, -hxcdf123478, -hxcdf123678, -hcxdf1234678) %>%
  drop_na()

# Impute for exposure values not detected
nondetects <- read_xpt("L28POC_B.xpt") %>%
  select(SEQN, PCB74=LBD074LC,
         PCB99=LBD099LC,
         PCB138=LBD138LC,
         PCB153=LBD153LC,
         PCB170=LBD170LC,
         PCB180=LBD180LC,
         PCB187=LBD187LC,
         PCB194=LBD194LC,
         PCB126=LBDPCBLC,
         PCB169=LBDHXCLC,
         PCB118=LBD118LC,
         hxcdd123678=LBDD03LC,# 1,2,3,6,7,8-hxcdd
         hpcdd1234678=LBDD05LC, # 1,2,3,4,6,7,8-hpcdd
         ocdd12346789=LBDD07LC , # 1,2,3,4,6,7,8,9-ocdd
         pcndf23478=LBDF03LC, # 2,3,4,7,8-pncdf 
         hxcdf123478=LBDF04LC, # 1,2,3,4,7,8-hxcdf 
         hxcdf123678=LBDF05LC, # 1,2,3,6,7,8-hxcdf
         hcxdf1234678=LBDF08LC # 1,2,3,4,6,7,8-hxcdf
  ) %>%
  drop_na()

round(colMeans(nondetects),2)

# NHANES already imputes using LOQ/sqrt(2)

# BMI
bmi_dat <- read_xpt("BMX_B.xpt") %>%
  select(SEQN, bmi=BMXBMI) %>%
  mutate(under=as.numeric(bmi < 25),
         # normal=bmi>=25 & bmi < 30,
         over = as.numeric(bmi >=30) )%>%
  select(-bmi) %>%
  drop_na()

# demographic
demo_dat <- read_xpt("DEMO_B.xpt") %>%
  select(SEQN, gender=RIAGENDR,
         age=RIDAGEYR,
         educ=DMDEDUC2,
         raceeth=RIDRETH1,
         ) %>%
  mutate(age2=age^2,
         female=as.numeric(gender==2),
         # lesshigh = educ<=2,
         high = as.numeric(educ==3),
         some_coll = as.numeric(educ==4),
         college=as.numeric(educ==5),
         # nh_white=raceeth==3,
         nh_black=as.numeric(raceeth==4),
         mexam=as.numeric(raceeth==1),
         othrace=as.numeric(raceeth==5)) %>%
  select(-educ, -raceeth, -gender) %>%
  drop_na()

# cotinine 
cot_dat <- read_xpt("L06_B.xpt") %>%
  select(SEQN, LBXCOT) %>%
  mutate(ln_cot = log(LBXCOT)) %>%
  select(-LBXCOT) %>%
  drop_na()

# blood cells
blood_dat <- read_xpt("L25_B.xpt") %>%
  select(SEQN, white_ct = LBXWBCSI, 
         lymph_pct = LBXLYPCT ,# pct lymphocytes
         mon_pct = LBXMOPCT , # pct monocytes
         neut_pct = LBXNEPCT,# pct (segmented) neutrophils
         eos_pct = LBXEOPCT, # pct Eosinophils  
         bas_pct = LBXBAPCT )  %>% # pct basophils
drop_na()

############################################################################
#                  Construct analytic sample                               #
############################################################################

dat <- inner_join(telo_dat, exp_dat) %>%
  inner_join(demo_dat) %>%
  inner_join(blood_dat) %>%
  inner_join(bmi_dat) %>%
  inner_join(cot_dat)
head(dat)

yy <- dat$ln_ltl
zz <- dat[,colnames(exp_dat)[-1]] %>%
  as.matrix() %>%
  t()

xx <- cbind(1,dat[,c(colnames(blood_dat)[-1], colnames(bmi_dat)[-1], 
             colnames(cot_dat)[-1], colnames(demo_dat)[-1])])
head(xx)

xx <- xx %>%
  as.matrix() %>%
  t()

save(dat, yy, zz, xx, file="nhanes_2001.Rdata")
