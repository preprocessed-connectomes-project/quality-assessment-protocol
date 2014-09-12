#!/usr/bin/env Rscript

library(plyr)

#' Combines the different QC measures in different files into one (like ABIDE)

#' Read in the data first
#+ read
qc              <- read.csv("abide_phenotypic_V1_0b_preprocessed1.csv")

#' Format subjects excluded from ABIDE paper as rater 4
#+ rater4
rater_4 <- qc$SUB_IN_SMP
qc$qc_rater_4 <- factor(rater_4, levels=c(0,1), labels=c("fail", "OK"))

#' Get the columns to extract
#+ cols
id_cols         <- c("subject", "SITE_ID") 
pheno_cols      <- c("DX_GROUP", "DSM_IV_TR", "AGE_AT_SCAN", "SEX", "HANDEDNESS_CATEGORY", "FIQ")
qc_cols         <- names(qc)[grep("^qc_", names(qc))]
qc_anat_cols    <- c(qc_cols[grep("^qc_rater", qc_cols)], qc_cols[grep("^qc_anat_rater", qc_cols)])
qc_func_cols    <- c(qc_cols[grep("^qc_rater", qc_cols)], qc_cols[grep("^qc_func_rater", qc_cols)])
anat_qc_cols    <- names(qc)[grep("^anat_", names(qc))]
func_qc_cols    <- names(qc)[grep("^func_", names(qc))]

#' Extract
#+ extract
qc_anat         <- subset(qc, select=c(id_cols, pheno_cols, anat_qc_cols, qc_anat_cols))
qc_func         <- subset(qc, select=c(id_cols, pheno_cols, func_qc_cols, qc_func_cols))

#' Rename
#+ rename
# step 1
#qc_anat         <- sub("anat_", "", names(qc_anat))
#qc_func         <- sub("func_", "", names(qc_func))
# step 2
qc_anat         <- rename(qc_anat, c(SITE_ID="site", DX_GROUP="dx", DSM_IV_TR="dsm", AGE_AT_SCAN="age", SEX="sex", HANDEDNESS_CATEGORY="handedness", FIQ="iq"))
qc_func         <- rename(qc_func, c(SITE_ID="site", DX_GROUP="dx", DSM_IV_TR="dsm", AGE_AT_SCAN="age", SEX="sex", HANDEDNESS_CATEGORY="handedness", FIQ="iq"))

#+ save
write.csv(qc_anat, file="../abide_anat.csv")
write.csv(qc_func, file="../abide_func.csv")


