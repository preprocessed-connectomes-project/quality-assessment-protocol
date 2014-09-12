#!/usr/bin/env Rscript

library(plyr)

#' Combines the different QC measures in different files into one (like ABIDE)

#' Read in the data first
#+ read
qc_anat         <- read.csv("FINAL_anat_qc.csv")
qc_func         <- read.csv("qc_filt_epi.csv")
qc_func_motion  <- read.csv("FINAL_func_qc.csv")

#' Select the desired columns
#+ select
qc_anat         <- subset(qc_anat, select=c("subject", "site", "cnr", "efc", "fber", "fwhm", "qi1", "snr"))
qc_func         <- subset(qc_func, select=c("subid", "session", "scan", "site.name", "efc", "fber", "fwhm", "gsr", "dvars", "quality"))
qc_func_motion  <- subset(qc_func_motion, select=c("Subject", "Scan", "site_x", "tr", "MeanFD_Jenkinson", "NumFD_greater_than_0.20", "PercentFD_greater_than_0.20"))

#' Make the columns a bit more consistent
#' subject, session for anat
#+ consistent-anat
tmp             <- strsplit(as.character(qc_anat$subject), "_") # check: all(sapply(tmp, length)==3)
qc_anat$subject <- as.character(sapply(tmp, function(x) x[1]))
qc_anat$session <- as.numeric(sapply(tmp, function(x) x[3]))

#' Have everything have a subject, session, and scan column for func
#' Note that we rename the scans for rockland sample to be 1-3
#+ consistent-func
# Func
qc_func         <- rename(qc_func, c(subid="subject"))
qc_func$subject <- sprintf("%07i", qc_func$subject)
qc_func$scan[qc_func$scan == 645] <- 1
qc_func$scan[qc_func$scan == 1400] <- 2
qc_func$scan[qc_func$scan == 2500] <- 3
# Func Motion
qc_func_motion  <- rename(qc_func_motion, c(Subject="subject", Scan="scan", site_x="site", MeanFD_Jenkinson="mean_fd", NumFD_greater_than_0.20="num_fd", PercentFD_greater_than_0.20="perc_fd"))
tmp             <- strsplit(as.character(qc_func_motion$subject), "_") # check: all(sapply(tmp, length)==3)
qc_func_motion$subject <- as.character(sapply(tmp, function(x) x[1]))
qc_func_motion$session <- as.numeric(sapply(tmp, function(x) x[3]))
tmp             <- as.character(qc_func_motion$scan)
tmp             <- sub("645_1", "1", tmp)
tmp             <- sub("1400_1", "2", tmp)
tmp             <- sub("2500_1", "3", tmp)
tmp2            <- strsplit(tmp, "_") # check: all(sapply(tmp2, length)==3)
qc_func_motion$scan <- as.numeric(sapply(tmp2, function(x) x[2]))

#' Merge the func together
#+ merge-func
tmp             <- merge(qc_func, qc_func_motion)
tmp             <- rename(tmp, c(site="orig_site", site.name="site"))
qc_func         <- tmp
# correct site
qc_anat2         <- ddply(qc_anat, .(site), function(x) {
  ind <- as.character(qc_func$subject) == as.character(x$subject[1])
  cat(sum(ind), as.character(x$site[1]), as.character(qc_func$site[ind][1]), "\n")
  if( sum(ind) != 0) {
    x$site <- as.character(qc_func$site[ind][1])
  }
  x
})
qc_anat2$site[qc_anat2$site == "MRN_1"] <- "MRN"
qc_anat2$site <- as.factor(qc_anat2$site)
qc_anat <- qc_anat2

#' Filter
#+ filter
qc_anat         <- subset(qc_anat, select=c("subject", "session", "site", "cnr", "efc", "fber", "fwhm", "qi1", "snr"))
qc_func         <- subset(qc_func, select=c("subject", "session", "scan", "site", "efc", "fber", "fwhm", "gsr", "dvars", "quality", "mean_fd", "num_fd", "perc_fd"))

#' There seem to be some duplicate values across sessions/scans
#' We'll fix it for the anatomical measures by removing the duplicates
#+ fix-duplicates
to_num <- function(x) as.numeric(as.character(x))
qc_anat <- ddply(qc_anat, .(site, subject), function(x) {
  if (nrow(x) > 1) {
    x2 <- as.matrix(x[,c("cnr", "efc", "fber", "fwhm", "qi1", "snr")])
    x2 <- apply(x2, 2, to_num)
    dmat <- as.matrix(dist(x2))
    nident <- sum(dmat[lower.tri(dmat)] == 0)
    if (nident > 1) {
      cat(as.character(x$site[1]), as.character(x$subject[1]), nident, "\n")
      bad_inds <- which(dmat[-1,1] == 0) + 1
      x <- x[-bad_inds,]
    }
  }
  x
})
# not doing anything about this for now
d_ply(qc_func, .(site, subject), function(x) {
  if (nrow(x) > 1) {
    x2 <- as.matrix(x[,c("efc", "fber", "fwhm", "gsr", "dvars", "quality", "mean_fd", "num_fd", "perc_fd")])
    x2 <- apply(x2, 2, to_num)
    dmat <- as.matrix(dist(x2))
    nident <- sum(dmat[lower.tri(dmat)] == 0)
    if (nident > 1) {
      cat(as.character(x$site[1]), as.character(x$subject[1]), nident, "\n")
      print(dmat)
    }
  }
})

#' Select removal of sites
#' IBA-TRT might have some duplicates so will remove
#' IPCAS 4 has duplicates in the functional data domain so will remove
#+ rm-sites
qc_anat <- subset(qc_anat, !(site %in% c("IBATRT", "IPCAS 4")))
qc_func <- subset(qc_func, !(site %in% c("IBATRT", "IPCAS 4")))

#' Add on
#+ addon
# Anat
inds            <- !(names(qc_anat) %in% c("subject", "session", "site"))
names(qc_anat)[inds] <- paste("anat", names(qc_anat)[inds], sep="_")
# Func
inds            <- !(names(qc_func) %in% c("subject", "session", "scan", "site"))
names(qc_func)[inds] <- paste("func", names(qc_func)[inds], sep="_")

#' Sort Columns
#+ sort
qc_anat <- qc_anat[with(qc_anat, order(site, subject, session)),]
qc_func <- qc_func[with(qc_func, order(site, subject, session, scan)),]

#' Save
#+ save
write.csv(qc_anat, file="../corr_anat.csv")
write.csv(qc_func, file="../corr_func.csv")


#"anat_cnr","anat_efc","anat_fber","anat_fwhm","anat_qi1","anat_snr",
#"func_efc","func_fber","func_fwhm","func_ghost_x","func_ghost_y",
#"func_dvars","func_outlier","func_quality",
#"func_mean_fd","func_num_fd","func_perc_fd",
#"func_gsr"
