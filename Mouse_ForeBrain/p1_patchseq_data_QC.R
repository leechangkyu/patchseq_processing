##################################################################################
#
# p1_patchseq_data_QC.R
#
#

# Load libraries
library(feather)
library(dplyr)
library(patchseqtools)  # devtools::install_github("AllenInstitute/patchseqtools")
library(patchSeqQC)     # devtools::install_github('PavlidisLab/patchSeqQC')
library(scrattch.io)    # devtools::install_github("AllenInstitute/scrattch"); scrattch::install_scrattch_deps(); scrattch::install_scrattch()
library(scrattch.hicat)
options(stringsAsFactors = FALSE)


##
# Read in the patch-seq data/metadata and perform the PatchseqQC
##
refMarkerFile  = "QC_markers.RData"
load(file.path(refFolder, refMarkerFile))

print("Format the reference and patch-seq data")
## -- NOTE: relevant reference data and type assignments are stored in refMarkerFile
tmp                 <- cpmQC # countsQC
rownames(tmp)       <- make.names(rownames(tmp))
facs_df             <- as.data.frame(t(tmp[allMarkers,])+1)
facs_df$sample_id   <- rownames(facs_df)
facs_df$major_type  <- as.character(classBr)     # Defined in previous section
facs_df$contam_type <- as.character(subclassF)   # Defined in previous section
facs_df             <- facs_df %>% filter(major_type != "Low Quality")

# from p0_patchseq_data_prep.R
if (patchseq_type=="Robj") {
   tmp           <- cpmPat
} else {
   tmp           <- cpm(datPat) #
}
#rownames(tmp)   <- make.names(rownames(tmp))
allMarkers       <- intersect(allMarkers,rownames(tmp))  # To account for differences in transcriptome
sdx              <- match(annoPat$sample_id, colnames(tmp))
pat_df           <- as.data.frame(t(tmp[allMarkers,sdx])+1)
pat_df$sample_id <- rownames(pat_df)


print("Define which type each patch-seq cell is assigned to, based on maximal marker expression.")  
nm          <- names(markers)
isOn        <- substr(nm,nchar(nm)-2,nchar(nm))=="_on"
useThese    <- nm[isOn&(!is.element(nm,paste0(nm,"_on")))]
subclassDat <- calcContamAllTypes(pat_df, markers[useThese])  # Identify subclass based on marker gene expression
subclass    <- colnames(subclassDat)[subclassDat %>% apply(1,which.max)]
subclass    <- gsub("_on","",subclass)

pat_df$contam_type <- subclass  
pat_df$major_type  <- as.character(classBr)[match(pat_df$contam_type,as.character(subclassF))]
pat_df$contam_type <- paste0(pat_df$contam_type,"_on")
pat_df             <- pat_df %>% filter(major_type != "Low Quality")


print("Remove genes not patch-seq transcriptome (ideally this wouldn't do anything)")
for (i in 1:length(markers)) markers[[i]] <- intersect(markers[[i]],allMarkers)
if ("Low Quality" %in% names(markers)) markers = markers[-match("Low Quality", names(markers))]  
print("Calculate the QC metrics")

qcMetrics <- calculatePatchSeqQCMetrics2(pat_df,facs_df,markers)
qcMetrics$quality_score[is.na(qcMetrics$quality_score)] = 0.0

print("Set NMS>0.4 flag and determine most contaminated type")
qcMetrics$Norm_Marker_Sum.0.4 <- c(TRUE,FALSE)[(qcMetrics$marker_sum_norm<0.40)+1]
cls               <- sort(setdiff(classBr,c("GABAergic","Glutamatergic", "Low Quality")))
contaminationType <- cls[apply(qcMetrics[,cls],1,which.max)]
qcMetrics$contaminationType   <- contaminationType


print("Update annotations")
anno <- read_feather(file.path(refFolder_pseq,"anno.feather"))

cn      <- c("quality_score","marker_sum_norm","Norm_Marker_Sum.0.4","contaminationType","contam_sum")
annoNew <- anno[match(pat_df$sample_id, anno$sample_id),]
for (i in 1:length(cn))      # Remove duplicate column names, if any
  annoNew <- annoNew[,!grepl(cn[i],colnames(annoNew))]

annoNew <- cbind(annoNew,qcMetrics[,cn])

annoNew$quality_score = pmax(annoNew$quality_score,0)
annoNew$marker_sum_norm = pmax(annoNew$marker_sum_norm,0)
annoNew$contam_sum      = pmax(annoNew$contam_sum,0)
annoNew <- annoNew %>% annotate_num("quality_score", scale = "linear")    
annoNew <- annoNew %>% annotate_num("marker_sum_norm", scale = "linear")  
annoNew <- annoNew %>% annotate_num("contam_sum", scale = "linear") 
annoNew <- annoNew %>% annotate_cat("Norm_Marker_Sum.0.4")
annoNew <- annoNew %>% annotate_cat("contaminationType")

if (file.exists(file.path(refFolder_pseq,"anno.feather"))){
  file.remove(file.path(refFolder_pseq,"anno.feather"))
}
write_feather(annoNew,file.path(refFolder_pseq,"anno.feather"))

#print(paste(nnoNew$quality_score"result are in ", refFolder_pseq))

