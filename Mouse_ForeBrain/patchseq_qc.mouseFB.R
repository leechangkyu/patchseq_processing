##########################################################################
# patchseq_mapping 
# !!! ver1
# !!! ssh hpc-login 
# !!! qsub -I -q celltypes -l walltime=100:00:00 -l mem=500G
# !!! singularity run docker://alleninst/mapping_on_hpc
# !!! R
#
# INPUT/params
# - species           : human/mouse
# - taxonomy_str      : 
# - ref_version_str   : reference id
# - mapping_version_str : mapping task version id
# - refFolder         : reference built from "taxonomy_str" taxonomy 
#                     0) "anno.feather/data.father"
#                     1) "reference.rda"
#                     2) "QC.markers.RData"
#                     3) "membership_information_reference.rda" (memb.ref & map.ref)
#                     4) "keep.KL.mapping.rda"  (saved for KL mapping)
#                     5) "norm.dat.rda" (needed for Seurat mapping)
#                     6) "dend.RData"/"dend.reference.RData"
# - fdir              : FACs folder
#                       1) dend_file : "VISp_dend_20180626/V1.dend.with.gen.cl.with.bp.40.rda"
#                       2) cluster_annotation file "cluster_anno_zy.tsv"
# - regions           : region of interest for mapping
# - facs.samp.dat.FN  : facs samp.date
# - patchseq_batch <- batch_date : batch id
# - patchseq_type     : patchseq type "Robj" / "feather"
# - patchseq_robj_dir : patchseq Robject folder 
# - lastmap_dir       : folder for the last mapping 
#                        1) lastmap_FN : (e.g. "mapping.df.with.bp.40.lastmap.csv")
# - Niter             : number of runs in tree_mapping (default=100)
#
# OUTPUT
# - refFolder_pseq    : reference folder expected to have following files for pathseq reference/mapping
# - mapping_dir       : patchseq mapping folder
#   (shiny.dir)
# - scmmviewerFolder  : single cell multimodal viewer folder
#
#########################################################################


library(dplyr)
library(scrattch.hicat)
library(scrattch.io)
library(scrattch.vis)
library(feather)
library(dendextend)
options(stringsAsFactors = F)


# shiny page
shiny_dir_mouse = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star"
shiny_dir_human = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/human"
my_script_rep   = "/allen/programs/celltypes/workgroups/rnaseqanalysis/changkyul/genomic_pipelines/patch-seq/mapping"

#
# INPUT/OUTPUT/Params
#
# input for patchseq data
patchseq_batch <- batch_date <- "20210818_BT014-RSC-286_mouse"
patchseq_dir        <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_Forebrain_20210818_collapsed40_cpm/"
#
# output
refFolder_pseq      <- getwd() #file.path(shiny_dir_mouse, paste(species, "patchseq", taxonomy_str, rundate, "again", sep="_"))

########################################################################
# NO CHANGE FROM HERE
rundate             = unlist(strsplit(patchseq_batch, "_"))[1]
patchseq_type       = "Robj"
patchseq_robj_dir   = "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Mouse/patchseq/R_Object/"


#---------- no change from here
refFolder           = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/reference_taxonomies/mouse/ForeBrain"

# taxonomy specific : anno.feather or  (annotation file & dendrogram)
regions             = c() #or mapping everything 
facs.samp.dat.FN    = ""

ref_version_str     = "20201204_mouse_ForeBrain"
mapping_version_str = "20201204_mouse_ForeBrain_v1"
Niter               = 100 
# if previous mapping available
lastmap_dir        = ""
lastmap_FN         = ""
script_rep         = "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_script_repository/"

##########################################################################
options(error=recover)
ref_version_str     = paste(species, taxonomy_str)
mapping_version_str = paste(rundate, species)
mapping_dir         = refFolder_pseq

if (!dir.exists(refFolder_pseq)) dir.create(refFolder_pseq)
setwd(refFolder_pseq)
print(refFolder_pseq)
mm=30 
if (!exists("reference")) {
   if (file.exists(file.path(refFolder,"references.rda"))) {
      load(file.path(refFolder,"references.rda"))
      reference = references[[as.character(mm)]]
   } else {
      load(file.path(refFolder,"reference.rda"))
   }
}
##########################################################################
# MAPPING
# prepare patchsests(file.path(refFolder, "cl.final.rda")) data

flag.OverWrite = FALSE
#print(" ==  prepare patchseq data")
source(file.path(my_script_rep, "p0_patchseq_data_anno_prep.R"))

print(" ==  patchseq data qc")
source(file.path(my_script_rep, "p1_patchseq_data_QC.R"))

 
