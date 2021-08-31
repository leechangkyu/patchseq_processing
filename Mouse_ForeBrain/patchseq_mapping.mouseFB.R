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
species             = "mouse"
taxonomy_str        = "ForeBrain"

flag_data_QC        = TRUE #TRUE # FALSE
flag_KL_div         = TRUE
flag_SHINY          = FALSE
mapping_methods     = c("Cor", "Seurat", "Tree")
script_rep          = "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_script_repository/"


# Human MTG
if (species=="human" && taxonomy_str =="MTG") {
   ref_version_str     = "CCN202107190"
   mapping_version_str = "standard_mapping_tool_2021"
 
   ref_version_str     = "CCN201908210"
   mapping_version_str = "2021"
   mapping_method      = mapping_methods[c(1,2,3)]
   flag_KL_div         = TRUE
   patchseq_batch <- batch_date <- "20210602_RSC-11-278_human"
   rundate             = unlist(strsplit(patchseq_batch, "_"))[1]
   Niter               = 100 

   # patchseq data
   patchseq_type       = "feather" # "Robj"
   patchseq_dir        = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/human/human_patchseq_MTG_20210604/"
   patchseq_robj_dir   = "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Human/patchseq/R_Object/"
   # output
   refFolder_pseq      = file.path(shiny_dir_human, paste(species, "patchseq", taxonomy_str, rundate, "u", sep="_"))


   #---------- no change from here
   refFolder           = "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/reference_taxonomies/human/MTG_great_ape_neuron"
   #refFolder            = "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/human/reference_taxonomies/MTG_great_ape_neuron/"
   # taxonomy specific : anno.feather or  (annotation file & dendrogram)
   fdir                = "//allen/programs/celltypes/workgroups/rnaseqanalysis/Nik/Analyses_for_great_apes_paper/Shiny_obj/human/"
   regions             = c() #or mapping everything 
   facs.samp.dat.FN    = "//allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Human/facs/R_Object/Archive/20180918_RSC-004-171_human_star2.0_samp.dat.Rdata"
   # if previous mapping available
   lastmap_dir         = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/human/human_patchseq_MTG_20210504/"
   lastmap_FN          = "mapping.df.lastmap.csv"
}

===========================
if (species=="mouse" && taxonomy_str =="ForeBrain") {
   mapping_method      = mapping_methods[c(1,2,3)]
   patchseq_batch <- batch_date <- "20210818_BT014-RSC-286_mouse"
   patchseq_dir        = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_ForeBrain_20210818_collapsed40_cpm/"

   # output
   refFolder_pseq      = getwd() #file.path(shiny_dir_mouse, paste(species, "patchseq", taxonomy_str, rundate, "again", sep="_"))

   rundate             = unlist(strsplit(patchseq_batch, "_"))[1]

   ref_version_str     = "20201204_mouse_ForeBrain"
   mapping_version_str = "20201204_mouse_ForeBrain_v1"
   Niter               = 100 

   # patchseq data
   patchseq_type       = "Robj"
   patchseq_robj_dir   = "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Mouse/patchseq/R_Object/"


   #---------- no change from here
   refFolder           = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/reference_taxonomies/mouse/ForeBrain"
   # taxonomy specific : anno.feather or  (annotation file & dendrogram)
   #fdir                = "/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/SmartSeq_cells/foreBrain/"
   fdir                = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_ForeBrain_20201203"
   regions             = c() #or mapping everything 
   facs.samp.dat.FN    = ""
   # if previous mapping available
   lastmap_dir        = ""
   lastmap_FN         = ""
   script_rep         = "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_script_repository/"
}


# Mouse V1
if (species=="mouse" && taxonomy_str =="V1-ALM") {
   mapping_method      = mapping_methods[c(1,2,3)]
   patchseq_batch <- batch_date <- "20210331_BT014-RSC-273_mouse"
   rundate             = unlist(strsplit(patchseq_batch, "_"))[1]
   ref_version_str     = "20180626_mouse_VISp_v93_bp40"
   mapping_version_str = "2018092l_mouse_conf70_marker1"
   Niter               = 100 
   # patchseq data
   patchseq_type       = "Robj"
   patchseq_dir        = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20210604_collapsed40_cpm/"
   patchseq_robj_dir   = "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Mouse/patchseq/R_Object/"

   # output
   refFolder_pseq      = file.path(shiny_dir_human, paste(species, "patchseq", taxonomy_str, rundate, "u", sep="_"))

   #---------- no change from here
   refFolder           = "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/reference_taxonomies/mouse/VISp"
   # taxonomy specific : anno.feather or  (annotation file & dendrogram)
   fdir                = "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/"
   regions             = c("VISp","TCx","FCx","MOp","TEa","HIPCA1") #or mapping everything 
   facs.samp.dat.FN    = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_facs_Rdata/20181009_RSC-11-176_mouse_star_samp.dat_updated_inj.Rdata"
   # if previous mapping available
   lastmap_dir        = "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20210403_collapsed40_cpm/"
   lastmap_FN         = "mapping.df.with.bp.40.lastmap.csv"
}


if (species=="mouse" && taxonomy_str =="CTX-HPF") {
   mapping_method      = mapping_methods[c(1,2,3)]
   patchseq_batch <- batch_date <- "20200318_BT014-RSC-252"
   rundate             = unlist(strsplit(patchseq_batch, "_"))[1]
   ref_version_str = "test_ref_v"
   mapping_version_str = "test_mapping_v"
   Niter             = 100 

   # patchseq data
   patchseq_type     = "Robj"
   patchseq_dir <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_CTX-HIP_20200318_collapsed40_cpm/"
   patchseq_robj_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Mouse/patchseq/R_Object/"
   
   #output
   refFolder_pseq      = file.path(shiny_dir_human, paste(species, "patchseq", taxonomy_str, rundate, "u", sep="_"))
   #ref.data:
   fdir              = "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/tomes/multi_method/CTX_HIP/"
   dend_Ver          = ""
   dend_FN           = "dend.AIT5.1.rda"
   cluster_anno_FN   = ""
   regions           = "" #or mapping everything 
   
   # if previous mapping available
   lastmapdir        = "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_CTX-HIP_20200313_collapsed40_cpm/"
   lastmap_FN        = "mapping.df.with.bp.40.lastmap.csv"
}

##########################################################################
options(error=recover)
ref_version_str     = paste(species, taxonomy_str)
mapping_version_str = paste(rundate, species)
mapping_dir         = refFolder_pseq

if (!dir.exists(refFolder_pseq)) dir.create(refFolder_pseq)
setwd(refFolder_pseq)
print(refFolder_pseq)
mm=30 #15,20,25,30
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
print(" ==  prepare patchseq data")
source(file.path(my_script_rep, "p0_patchseq_data_anno_prep.R"))

print(" ==  patchseq data qc")
if (flag_data_QC) source(file.path(my_script_rep, "p1_patchseq_data_QC.R"))

#browser()
# run mapping
print(" == Mapping")
if ("Cor"    %in% mapping_method) source(file.path(my_script_rep, "p2_patchseq_data_COR_mapping.R"))
if ("Tree"   %in% mapping_method) source(file.path(my_script_rep, "p3_patchseq_data_tree_mapping.R"))
if ("Seurat" %in% mapping_method) source(file.path(my_script_rep, "p4_patchseq_data_Seurat_mapping.R"))

# run KL_div for confidence of mapping
if (flag_KL_div) source(file.path(my_script_rep, "p5_patchseq_data_tree_KL.R"))

print(" == mapping done : Compare with previous mappings!!!")
#source(file.path(my_script_rep,"debug.compare.R"))
print(" Comparison DONE! ")

if (flag_SHINY) {
   print(" == Updating anno with class/membership/mapping for SHINY")
   source(file.path(my_script_rep, "prepare_anno_dend_data.R"))
   source(file.path(my_script_rep, "update_anno_w_mapping.R"))
   
   print(" == creating shiny links")
   setwd(mapping_dir)
   source(file.path(my_script_rep,"create_shiny_links_per_cell.R"))
   
   print(" == moving feathers to SC-multimodal viewer folder") # Jeremy's script
   source(file.path(my_script_rep,"move2SC_MMviewer.R"))
}
print(" == Shiny DONE!!!")
print(" == Compare with previous mappings!!!")
#source(file.path(my_script_rep,"debug.compare.R"))

 
