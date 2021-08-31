# patchseq_processing


Human_MTG       = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/human/human_patchseq_MTG_20210818"
Mouse_CTX_HIP   = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_CTX-HIP_20210818_collapsed40_cpm"
Mouse_ForeBrain = "/allen/programs/celltypes/workgroups/rnaseqanalysis/changkyul/genomic_pipelines/patch-seq/mapping"
Mouse_V1_ALM    = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20210818_collapsed40_cpm"

system(paste("cp ", file.path(Human_MTG, "build_feather_files_server.R"), " Human_MTG/."))
system(paste("cp ", file.path(Mouse_CTX_HIP, "build_feather_files_server.R"), " Mouse_CTX-HIP/."))
system(paste("cp ", file.path(Mouse_V1_ALM, "build_feather_files_server.R"), " Mouse_V1-ALM/."))
system(paste("cp ", file.path(Mouse_ForeBrain, "patchseq_mapping.mouseFB.R"), " Mouse_ForeBrain/."))
system(paste("cp ", file.path(Mouse_ForeBrain, "patchseq_qc.mouseFB.R"), " Mouse_ForeBrain/."))
