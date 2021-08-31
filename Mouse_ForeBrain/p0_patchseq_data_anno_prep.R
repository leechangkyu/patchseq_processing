########################################################################
# - p0_prep_patchseq_data
#
# Read in patchserq data and set up 
#     datPat, samp.dat, query.dat.norm (query.dat.cells)
########################################################################

##
#  directory set up
##
if (!patchseq_type %in% c("Robj", "feather")) {
   print(paste0( "Available patchseq_type : 'Robj', 'feather'
                  you entered", patchseq_type))
   stop()
} else {
   if (patchseq_type == "Robj") {
      print("- Robj type : samp & cpm ready")
      load(file.path(patchseq_robj_dir, paste0(patchseq_batch, "_patchseq_star2.0_samp.dat.Rdata")))
      load(file.path(patchseq_robj_dir, paste0(patchseq_batch, "_patchseq_star2.0_cpm.Rdata")))

      annoPat <- Samp.datp <- samp.dat
      annoPat$sample_id <- Samp.datp$sample_id <-  as.character(annoPat$patched_cell_container)

      query.dat = cpmR[,as.character(samp.dat$exp_component_name)]
      colnames(query.dat)=as.character(samp.dat$patched_cell_container)
      datPat = log2(as.matrix(query.dat+1))
      cpmPat = query.dat
   }
   if (patchseq_type == "feather") {
      print("- feather type : anno & data ready ")
      Samp.datp   = read_feather(file.path(patchseq_dir,"anno.feather"))
      Expr.datp   = read_feather(file.path(patchseq_dir,"data.feather"))   # FPKM
      common.sample_id = intersect(Expr.datp$sample_id,Samp.datp$sample_id)
      Samp.datp   = Samp.datp[match(common.sample_id, Samp.datp$sample_id),]
      Expr.datp   = Expr.datp[match(common.sample_id, Expr.datp$sample_id),]
      #Samp.datp$subclass_label    <- "cell"      # I don't know why we need this or the next line
      #Samp.datp$dendcluster_color <- Samp.datp$cluster_color
      samp.dat <- annoPat <- Samp.datp
      datPat   = as.matrix(Expr.datp[,names(Expr.datp)!="sample_id"])
      rownames(datPat)= annoPat$sample_id
      datPat   = logCPM(t(datPat))
   }
}

#load(file=file.path(refFolder, "keep_ref.KL.mapping.rda") )
load(file=file.path(refFolder, "reference.rda"))
ref.genename = rownames(reference$cl.dat)

print("select reference data genes in patchseq")
#common.genes = intersect(select.markers$markers, intersect(ref.genename, rownames(datPat)))
common.genes = intersect(ref.genename, rownames(datPat))
gdx=match(common.genes, rownames(datPat))
print(paste(length(gdx), " selected genes out of ", length(ref.genename)))

print("filter by region of interest")
if (length(regions) > 0) {
   keepcells = which(samp.dat$Region %in% regions) #| as.numeric(gsub("RSC-","",samp.dat$batch_vendor_name))>=244) #RSC-244 is the latest batch to have
} else {
   keepcells = 1:nrow(samp.dat)
}

#
# query.dat.norm & query.dat.cell : 
#

print("      Make query data set ready with gene/region     ")
query.dat.norm=as.matrix(datPat[gdx, keepcells, drop = F])
rownames(query.dat.norm) = rownames(datPat)[gdx]
query.dat.cells <- colnames(query.dat.norm) <- colnames(datPat)[keepcells]

print(" query data with all genes : keepcells are set in p0_patchseq_data_prep.R")
query.allgene.norm=as.matrix(datPat[, keepcells])
query.allgene.cells <- colnames(query.allgene.norm) <- colnames(datPat)[keepcells]

print(" # == patchseq data is ready ")
print(" # (query.dat.norm, query.dat.cells) ")
save(query.dat.norm, query.dat.cells, file=file.path(refFolder_pseq, "query.dat.norm.cells.rda"))
save(query.allgene.norm, query.allgene.cells, file=file.path(refFolder_pseq, "query.allgene.norm.cells.rda"))
print(" # (query.allgene.norm, query.allgene.cells) ")


print("Create anno.feather file. ")
if (file.exists(file.path(refFolder_pseq,"anno.feather"))){
  file.remove(file.path(refFolder_pseq,"anno.feather"))
}
write_feather(annoPat,file.path(refFolder_pseq,"anno.feather"))

