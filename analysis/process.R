# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,Rmd,R
#     text_representation:
#       extension: .R
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.7.1
#   kernelspec:
#     display_name: R
#     language: R
#     name: ir
# ---

# # 1. Pre-reqs
# We may need to install the rrscale and memnorm packages
# +
#install.packages('../r_packages/rrscale_rpkg/rrscale_1.0.3.tar.gz', repos=NULL,type="source")
#install.packages('../r_packages/mema_norm_rpkg/memanorm_0.0.0.9004.tar.gz', repos=NULL,type="source")
# -

library('rrscale')
library('memanorm')
source('processing_scripts/util.R')

# and then we'll set up a directory for the processed data
# +
# identifying tag for the data
tag = "MCF10A"

# directory to put processed data
ddir ='processed_data/'%+%tag%+%"/"
dir.create(ddir,showWarning=FALSE,recursive=TRUE)
# -

# Now let's make a list containing the feature matrices 
source('processing_scripts/make_data_matrices.R') # helper script
tsv_dir = 'raw_data/'%+%tag%+%'/'
data = make_data_matrices(tsv_dir)

# here each feature matrix is has a row for each well and a column for each spot. We can look at the dimension of the first few feature matrices

head(lapply(data,dim))

# we'll make the list of features we want to RR transform

to_rr = grep("_PA_|_CP_",names(data),value=TRUE)
to_rr = to_rr[!grepl("_SE",to_rr)]
unum = sapply(to_rr,function(x)length(unique(unlist(data[[x]]))))
to_rr = names(unum[unum>20])
sample(to_rr,5)

# and then save the feature matrices and list of features

dl = list(data=data,rr_feats = to_rr)
saveRDS(dl,ddir%+%tag%+%'.rds')

# now we can apply the RR transformation to the matrices (this can take a while)

# # 2. RRScale

source('processing_scripts/estimate.R',chdir=TRUE)
scaledY_file = ddir%+%"scaledY_"%+%tag%+%".rds"
if(file.exists(scaledY_file)){
    cat("Reading from cache.\n")
    scaledY = readRDS(scaledY_file)
} else {
    cat("No cache, re-running.\n")
    scaledY = scale_data(dl)
    saveRDS(scaledY,file=scaledY_file)
}

# Finally, we'll print the session info for reproducibility

# # 3. Clean up spatial effects

# +
base_fmats = dl$data
rr_fmats = lapply(scaledY,"[[","RR")

ecmps = factor(unlist(base_fmats$ECMp[1,]))
ligands = factor(unlist(base_fmats$Ligand[,1]))
# -

# find which features are measured in all batches for normalization (the "dapi" features)

all_missing = sapply(rr_fmats,function(y)apply(y,1,function(x)all(!is.finite(x))))
all_batches = apply(!all_missing,2,all)
all_batches = names(all_batches[all_batches])
all_batches = all_batches[!grepl("Loess",all_batches)]
all_batches

#devtools::load_all('../r_packages/mema_norm_rpkg/memanorm')
d_none = 0
unadj_norm_out = remove_spatial(f_mats = rr_fmats,
                  adj_idx = all_batches,
                  row_rep = ligands,
                  col_rep = ecmps,
                  d=d_none,
                  verbose=FALSE)
unadj_fmats=unadj_norm_out$f_adj
names(unadj_fmats) <- names(unadj_fmats)%+%"_RR"
saveRDS(unadj_fmats,ddir%+%"unadjusted_rr_"%+%d_none%+%"_"%+%tag%+%".rds")

library('memanorm')
d_adj=15
norm_out = remove_spatial(f_mats = rr_fmats,
                  adj_idx = all_batches,
                  row_rep = ligands,
                  col_rep = ecmps,
                  d=d_adj,
                  verbose=FALSE)
adj_fmats=norm_out$f_adj
names(adj_fmats) <- names(adj_fmats)%+%"_RR_ADJ"
saveRDS(adj_fmats,ddir%+%"adjusted_rr_"%+%d_adj%+%"_"%+%tag%+%".rds")

# unadjusted features (same exact procedure but not removing spatial)

# # 4. make .csv data matrices

dmsn = dim(base_fmats[[1]])
fmat_col = array(rep(1:dmsn[2],each=dmsn[1]),dmsn)
fmat_row = t(array(rep(1:dmsn[1],each=dmsn[2]),rev(dmsn)))
fmat_dims = list("fcol"=fmat_col,"frow"=fmat_row)

library('data.table')
all_fmats = c(base_fmats,unadj_fmats,adj_fmats,fmat_dims)
all_fmats = lapply(all_fmats,data.table)

ul_all_fmats = lapply(all_fmats,unlist)
d = data.frame(ul_all_fmats)

dir.create(ddir%+%"csv/",showWarnings=FALSE,recursive=TRUE)
write.csv(x=d,file=ddir%+%"csv/"%+%tag%+%"_"%+%d_adj%+%"_df.csv")
write.csv(names(adj_fmats),file=ddir%+%"csv/"%+%tag%+%"_"%+%d_adj%+%"_features.csv")
