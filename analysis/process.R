# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,Rmd,R
#     text_representation:
#       extension: .R
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.4
#   kernelspec:
#     display_name: R
#     language: R
#     name: ir
# ---

# # 1. Pre-reqs
# We may need to install the rrscale and memnorm packages

# +
#install.packages('../r_packages/rrscale_rpkg/rrscale_1.0.3.tar.gz', repos=NULL,type="source")
#install.packages('../r_packages/mema_norm_rpkg/memanorm_0.0.0.9005.tar.gz', repos=NULL,type="source")
# -

library('rrscale')
devtools::load_all('../r_packages/mema_norm_rpkg/memanorm/')
source('processing_scripts/util.R')

# # and then we'll set up a directory for the processed data
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

# we'll make the list of features we want to RR transform

to_rr = grep("_PA_|_CP_",names(data),value=TRUE)
to_rr = to_rr[!grepl("_SE",to_rr)]
unum = sapply(to_rr,function(x)length(unique(unlist(data[[x]]))))
to_rr = names(unum[unum>20])
sample(to_rr,5)

# and then save the feature matrices and list of features

dl = list(data=data,rr_feats = to_rr)

fn = ddir%+%tag%+%'.rds'
fn

if(!file.exists(fn))
    saveRDS(dl,fn)

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

# # 3. Clean up spatial effects

# +
base_fmats = dl$data
rr_fmats_unproc = lapply(scaledY,"[[","RR")

ecmps = factor(unlist(base_fmats$ECMp[1,]))
ligands = factor(unlist(base_fmats$Ligand[,1]))
# -

# find which features are measured in all batches for normalization (the "dapi" features)

rr_fmats_im = lapplyp(rr_fmats_unproc, memanorm::row_impute_fn, ecmps, verbose = TRUE, .name = "Imputing")

rr_fmats = rr_fmats_im

all_batches = names(rr_fmats)

d_none = 0
unadj_norm_out = remove_spatial(f_mats = rr_fmats,
                  adj_idx = all_batches,
                  row_rep = ligands,
                  col_rep = ecmps,
                  d=d_none,
                  verbose=FALSE,
                  row_center=FALSE,
                  row_impute=FALSE)
unadj_fmats=unadj_norm_out$f_adj
names(unadj_fmats) <- names(unadj_fmats)%+%"_RR"

unadj_fmats0 = rr_fmats_unproc
names(unadj_fmats0) <- names(unadj_fmats0)%+%"_RR"

adj_fmat_l = list()

d_adj_list = c(0,1,2,3,5,7,10,15,25)
for(d_adj in d_adj_list){
    cat(paste0("Adjusting for spatial effects d=",d_adj,"\n"))
    flush.console()
    norm_out = remove_spatial(f_mats = rr_fmats,
                      adj_idx = all_batches,
                      row_rep = ligands,
                      col_rep = ecmps,
                      d=d_adj,
                      verbose=FALSE,
                      full.results=TRUE,
                      row_center=FALSE,
                      row_impute=FALSE)
    adj_fmats=norm_out$f_adj
    names(adj_fmats) <- names(adj_fmats)%+%"_RR_ADJ"
    adj_fmats_spatial = adj_fmats
    
    adj_fmat_l[[length(adj_fmat_l)+1]] = adj_fmats
}


# # plot $\hat{S}$

library('ggplot2')

norm_out = remove_spatial(f_mats = rr_fmats,
                      adj_idx = all_batches,
                      row_rep = ligands,
                      col_rep = ecmps,
                      d=15,
                      verbose=FALSE,
                      row_impute=FALSE,
                      full.results=TRUE)

k=15

sh = norm_out$s_hat[1:k,,drop=FALSE]
sh = data.frame(t(sh))
colnames(sh) = paste0("S",1:k)
dim(sh)
rw = base_fmats[['PrintRow']]
cl = base_fmats[['PrintColumn']]
sh$rw = rw[1,]
sh$cl = cl[1,]
msh = melt(sh,id.vars=c('rw','cl'))
options(repr.plot.width = 15, repr.plot.height = 8)
plt = ggplot(data=msh,mapping=aes(x=cl,y=rw,fill=value))+
    facet_wrap(~variable,nrow=3,scale='free')+
    geom_tile()+
    theme_classic()+
    scale_fill_gradient2()+
    theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())+
            theme(axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank())
plt
ggsave(plot=plt,file="S.pdf",width=5,height=5)

# # 4. make .csv data matrices

library('data.table')

to_csv = function(adj_fmats,d_adj,tg=''){
    cat(d_adj,"\n")
    flush.console()
    dmsn = dim(base_fmats[[1]])
    fmat_col = array(rep(1:dmsn[2],each=dmsn[1]),dmsn)
    fmat_row = t(array(rep(1:dmsn[1],each=dmsn[2]),rev(dmsn)))
    fmat_dims = list("fcol"=fmat_col,"frow"=fmat_row)
    
    all_fmats = c(base_fmats,unadj_fmats,adj_fmats,fmat_dims)
    all_fmats = lapply(all_fmats,data.table)
    
    ul_all_fmats = lapply(all_fmats,unlist)
    d = data.frame(ul_all_fmats)
    
    dir.create(ddir%+%"csv/",showWarnings=FALSE,recursive=TRUE)
    write.csv(x=d,file=ddir%+%"csv/"%+%tag%+%"_"%+%d_adj%+%"_df"%+%tg%+%".csv")
    write.csv(names(adj_fmats),file=ddir%+%"csv/"%+%tag%+%"_"%+%d_adj%+%"_features"%+%tg%+%".csv")
    return(d)
}

stopifnot(length(d_adj_list)==length(adj_fmat_l))

toss = lapply(1:length(d_adj_list),function(i)to_csv(adj_fmat_l[[i]],d_adj_list[[i]],tg=''))
