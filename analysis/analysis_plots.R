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

# # Data Pre-processing
#
# First we'll load some packages and utility scripts

source('./processing_scripts/util.R')
tag="MCF10A"
library('ggplot2')
library('data.table')
library('cowplot')
options(warn=-1) # turn off warnings

# and then read in the data

tg = ""

dm=15
ddir = './processed_data/'%+%tag%+%"/csv/"
dir(ddir)
d = read.csv(ddir%+%tag%+%"_"%+%dm%+%"_df"%+%tg%+%".csv",row.names=1)
d = data.table(d)
d$PlateWell <- d$Plate%+%"."%+%d$Well

# We also add in the feature "PlateWell" which is the interaction of Plate and Well. The variable `dm` is the number of dimensions of unwanted variation to remove. For this example we set it to 15. In this representation of the data there are 128064 rows (192 wells x 667 spots) and 486 columns (features)

# We can also plot the design layout

seld = d[Barcode == unique(Barcode)[1],]

plt = ggplot(data=seld,mapping=aes(y=ArrayRow,x=ArrayColumn,color=ECMp))+
    geom_point(size=5)
plt = plt + theme_bw()
plt = plt + guides(color=FALSE) + coord_fixed()
plt = plt + labs(x=NULL,y=NULL)
plt = plt + theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
plt = plt + theme(plot.margin=unit(c(0,0,0,0),"mm")) 
plt
ggsave(plot=plt,file="design.png",width=4,height=7.1)

# The first several features are meta-data concerning the spots
head(colnames(d),n=20)

# and for each of the measured features there are three versions:
# 1. the raw version, 
# 2. the RR transformed version,
# 3. the adjusted RR version

colnames(d)[c(99,282,398)] # 3 versions of Nuclei Area

numeric_cols = colnames(d)[71:ncol(d)-3][-(180:184)]
numeric_cols = numeric_cols[which(sapply(numeric_cols,function(s)grepl("_RR$|_RR_ADJ$",s)))]
texture_cols = numeric_cols[grepl("Texture",numeric_cols)]
non_texture = setdiff(numeric_cols,texture_cols)                                         

d_c = as.data.frame(d)
d_c = data.table(d_c)
d_c = d_c[ , (numeric_cols) := lapply(.SD, function(x)x-mean(x,na.rm=TRUE)), .SDcols = numeric_cols,by=.(PlateWell)]

# In addition to the format of the .csv file it will be useful to have a list of feature matrices. For each matrix a row is a well and each column a spot.

source('./processing_scripts/make_fmat.R')
fmats = lapply(names(d),function(v)make_fmat(d,v))
names(fmats) <- names(d)

# For this data we have three staining batches we'll call 1, 2 and 3. We want to take note of which features are measured in which batch, here "1", "2", or "3" means the feature is measured in those batches alone, "13" means the feature is measured in batch 1 and 3, and "123" means it is measured in all batches. Those features measured in all batches were used to adjust the other features. 

staining_set = substr(apply(fmats$StainingSet,1,unique),3,3)
missing_rows = sapply(fmats,function(X)apply(X,1,function(x)all(is.na(x))))
feature_batch = apply(missing_rows,2,function(x)paste0(sort(unique(staining_set[!x])),collapse=""))
head(sample(feature_batch))

# For this dataset features are measured in either batches 1, 2, or 3 alone, batches 1 and 3, or all three batches. Note here a missing "feature batch" indicates meta-data.

# # Plots in Paper

# ## Fig 2

# +
cmb = cbind(fmats[["Plate"]][,1],fmats[["Well"]][,1])
pws = apply(cmb,1,function(x)paste0(x,collapse="."))

source('plot_scripts/plot_fns.R')

library('viridis')

warn = getOption("warn")
options(warn=-1)


sel_feats = c('Cytoplasm_CP_Intensity_MedianIntensity_Actin_RR',
                         'Nuclei_CP_AreaShape_Area_RR',
                         'Cytoplasm_CP_Intensity_MedianIntensity_MitoTracker_RR',
                        'Nuclei_PA_Cycle_DNA4NProportion_RR')
short_names = c("Cyto. Med. Actin",
               "Nuc. Area",
               "Cyto. Med. MT",
               "Nuc. DNA4 Prop."
              )
sel_pw = "LI8X00525.A01"

wp1=well_plot(feats=sel_feats,
    fnames = short_names,
    pw=sel_pw,
    facet="~Feature",nr=1)+scale_fill_gradient2(limits=c(-4,4),na.value='white')

sel_feats = paste0(sel_feats,"_ADJ")
wpa1=well_plot(feats=sel_feats,
    fnames = short_names,
    pw=sel_pw,
    facet="~Feature",nr=1)+scale_fill_gradient2(limits=c(-4,4),na.value='white')


ft = 'Nuclei_PA_Cycle_DNA2NProportion_RR'
ft = 'Cytoplasm_CP_Intensity_MedianIntensity_MitoTracker_RR'
scc = fmats[[ft]] #fmat

# find wells that are highly cor with other wells generally
# we will search among these
cor_mat = cor(t(scc),use="pairwise.complete.obs")^2 # squared cors btwn wells

sel_feats=ft
short_names = 'Cyto. Med. MT'
sel_pw = c("LI8X00523.B02","LI8X00524.B02","LI8X00525.B01","LI8X00525.B02")#pws[order(cor_mat[which(pws=='LI8X00523.B02'),],decreasing=TRUE)[1:4]]

wpw3=well_plot(feats=sel_feats,
fnames = short_names,
pw=sel_pw,
facet="~Well",
scale_name=short_names,nr=1)+scale_fill_gradient2(limits=c(-4,4),na.value='white')

sel_feats = paste0(sel_feats,"_ADJ")
wpwa3=well_plot(feats=sel_feats,
fnames = short_names,
pw=sel_pw,
facet="~Well",
scale_name=short_names,nr=1)+scale_fill_gradient2(limits=c(-4,4),na.value='white')

lbls = c('(A) Spatial Sig. Across Features',
         '(B) Spatial Sig. Across Wells/Plates',
         '(A) Spatial Sig. Across Features (Adj.)',
         '(B) Spatial Sig. Across Wells/Plates (Adj.)'
        )

out = plot_grid(wp1+theme(legend.position='bottom'),
                wpw3+theme(legend.position='bottom'),
                labels=lbls[1:2],nrow=2,
                hjust=0,vjust=1.1,
                scale=.925
               )
ggsave(plot=out,
       file=paste0('./',tg,'patterns.pdf'),width=6,height=8)
# -

out

# ## Fig 4

# +
out2 = plot_grid(wpa1+theme(legend.position='bottom'),
                wpwa3+theme(legend.position='bottom'),
                labels=lbls[3:4],nrow=2,
                hjust=0,vjust=1.1,
                scale=.925)

ggsave(plot=out2,
       file=paste0('./',tg,'patterns_clean.pdf'),width=6,height=8)
# -

out2

# # Supplementary figure 2

ft = 'Nuclei_CP_AreaShape_Compactness_RR'
F = (scale((fmats[[ft]]),center=FALSE,scale=FALSE))
FA = (scale((fmats[[ft%+%"_ADJ"]]),center=FALSE,scale=FALSE))
svd1 = memanorm::average_svd_na(list(F))
svd1_adj = memanorm::average_svd_na(list(FA))

udf = function(svd,hl=NULL){
    udf = svd$u[,1:3]
    colnames(udf) = paste0("PC",1:ncol(udf))
    udf = data.frame(udf)
    udf$sset = unlist(fmats[['StainingSet']][,1])
    udf$bc = unlist(fmats[['Barcode']][,1])
    udf$wl = unlist(fmats[['WellLetter']][,1])
    udf$lig = unlist(fmats[['Ligand']][,1])
    #if(!is.null(hl)){
    #    udf$lig[!(udf$lig%in%hl)] <- "Other"
    #    udf$lig <- factor(udf$lig,levels=c("Other",hl))
    #}
    return(udf)
}

udf1 = udf(svd1)
udf1_adj = udf(svd1_adj)

p1 = ggplot(data=udf1,mapping=aes(x=PC1,y=PC2,color=sset))+geom_point(size=5)+
    scale_shape_manual(values=c(16,17,1:14,48:57,65:95))+
    labs(color='Staining Batch')+
    ggtitle("Unadjusted")+theme_bw()

p2 = ggplot(data=udf1,mapping=aes(x=PC1,y=PC2,color=wl,shape=wl))+geom_point(size=5)+
    scale_shape_manual(values=c(16,17,1:14,48:57,65:95))+
    labs(color='Well Position',shape='Well Position')+
    ggtitle("Unadjusted")+theme_bw()

p3 = ggplot(data=udf1_adj,mapping=aes(x=PC1,y=PC2,color=sset))+geom_point(size=5)+
    labs(color='Staining Batch')+
    ggtitle("Adjusted")+theme_bw()

p4 = ggplot(data=udf1_adj,mapping=aes(x=PC1,y=PC2,color=wl,shape=wl))+geom_point(size=5)+
    scale_shape_manual(values=c(16,17,1:14,48:57,65:95))+
    ggtitle("Adjusted")+theme_bw()+
    labs(color='Well Position',shape='Well Position')

leg1 = get_legend(p1)
leg2 = get_legend(p2)

p1 = p1 + theme_classic()+guides(color='none',shape='none')
p2 = p2 + theme_classic()+guides(color='none',shape='none')
p3 = p3 + theme_classic()+guides(color='none',shape='none')
p4 = p4 + theme_classic()+guides(color='none',shape='none')

options(repr.plot.width = 12, repr.plot.height = 6)
pg = plot_grid(p1+guides(color='none'),p3,rel_widths=c(1,1),labels=c('(A)','(B)'))

pg2 = plot_grid(p2+guides(color='none'),p4,rel_widths=c(1,1),labels=c('(C)','(D)'))

plt = plot_grid(pg,leg1,pg2,leg2,rel_widths=c(1,.25,1,.25),nrow=2)
plt

ggsave(plot=plt,file=paste0(tg,"batch_svd.pdf"),width=12,height=7)

# # Figure 5

ft = "Nuclei_PA_Cycle_DNA2NProportion_RR"

F = t(scale(t(fmats[[ft]]),center=TRUE,scale=FALSE))
FA = t(scale(t(fmats[[ft%+%"_ADJ"]]),center=TRUE,scale=FALSE))

svd1 = memanorm::average_svd_na(list(F))
svd1_adj = memanorm::average_svd_na(list(FA))

source('plot_scripts/plot_fns.R')
hl="THBS1_1"
options(repr.plot.width = 12, repr.plot.height = 5)
sctr_thbs1=asvd_compare(svd1,svd1_adj,hl=hl,fmats=fmats)

sctr_row=asvd_compare(svd1,svd1_adj,fn=asvd_scatter_layout,hl="PrintRow",fmats=fmats)

ecmp_sctr = asvd_scatter(svd1,fmats=fmats,k=2)[[1]]
ecmp_sctr_adj = asvd_scatter(svd1_adj,fmats=fmats,k=2)[[1]]

scatter1 = plot_grid(sctr_thbs1,sctr_row,labels=c('(A)','(B)'))
scatter1

ggsave(plot=scatter1, filename=paste0(tg,"compact_scatter.pdf"),width=10,height=5)

# ## Figure 6

ss_cca = function(ft,K=1){
    cat(ft,"\n")
    flush.console()
    F = t(scale(t(fmats[[ft]]),center=TRUE,scale=FALSE))
    FA = t(scale(t(fmats[[ft%+%"_ADJ"]]),center=TRUE,scale=FALSE))
    svd1 = memanorm::average_svd_na(list(F))
    svd1_adj = memanorm::average_svd_na(list(FA))
    B = model.matrix(~-1+unlist(fmats$ECMp[1,]))
    cca1 = mean(cancor(svd1$v[,1:K],B)$cor^2)
    cca1_adj = mean(cancor(svd1_adj$v[,1:K],B)$cor^2)
    return(list(cca1=cca1,cca1_adj=cca1_adj,diff=cca1_adj-cca1))
}

fts=sapply(numeric_cols,function(s)gsub("_ADJ","",s))
ss_ccas = lapply(fts,ss_cca)

# +
options(repr.plot.width = 12, repr.plot.height = 6)
kp_ccas = ss_ccas[!duplicated(fts)]
cca_df = data.frame(diff = sapply(kp_ccas,"[[","diff"),
                   unadj=sapply(kp_ccas,"[[","cca1"),
                   adj=sapply(kp_ccas,"[[","cca1_adj"))
cca_df$feature <- names(kp_ccas)
cca_df$feature = gsub("_RR|_CP|_PA","",cca_df$feature)
cca_df$feature = gsub("_AreaShape|_Texture|_Gated|IntegratedIntensity_","",cca_df$feature)
cca_df$feature = gsub("Intensity_","",cca_df$feature)
cca_df$feature = gsub("Cytoplasm","Cyto",cca_df$feature)
cca_df$feature = gsub("Nuclei","Nuc",cca_df$feature)
cca_df$feature = gsub("Positive","Pos",cca_df$feature)
cca_df$feature = gsub("Proportion","Prop",cca_df$feature)
cca_df$feature = gsub("Negative","Neg",cca_df$feature)
cca_df$feature = gsub("Luminal","Lum",cca_df$feature)
cca_df$feature = gsub("Angular","Ang",cca_df$feature)
cca_df$feature = gsub("Second","Sec",cca_df$feature)
cca_df$feature = gsub("Basal","Bas",cca_df$feature)
cca_df$feature = gsub("Difference","Diff",cca_df$feature)
cca_df$feature = gsub("Moment","Mom",cca_df$feature)
cca_df$feature = gsub("Inverse","Inv",cca_df$feature)
cca_df$feature = gsub("Fibrillarin","Fib",cca_df$feature)
cca_df$feature = gsub("Variance","Var",cca_df$feature)
cca_df$feature = gsub("Correlation","Cor",cca_df$feature)
cca_df$feature = gsub("Entropy","Ent",cca_df$feature)
cca_df$feature = gsub("Median","Med",cca_df$feature)
cca_df$feature = gsub("Average","Avg",cca_df$feature)
cca_df$feature <- factor(cca_df$feature,levels=cca_df$feature[order(cca_df$diff,decreasing=FALSE)])

mcca = melt(cca_df)
mcca = mcca[mcca$variable%in%c('adj','unadj'),]
colnames(mcca)[2] <- "Adjustment"
levels(mcca$Adjustment) <- c('','Unadjusted','Adjusted')
plt = ggplot(data=mcca,mapping=aes(y=value,x=feature,color=Adjustment))+geom_point()
plt = plt + theme_bw()
plt = plt+theme(plot.margin = unit(c(.1,.1,.1,.1), "cm"),
               axis.text.x = element_text(angle = 90, vjust = 1, hjust=1,size=8))+
    labs(y="Squared Cannonical Correlation\n between ECMp and PC1",x="Feature")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
plt = plt + scale_colour_manual(values=cbPalette)+theme(legend.position='bottom')
ggsave(plot=plt,file="cca.pdf",width=10,height=5)
# -

plt

ggsave(plot=plt, filename=paste0(tg,"cca.pdf"),width=10,height=5)

# ## 
