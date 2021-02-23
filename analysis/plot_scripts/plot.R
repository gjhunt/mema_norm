## read in data
source('../processing_scripts/util.R')
tag="MCF10A"
library('ggplot2')
library('data.table')

# read in csv df
ddir = '../data/'%+%tag%+%"/"
d = read.csv(ddir%+%tag%+%"_df.csv",row.names=1)
d = data.table(d)

# produce list of feature matrices
source('make_fmat.R')
fmats = lapply(names(d),function(v)make_fmat(d,v))
names(fmats) <- names(d)

# determine which features are measured in which batches
staining_set = substr(apply(fmats$StainingSet,1,unique),3,3)
missing_rows = sapply(fmats,function(X)apply(X,1,function(x)all(is.na(x))))
feature_batch = apply(missing_rows,2,function(x)paste0(sort(unique(staining_set[!x])),collapse=""))

print(head(sample(feature_batch)))
print(unique(feature_batch))

## asvd layout plots
library('memanorm')

vdf = function(asvd){
    V = data.table(asvd$v)
    colnames(V) = paste0("PC",1:ncol(V))
    V$PrintRow <- unlist(fmats$PrintRow[1,])
    V$PrintColumn <- unlist(fmats$PrintColumn[1,])
    V$ECMp <- unlist(fmats$ECMp[1,])
    return(V)
}

asvd_layout = function(asvd,k=5){
    V = vdf(asvd)
    mV = melt(V,id.vars=c('PrintRow','PrintColumn','ECMp'))
    ggplot(data=mV[variable%in%paste0("PC",1:k),],mapping=aes(x=PrintColumn,y=PrintRow,fill=value))+geom_tile()+scale_fill_gradient2()+facet_wrap(~variable)+coord_fixed()
}


# asvd pairwise scatter plots
asvd_scatter = function(asvd,k=3,hl=NULL){
    V = vdf(asvd)
    if(!is.null(hl))
        V$ECMp[V$ECMp!=hl] <- "Other"
    sfn = function(k){
        k1=k[1];k2=k[2]
        ggplot(data=V,mapping=aes_string(x="PC"%+%k1,y="PC"%+%k2,color="ECMp"))+geom_point()
    }
    pl = apply(combn(k,2),2,sfn)
    return(pl)
}

## plots
library('cowplot')
ss1_rr = grep("_RR$",names(which(feature_batch=="1")),value=TRUE)
ss1_rr_adj = grep("_RR_ADJ$",names(which(feature_batch=="1")),value=TRUE)
asvd1_rr = memanorm::average_svd_na(fmats[ss1_rr])
asvd1_rr_adj = memanorm::average_svd_na(fmats[ss1_rr_adj])

plot_grid(asvd_layout(asvd1_rr,k=10),
          asvd_layout(asvd1_rr_adj,k=10))


hl="THBS1_1"
plot_grid(plot_grid(plotlist=asvd_scatter(asvd1_rr,hl=hl)),
          plot_grid(plotlist=asvd_scatter(asvd1_rr_adj,hl=hl)))

# plots for group 123
ss123_rr = grep("_RR$",names(which(feature_batch=="123")),value=TRUE)
ss123_rr_adj = grep("_RR_ADJ$",names(which(feature_batch=="123")),value=TRUE)
asvd123_rr = memanorm::average_svd_na(fmats[ss123_rr])
asvd123_rr_adj = memanorm::average_svd_na(fmats[ss123_rr_adj])

plot_grid(asvd_layout(asvd123_rr,k=10),
          asvd_layout(asvd123_rr_adj,k=10))

hl="THBS1_1"
plot_grid(plot_grid(plotlist=asvd_scatter(asvd123_rr,hl=hl)),
          plot_grid(plotlist=asvd_scatter(asvd123_rr_adj,hl=hl)))

###

View(sapply(unique(staining_set),function(ss)which(staining_set==ss)))
batch = sapply(apply(missing_rows,2,which),paste0,collapse="_")


#feat = "Nuclei_CP_AreaShape_Extent_RR_ADJ"
feat = 'Cytoplasm_CP_Intensity_MedianIntensity_Actin_RR_ADJ'
fr=113
dta = d[frow==fr,]
p=ggplot(data=dta,mapping=aes_string(x="PrintColumn",y="PrintRow",fill=feat))
p = p+geom_tile()+scale_fill_gradient2()+facet_grid(WellLetter~WellNumber)
p=p+theme(legend.position="bottom")+coord_fixed()+scale_y_reverse()
p

library('ggplot2')
plate_plot = function(feat, plte){
    p = ggplot(data=d[Plate==plte,],mapping=aes_string(x="PrintColumn",y="PrintRow",fill=feat))
    p = p+geom_tile()+scale_fill_gradient2()+facet_grid(WellLetter~WellNumber)
    p=p+ggtitle(plte)+theme(legend.position="bottom")+coord_fixed()
    return(p)
}

ft = "Nuclei_CP_AreaShape_Area_RR"
plts = lapply(unique(d$Plate),function(plte)plate_plot(ft,plte))

library('cowplot')
plot_grid(plotlist=plts,nrow=3)

pp(plts,fn="test.pdf")

(as.integer(factor(d$Plate))-1)*8 + as.integer(factor(d$WellIndex

ggplot(data=d[Plate==plte,],mapping=aes(x=PrintRow,y=PrintColumn,fill=Nuclei_CP_AreaShape_Area_RR))+geom_tile()+facet_grid(WellLetter~WellNumber)+ggtitle(plte)+scale_fill_gradient2()

