library('ggplotify')

library('GGally')

well_scatter=function(feats,fnames,pw,facet){  
  stopifnot(length(feats)==length(fnames))
  dta = d[PlateWell%in%pw,]
  sdta = dta[,c(feats,'PrintRow','PrintColumn','PlateWell'),with=FALSE]
  colnames(sdta) = c(fnames,'PrintRow','PrintColumn','Well')
  p=ggpairs(data=sdta,columns=1:length(feats))
  return(p)
}

well_scatter_byw=function(feats,fnames,pw,facet){  
  stopifnot(length(feats)==length(fnames))
  dta = d[PlateWell%in%pw,]
  sdta = dta[,c(feats,'PrintRow','PrintColumn','PlateWell'),with=FALSE]
  colnames(sdta) = c(fnames,'PrintRow','PrintColumn','Well')
  ssdta = sdta%>% spread(Well,fnames)
  p=ggpairs(ssdta,columns=3:(length(pw)+2))
  return(p)
}


# a useful function for plotting wells individually
well_plot=function(feats,fnames,pw,facet,scale_name="value",nr=2){
  stopifnot(length(feats)==length(fnames))
  dta = d_c[PlateWell%in%pw,]
  sdta = dta[,c(feats,'PrintRow','PrintColumn','PlateWell'),with=FALSE]
  colnames(sdta) = c(fnames,'PrintRow','PrintColumn','Well')
  mdta=melt(sdta,id.vars=c('PrintRow','PrintColumn','Well'))
  colnames(mdta) = c('PrintRow','PrintColumn','Well','Feature','value')
  p=ggplot(data=mdta,mapping=aes_string(x="PrintColumn",y="PrintRow",fill="value"))
  p = p+theme_classic()
  p = p+geom_tile()+scale_fill_gradient2(name=scale_name,na.value='white')
  p = p +coord_fixed()+scale_y_reverse()
  p=p+theme(legend.position="bottom")
  p=p+facet_wrap(as.formula(facet),nrow=nr)
  p=p+theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
  p=p+theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())+
    theme(panel.background=element_rect(fill="white", colour="white"))
        
  return(p)
}

## asvd scatter plots
library('memanorm')

vdf = function(asvd,fmats,scale=TRUE){
  if(scale)
      V = data.table(t(t(asvd$v)*asvd$dv))
  else
      V = data.table(asvd$v)   
  colnames(V) = paste0("PC",1:ncol(V))
  V$PrintRow <- unlist(fmats$PrintRow[1,])
  V$PrintColumn <- unlist(fmats$PrintColumn[1,])
  V$ECMp <- unlist(fmats$ECMp[1,])
  return(list(v=V,dv=asvd$dv))
}

library('scales')
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# asvd pairwise scatter plots
asvd_scatter = function(asvd,fmats,k=3,hl=NULL,ttl=NULL,legnd=TRUE){
  vdf.out = vdf(asvd,fmats)
  V = vdf.out$v
  pcts = vdf.out$dv^2
  pcts = pcts/sum(pcts)*100
  if(!is.null(hl)){
    V$ECMp[V$ECMp!=hl] <- "Other"
    V$ECMp <- ordered(V$ECMp,levels=c("Other",hl))
  }
  sfn = function(k){
    k1=k[1];k2=k[2]
    xl="PC"%+%k1%+%" ("%+%round(pcts[k1],3)%+%"%)" 
    yl="PC"%+%k2%+%" ("%+%round(pcts[k2],3)%+%"%)" 
      plt=ggplot(data=V,mapping=aes_string(x="PC"%+%k1,y="PC"%+%k2,color="ECMp",shape="ECMp"))+geom_point()+theme_classic()+scale_color_manual(values=gg_color_hue(length(unique(V$ECMp))))+scale_shape_manual(values=c(16,17,1:14,48:57,65:90))+ggtitle(ttl)+labs(x=xl,y=yl)#+coord_fixed()
    if(!is.null(hl))
      plt = plt + geom_point(data=V[V$ECMp==hl],size=3)
    if(!legnd)
      plt = plt + guides(color=FALSE,shape=FALSE)
    return(plt)
  }
  pl = apply(combn(k,2),2,sfn)
  return(pl)
}

asvd_scatter_layout = function(asvd,fmats,k=3,hl=NULL,ttl=NULL,legnd=TRUE){
  vdf.out = vdf(asvd,fmats)
  V = vdf.out$v
  pcts = vdf.out$dv^2
  pcts = pcts/sum(pcts)*100
  sfn = function(k){
    k1=k[1];k2=k[2]
    xl="PC"%+%k1%+%" ("%+%round(pcts[k1],1)%+%"%)" 
    yl="PC"%+%k2%+%" ("%+%round(pcts[k2],1)%+%"%)"
    plt=ggplot(data=V,mapping=aes_string(x="PC"%+%k1,y="PC"%+%k2,color=hl))+geom_point()+theme_classic()+ggtitle(ttl)+labs(x=xl,y=yl)
      #+coord_cartesian(xlim=c(-.15,.15),ylim=c(-.15,.15))
    return(plt)
  }
    pl = apply(combn(k,2),2,sfn)
  return(pl)
}

asvd_compare=function(asvd,asvd_adj,fmats,rw_leg = .15, hl=NULL,fn = asvd_scatter){
  unadj = lapply(fn(asvd,fmats,hl=hl,k=2,ttl="Unadjusted"),function(x)x+guides(shape=FALSE,color=FALSE))
  adj = lapply(fn(asvd_adj,fmats,hl=hl,k=2,ttl="   Adjusted"),function(x)x+guides(shape=FALSE,color=FALSE))
  leg = as.ggplot(get_legend(fn(asvd,fmats,hl=hl,k=2,ttl="Unadjusted")[[1]]))
  plot_grid(plotlist=list(plot_grid(plotlist=c(unadj,adj),nrow=2),leg),nrow=1,rel_widths=c(1,rw_leg))
               #,labels=c('(i)','(ii)'),nrow=1,rel_widths=c(1,1,.35),hjust=-1.5)
}
## asvd layout
asvd_layout = function(asvd,fmats,k=5,ttl=NULL,nrow=NULL,scale=TRUE){
  V = vdf(asvd,fmats,scale=scale)$v
  colnames(V) <- gsub("PC","S",colnames(V))
  mV = melt(V,id.vars=c('PrintRow','PrintColumn','ECMp'))
  p=ggplot(data=mV[variable%in%paste0("S",1:k),],mapping=aes(x=PrintColumn,y=PrintRow,fill=value))+geom_tile()+scale_fill_gradient2()+facet_wrap(~variable,nrow=nrow)+coord_fixed()
  p=p+theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
  p=p+theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())
  p = p +theme_classic()+ggtitle(ttl)
  return(p)
}

bcss = function(SVD,fmats){
    V = data.table(vdf(SVD,fmats)$v[,c('PC1','PC2','PrintRow','PrintColumn','ECMp')])
    W = V[,.(WSS = (PC1-mean(PC1))^2+(PC2-mean(PC2))^2),by=ECMp]
    W = W[,.(WSS=sum(WSS)),by=.(ECMp)]
    T = V[,.(TSS = (PC1-mean(PC1))^2+(PC2-mean(PC2))^2)]
    T = T[,.(TSS=sum(TSS))]
    T = unlist(T)
    WW = W[,.(WSS=sum(WSS))]
    B = unlist(T) - unlist(WW)
    W$WT = W$WSS/T
    return(W)
}

bcss_compare = function(SVD,SVD_adj){
    b_un = bcss(SVD)
    b_adj = bcss(SVD_adj)
    bc = merge(b_un,b_adj,by='ECMp')
    colnames(bc) <- c('ECMp','WSS','Unadjusted','WSS_adj','Adjusted')
    bc$d = bc$Adjusted-bc$Unadjusted
    mbc = melt(bc[,.(ECMp,Unadjusted,Adjusted)],id.vars='ECMp')
    mbc$ECMp <- factor(mbc$ECMp,levels=bc$ECMp[order(bc$d,decreasing=FALSE)])
    bt = bc[,.(WT=sum(Unadjusted),WT_adj=sum(Adjusted))]
    return(list(bc=bc,mbc=mbc,bt=bt))
}
