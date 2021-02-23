```{r}
feat = 'Cytoplasm_CP_Intensity_MedianIntensity_Actin_RR'
fts=c(feat,feat%+%"_ADJ")
dta = d[,c(fts,"ECMp","PlateWell"),with=FALSE]
colnames(dta)[1:2] = c('Unadjusted','Adjusted')
dta = dta[,.(Unadjusted=Unadjusted-mean(Unadjusted),Adjusted=Adjusted-mean(Adjusted),ECMp=ECMp),by=.(PlateWell)]
mdta = melt(dta,measure.vars=c('Unadjusted','Adjusted'))
colnames(mdta) = c('PlateWell','ECMp','Method','value')
vs = mdta[,.(var=var(value,na.rm=TRUE)),by=.(Method,ECMp)]
medians = mdta[,.(mdn=median(value,na.rm=TRUE)),by=.(Method,ECMp)]
ordering = medians[,.(mdn=abs(median(mdn))),by=ECMp]
ecmp_selection = ordering[order(abs(ordering$mdn),decreasing=TRUE),]$ECMp[1:12]
#set.seed(4321)
#ecmp_selection = sample(unique(dta$ECMp))[1:12]
ggplot(data=mdta[ECMp%in%ecmp_selection,],mapping=aes(x=value,fill=Method))+geom_density(alpha=.5)+facet_wrap(~ECMp,scales='free')+theme(legend.position="bottom")+labs(x="Actin Intensity")+theme_classic()+coord_cartesian(xlim=c(-2,2))
```