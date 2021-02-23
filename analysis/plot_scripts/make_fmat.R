library('tidyr')

make_fmat = function(d,v){
    sd = d[,.(Plate,PrintSpot,Well,v=get(v))]
    sd = sd %>% spread(PrintSpot,"v")
    rns = apply(sd[,1:2],1,function(x)paste(x,collapse="_"))
    sd = sd[,-(1:2)]
    rownames(sd) = rns
    return(sd)
}


make_df = function(L){
    ul_all_fmats = lapply(all_fmats,unlist)
    df = data.frame(ul_all_fmats)
    return(df)
}
