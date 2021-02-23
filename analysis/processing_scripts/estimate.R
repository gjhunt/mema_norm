source('util.R')

scale_data = function(data){

    Ys = data$data[data$rr_feats]
    
    rrscale_it <- function (ii) {
        nm=data$rr_feats[ii]
        cat(nm,"\n")
        cat(paste0(ii,"/",length(dl$rr_feats),"\n"))
        T0 = Sys.time()
        cat("Start: " %+% T0 %+%"\n")
        flush.console()
        Y = Ys[[nm]]
        if(all(Y>=0,na.rm=TRUE)){
            tl = list(box_cox_negative=box_cox_negative)
            ll=list(box_cox_negative = c(-100,100))
        } else {
            tl = list(asinh=asinh)
            ll = list(asinh=list(0,100))
        }
        ctrl = list("algorithm"="NLOPT_GN_CRS2_LM",
                                'maxeval'=500,
                                'xtol_rel'=1E-5,
                                'print_level'=1)
        ret=tryCatch({
            rrscale::rrscale(Y=Y,
                             trans_list=tl,lims_list=ll,
                             verbose=1,zeros=1E-10,log_dir='.rrscale/'%+%nm%+%'/',
                             opt_method='nloptr',opt_control=ctrl)
        },error=function(e){
            cat(e)
            cat("Failed....")
            return(NULL)
        })
        T1 = Sys.time()
        cat("Stop: " %+% T1 %+%"\n")        
        print(T1-T0)
        flush.console()
        return(ret)
    }

    cat("Scaling.\n")
    scaledY = lapply(1:length(data$rr_feats),rrscale_it)
    names(scaledY) = data$rr_feats
    names(scaledY) = names(Ys)
    
    return(scaledY)
}
