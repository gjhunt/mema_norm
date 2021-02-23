library('dplyr')
library('data.table')
library('tidyr')
library('readr')

make_data_matrices = function(data_dir){

    # Read in data
    files = dir(data_dir, full.names=TRUE)
    read_files = lapply(files,read.csv,sep="\t",header=TRUE,stringsAsFactors=FALSE)

    # Bind together rows of plates
    d = bind_rows(read_files)
    d$PrintRow = floor((d$PrintSpot-1)/20)+1
    d$PrintColumn = (d$PrintSpot-1)%%20+1

    #Remove punctuation from ECMP names
    d$ECMp = gsub('[[:punct:]]+','_',d$ECMp)

    meta_cols = colnames(d)
    meta_cols = meta_cols[!grepl("PA|CP|_SE|Row|Col",meta_cols)]

    for(j in meta_cols){
        d[,j] <- as.factor(d[,j])
    }
    d = data.table(d,stringsAsFactors=FALSE)

    make = function(to_rmv){
        d = d[!grepl(to_rmv,d$ECMp),]

        # Make well by spot matrices
        make_v = function(v){
            sd = d[,.(Barcode,PrintSpot,Well,v=get(v))]
            sd = sd %>% spread(PrintSpot,"v")
            rns = apply(sd[,1:2],1,function(x)paste(x,collapse="_"))
            sd = sd[,-(1:2)]
            rownames(sd) = rns
            return(sd)
        }

        v_names = colnames(d)
        d_list = lapply(v_names,make_v)
        names(d_list) = v_names
        return(d_list)
    }

    data = make("fiducial|Fiducial|gelatin|blank|air|PBS|ELN|NID")

    # impute factors with missing values
    data$Ligand <- uImpute(data$Ligand,1)
    data$Well <- uImpute(data$Well,1)
    data$StainingSet <- uImpute(data$StainingSet,1)
    data$Plate <- uImpute(data$Barcode,1)
    
    data$ECMp <- t(uImpute(data$ECMp,2))
    data$PrintRow <- t(uImpute(data$PrintRow,2))
    data$PrintColumn <- t(uImpute(data$PrintCol,2))
    data$PrintSpot <- t(uImpute(data$PrintSpot,2))

    # add some useful factors
    data[["WellLetter"]] <- data.table(apply(data[["Well"]],c(1,2),function(x)gsub('[[:digit:]]+', '', x)))
    data[["WellNumber"]] <- data.table(apply(data[["Well"]],c(1,2),parse_number))
    
    return(data)
}

uImpute=function(mtx,margin){
    mtx = t(apply(mtx,margin,function(x)rep(unique(x[!is.na(x)]),length(x))))
    return(data.table(mtx))
}
