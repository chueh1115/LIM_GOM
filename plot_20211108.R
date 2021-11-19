flowname
compartment<-c("SED","BAC","MEI","MAC","MEG","FIS",
               "POC_w","EXP_S","DIC_W","EXP_B")
n<-length(compartment)
flow<-matrix(0,n,n)
flow
colnames(flow)=rownames(flow)=compartment
data<-data.matrix(value$mean)
data
# rownames(data)=flowname
# typeof(data)

# convert to dataframe
data <- data.frame(data)
# extract flow directiion for for loop filling
data$before <- gsub("->.*", "", flowname)
data$after <- gsub(".*->", "", flowname)


value

d_expand <- expand.grid(before = compartment, after = compartment)
d_expand[1,1:2] == data[1, 2:3]

d_expand$before[1] == 
data[1,]
d_expand$data <- NA
data$data[match(d_expand$before, data$before) & match(d_expand$after, data$after)]


data

flowname

b %in% data$before & a %in% data$after

# i as row
for(i in 1:nrow(flow)) {
  # j as col
  for(j in 1:ncol(flow)){
    b <- colnames(flow)[j]
    a <- row.names(flow)[i]
    # matching
    if(b %in% data$before && a %in% data$after){
      flow[i,j] <- subset(data, 
                          before == b & after == a)$data
    } else{
      flow[i,j] <- 0
    }  
  }
}

subset(data, before == colnames(flow)[1] & after == row.names(flow)[2])$data

subset(data, before == 1 & after == )
nrow(data)
