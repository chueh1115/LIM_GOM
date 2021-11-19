rm(list=ls())

setwd("C:/Users/user/Downloads/labWei/Tung_thesis/LIM_gulf of mexico/")

#-- Load the LIM (and limSolve) package 
library(LIM)
library(splus2R)

##   1)   MODEL SETUP
#-- Define directory that contains the input file
DataDir <- "C:/Users/user/Downloads/labWei/Tung_thesis/LIM_gulf of mexico/"
#-- Read the ascii files
File<- paste(DataDir,"gom_S42.input",sep="")  
LIM<- Setup(file=File) 
##   2)  Parsimonious method
# Find the solution range of each flow
flowSol <- Xranges(LIM)
# Find the parsimonious solution of each flow
pars <- Lsei(LIM, parsimonious = TRUE)
# Print your ranges and parsimonious solution on screen
SSA<-data.frame(flowSol, parsimonious=pars$X)
SSA
# plot parsimonious result
plotweb(Flowmatrix(LIM),main="gom_MT1",
        sub="(mgC/m2/d)",val=F, val.size = 0.6,
        lab.size=0.8)


##   3)  Likelihood method
#-- SSA is not based on ecological theory
xranges<-Xranges(LIM)
x0 <- lsei(E=LIM$A,
           F=LIM$B,
           A=diag(LIM$NUnknowns),
           B=rowMeans(xranges),
           G=LIM$G,
           H=LIM$H)$X
iter<-100
jumpsize<-10000
xs <- xsample(E    = LIM$A,
              F    = LIM$B,
              G    = LIM$G,
              H    = LIM$H,
              jmp  = (xranges[,2] - xranges[,1])/jumpsize,
              x0   = x0,
              iter = iter)
nameoutput <- "gom_S42_1000.Rdata" #100=iteration
save(xs, LIM, file=nameoutput)

#check: 
#1) if the number of iterations is high enough to produces convergence of mean and sd
#2) if the range of sampled values cover the range of possible solutions.


step <- 3 #step size you want to sample on iteration
#if iteration <1000 use 1; if >1000, use iteration/1000 
flownr  <- 9 #which flow you want to examine
#advice: examine a couple of flows(3 or so) of different magnitudes
load(file=nameoutput)
nsample <- iter/step
means   <- vector("numeric", nsample)
sdevs   <- vector("numeric", nsample)

for (i in 1:nsample){
  random   <- sort(runif(iter), index.return=TRUE)
  sdevs[i] <-   sd(xs$X[random$ix[1:(i*step)], flownr])
  means[i] <- mean(xs$X[random$ix[1:(i*step)], flownr])
}
plot(means)
meanvalues <- cbind(LIM$Unknowns, colMeans(xs$X))
standarddev <- cbind(LIM$Unknowns, sqrt(diag(var(xs$X))))
# need a measure of when mean and sd not significantly fluctuate anymore
# calculate fluctuation reduction
# Error margin of 2% of average stand. dev.
stdevofflow <- as.numeric(standarddev[flownr, 2])
standarddev[flownr,2]
errormargin <- stdevofflow * 0.02
plot(sdevs)
abline(h = stdevofflow + errormargin/2)
abline(h = stdevofflow - errormargin/2)

LA<-data.frame(flow=LIM$Unknowns, 
                 mean=colMeans(xs$X),
                 sd=sqrt(diag(var(xs$X))))

#check jumpsize before increase the number of iteration 
#check jumpsize=coverage of the whole range solution
samplerange <- data.frame(xranges)
samplerange$samplemin <- apply(xs$X, 2, min)
#apply(X, margin=1(rows) or 2(column), function)
samplerange$samplemax <- apply(xs$X, 2, max)
samplerange$percCover <- ((samplerange$samplemax - samplerange$samplemin) / (samplerange$max - samplerange$min) * 100)
print(samplerange, digits = 2)

# Mean percentage covered range
mean(samplerange$percCover, na.rm = T)
#increase jumpsize= decrease jumpsize parameter
#take a longer time to run

#segment plot: compare results from SSA and LA method
name<-LIM$Unknowns
windows()
dotchart(x=pars$X,col = 1,
         pch=16,xlim = c(0,400))
points(x=LA$mean,1:27,col=10,pch=18)
segments(LA$sd,1:27,LA$sd,1:27)
legend("right",pch=c(16,18,NA),lty=c(NA,NA,1),col=c(1,10,1),
       legend = c("Parsimonious","Montecarlo mean","sd"))



#caluculation
biosed<-data.frame(flow=c("SED-BAC",
                          "SED-MEI",
                          "SED-MAC",
                          "SED-MEG",
                          "SED-FIS"),
                   value=c(LA$mean[3]-LA$mean[8],
                           LA$mean[4]-LA$mean[13],
                           LA$mean[5]-LA$mean[17],
                           LA$mean[6]-LA$mean[21],
                           LA$mean[7]-LA$mean[25]))
sedpar<-data.frame(use=c(LA$flow[2:7]),
                   perc=c(LA$mean[2]/sum(LA$mean[2:7]),
                          LA$mean[3]/sum(LA$mean[2:7]),
                          LA$mean[4]/sum(LA$mean[2:7]),
                          LA$mean[5]/sum(LA$mean[2:7]),
                          LA$mean[6]/sum(LA$mean[2:7]),
                          LA$mean[7]/sum(LA$mean[2:7])))
bacpar<-data.frame(use=c(LA$flow[8:12]),
                   perc=c(LA$mean[8]/sum(LA$mean[8:12]),
                          LA$mean[9]/sum(LA$mean[8:12]),
                          LA$mean[10]/sum(LA$mean[8:12]),
                          LA$mean[11]/sum(LA$mean[8:12]),
                          LA$mean[12]/sum(LA$mean[8:12])))
meipar<-data.frame(use=c(LA$flow[13:16]),
                   perc=c(LA$mean[13]/sum(LA$mean[13:16]),
                          LA$mean[14]/sum(LA$mean[13:16]),
                          LA$mean[15]/sum(LA$mean[13:16]),
                          LA$mean[16]/sum(LA$mean[13:16])))
macpar<-data.frame(use=c(LA$flow[13:16]),
                   perc=c(LA$mean[13]/sum(LA$mean[13:16]),
                          LA$mean[14]/sum(LA$mean[13:16]),
                          LA$mean[15]/sum(LA$mean[13:16]),
                          LA$mean[16]/sum(LA$mean[13:16])))

# Sankey Diagram
library(networkD3)
nodes = data.frame("name" = 
                     c("SED", # Node 0
                       "BAC", # Node 1
                       "MEI", # Node 2
                       "MAC", # Node 3
                       "MEG", # Node 4
                       "FIS", # Node 5
                       "SED", # Node 6
                       "BAC", # Node 7
                       "MEI", # Node 8
                       "MAC", # Node 9
                       "MEG", # Node 10
                       "FIS", # Node 11
                       "POC_W",# Node 12
                       "DIC_W",# Node 13
                       "EXP_B",# Node 14
                       "EXP_S" # Node 15
                       ))
# Each row represents a link. The first number
# represents the node being conntected from. 
valuelinks = as.data.frame(matrix(c(
  0, 7, LA$mean[3], 0, 8, LA$mean[4], 0, 9, LA$mean[5], 0, 10, LA$mean[6], 0, 11, LA$mean[7],
  1, 6, LA$mean[8], 1, 8, LA$mean[9], 1, 9, LA$mean[10], 1, 10, LA$mean[11],
  2, 6, LA$mean[13], 2, 9, LA$mean[14], 2, 10, LA$mean[15],
  3, 6, LA$mean[17], 3, 10, LA$mean[18], 3, 11, LA$mean[19],
  4, 6, LA$mean[21], 4, 11, LA$mean[22],
  5, 6, LA$mean[25],
  12,0, LA$mean[1], 
  1,13, LA$mean[12], 2,13,LA$mean[16], 3,13, LA$mean[20], 4,13, LA$mean[23],5,13, LA$mean[26],  
  4,14,  LA$mean[24], 5,14, LA$mean[27],
  0,15, LA$mean[2]),# The third number is the value of the node
  byrow = TRUE, ncol = 3))
names(valuelinks) = c("source", "target", "value")

valuelinks$group = c("type_0", "type_0","type_0","type_0","type_0",
                "type_1", "type_1", "type_1","type_1", 
                "type_2","type_2","type_2",
                "type_3","type_3","type_3",
                "type_4","type_4",
                "type_5",
                "type_6",
                "type_1","type_2","type_3","type_4","type_5",
                "type_4","type_5",
                "type_0")
## Create custom color list using d3 for each node
library(RColorBrewer)
display.brewer.all()
brewer.pal(n = 12, name = "Spectral")
node_color <- 'd3.scaleOrdinal().domain(["SED","BAC","MEI","MAC","MEG","FIS",
                                        "POC_W","DIC_W","EXP_B","EXP_S",
                                        "type_0", "type_1", "type_2","type_3", "type_4","type_5","type_6"])
                                .range(["#9E0142","#D53E4F","#F46D43","#FDAE61","#FEE08B","#E6F598",
                                        "#ABDDA4","#66C2A5","#3288BD","#5E4FA2",
                                        "#9E0142","#D53E4F","#F46D43","#FDAE61","#FEE08B","#E6F598", "#ABDDA4"])'
## Draw Sankey Diagram (value)
p = sankeyNetwork(Links = valuelinks, Nodes = nodes,
                  Source = "source", Target = "target",
                  Value = "value", NodeID = "name",
                  fontSize = 16, nodeWidth = 40,
                  colourScale = node_color,
                  LinkGroup="group",unit="mgC/m2/d")%>%saveNetwork("MT_1_value.html")

perclinks = as.data.frame(matrix(c(
  0, 7, LA$mean[3]/sum(LA$mean[2:7]), 0, 8, LA$mean[4]/sum(LA$mean[2:7]), 0, 9, LA$mean[5]/sum(LA$mean[2:7]), 0, 10, LA$mean[6]/sum(LA$mean[2:7]), 0, 11, LA$mean[7]/sum(LA$mean[2:7]),
  1, 6, LA$mean[8]/sum(LA$mean[8:12]), 1, 8, LA$mean[9]/sum(LA$mean[8:12]), 1, 9, LA$mean[10]/sum(LA$mean[8:12]), 1, 10, LA$mean[11]/sum(LA$mean[8:12]),
  2, 6, LA$mean[13]/sum(LA$mean[13:16]), 2, 9, LA$mean[14]/sum(LA$mean[13:16]), 2, 10, LA$mean[15]/sum(LA$mean[13:16]),
  3, 6, LA$mean[17]/sum(LA$mean[17:20]), 3, 10, LA$mean[18]/sum(LA$mean[17:20]), 3, 11, LA$mean[19]/sum(LA$mean[17:20]),
  4, 6, LA$mean[21]/sum(LA$mean[21:24]), 4, 11, LA$mean[22]/sum(LA$mean[21:24]),
  5, 6, LA$mean[25]/sum(LA$mean[25:27]),
  12,0, LA$mean[1]/LA$mean[1], 
  1,13, LA$mean[12]/sum(LA$mean[8:12]), 2, 13,LA$mean[16]/sum(LA$mean[13:16]), 3, 13, LA$mean[20]/sum(LA$mean[17:20]), 4,13,  LA$mean[23]/sum(LA$mean[21:24]), 5,13, LA$mean[26]/sum(LA$mean[25:27]),  
  4, 14,LA$mean[24]/sum(LA$mean[21:24]),5, 14, LA$mean[27]/sum(LA$mean[25:27]),
  0,15, LA$mean[2]/sum(LA$mean[2:7])),# The third number is the value of the node
  byrow = TRUE, ncol = 3))
perclinks$V3<-perclinks$V3*100
names(perclinks) = c("source", "target", "value")

perclinks$group = c("type_0", "type_0","type_0","type_0","type_0",
                     "type_1", "type_1", "type_1","type_1", 
                     "type_2","type_2","type_2",
                     "type_3","type_3","type_3",
                     "type_4","type_4",
                     "type_5",
                     "type_6",
                     "type_1","type_2","type_3","type_4","type_5",
                     "type_4","type_5",
                     "type_0")
## Draw Sankey Diagram (%)
p = sankeyNetwork(Links = perclinks, Nodes = nodes,
                  Source = "source", Target = "target",
                  Value = "value", NodeID = "name",
                  fontSize = 16, nodeWidth = 40,
                  colourScale = node_color,
                  LinkGroup="group",unit="%")%>%saveNetwork("MT_1_percentage.html")

library(wesanderson)
  ?wes_palette()
wes_palette("Royal1")

#Extracts and sum the flows which formed the variables 
#Varranges():find solution ranges of the defined variable
varrange<-Varranges(LIM)
#Need to rethink again!
#create new df to restore data
# parsval<- data.frame(varrange)
# parsval$parsimonious<-numeric(LIM$NVariables)
# vareq<-LIM$VarA #contain matrix which define the variable equations
# parvec<-LIM$Parameters$val#vecotr with all parameter values in the right order
# flowvec<-pars$X#vector with all parsimonious flow values(also levels)
# #loop which calculate the parsimonious solution of the first variable
# #updates the solution vector
# #calculate the solution of the next variables
# #updates the solution vector
# for(i in 1:LIM$NVariables){
#   #refresh parsimonious solution vector
#   varvec<-parsval$parsimonious
#   #get subset of the equation matrix which all have the same equation nr
#   subset<-vareq[vareq$nr==i,]
# 
#   #take the parameter, variable or flow nr
#   #from the subset to use as index to find the corresponding
#   #parsimonious values and multiply them by 'val', which is the
#   #coefficient (1 or -1), and sum them
#   sum<-
#     sum(parvec[subset$par1]*subset$val,na.rm = T)+
#     sum(parvec[subset$par2]*subset$val,na.rm = T)+
#     sum(parvec[subset$par3]*subset$val,na.rm = T)+
#     sum(parvec[subset$par4]*subset$val,na.rm = T)+
#     sum(varvec[subset$var]*subset$val, na.rm = T)+
#     sum(flowvec[subset$flow]*subset$val, na.rm = T)
#   #replace 0 in df to actual sum value
#   parsval$parsimonious[i]=sum
# 
# }
# parsval