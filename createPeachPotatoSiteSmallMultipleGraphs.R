# create code to look at the January February temperature relationship between sites around the country and
# first flight of peach potato aphid
require(raster)

# emergeStack <- stack()
tempStack <- stack()

outDir <- "/project/gisclim/Projects/Active/PlantPest/peachPotWin/emergence1km/"
outDir <- "/project/gisclim/Projects/Active/PlantPest/peachPotWin/emergence1km/jfTemp/"

for (year in seq(1960,2017)) {
  outEmerge <- paste(outDir,"/emerge_",toString(year),'.tif',sep='')
  outTemp <-  paste(outDir,"/jf_",toString(year),'.tif',sep='')
  
  #emergeStack <- stack(emergeStack,raster(outEmerge))
  tempStack <- stack(tempStack,raster(outTemp))
}

sRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

getTmeanDoyCor <- function(tmean,doy){
  
  fit = lm(tmean ~ doy, na.action=na.exclude)
  lmResult <- lm(formula = doy ~ tmean, na.action=na.exclude)
  p1 = coef(lmResult)["(Intercept)"]
  p2 = coef(lmResult)[[2]]
  sumFit <- summary(fit)
  rsq = sumFit[8][[1]]
  rsq2 <- paste(round(rsq, 3))
  eq2 = paste0("y = ", round(p2, 3), "x +", round(p1, 1), "\n  RSq=", rsq2)
  
  print(sumFit)
  
  # return the intercept, correlation and slope
  return(c(p1,p2,round(rsq, 3)))
  
}


sites <- read.csv("/net/home/h05/hadnk/public_html/peachPot/SiteLookup1.csv")
sites$site2=""


for (row in 1:nrow(sites)) {
  inSite <- strsplit(as.character(sites$Sitename), "ST")[[row]][1]
  sites[row, "site2"] =inSite
}



peachStats <- read.csv("/project/gisclim/Projects/Active/PlantPest/peachPotRoth.csv")

peachStats <- read.csv("/project/gisclim/Projects/Active/PlantPest/peachPotRoth_min1StarCross.csv")

uniqueId <- unique(peachStats$id)
d1 <- combn(uniqueId,2)


uniqueId[length(uniqueId)+1]=999

df <- expand.grid(uniqueId,uniqueId)





# read in the DePreSys ensemble JanFeb Anomaly values, these are based on a 4.57C JF UK Climatology
df1 <- read.table("/net/home/h05/hadnk/public_html/peachPot/DePreSys_UK_Jan_Feb_anoms.txt")
names(df1) <- "val"
df1$typeStr="mod"

outDf <- data.frame(id=uniqueId,lat=rep(0,length(uniqueId)),long=rep(0,length(uniqueId)),numYears=rep(0,length(uniqueId)),siteName=rep("",length(uniqueId)),intercept=rep(0,length(uniqueId)),slope=rep(0,length(uniqueId)),rsq1st=rep(0,length(uniqueId)),rsq5th=rep(0,length(uniqueId)),rsq10th=rep(0,length(uniqueId)),rsq25th=rep(0,length(uniqueId)),rsq50th=rep(0,length(uniqueId)),rsq1st_5th=rep(0,length(uniqueId)),stringsAsFactors = FALSE)

myBigList = list()

cnt=1
stYear=1950


# loop through the sites
for (row in 1:nrow(sites)) {
  #for (row in seq(1,4)) {
  
  x <- sites[row, "X"]
  y <- sites[row, "Y"]
  lat <- sites[row, "Lat"]
  long <- sites[row, "Long"]
  
  idIn <- sites[row, "id"]
  siteName <- sites[row, "Sitename"]
  siteExists = (idIn %in% uniqueId)
  
  if (siteExists) {
    peachSite <- subset(peachStats,peachStats$id==idIn)	
    point <- cbind(x,y)
    jfTemp <- as.numeric(extract(tempStack,point))
    #ids <- rep(id,58)
    Season=seq(1960,2017)
    fjDF <- data.frame(Season,jfTemp)
    
    mergeDf1 <- na.omit(merge(fjDF,peachSite))
    mergeDf <- subset(mergeDf1,mergeDf1$Season>stYear)
    
    myBigList[[cnt]] <- mergeDf
    cnt=cnt+1
    
    numYears <- nrow(mergeDf)
    eqSite<- getTmeanDoyCor(unlist(mergeDf["jfTemp"]),unlist(mergeDf["X1st.flight"]))
    eqSite1<- getTmeanDoyCor(unlist(mergeDf["jfTemp"]),unlist(mergeDf["X5..flight"]))
    eqSite2<- getTmeanDoyCor(unlist(mergeDf["jfTemp"]),unlist(mergeDf["X10..flight"]))
    eqSite3<- getTmeanDoyCor(unlist(mergeDf["jfTemp"]),unlist(mergeDf["X25..flight"]))
    eqSite4<- getTmeanDoyCor(unlist(mergeDf["jfTemp"]),unlist(mergeDf["X50..flight"]))
    
    eqRsq <- getTmeanDoyCor(unlist(mergeDf["X1st.flight"]),unlist(mergeDf["X5..flight"]))
    
    # get row with id we are looking at 
    row=which(outDf[,1] == idIn)
    
    outDf[row,"siteName"]=as.character(siteName)
    outDf[row,"lat"]=lat
    outDf[row,"long"]=long
    outDf[row,"intercept"]=eqSite[1]
    outDf[row,"slope"]=eqSite[2]
    outDf[row,"numYears"]=numYears
    outDf[row,"rsq1st"]=eqSite[3]
    outDf[row,"rsq5th"]=eqSite1[3]
    outDf[row,"rsq10th"]=eqSite2[3]
    outDf[row,"rsq25th"]=eqSite3[3]
    outDf[row,"rsq50th"]=eqSite4[3]
    outDf[row,"rsq1st_5th"]=eqRsq[3]
    
    
  }
  
}

row=length(uniqueId)

allJf=do.call(rbind, myBigList)

eqSite<- getTmeanDoyCor(unlist(allJf["jfTemp"]),unlist(allJf["X1st.flight"]))
eqSite1<- getTmeanDoyCor(unlist(allJf["jfTemp"]),unlist(allJf["X5..flight"]))
eqSite2<- getTmeanDoyCor(unlist(allJf["jfTemp"]),unlist(allJf["X10..flight"]))
eqSite3<- getTmeanDoyCor(unlist(allJf["jfTemp"]),unlist(allJf["X25..flight"]))
eqSite4<- getTmeanDoyCor(unlist(allJf["jfTemp"]),unlist(allJf["X50..flight"]))
eqRsq <- getTmeanDoyCor(unlist(allJf["X1st.flight"]),unlist(allJf["X5..flight"]))

outDf[row,"siteName"]="ALL"
outDf[row,"lat"]=0
outDf[row,"long"]=0
outDf[row,"intercept"]=eqSite[1]
outDf[row,"slope"]=eqSite[2]
outDf[row,"numYears"]=numYears
outDf[row,"rsq1st"]=eqSite[3]
outDf[row,"rsq5th"]=eqSite1[3]
outDf[row,"rsq10th"]=eqSite2[3]
outDf[row,"rsq25th"]=eqSite3[3]
outDf[row,"rsq50th"]=eqSite4[3]
outDf[row,"rsq1st_5th"]=eqRsq[3]

write.csv(outDf,file=paste0(outDir,"siteRsq_StarcrossAdj.csv"),row.names=F,quote=F)

outDfStat <- outDf[-5,] 
outDfStat1 <- outDfStat[-nrow(outDfStat),] 

require(ggplot2)

legPal <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99")
legPal <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')

legPal <- c('#e6194b', '#3cb44b', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')



require(RColorBrewer)
myColors <- brewer.pal(11,"Paired"); 
myColors <- legPal[1:11]

names(myColors) <- levels(droplevels(allJf$Trap))


#scale_color_manual(values = c("foo" = "#999999", "bar" = "#E69F00"))

xmin <- min(allJf$jfTemp)
xmax <- max(allJf$jfTemp)
ymin <- min(allJf["X5..flight"])
ymax <- max(allJf["X5..flight"])

print(paste(xmin,xmax,ymin,ymax))



#for (i in seq(1,nrow(d1)) {
for (i in seq(2,2)) {    
  
  pairId <- d1[,i]
  
  selJF <- subset(allJf,id %in% pairId)
  
  trapNames <- as.character(unique(droplevels(selJF$Trap)))
  
  colsToUse <- myColors[which(attr(myColors,"names") %in% trapNames)]
  
  
  
  
  
}

sp







