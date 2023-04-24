library(dplyr) 
library(vegan) 
library(glmmTMB)
library(splines)
library(emmeans)
library(ggplot2)
library(plyr)
library(gridExtra) 
library(lme4)
library(readr)
library(ggpubr)

########### Community Composition ###########
copyFishData <- read.csv("copyFishData.csv", fileEncoding="UTF-8-BOM", check.names = F)
rownames(copyFishData) <- copyFishData$Code
copyFishData <- copyFishData %>% dplyr::select(-Code)

#create subset of data with only sets 1 and no May data
data2 <- subset(copyFishData, substr(rownames(copyFishData), start = 8, stop = 8) == "1" & substr(rownames(copyFishData), start = 3, stop = 4) != "05")
data2 <- subset(data2, row.names(data2) != "100604B1")
  
#compares community composition between years for each site (all data) -- no significant changes
getNMDSsigYear <- function(data, SiteName){
  set.seed(1)
  if(SiteName == "all"){
    subData = data
  }
  else{
    subData = subset(data, substr(rownames(data), start = 7, stop = 7) ==  substr(SiteName, start=1, stop=1))
  }
  info <- data.frame(Year = gsub(" ", "", paste("20", substr(rownames(subData), start=1, stop=2))))
  fit <- adonis2(subData ~ Year, data = info, permutations = 999, method = "bray")
  fit
}
getNMDSsigYear(data2, "Slaughterhouse") #P-value = 0.508
getNMDSsigYear(data2, "Old Orchard") #P-value = 0.778
getNMDSsigYear(data2, "Reed Point") #P-value = 0.688
getNMDSsigYear(data2, "Barnett") #P-value = 0.604
getNMDSsigYear(data2, "Dockrill") #P-value = 0.056
getNMDSsigYear(data2, "all") #P-value = 0.159

#Creates Graphs for each site
subGraph = function(dataFrame, SiteName){
  if(SiteName == "Total"){
    subData = dataFrame
  }else{
    subData = subset(dataFrame, substr(rownames(dataFrame), start = 7, stop = 7) ==  substr(SiteName, start=1, stop=1))
  }
  #Create NMDS graph
  set.seed(1)
  NMDS <- metaMDS(subData, distance = "bray", k = 2, try = 20, trymax = 500)
  
  data.scores <- as.data.frame(scores(NMDS)$sites)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
  data.scores$site <- substr(rownames(data.scores), start = 7, stop = 7)  # create a column of site names, from the rownames of data.scores
  data.scores$Year <- gsub(" ", "", paste("20", substr(rownames(data.scores), start=1, stop=2)))
  
  Y10 <- data.scores[data.scores$Year == "2010", ][chull(data.scores[data.scores$Year == 
                                                                       "2010", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
  Y22 <- data.scores[data.scores$Year == "2022", ][chull(data.scores[data.scores$Year == 
                                                                       "2022", c("NMDS1", "NMDS2")]), ]  # hull values for grp B
  hull.data <- rbind(Y10, Y22)
  gg <- ggplot() + 
    geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=Year,group=Year),alpha=0.30) +
    geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=Year,colour=Year),size=3) + # add the point markers
    scale_colour_manual(values=c("2010" = "red", "2022" = "blue")) +
    coord_equal() +
    labs(title = SiteName) +
    theme_bw() +
    theme(axis.text.x = element_blank(),  # remove x-axis text
          axis.text.y = element_blank(), # remove y-axis text
          axis.ticks = element_blank(),  # remove axis ticks
          axis.title.x = element_text(size=12), # remove x-axis labels
          axis.title.y = element_text(size=12), # remove y-axis labels
          panel.background = element_blank(), 
          panel.grid.major = element_blank(),  #remove major-grid labels
          panel.grid.minor = element_blank(),  #remove minor-grid labels
          plot.background = element_blank(),
          plot.title = element_text(hjust = 0.5),
          text = element_text(size = 12),
          aspect.ratio = 1)
  gg
}
subGraph(data2, "Total")
#Creates grid arrange with all 5 sites based on a given data frame
createNMDSgraph <- function(dataFrame, graphTitle){
  s1 = subGraph(dataFrame, "Barnett")
  s2 = subGraph(dataFrame, "Dockrill")
  s3 = subGraph(dataFrame, "Old Orchard")
  s4 = subGraph(dataFrame, "Reed Point")
  s5 = subGraph(dataFrame, "Slaughterhouse")
  #s6 = subGraph(dataFrame, "Total")
  
  get_legend<-function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
  s1 <- s1 + theme(legend.text = element_text(size=12), legend.title = element_text(size=15), legend.key.size = unit(1, 'cm'))
  legend <- get_legend(s1)
  s1 <- s1 + theme(legend.position="none")
  s2 <- s2 + theme(legend.position="none")
  s3 <- s3 + theme(legend.position="none")
  s4 <- s4 + theme(legend.position="none")
  s5 <- s5 + theme(legend.position="none")
  #s6 <- s6 + theme(legend.position="none")
  
  title1 = text_grob(graphTitle, size = 15)
  nmdsGraph <- grid.arrange(s1, s2, s3, s4, s5, legend, nrow=2, top = title1, ncol=3)
  nmdsGraph
}

nmdsGraph1 <- createNMDSgraph(data2, "")

########### Abundance, Diversity, and Richness Models #############
modelData <- read.csv("modelData.csv", fileEncoding="UTF-8-BOM", check.names = F)
rownames(modelData) <- modelData$Code
modelData <- modelData %>% dplyr::select(-Code)
modelData$Site <- as.factor(modelData$Site)
modelData$Year <- as.factor(modelData$Year)

quasi_table <- function(model) {
  ctab=coef(summary(model))
  phi <- sum(residuals(model, type="pearson")^2)/df.residual(model) 
  qctab <- within(as.data.frame(ctab$cond),
                  {   `Std. Error` <- `Std. Error`*sqrt(phi)
                  `z value` <- Estimate/`Std. Error`
                  `Pr(>|z|)` <- 2*pnorm(abs(`z value`), lower.tail=FALSE)
                  })
  return(list(qctab = qctab, phi = phi))
}

#Richness
model1 = glmmTMB(Richness ~ Year + bs(DSJ, degree=3) + (1|Site), data=modelData, family="poisson")
summary(model1)
#Abundance
model2 <- glmmTMB(Abundance ~ Year + bs(DSJ) + (1|Site), family = "poisson", data = modelData)
quasi_table(model2)
#Diversity
model3 <- lmer(H ~ Year + bs(DSJ) + (1|Site), data=modelData)
coefs <- data.frame(coef(summary(model3)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

newdata =  expand.grid(Year=c(2010, 2022), DSJ=seq(160, 235, 2),  Site=NA) # 156 is June 5; 235 is Aug 23
pred.newdata = predict(model1, newdata=newdata, type='response', se=T, allow.new.levels=TRUE)

plot(newdata$DSJ[newdata$Year==2022], pred.newdata$fit[newdata$Year==2022], col='#0000FF', 
     ylim=c(0,8),pch=16,type='l', lwd=3, axes=F, ylab="", xlab="")
polygon(c(newdata$DSJ[newdata$Year==2022],rev(newdata$DSJ[newdata$Year==2022])),
        c(pred.newdata$fit[newdata$Year==2022]+1.96*pred.newdata$se.fit[newdata$Year==2022],
          rev(pred.newdata$fit[newdata$Year==2022]-1.96*pred.newdata$se.fit[newdata$Year==2022])),
        col=rgb(0/256,191/256,196/256, .3), border=F)
points(newdata$DSJ[newdata$Year==2010], 	pred.newdata$fit[newdata$Year==2010],col='#FF0000', pch=16, type='l', cex=.3, lwd=3)
polygon(c(newdata$DSJ[newdata$Year==2010],rev(newdata$DSJ[newdata$Year==2010])),
        c(pred.newdata$fit[newdata$Year==2010]+1.96*pred.newdata$se.fit[newdata$Year==2010],
          rev(pred.newdata$fit[newdata$Year==2010]-1.96*pred.newdata$se.fit[newdata$Year==2010])),
        col=rgb(248/256,118/256,109/256, .3), border=F)
box()
text(seq(160, 235, 14), rep(-.6, length(seq(160, 235, 14))), srt=45, adj=1, xpd=NA,
     labels=as.Date(seq(160, 235, 14), origin="01/01/22", "%m/%d/%y"))
axis(2, las=2)
axis(1, at=seq(160, 235, 14), labels=F)
mtext(side=2, line=3, "Fish Species Richness")
legend('bottomright', bty='n', legend=c('2010', 2022), fill=c(rgb(248/256,118/256,109/256, .5), rgb(0/256,191/256,196/256, .5)), border=F, pt.cex=3)

newdata2 =  expand.grid(Year=c(2010, 2022), DSJ=197,  Site=NA) # 156 is June 5; 235 is Aug 23
pred.newdata2 = predict(model1, newdata=newdata2, type='response', se=T, allow.new.levels=TRUE)

#2010 95% confidence interval
3.3 + 1.96*pred.newdata2$se.fit[newdata2$Year==2010]
3.3 - 1.96*pred.newdata2$se.fit[newdata2$Year==2010]
#2022 95% confidence interval
5.5 + 1.96*pred.newdata2$se.fit[newdata2$Year==2022]
5.5 - 1.96*pred.newdata2$se.fit[newdata2$Year==2022]

########### Length Violin Graphs ###########
lengthData <- read.csv("FishLengthData.csv", fileEncoding="UTF-8-BOM", check.names = F)[1:250,]

createDataFrame <- function(){
  species <- list()
  dates <- list()
  sets <- list()
  sites <- list()
  times <- list()
  numbers <- list()
  lengths <- list()
  for(index in 1:nrow(lengthData)){
    number = ""
    if(lengthData$Length[index] != ""){
      for(char in 1:nchar(lengthData$Length[index])){
        ascii = as.numeric(charToRaw(substr(lengthData$Length[index], char, char)))
        if(ascii >= 48 && ascii <= 57){
          number = paste(number, substr(lengthData$Length[index], char, char), sep = "")
          if(char == nchar(lengthData$Length[index])){
            numbers <- append(numbers, lengthData$Number[index])
            species <- append(species, lengthData$Species[index])
            lengths <- append(lengths, as.numeric(number))
            sites <- append(sites, lengthData$Site[index])
            sets <- append(sets, lengthData$Set[index])
            times <- append(times, lengthData$Time[index])
            dates <- append(dates, lengthData$Date[index])
          }
        }else{
          if(number != ""){
            numbers <- append(numbers, lengthData$Number[index])
            species <- append(species, lengthData$Species[index])
            lengths <- append(lengths, as.numeric(number))
            sites <- append(sites, lengthData$Site[index])
            sets <- append(sets, lengthData$Set[index])
            times <- append(times, lengthData$Time[index])
            dates <- append(dates, lengthData$Date[index])
            number = ""
          }
        }
      }
    }else{
      numbers <- append(numbers, lengthData$Number[index])
      species <- append(species, lengthData$Species[index])
      lengths <- append(lengths, NA)
      sites <- append(sites, lengthData$Site[index])
      sets <- append(sets, lengthData$Set[index])
      times <- append(times, lengthData$Time[index])
      dates <- append(dates, lengthData$Date[index])
    }
  }
  
  totalLength <- data.frame(unlist(dates), unlist(sites), unlist(sets), unlist(times), unlist(numbers), unlist(species), unlist(lengths))
  names(totalLength) <- c("Date", "Site", "Set", "Time", "Number", "Species", "Length")
  totalLength
}
totalLength <- createDataFrame()
totalLength$Date <- as.Date(totalLength$Date, "%Y-%m-%d")
totalLength$Date <- format(as.Date(totalLength$Date), "%m/%d")
totalLength$Length <- as.numeric(totalLength$Length)
totalLength <- subset(totalLength, Species == "Arrow Goby" | Species == "Staghorn Sculpin" | Species == "Shiner Perch" | Species == "Three-Spined Stickleback")


uniqueDF <- totalLength %>% distinct(Date, Site, Set, Species, .keep_all=TRUE)
allIndividuals <- data.frame()
for(row in 1:nrow(uniqueDF)){
  subset <- subset(totalLength, Date == uniqueDF$Date[row] & Site == uniqueDF$Site[row] & Set == uniqueDF$Set[row] & Time == uniqueDF$Time[row] & Number == uniqueDF$Number[row] & Species == uniqueDF$Species[row])
  numRow = nrow(subset)
  numIndividuals = as.integer(uniqueDF$Number[row])
  
  species <- list()
  dates <- list()
  lengths <- list()
  for(copy in 1:round(numIndividuals/numRow)){
    species <- append(species, subset$Species)
    lengths <- append(lengths, subset$Length)
    dates <- append(dates, subset$Date)
  }
  dfAppend <- data.frame(unlist(dates), unlist(species), unlist(lengths))
  allIndividuals <- rbind(allIndividuals, dfAppend)
}
names(allIndividuals) <- c("Date", "Species", "Length")

lengthSubset <- function(speciesName){
  subset(allIndividuals, Species == speciesName)
}

goby <- lengthSubset("Arrow Goby")
perch <- lengthSubset("Shiner Perch")
staghorn <- lengthSubset("Staghorn Sculpin")
stickleback <- lengthSubset("Three-Spined Stickleback")

violinGraph <- function(df, title){
  p <- ggplot(df, aes(x=as.factor(Date), y=Length)) + 
    geom_violin() + 
    #xlab("Date") + 
    #ylab("Length (in cm)") +
    xlab("") + 
    ylab("") +
    ggtitle(title)+
    theme_bw() +
    theme(axis.ticks = element_blank(),  # remove axis ticks
          axis.title.x = element_text(size=12), # remove x-axis labels
          axis.title.y = element_text(size=12), # remove y-axis labels
          panel.background = element_blank(), 
          panel.grid.major = element_blank(),  #remove major-grid labels
          panel.grid.minor = element_blank(),  #remove minor-grid labels
          plot.background = element_blank(),
          plot.title = element_text(hjust = 0.5),
          text = element_text(size = 12),
          #axis.text.x = element_text(angle = 40, vjust = 0.5, hjust = 0.5), 
          aspect.ratio = 0.5)
  p
}
s1 <- violinGraph(goby, "Arrow Goby")
s2 <- violinGraph(perch, "Shiner Perch")
s3 <- violinGraph(staghorn, "Staghorn Sculpin")
s4 <- violinGraph(stickleback, "Three-Spined Stickleback")
grid.arrange(s1, s2, s3, s4, nrow=2, ncol=2)


########### Beach Profiles ########
beachProfiles <- read.csv("BeachProfiles.csv", fileEncoding="UTF-8-BOM", check.names = F)
yearSiteSubset <- function(year, site){
  subset <- subset(beachProfiles, Year == year & Site == site)
  subset$Distance = as.numeric(subset$Distance)
  subset$Height = as.numeric(subset$Height)
  if(year == 2022){
    max <- max(subset$Distance)
    subset$Distance <- max - subset$Distance
  }
  profiles <- subset %>% ggplot(aes(x=as.numeric(Distance), y=as.numeric(Height))) +
    geom_line() + 
    xlab("Horizontal Distance (m)") + 
    ylab("Height (m)") +
    ggtitle(site) +
    theme_bw() +
    theme(axis.ticks = element_blank(),  # remove axis ticks
          axis.title.x = element_text(size=10), # remove x-axis labels
          axis.title.y = element_text(size=10), # remove y-axis labels
          panel.background = element_blank(), 
          panel.grid.major = element_blank(),  #remove major-grid labels
          panel.grid.minor = element_blank(),  #remove minor-grid labels
          plot.background = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          text = element_text(size = 12))
  if(site == "Slaughterhouse" & year == 2010){
    profiles <- profiles + ylim(0,4)
  }
  profiles
}
B10 <- yearSiteSubset(2010, "Barnett")
D10 <- yearSiteSubset(2010, "Dockrill")
O10 <- yearSiteSubset(2010, "Old Orchard")
R10 <- yearSiteSubset(2010, "Reed Point")
S10 <- yearSiteSubset(2010, "Slaughterhouse")

B22 <- yearSiteSubset(2022, "Barnett")
D22 <- yearSiteSubset(2022, "Dockrill")
O22 <- yearSiteSubset(2022, "Old Orchard")
R22 <- yearSiteSubset(2022, "Reed Point")
S22 <- yearSiteSubset(2022, "Slaughterhouse")

profileGraph10 <- grid.arrange(B10, D10, O10, R10, S10, nrow=3, ncol=2)
profileGraph10 

profileGraph22 <- grid.arrange(B22, D22, O22, R22, S22, nrow=3, ncol=2)
profileGraph22 

########### Site Richness and Diversity Boxplots ########
Year = gsub(" ", "", paste("20", substr(rownames(data2), start=1, stop=2)))
Site = substr(rownames(data2), start = 7, stop = 7)
H = diversity(data2, "shannon")
richness = rowSums(data2!=0)
YearSite = paste(Site, Year)
d <- data.frame(Year, Site, YearSite, H, richness)

# Average diversity at each site (2010 and 2022)
tapply(d$H, d$YearSite, mean)
# Average richness at each site (2010 and 2022)
tapply(d$richness, d$YearSite, mean)

getSWDgraph = function(data, metric){
  if(metric == "Richness"){
    ggplot(d, aes(x=Site, y=richness, fill=Year)) + 
      geom_boxplot() +
      #scale_fill_manual(values=c("#FD3A37", "#0E47ED")) +
      theme(panel.background = element_blank()) +
      ylab("Mean Number of Species per Set")
  }else if(metric == "Diversity"){
    ggplot(data, aes(x=Site, y=H, fill= Year)) + 
      geom_boxplot() +
      theme(panel.background = element_blank()) +
      ylab("Diversity")
  }else{
    print("Did not work. Use either Diversity or Richness")
  }
}
getSWDgraph(d, "Richness")
getSWDgraph(d, "Diversity")


########### Set Comparison ##########
setData <- read.csv("SetData.csv", fileEncoding="UTF-8-BOM", check.names = F)
rownames(setData) <- setData$Code
setData <- setData %>% dplyr::select(-Code)

DAR <- data.frame(H = diversity(data1), Abundance = rowSums(data1), Richness = rowSums(data1 != 0))
set1 <- subset(DAR, substr(rownames(DAR), start = 8, stop = 8) == "1")
set2 <- subset(DAR, substr(rownames(DAR), start = 8, stop = 8) == "2")

# COMPOSITION
getNMDSsigSite <- function(data){
  set.seed(1)
  info <- data.frame(Site = as.factor(substr(rownames(data), start=8, stop=8)))
  fit <- adonis2(data ~ Site, data = info, permutations = 999, method = "bray")
  fit
}
getNMDSsigSite(setData) #p-value = 0.95, R2 = 0.00859 -- so no difference in community composition

# DIVERSITY
wilcox.test(set1$H, set2$H) #p-value = 0.3808 -- so not significantly different

# ABUNDANCE
wilcox.test(set1$Abundance, set2$Abundance) #p-value = 0.9342 -- so abundance not significantly different

# RICHNESS
wilcox.test(set1$Richness, set2$Richness) #p-value = 0.09605 -- so richness not significantly different

# GRAPH - Composition set 1 vs 2
subGraph = function(){
  set.seed(1)
  NMDS <- metaMDS(setData, distance = "bray", k = 2, try = 20, trymax = 500)
  
  data.scores <- as.data.frame(scores(NMDS)$sites)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
  data.scores$set <- substr(rownames(data.scores), start = 8, stop = 8)  # create a column for set 1 and 2, from the rownames of data.scores
  
  set1 <- data.scores[data.scores$set == "1", ][chull(data.scores[data.scores$set == "1", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
  set2 <- data.scores[data.scores$set == "2", ][chull(data.scores[data.scores$set == "2", c("NMDS1", "NMDS2")]), ]  # hull values for grp B
  hull.data <- rbind(set1, set2)
  gg <- ggplot() + 
    geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=set,group=set),alpha=0.30) +
    geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=set,colour=set),size=3) + # add the point markers
    scale_colour_manual(values=c("1" = "red", "2" = "blue")) +
    coord_equal() +
    theme_bw() +
    theme(axis.text.x = element_blank(),  # remove x-axis text
          axis.text.y = element_blank(), # remove y-axis text
          axis.ticks = element_blank(),  # remove axis ticks
          axis.title.x = element_text(size=12), # remove x-axis labels
          axis.title.y = element_text(size=12), # remove y-axis labels
          panel.background = element_blank(), 
          panel.grid.major = element_blank(),  #remove major-grid labels
          panel.grid.minor = element_blank(),  #remove minor-grid labels
          plot.background = element_blank(),
          plot.title = element_text(hjust = 0.5),
          text = element_text(size = 12),
          aspect.ratio = 1)
  gg
}
subGraph()



