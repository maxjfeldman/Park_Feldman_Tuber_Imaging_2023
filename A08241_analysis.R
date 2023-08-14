## This is a R script to perform analysis presented in Park, Feldman, et al., 2021

##########################################################################################
## Import Libraries
##########################################################################################

library(ggplot2)
library(reshape2)
library(dplyr)
library(lme4)
library(plotrix)
library(corrplot)
library(ggpmisc)
library(ggpubr)
library(GGally)
library(grDevices)
library(farver)


##########################################################################################
## Define functions
##########################################################################################

## Heritability fxn
get_h2<-function(data){
  variance.out<-c()
  H2<-c()
  e2<-c()
  # For each treatment.phenotype calculate variance
  for(i in 3:length(colnames(data))){
    print(i)
    # Use only RILs with all measurements for each treatment.phenotype
    cc.data<-data[complete.cases(data[,i]),c(1:2,i)]
    # Build linear model each cofactor is a random effect
    model<-lmer(cc.data[,3]~(1|clone), data=cc.data, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
    # Extract variance from model object, save individual components in vector 're' and residual variance as a scalar named 'res'
    re<-as.numeric(VarCorr(model))
    res<-attr(VarCorr(model), "sc")^2
    # Extract individual components (order will remain the same)
    geno.var<-re[1]
    # Total variance is sum of all variances
    tot.var<-sum(re, res)
    # Get proportion of variance
    h<-geno.var/tot.var
    e<-res/tot.var
    # Append variables to a vector of variables
    H2<-c(H2,h)
    e2<-c(e2,e)
  }
  variance<-rbind(H2, e2)
  colnames(variance)<-colnames(data)[3:length(data)]
  rownames(variance)<-c('Genotype', 'Error')
  return(variance)
}  


## Get fx estimates
get_fx<-function(data){
  variance.out<-c()
  genotype<-c()
  replicate<-c()
  side<-c()
  residual<-c()
  # For each treatment.phenotype calculate variance
  for(i in 4:length(colnames(data))){
    print(i)
    # Use only RILs with all measurements for each treatment.phenotype
    cc.data<-data[complete.cases(data[,i]),c(1:3,i)]
    # Build linear model each cofactor is a random effect
    model<-lmer(cc.data[,4]~(1|clone) + (1|rep) + (1|side), data=cc.data, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
    # Extract variance from model object, save individual components in vector 're' and residual variance as a scalar named 'res'
    re<-as.data.frame(VarCorr(model))
    # Extract individual components (order will remain the same)
    geno.var<-re[re$grp == 'clone', 'sdcor']
    rep.var<-re[re$grp == 'rep', 'sdcor']
    side.var<-re[re$grp == 'side', 'sdcor']
    res.var<-re[re$grp == 'Residual', 'sdcor']
    # Total variance is sum of all variances
    tot.var<-sum(re[,'sdcor'])
    # Get proportion of variance
    geno<-geno.var/tot.var
    rep<-rep.var/tot.var
    s<-side.var/tot.var
    e<-res.var/tot.var
    # Append variables to a vector of variables
    genotype<-c(genotype, geno)
    replicate<-c(replicate, rep)
    side<-c(side, s)
    residual<-c(residual,e)
  }
  variance.out<-rbind(genotype, replicate, side, residual)
  colnames(variance.out)<-colnames(data)[4:length(data)]
  rownames(variance.out)<-c('Genotype', 'Replicate', 'Side', 'Error')
  return(variance.out)
}  



## Get fx estimates
get_fx_scan<-function(data){
  variance.out<-c()
  genotype<-c()
  replicate<-c()
  residual<-c()
  # For each treatment.phenotype calculate variance
  for(i in 4:length(colnames(data))){
    print(i)
    # Use only RILs with all measurements for each treatment.phenotype
    cc.data<-data[complete.cases(data[,i]),c(1:2,i)]
    # Build linear model each cofactor is a random effect
    model<-lmer(cc.data[,3]~(1|clone) + (1|rep), data=cc.data, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
    # Extract variance from model object, save individual components in vector 're' and residual variance as a scalar named 'res'
    re<-as.data.frame(VarCorr(model))
    # Extract individual components (order will remain the same)
    geno.var<-re[re$grp == 'clone', 'sdcor']
    rep.var<-re[re$grp == 'rep', 'sdcor']
    res.var<-re[re$grp == 'Residual', 'sdcor']
    # Total variance is sum of all variances
    tot.var<-sum(re[,'sdcor'])
    # Get proportion of variance
    geno<-geno.var/tot.var
    rep<-rep.var/tot.var
    e<-res.var/tot.var
    # Append variables to a vector of variables
    genotype<-c(genotype, geno)
    replicate<-c(replicate, rep)
    residual<-c(residual,e)
  }
  variance.out<-rbind(genotype, replicate, residual)
  colnames(variance.out)<-colnames(data)[4:length(data)]
  rownames(variance.out)<-c('Genotype', 'Replicate', 'Error')
  return(variance.out)
}  




##########################################################################################
## Load data
##########################################################################################

## Define directories used and read in data
setwd("/Users/max.feldman/Documents/data/analysis_test")
base.dir<-getwd()

## Black background
std.filepath<-paste(base.dir, "/A08241_potato_measurements_std_shape_ColorCorrected.csv", sep="")
std<-read.csv(std.filepath)
std$red_ave<-gsub("\\[", "", std$red_ave)
std$red_ave<-gsub("\\]", "", std$red_ave)
std$red_ave<-as.numeric(as.character(std$red_ave))

std$green_ave<-gsub("\\[", "", std$green_ave)
std$green_ave<-gsub("\\]", "", std$green_ave)
std$green_ave<-as.numeric(as.character(std$green_ave))

std$blue_ave<-gsub("\\[", "", std$blue_ave)
std$blue_ave<-gsub("\\]", "", std$blue_ave)
std$blue_ave<-as.numeric(as.character(std$blue_ave))

std$red_sd<-gsub("\\[", "", std$red_sd)
std$red_sd<-gsub("\\]", "", std$red_sd)
std$red_sd<-as.numeric(as.character(std$red_sd))

std$green_sd<-gsub("\\[", "", std$green_sd)
std$green_sd<-gsub("\\]", "", std$green_sd)
std$green_sd<-as.numeric(as.character(std$green_sd))

std$blue_sd<-gsub("\\[", "", std$blue_sd)
std$blue_sd<-gsub("\\]", "", std$blue_sd)
std$blue_sd<-as.numeric(as.character(std$blue_sd))
std_color<-std[,c(2:21)]

## Extract shape values to seperate data.frame
shape<-std[,c(2:21,790:889)]

## Get standard values from black background
#std<-std[,c(1:789)]

## Lightbox background
box.filepath<-paste(base.dir, "/A08241_potato_measurements_box_shape_ColorCorrected.csv", sep="")
box<-read.csv(box.filepath)
box$red_ave<-gsub("\\[", "", box$red_ave)
box$red_ave<-gsub("\\]", "", box$red_ave)
box$red_ave<-as.numeric(as.character(box$red_ave))

box$green_ave<-gsub("\\[", "", box$green_ave)
box$green_ave<-gsub("\\]", "", box$green_ave)
box$green_ave<-as.numeric(as.character(box$green_ave))

box$blue_ave<-gsub("\\[", "", box$blue_ave)
box$blue_ave<-gsub("\\]", "", box$blue_ave)
box$blue_ave<-as.numeric(as.character(box$blue_ave))

box$red_sd<-gsub("\\[", "", box$red_sd)
box$red_sd<-gsub("\\]", "", box$red_sd)
box$red_sd<-as.numeric(as.character(box$red_sd))

box$green_sd<-gsub("\\[", "", box$green_sd)
box$green_sd<-gsub("\\]", "", box$green_sd)
box$green_sd<-as.numeric(as.character(box$green_sd))

box$blue_sd<-gsub("\\[", "", box$blue_sd)
box$blue_sd<-gsub("\\]", "", box$blue_sd)
box$blue_sd<-as.numeric(as.character(box$blue_sd))


## Combine datasets
both<-rbind(std, box)

## Load in scanner measurements
scan.filepath<-paste(base.dir, "/A08241_potato_measurements_scanner.csv", sep="")
scan<-read.csv(scan.filepath)

scan$red_ave<-gsub("\\[", "", scan$red_ave)
scan$red_ave<-gsub("\\]", "", scan$red_ave)
scan$red_ave<-as.numeric(as.character(scan$red_ave))

scan$green_ave<-gsub("\\[", "", scan$green_ave)
scan$green_ave<-gsub("\\]", "", scan$green_ave)
scan$green_ave<-as.numeric(as.character(scan$green_ave))

scan$blue_ave<-gsub("\\[", "", scan$blue_ave)
scan$blue_ave<-gsub("\\]", "", scan$blue_ave)
scan$blue_ave<-as.numeric(as.character(scan$blue_ave))

scan$red_sd<-gsub("\\[", "", scan$red_sd)
scan$red_sd<-gsub("\\]", "", scan$red_sd)
scan$red_sd<-as.numeric(as.character(scan$red_sd))

scan$green_sd<-gsub("\\[", "", scan$green_sd)
scan$green_sd<-gsub("\\]", "", scan$green_sd)
scan$green_sd<-as.numeric(as.character(scan$green_sd))

scan$blue_sd<-gsub("\\[", "", scan$blue_sd)
scan$blue_sd<-gsub("\\]", "", scan$blue_sd)
scan$blue_sd<-as.numeric(as.character(scan$blue_sd))

scan<-scan[,c(2:19)]

## Load in ground truth data
gt.filepath<-paste(base.dir, "/ground_truth_data.csv", sep="")
gt<-read.csv(gt.filepath)

## Lets make an directory to output figures for the manuscript
fig.filepath<-paste(base.dir, "/figures", sep="")
dir.create(fig.filepath)

## Lets make an directory to output tables for the manuscript
tab.filepath<-paste(base.dir, "/tables", sep="")
dir.create(tab.filepath)

##########################################################################################
## Analysis of size standards
##########################################################################################

## Get only values for the size marker (poker chip) from the top-down imaging config.

## Get only poker chip contours
marker.top_down<-both[both$tuber == 'marker',]
marker.scan<-scan[scan$tuber == 'marker',]

# Remove columns that are not helpful 
marker.top_down<-marker.top_down[,-which(names(marker.top_down) %in% c("X","img_name"))]
marker.scan<-marker.scan[,-which(names(marker.scan) %in% c("X","img_name"))]

## Lets only keep the images from side 1 because the scanner images only contain a single side
marker.top_down.s1<-marker.top_down[marker.top_down$side == '1',]
marker.top_down.s1<-marker.top_down.s1[,-which(names(marker.top_down.s1) %in% c("side"))]
marker.top_down.s1<-marker.top_down.s1[,c(1:18)]

## Lets add a categorical variable named 'light' to the scan data.frame 
## This identifies the collection platform
marker.scan$light<-rep("scanner", nrow(marker.scan))
marker.scan<-marker.scan[,names(marker.top_down.s1)]

## Lets make sure the column order are the same between the top-down and scan data.frames
## Then merge them using rbind
marker.all<-rbind(marker.top_down.s1, marker.scan)

## Lets calculate length and width in mm (poker chip size standard has diameter of 37 mm)
marker.all$mm_per_px<-37/marker.all$length
marker.all$length_mm<-marker.all$length * marker.all$mm_per_px
marker.all$width_mm<-marker.all$width * marker.all$mm_per_px
#marker.all$pct_area<-marker.all$area/mean(marker.all$area)
marker.all$rel_area<-marker.all$area * (marker.all$mm_per_px)^2
marker.all$light<-as.character(marker.all$light)
marker.all[marker.all$light == "std", "light"]<-c("Black background")
marker.all[marker.all$light == "box", "light"]<-c("Illumination box")
marker.all[marker.all$light == "scanner", "light"]<-c("Flatbed scanner")




## Subset data.frame and convert to long form
marker.all<-marker.all[,c("clone", "rep", "light", "tuber", "area", "rel_area", "length", "length_mm", "width", "width_mm", "ratio")]


id.vars<-colnames(marker.all)[c(1:4)]
measure.vars<-colnames(marker.all)[c(5:(ncol(marker.all)))]

marker.all_long<-melt(marker.all,
                  # ID variables - all the variables to keep but not split apart on
                  id.vars=id.vars,
                  # The source columns
                  measure.vars=measure.vars,
                  # Name of the destination column that will identify the original
                  # column that the measurement came from
                  variable.name="Trait",
                  value.name="Value"
)

colnames(marker.all_long)[3]<-c("Configuration")


## Here we can make a plot of all quantities
p<-ggplot(marker.all_long, aes(x=Configuration, y=Value, col=Configuration)) + geom_jitter(size = 0.1, width=0.02, height=0) + facet_wrap(~Trait, scales="free_y", ncol=2) + theme_bw()
x<-p + stat_summary(fun.y="mean", geom="point", size=6, shape=4, col="black") + scale_color_manual(values = c("red", "blue", "black"))
q<-x +  stat_summary(fun.data=mean_se, geom = "errorbar", width=0.2) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(q)

## Lets make Fig. 1A
rel_area<-marker.all_long[marker.all_long$Trait == 'rel_area',]

p<-ggplot(rel_area, aes(x=Configuration, y=Value, col=Configuration)) + geom_jitter(size = 0.1, width=0.02, height=0) + theme_bw()
x<-p + stat_summary(fun.y="mean", geom="point", size=6, shape=4, col="black") + scale_color_manual(values = c("red", "blue", "black"))
q<-x +  stat_summary(fun.data=mean_se, geom = "errorbar", width=0.2) + theme(text = element_text(size=15), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") + ylab("Area (mm2)") + xlab("")

fig.S1a_path<-paste(fig.filepath, "/Fig_S1a.pdf", sep="")
pdf(fig.S1a_path, height=6, width=6)
print(q)
dev.off()

## Fig. 1B
width_mm<-marker.all_long[marker.all_long$Trait == 'width_mm',]

p<-ggplot(width_mm, aes(x=Configuration, y=Value, col=Configuration)) + geom_jitter(size = 0.1, width=0.02, height=0) + theme_bw()
x<-p + stat_summary(fun.y="mean", geom="point", size=6, shape=4, col="black") + scale_color_manual(values = c("red", "blue", "black"))
q<-x +  stat_summary(fun.data=mean_se, geom = "errorbar", width=0.2) + theme(text = element_text(size=15), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")  + ylab("Distance (mm)") + xlab("")

fig.S1b_path<-paste(fig.filepath, "/Fig_S1b.pdf", sep="")
pdf(fig.S1b_path, height=6, width=6)
print(q)
dev.off()


## Fig. 1C
r_area.blk_bkgd<-marker.all_long[marker.all_long$Trait == 'rel_area' & marker.all_long$Configuration == "Black background",]
r_area.ill_box<-marker.all_long[marker.all_long$Trait == 'rel_area' & marker.all_long$Configuration == "Illumination box",]
r_area.scan<-marker.all_long[marker.all_long$Trait == 'rel_area' & marker.all_long$Configuration == "Flatbed scanner",]

r_area.blk.ave<-mean(r_area.blk_bkgd$Value)
r_area.box.ave<-mean(r_area.ill_box$Value)
r_area.scan.ave<-mean(r_area.scan$Value)

r_area.blk.sd<-sd(r_area.blk_bkgd$Value)
r_area.box.sd<-sd(r_area.ill_box$Value)
r_area.scan.sd<-sd(r_area.scan$Value)

r_area.blk_bkgd$Pct_Error<-r_area.blk_bkgd$Value / r_area.blk.ave
r_area.ill_box$Pct_Error<-r_area.ill_box$Value / r_area.box.ave
r_area.scan$Pct_Error<-r_area.scan$Value / r_area.scan.ave

r_area.all<-rbind(r_area.blk_bkgd,r_area.ill_box,r_area.scan)
r_area.all$Pct_Error<-r_area.all$Pct_Error*100
r_area.all$Pct_Error<-r_area.all$Pct_Error-100

p<-ggplot(r_area.all, aes(x=Configuration, y=Pct_Error, col=Configuration)) + geom_jitter(size = 0.1, width=0.02, height=0) + theme_bw()
x<-p + stat_summary(fun.y="mean", geom="point", size=6, shape=4, col="black") + scale_color_manual(values = c("red", "blue", "black"))
q<-x +  stat_summary(fun.data=mean_se, geom = "errorbar", width=0.2) + theme(text = element_text(size=15), axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("Difference from average (% area)") + xlab("")


fig.S1c_path<-paste(fig.filepath, "/Fig_S1c.pdf", sep="")
pdf(fig.S1c_path, height=6, width=8)
print(q)
dev.off()


## Do the same for distance measurment (width_mm)
r_dist.blk_bkgd<-marker.all_long[marker.all_long$Trait == 'width_mm' & marker.all_long$Configuration == "Black background",]
r_dist.ill_box<-marker.all_long[marker.all_long$Trait == 'width_mm' & marker.all_long$Configuration == "Illumination box",]
r_dist.scan<-marker.all_long[marker.all_long$Trait == 'width_mm' & marker.all_long$Configuration == "Flatbed scanner",]

r_dist.blk.ave<-mean(r_dist.blk_bkgd$Value)
r_dist.box.ave<-mean(r_dist.ill_box$Value)
r_dist.scan.ave<-mean(r_dist.scan$Value)

## calculate standard deviation
r_dist.blk.sd<-sd(r_dist.blk_bkgd$Value)
r_dist.box.sd<-sd(r_dist.ill_box$Value)
r_dist.scan.sd<-sd(r_dist.scan$Value)

r_dist.blk_bkgd$Pct_Error<-r_dist.blk_bkgd$Value / r_dist.blk.ave
r_dist.ill_box$Pct_Error<-r_dist.ill_box$Value / r_dist.box.ave
r_dist.scan$Pct_Error<-r_dist.scan$Value / r_dist.scan.ave

r_dist.all<-rbind(r_dist.blk_bkgd,r_dist.ill_box,r_dist.scan)
r_dist.all$Pct_Error<-r_dist.all$Pct_Error*100
r_dist.all$Pct_Error<-r_dist.all$Pct_Error-100


p<-ggplot(r_dist.all, aes(x=Configuration, y=Pct_Error, col=Configuration)) + geom_jitter(size = 0.1, width=0.02, height=0) + theme_bw()
x<-p + stat_summary(fun.y="mean", geom="point", size=6, shape=4, col="black") + scale_color_manual(values = c("red", "blue", "black"))
q<-x +  stat_summary(fun.data=mean_se, geom = "errorbar", width=0.2) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Imaging configuration") + ylab("Difference from average (%)")


##########################################################################################
## Analysis of tuber size on black background
##########################################################################################

## Convert values to mm
marker.std<-std[std$tuber == 'marker',]

imgs<-unique(marker.std$img_name)
std$px_per_mm<-c(NA)
for (i in imgs){
  px_per_mm<-37/marker.std[marker.std$img_name == i,'length']
  std[std$img_name == i, 'px_per_mm']<-px_per_mm
}


## Get px_per_mm for each img
std$area_mm<-std$area * (std$px_per_mm)^2
std$length_mm<-std$length * std$px_per_mm
std$width_mm<-std$width * std$px_per_mm
std$perimeter_mm<-std$perimeter * std$px_per_mm

std_size<-std[,c("img_name", "clone", "rep", "side", "tuber", "area", "perimeter", "length", "width", "ratio", "eccentricity", "red_ave", "green_ave", "blue_ave", "red_sd", "green_sd", "blue_sd", "area_mm", "length_mm", "width_mm","perimeter_mm")]
std_size<-std_size[std_size$tuber != 'marker',]

## Lets get the measurements only from side 1
std_size.one<-std_size[std_size$side == 1,]

std_size.one<-std_size.one[,c(2,5:ncol(std_size.one))]
colnames(std_size.one)[1:2]<-c('clone', 'tuber')

## Merge the black background and ground truth measurements
tuber_size<-merge(std_size.one, gt, by=c("clone", "tuber"))

## Calculate correlation between various measurements
cor(tuber_size$area_mm, tuber_size$weight, use="complete.obs")

cor(tuber_size$length_mm, tuber_size$caliper_length, use="complete.obs")

cor(tuber_size$width_mm, tuber_size$caliper_width, use="complete.obs")

## Lets calculate % variance for each factor (clone, replicate, tuber size)
# Build linear model each cofactor is a random effect
model<-lmer(std_size$area_mm~(1|clone)+(1|rep)+(1|side), data=std_size)

re<-as.numeric(VarCorr(model))
res<-attr(VarCorr(model), "sc")^2
tot.var<-sum(re, res)

clone.var<-re[1]/tot.var
rep.var<-re[2]/tot.var
side.var<-re[3]/tot.var


## Lets calculate some summary statistics on the size measurements 
## Average
tuber_size.ave<-aggregate(.~ clone, data=tuber_size[,c(1,3:22)], mean)
colnames(tuber_size.ave)[2:ncol(tuber_size.ave)]<-paste(colnames(tuber_size.ave[2:ncol(tuber_size.ave)]), "ave", sep=".")

## Standard deviation
tuber_size.sd<-aggregate(.~ clone, data=tuber_size[,c(1,3:22)], sd)
colnames(tuber_size.sd)[2:ncol(tuber_size.sd)]<-paste(colnames(tuber_size.sd[2:ncol(tuber_size.sd)]), "sd", sep=".")

## Lets merge these parameters into the same data.frame
tuber_size.par<-merge(tuber_size.ave, tuber_size.sd, by=c("clone"))
tuber_size.par<-tuber_size.par[complete.cases(tuber_size.par),]


## Is there correlation between weight and variance?
cor(tuber_size.par$weight.ave, tuber_size.par$weight.sd)

## What is the mean size?
clone.ave<-mean(tuber_size.par$weight.ave)
## What is the standard devatiation of the mean? 
clone.ave.sd<-mean(tuber_size.par$weight.sd)

## How about min and max
clone.ave.min<-min(tuber_size.par$weight.ave)
clone.ave.max<-max(tuber_size.par$weight.ave)

## What about the co-efficient of variation?
clone.ave.sd/clone.ave

## Lets cacluation broad sense H2 for size measurements
tuber_size_for_h2<-tuber_size[,c(1:22)]

h2_tuber.size<-get_h2(tuber_size_for_h2)

## Lets do the same for tuber size variance
tuber_sd_for_h2<-tuber_size[,c(1:22)]
#tuber_sd_for_h2<-aggregate(.~ clone + replicate, data=tuber_sd_for_h2, sd)
tuber_sd_for_h2<-aggregate(.~ clone, data=tuber_sd_for_h2, sd)
tuber_sd_for_h2<-tuber_sd_for_h2[,-c(2,3)]

h2_tuber.size.sd<-get_h2(tuber_sd_for_h2)

## Make plots of correlation between machine vision and ground truth measurements
#p<-ggplot(tuber_size, aes(x=caliper_length, y=length_mm)) + geom_point(size=0.5, color=brown"") + theme_bw() + xlab("Caliper length (mm)") + ylab("Machine vision length (mm)") + theme(text = element_text(size=15),legend.position = "none") + geom_smooth(method = "lm", se = FALSE, color="red",linetype="dashed") + stat_cor(aes(size = 2, label =  ..rr.label..))
p<-ggplot(tuber_size, aes(x=caliper_length, y=length_mm)) + geom_point(size=0.5, color="grey30") + theme_bw() + xlab("Caliper length (mm)") + ylab("Machine vision length (mm)") + theme(text = element_text(size=15),legend.position = "none") + geom_smooth(method = "lm", se = FALSE, color="black",linetype="dashed") 

fig.1b_path<-paste(fig.filepath, "/Fig_1b.pdf", sep="")
pdf(fig.1b_path, height=4, width=4)
print(p)
dev.off()

fig.1c_path<-paste(fig.filepath, "/Fig_1c.pdf", sep="")

p<-ggplot(tuber_size, aes(x=caliper_width, y=width_mm)) + geom_point(size=0.5, color="grey50") + theme_bw() + xlab("Caliper width (mm)") + ylab("Machine vision width (mm)")  + theme(text = element_text(size=15),legend.position = "none")+ geom_smooth(method = "lm", se = FALSE, color="black",linetype="dashed") 


pdf(fig.1c_path, height=4, width=4)
print(p)
dev.off()

fig.1a_path<-paste(fig.filepath, "/Fig_1a.pdf", sep="")
p<-ggplot(tuber_size, aes(x=weight, y=area_mm)) + geom_point(size=0.5, color=c("grey10")) + theme_bw() + xlab("Weight (oz)") + ylab("Area (mm2)")  + theme(text = element_text(size=15),legend.position = "none") + geom_smooth(method = "lm", se = FALSE, color="black",linetype="dashed") 

pdf(fig.1a_path, height=4, width=4)
print(p)
dev.off()


tuber_size<-tuber_size[complete.cases(tuber_size),]
tuber_size$clone = reorder(tuber_size$clone, tuber_size$weight, median)
p<-ggplot(tuber_size, aes(x=as.factor(clone), y=weight), color="sienna4") + geom_boxplot() + theme_bw() + xlab("Clone") + ylab("Tuber size (oz)")  + theme(text = element_text(size=15),legend.position = "none",axis.text.x = element_text(angle = 90, hjust = 1, size=6))

fig.1d_path<-paste(fig.filepath, "/Fig_1d.pdf", sep="")

pdf(fig.1d_path, height=3, width=14)
print(p)
dev.off()



std_size.model<-std_size[,c(2:4,18:ncol(std_size))]
size.variance_output<-get_fx(std_size.model)


## Lets look at correlation of measurements between each tuber side
side1<-std_size[std_size$side == 1,c(2:5,18:21)]
side2<-std_size[std_size$side == 2,c(2:5,18:21)]
colnames(side1)[c(5:ncol(side1))]<-paste(colnames(side1)[c(5:ncol(side1))], "1", sep="_")
colnames(side2)[c(5:ncol(side2))]<-paste(colnames(side2)[c(5:ncol(side2))], "2", sep="_")

stds_size_compareSides<-merge(side1[,-c(3)], side2[,-c(3)], by=c('clone', 'rep', 'tuber'), all=T)
stds_size_compareSides<-stds_size_compareSides[complete.cases(stds_size_compareSides),]

area.cor<-cor(stds_size_compareSides$area_mm_1, stds_size_compareSides$area_mm_2)
length.cor<-cor(stds_size_compareSides$length_mm_1, stds_size_compareSides$length_mm_2)
width.cor<-cor(stds_size_compareSides$width_mm_1, stds_size_compareSides$width_mm_2)
perimeter.cor<-cor(stds_size_compareSides$perimeter_mm_1, stds_size_compareSides$perimeter_mm_2)

p<-ggplot(stds_size_compareSides, aes(x=area_mm_1, y=area_mm_2)) + geom_point() + theme_bw() + xlab("Area (mm2): Tuber side 1") +  ylab("Area (mm2): Tuber side 2")
p<-ggplot(stds_size_compareSides, aes(x=length_mm_1, y=length_mm_2)) + geom_point() + theme_bw() + xlab("Length (mm): Tuber side 1") +  ylab("Length (mm): Tuber side 2")



## Lets look at standard error for different sample sizes

## First step is to identify samples with 10 tubers

complete.clones<-names(table(tuber_size$clone)[table(tuber_size$clone) > 9])
tuber_size_complete<-tuber_size[tuber_size$clone %in% complete.clones,]

tuber_size_complete.z<-apply(tuber_size_complete[,c(3:22)], 2, scale)
tuber_size_complete.z<-cbind(tuber_size_complete[,c(1:2)], tuber_size_complete.z)

std_error_size<-c()
for(i in 3:10){
  clones<-unique(tuber_size_complete$clone)
  for (c in clones){
    clone.dat<-tuber_size_complete.z[tuber_size_complete.z$clone == c,]
    tubers<-sample(c(1:10), i, replace=F)
    small.dat<-clone.dat[clone.dat$tuber %in% tubers, ]
    weight.sd.err<-std.error(small.dat$weight)
    area.sd.err<-std.error(small.dat$area_mm)
    length.sd.err<-std.error(small.dat$length_mm)
    width.sd.err<-std.error(small.dat$width_mm)
    perimeter.sd.err<-std.error(small.dat$perimeter_mm)
    caliper_L.sd.err<-std.error(small.dat$caliper_length)
    caliper_W.sd.err<-std.error(small.dat$caliper_width)
    temp<-c(i,c,weight.sd.err,area.sd.err,length.sd.err, width.sd.err, perimeter.sd.err,caliper_L.sd.err,caliper_W.sd.err)
    std_error_size<-rbind(std_error_size, temp)
  }
}

std_error_size<-as.data.frame(std_error_size)
colnames(std_error_size)<-c("Replicates", "Clone", "Weight","Area","Length","Width","Perimeter","Caliper Length", "Caliper Width")

for(i in 1:ncol(std_error_size)){
  std_error_size[,i]<-as.numeric(as.character(std_error_size[,i]))
}

std_error_size.ag<-aggregate(.~Replicates, data=std_error_size[,-c(2)], mean)

id.vars<-colnames(std_error_size.ag)[c(1)]
measure.vars<-colnames(std_error_size.ag)[c(2:(ncol(std_error_size.ag)))]
std_error_size.ag.long<-melt(std_error_size.ag,
                    # ID variables - all the variables to keep but not split apart on
                    id.vars=id.vars,
                    # The source columns
                    measure.vars=measure.vars,
                    # Name of the destination column that will identify the original
                    # column that the measurement came from
                    variable.name="Trait",
                    value.name="SE"
)


p<-ggplot(std_error_size.ag.long, aes(x=Replicates, y=SE, color=Trait)) + geom_line() + theme_bw()
q<-p+geom_point() + xlab("Number of tubers") + ylab("Standard Error")


#fig.S3a_path<-paste(fig.filepath, "/Fig_S3a.pdf", sep="")
#pdf(fig.1d_path, height=8, width=5)
#print(p)
#dev.off()

##########################################################################################
## Analysis of size standard shape on black background
##########################################################################################

marker_shape<-marker.top_down[marker.top_down$light == 'std' & marker.top_down$side == 1,]

marker_shape<-marker_shape[,c(1:2,5,8:13)]

## Lets calculate length and width in mm
marker_shape$mm_per_px<-37/marker_shape$length
marker_shape$length_mm<-marker_shape$length * marker_shape$mm_per_px
marker_shape$width_mm<-marker_shape$width * marker_shape$mm_per_px
#marker_shape$pct_area<-marker_shape$area/mean(marker_shape$area)
marker_shape$rel_area<-marker_shape$area * marker_shape$mm_per_px

mean(marker_shape$ratio)
sd(marker_shape$ratio)
mean(marker_shape$eccentricity)
sd(marker_shape$eccentricity)

##########################################################################################
## Analysis of tuber  shape on black background
##########################################################################################

tuber_shape<-tuber_size


## Calculate L/W ratio as measured by calipers
tuber_shape$caliper_ratio<-tuber_shape$caliper_length/tuber_shape$caliper_width

## Assess correlation between different measurement types
cor(tuber_shape$sva, tuber_shape$caliper_ratio, use="complete.obs")

cor(tuber_shape$ratio, tuber_shape$caliper_ratio, use="complete.obs")

cor(tuber_shape$ratio, tuber_shape$eccentricity, use="complete.obs")

## Make plots of correlation

## Figure 2
p<-ggplot(tuber_shape, aes(x=caliper_ratio, y=ratio)) + geom_point(size=0.5, color=c("grey10")) + theme_bw() + xlab("L/W ratio (caliper)") + ylab("L/W ratio (machine vision)")  + theme(text = element_text(size=15),legend.position = "none") + geom_smooth(method = "lm", se = FALSE, color="black",linetype="dashed") 

fig.2a_path<-paste(fig.filepath, "/Fig_2a.pdf", sep="")

pdf(fig.2a_path, height=4, width=4)
print(p)
dev.off()

p<-ggplot(tuber_shape, aes(x=sva, y=caliper_ratio)) + geom_point(size=0.5, color=c("grey30")) + theme_bw() + xlab("SVA") + ylab("L/W ratio (caliper)")  + theme(text = element_text(size=15),legend.position = "none") + geom_smooth(method = "lm", se = FALSE, color="black",linetype="dashed") 
fig.2b_path<-paste(fig.filepath, "/Fig_2b.pdf", sep="")

pdf(fig.2b_path, height=4, width=4)
print(p)
dev.off()


p<-ggplot(tuber_shape, aes(x=caliper_ratio, y=eccentricity)) + geom_point(size=0.5, color=c("grey50")) + theme_bw() + xlab("L/W ratio (caliper)") + ylab("Eccentricity")  + theme(text = element_text(size=15),legend.position = "none") + geom_smooth(method = "lm", se = FALSE, color="black",linetype="dashed") 

fig.2c_path<-paste(fig.filepath, "/Fig_2c.pdf", sep="")

pdf(fig.2c_path, height=4, width=4)
print(p)
dev.off()


## Plot L/W ratio by clone
tuber_shape<-tuber_shape[complete.cases(tuber_shape),]
tuber_shape$clone = reorder(tuber_shape$clone, tuber_shape$caliper_ratio, median)
p<-ggplot(tuber_shape, aes(x=as.factor(clone), y=caliper_ratio)) + geom_boxplot() + theme_bw()  + theme(axis.text.x = element_text(angle = 90)) + ylab("L/W Ratio (caliper)") + xlab("Clone") + theme(text = element_text(size=15),legend.position = "none",axis.text.x = element_text(angle = 90, hjust = 1, size=6))
p


fig.2d_path<-paste(fig.filepath, "/Fig_2d.pdf", sep="")
pdf(fig.2d_path, height=3, width=14)
print(p)
dev.off()



##########################################################################################
## Calculate tuber biomass profiles
##########################################################################################

shape<-shape[shape$tuber != "marker",]
shape.mdl<-shape


## Lets see how repeatable the measurements of shape are by looking at both tuber sides
shape.all<-shape

shape.cols<-shape.all[,c(21:ncol(shape.all))]
shape.cols_nz<-shape.cols[ , which(apply(shape.cols, 2, var) != 0)]
shape.pca<-prcomp(shape.cols_nz, scale = TRUE, center = TRUE)
## Get % variance explained by first 2 PCs
pct.explained<-summary(shape.pca)$importance[2,1:2] * 100
pct_variance<-shape.pca$sdev^2/sum(shape.pca$sdev^2)
pct_variance<-pct_variance[1:15]
PCs<-c("PC1", "PC2", "PC3", "PC4","PC5","PC6", "PC7", "PC8", "PC9","PC10","PC11", "PC12", "PC13", "PC14","PC15")
skree<-cbind(PCs, pct_variance)
skree<-as.data.frame(skree)
colnames(skree)<-c("PC", "Variance_Explained")
skree$PC<-factor(skree$PC, levels=c("PC1", "PC2", "PC3", "PC4","PC5","PC6", "PC7", "PC8", "PC9","PC10","PC11", "PC12", "PC13", "PC14","PC15"))
skree$Variance_Explained<-as.numeric(as.character(skree$Variance_Explained))
skree$Variance_Explained<-skree$Variance_Explained*100
p<-ggplot(skree, aes(x=PC, y=Variance_Explained)) + geom_point(color="red") + ylab("% of Variance explained") + xlab("PC") + theme_bw()
#pdf("scree_plot_shape_PCA.pdf", height=4, width=6)
print(p)
#dev.off()

PC13<-as.data.frame(shape.pca$x[,1:3])
PC13$clone<-shape.all$clone
PC13$tuber<-shape.all$tuber
PC13$rep<-shape.all$rep
PC13$side<-shape.all$side

## Lets merge PC and biomass profile data
shape_PCA<-merge(shape.all, PC13, by=c('clone', 'tuber', 'rep', 'side'), all=T)
#write.csv(shape_PCA, file="biomass_profile_PCA_data_NM.csv", quote=F, row.names = F)

#p<-ggplot(PC16, aes(x=PC1, y=PC2)) + geom_point() + theme_bw() + theme(legend.position="none")  
p<-ggplot(PC13, aes(x=PC1, y=PC2)) + geom_point() + theme_bw() + theme(legend.position="none")  
q<- p + xlab(paste("PC1 (", pct.explained[1], "%)", sep="")) + ylab(paste("PC2 (", pct.explained[2], "%)", sep=""))
q<-q+ggtitle("PCA of Tuber biomass profile")

## Merge PC values and other measurements of shape
tuber_shape.all<-merge(shape.all[,c(2,3,4,6,13,14)], PC13, by=c('clone', 'tuber', 'rep', 'side'), all=T)

shape_side1<-tuber_shape.all[tuber_shape.all$side == 1,]
shape_side2<-tuber_shape.all[tuber_shape.all$side == 2,]

colnames(shape_side1)[c(5:ncol(shape_side1))]<-paste(colnames(shape_side1)[c(5:ncol(shape_side1))], "1", sep="_")
colnames(shape_side2)[c(5:ncol(shape_side2))]<-paste(colnames(shape_side2)[c(5:ncol(shape_side2))], "2", sep="_")

shape_compareSides<-merge(shape_side1[,-c(4)], shape_side2[,-c(4)], by=c('clone', 'rep', 'tuber'), all=T)
shape_compareSides<-shape_compareSides[complete.cases(shape_compareSides),]

cor(shape_compareSides$ratio_1, shape_compareSides$ratio_2)
cor(shape_compareSides$eccentricity_1, shape_compareSides$eccentricity_2)
cor(shape_compareSides$PC1_1, shape_compareSides$PC1_2)
cor(shape_compareSides$PC2_1, shape_compareSides$PC2_2)
cor(shape_compareSides$PC3_1, shape_compareSides$PC3_2)


## Keep only one side of the tuber
shape<-shape[shape$side == 1,]

shape.cols<-shape[,c(21:ncol(shape))]
shape.cols_nz<-shape.cols[ , which(apply(shape.cols, 2, var) != 0)]
shape.pca<-prcomp(shape.cols_nz, scale = TRUE, center = TRUE)
## Get % variance explained by first 2 PCs
pct.explained<-summary(shape.pca)$importance[2,1:2] * 100


pct_variance<-shape.pca$sdev^2/sum(shape.pca$sdev^2)
pct_variance<-pct_variance[1:15]
PCs<-c("PC1", "PC2", "PC3", "PC4","PC5","PC6", "PC7", "PC8", "PC9","PC10","PC11", "PC12", "PC13", "PC14","PC15")
skree<-cbind(PCs, pct_variance)
skree<-as.data.frame(skree)
colnames(skree)<-c("PC", "Variance_Explained")
skree$PC<-factor(skree$PC, levels=c("PC1", "PC2", "PC3", "PC4","PC5","PC6", "PC7", "PC8", "PC9","PC10","PC11", "PC12", "PC13", "PC14","PC15"))
skree$Variance_Explained<-as.numeric(as.character(skree$Variance_Explained))
skree$Variance_Explained<-skree$Variance_Explained*100
p<-ggplot(skree, aes(x=PC, y=Variance_Explained)) + geom_point(color="red") + ylab("% of Variance explained") + xlab("PC") + theme_bw()
#pdf("scree_plot_shape_PCA.pdf", height=4, width=6)
print(p)
#dev.off()

PC16<-as.data.frame(shape.pca$x[,1:6])
PC16$clone<-shape$clone
PC16$tuber<-shape$tuber
PC16$replicate<-shape$rep

## Merge PC values and other measurements of shape
tuber_shape<-merge(tuber_shape, PC16, by=c('clone', 'tuber', 'replicate'), all=T)

p<-ggplot(PC16, aes(x=PC1, y=PC2)) + geom_point() + theme_bw() + theme(legend.position="none")  
q<- p + xlab(paste("PC1 (", pct.explained[1], "%)", sep="")) + ylab(paste("PC2 (", pct.explained[2], "%)", sep=""))
q<-q+ggtitle("PCA of Tuber biomass profile")


## Get correlation between L/W ratio and PC values
cor(tuber_shape$ratio, tuber_shape$PC1, use="complete.obs")
cor(tuber_shape$sva, tuber_shape$PC1, use="complete.obs")

p<-ggplot(tuber_shape, aes(x=ratio, y=PC1)) + geom_point(size=0.5, color=c("grey10")) + theme_bw() + xlab("L/W ratio (machine vision)") + ylab("PC1shape")  + theme(text = element_text(size=15),legend.position = "none") + geom_smooth(method = "lm", se = FALSE, color="black",linetype="dashed") 

fig.3a_path<-paste(fig.filepath, "/Fig_3a.pdf", sep="")

pdf(fig.3a_path, height=4, width=4)
print(p)
dev.off()


cor(tuber_shape$ratio, tuber_shape$PC2, use="complete.obs")

p<-ggplot(tuber_shape, aes(x=ratio, y=PC2)) + geom_point(size=0.5, color=c("grey30")) + theme_bw() + xlab("L/W ratio (machine vision)") + ylab("PC2shape")  + theme(text = element_text(size=15),legend.position = "none") + geom_smooth(method = "lm", se = FALSE, color="black",linetype="dashed") 

fig.3b_path<-paste(fig.filepath, "/Fig_3b.pdf", sep="")

pdf(fig.3b_path, height=4, width=4)
print(p)
dev.off()

cor(tuber_shape$ratio, tuber_shape$PC3, use="complete.obs")

p<-ggplot(tuber_shape, aes(x=ratio, y=PC3)) + geom_point(size=0.5, color=c("grey50")) + theme_bw() + xlab("L/W ratio (machine vision)") + ylab("PC3shape")  + theme(text = element_text(size=15),legend.position = "none") + geom_smooth(method = "lm", se = FALSE, color="black",linetype="dashed") 


fig.3c_path<-paste(fig.filepath, "/Fig_3c.pdf", sep="")

pdf(fig.3c_path, height=4, width=4)
print(p)
dev.off()

cor(tuber_shape$ratio, tuber_shape$PC4, use="complete.obs")

cor(tuber_shape$ratio, tuber_shape$PC5, use="complete.obs")

cor(tuber_shape$ratio, tuber_shape$PC6, use="complete.obs")


## Lets plot boxplot of PC1shape
tuber_shape<-tuber_shape[complete.cases(tuber_shape),]
tuber_shape$clone = reorder(tuber_shape$clone, tuber_shape$PC1, median)
p<-ggplot(tuber_shape, aes(x=as.factor(clone), y=PC1)) + geom_boxplot() + theme_bw()  + theme(axis.text.x = element_text(angle = 90)) + ylab("PC1shape") + xlab("Clone") + theme(text = element_text(size=15),legend.position = "none",axis.text.x = element_text(angle = 90, hjust = 1, size=6))
p


fig.3d_path<-paste(fig.filepath, "/Fig_3d.pdf", sep="")
pdf(fig.3d_path, height=3, width=14)
print(p)
dev.off()


tuber_shape<-tuber_shape[complete.cases(tuber_shape),]
tuber_shape$clone = reorder(tuber_shape$clone, tuber_shape$PC2, median)
p<-ggplot(tuber_shape, aes(x=as.factor(clone), y=PC2)) + geom_boxplot() + theme_bw()  + theme(axis.text.x = element_text(angle = 90)) + ylab("PC2shape") + xlab("Clone") + theme(text = element_text(size=15),legend.position = "none",axis.text.x = element_text(angle = 90, hjust = 1, size=6))
p


## Lets plot PC1 vs PC2
p<-ggplot(data=tuber_shape, aes(x=PC1, y=PC2)) + geom_point() + theme_bw() + ylim(c(-15, 25)) 


## Get average of PCs on a per clone basis
tuber_shape.ave<-aggregate(.~ clone, data=tuber_shape[,-c(2,3)], mean)
colnames(tuber_shape.ave)[2:ncol(tuber_shape.ave)]<-paste(colnames(tuber_shape.ave)[2:ncol(tuber_shape.ave)], ".ave", sep="")

tuber_shape.sd<-aggregate(.~ clone, data=tuber_shape[-c(2,3)], sd)
colnames(tuber_shape.sd)[2:ncol(tuber_shape.sd)]<-paste(colnames(tuber_shape.sd)[2:ncol(tuber_shape.sd)], ".sd", sep="")


tuber_shape.par<-merge(tuber_shape.ave, tuber_shape.sd, by=c("clone"))
tuber_shape.par<-tuber_shape.par[complete.cases(tuber_shape.par),]

cor(tuber_shape.par$weight.ave, tuber_shape.par$weight.sd)

clone.ratio<-mean(tuber_shape.ave$ratio)
clone.ratio.sd<-mean(tuber_shape.sd$ratio)

clone.ratio.max<-max(tuber_shape.ave$ratio)
tuber_shape.ave[tuber_shape.ave$ratio == clone.ratio.max, "clone"]
clone.ratio.min<-min(tuber_shape.ave$ratio)
tuber_shape.ave[tuber_shape.ave$ratio == clone.ratio.min, "clone"]

clone.ratio.sd/clone.ratio

cor(tuber_shape.par$ratio.ave, tuber_shape.par$ratio.sd)


clone.PC1<-mean(tuber_shape.ave$PC1, na.rm=TRUE)
clone.PC1.sd<-mean(tuber_shape.sd$PC1, na.rm=TRUE)

clone.PC1.max<-max(tuber_shape.ave$PC1, na.rm=TRUE)
clone.PC1.min<-min(tuber_shape.ave$PC1, na.rm=TRUE)

clone.PC1.sd/clone.PC1

cor(tuber_shape.sd$PC1, tuber_shape.ave$PC1, use="complete.obs")




clone.PC2<-mean(tuber_shape.ave$PC2, na.rm=TRUE)
clone.PC2.sd<-mean(tuber_shape.sd$PC2, na.rm=TRUE)

clone.PC2.max<-max(tuber_shape.ave$PC2, na.rm=TRUE)
clone.PC2.min<-min(tuber_shape.ave$PC2, na.rm=TRUE)

clone.PC2.sd/clone.PC2

cor(tuber_shape.sd$PC2, tuber_shape.ave$PC2, use="complete.obs")


tuber_shape_for_h2<-tuber_shape[,c(1:2,8,9,16:ncol(tuber_shape))]

h2_tuber.shape<-get_h2(tuber_shape_for_h2)

## Lets do the same for tuber shape variance
tuber_sd_for_h2<-tuber_shape[,c(1:3,8,9,16:ncol(tuber_shape))]
tuber_sd_for_h2<-aggregate(.~ clone + replicate, data=tuber_sd_for_h2, sd)
tuber_sd_for_h2<-tuber_sd_for_h2[,-c(3)]

h2_tuber.shape.sd<-get_h2(tuber_sd_for_h2)

## Lets reformat for H2 estimates
PC16<-PC16[,c(7,8,9,1:5)]


PC16.ave<-tuber_shape.ave[,c(1,23:ncol(tuber_shape.ave))]
colnames(PC16.ave)[2:ncol(PC16.ave)]<-gsub(".ave", "", colnames(PC16.ave)[2:ncol(PC16.ave)])

id.vars<-colnames(PC16.ave)[c(1)]
measure.vars<-colnames(PC16.ave)[c(2:(ncol(PC16.ave)))]
PC16.ave.long<-melt(PC16.ave,
                    # ID variables - all the variables to keep but not split apart on
                    id.vars=id.vars,
                    # The source columns
                    measure.vars=measure.vars,
                    # Name of the destination column that will identify the original
                    # column that the measurement came from
                    variable.name="Trait",
                    value.name="Value"
)


PC16.sd<-tuber_shape.sd[,c(1,23:ncol(tuber_shape.sd))]

id.vars<-colnames(PC16.sd)[1]
measure.vars<-colnames(PC16.sd)[c(2:(ncol(PC16.sd)))]
PC16.sd.long<-melt(PC16.sd,
                       # ID variables - all the variables to keep but not split apart on
                       id.vars=id.vars,
                       # The source columns
                       measure.vars=measure.vars,
                       # Name of the destination column that will identify the original
                       # column that the measurement came from
                       variable.name="Trait",
                       value.name="Value"
)



PC16.ag<-merge(PC16.ave, PC16.sd, by=c("clone"))

p<-ggplot(PC16.ag, aes(x=PC1, y=PC2)) + geom_errorbarh(aes(xmin=PC1-PC1.sd, xmax=PC1+PC1.sd)) + geom_errorbar(aes(ymin=PC2-PC2.sd, ymax=PC2+PC2.sd)) + theme_bw() + geom_point(size=1.5, shape=19, color="tan4") + xlab("PC1") + ylab("PC2") + ggtitle("Tuber biomass profile (genotype mean)")  
tempPC1<-PC16.ag[PC16.ag$clone %in% c('2', '189', '170'),]
tempPC1$clone<-factor(tempPC1$clone, levels=c('2', '189', '170'))

q<-p + geom_point(data = tempPC1, aes(x=PC1, y=PC2, fill=as.factor(clone)), shape=21, size = 6) + scale_fill_manual(values=c("firebrick4", "red", "pink"),guide="none")
x<-q + geom_text(data = tempPC1, aes(x=PC1-1, y=PC2-1, label = as.character(clone), color=as.factor(clone)), size = 4) + scale_color_manual(values=c("firebrick4", "red", "pink"),guide="none")
skree$Variance_Explained<-round(skree$Variance_Explained, 1)
z<-x  + xlab(paste("PC1 (", skree[1,2], "%)", sep="")) + ylab(paste("PC2 (", skree[2,2], "%)", sep=""))

#x<-q + geom_text(data = tempPC1, aes(x=PC1-1, y=PC2-1, label = as.character(clone)), color="red3", size = 4)

fig.4_path<-paste(fig.filepath, "/Fig_4.pdf", sep="")
pdf(fig.4_path, height=6, width=6)
print(z)
dev.off()


## This can be a supplmentary figure if needed
tempPC2<-PC16.ag[PC16.ag$clone %in% c('72', '146', '174'),]
p<-x + geom_point(data = tempPC2, aes(x=PC1, y=PC2), color="blue2", size = 6) 
y<-p + geom_text(data = tempPC2, aes(x=PC1-1, y=PC2-1, label = as.character(clone)), color="blue2", size = 4)


#q<-p + geom_point(data = tempPC1, aes(x=PC1, y=PC2), color="red3", size = 6) 
#x<-q + geom_text(data = tempPC1, aes(x=PC1-1, y=PC2-1, label = as.character(clone)), color="red3", size = 4)
#tempPC2<-PC16.ag[PC16.ag$clone %in% c('72', '146', '174'),]
#p<-x + geom_point(data = tempPC2, aes(x=PC1, y=PC2), color="blue2", size = 6) 
#y<-p + geom_text(data = tempPC2, aes(x=PC1-1, y=PC2-1, label = as.character(clone)), color="blue2", size = 4)




#fig.4_path<-paste(fig.filepath, "/Fig_4.pdf", sep="")
#pdf(fig.4_path, height=6, width=6)
#print(p)
#dev.off()
#pdf("PCA_genotype_means.pdf", height=6, width=6)
print(y)
#dev.off()



##### Lets assess fx 
shape.mdl.cols<-shape.mdl[,c(21:ncol(shape.mdl))]
shape.mdl.cols_nz<-shape.mdl.cols[ , which(apply(shape.mdl.cols, 2, var) != 0)]
shape.mdl.pca<-prcomp(shape.mdl.cols_nz, scale = TRUE, center = TRUE)

shape.mdl.PC15<-as.data.frame(shape.mdl.pca$x[,1:5])
shape.mdl.PC15$clone<-shape.mdl$clone
shape.mdl.PC15$rep<-shape.mdl$rep
shape.mdl.PC15$side<-shape.mdl$side
shape.mdl.PC15$tuber<-shape.mdl$tuber

## Merge PC values and other measurements of shape

tuber_shape.mdl<-merge(shape.mdl[,c(2:4,6,13,14)], shape.mdl.PC15, by=c('clone', 'rep', 'side', 'tuber'), all=T)

tuber_shape.model<-tuber_shape.mdl[,c(1:3,5:ncol(tuber_shape.mdl))]
shape.variance_output<-get_fx(tuber_shape.model)


## Lets look at standard error for different sample sizes

## First step is to identify samples with 10 tubers

complete.clones<-names(table(tuber_shape$clone)[table(tuber_shape$clone) > 9])
tuber_shape_complete<-tuber_shape[tuber_shape$clone %in% complete.clones,]

tuber_shape_complete.z<-apply(tuber_shape_complete[,c(4:30)], 2, scale)
tuber_shape_complete.z<-cbind(tuber_shape_complete[,c(1:2)], tuber_shape_complete.z)

std_error_shape<-c()
for(i in 3:10){
  clones<-unique(tuber_shape_complete$clone)
  for (c in clones){
    clone.dat<-tuber_shape_complete.z[tuber_shape_complete.z$clone == c,]
    tubers<-sample(c(1:10), i, replace=F)
    small.dat<-clone.dat[clone.dat$tuber %in% tubers, ]
    ratio.sd.err<-std.error(small.dat$ratio)
    eccentricity.sd.err<-std.error(small.dat$eccentricity)
    caliper_ratio.sd.err<-std.error(small.dat$caliper_ratio)
    sva.sd.err<-std.error(small.dat$sva)
    PC1.sd.err<-std.error(small.dat$PC1)
    PC2.sd.err<-std.error(small.dat$PC2)
    PC3.sd.err<-std.error(small.dat$PC3)
    temp<-c(i,c,ratio.sd.err,eccentricity.sd.err,caliper_ratio.sd.err, sva.sd.err, PC1.sd.err,PC2.sd.err,PC3.sd.err)
    std_error_shape<-rbind(std_error_shape, temp)
  }
}

std_error_shape<-as.data.frame(std_error_shape)
colnames(std_error_shape)<-c("Replicates","Clone", "L/W ratio (MV)", "Eccentricity","L/W ratio (caliper)","SVA","PC1.shape","PC2.shape","PC3.shape")

for(i in 1:ncol(std_error_shape)){
  std_error_shape[,i]<-as.numeric(as.character(std_error_shape[,i]))
}

std_error_shape.ag<-aggregate(.~Replicates, data=std_error_shape[,-c(2)], mean)

id.vars<-colnames(std_error_shape.ag)[c(1)]
measure.vars<-colnames(std_error_shape.ag)[c(2:(ncol(std_error_shape.ag)))]
std_error_shape.ag.long<-melt(std_error_shape.ag,
                             # ID variables - all the variables to keep but not split apart on
                             id.vars=id.vars,
                             # The source columns
                             measure.vars=measure.vars,
                             # Name of the destination column that will identify the original
                             # column that the measurement came from
                             variable.name="Trait",
                             value.name="SE"
)


p<-ggplot(std_error_shape.ag.long, aes(x=Replicates, y=SE, color=Trait)) + geom_line() + theme_bw()
q<-p+geom_point() + xlab("Number of tubers") + ylab("Standard Error")



## Lets look at the eigenvectors from the Biomass profiles
## Load in ground truth data
BPE.filepath<-paste(base.dir, "/BiomassProfileEigenVectors.csv", sep="")
biomass_profile.eigenvector<-read.csv(BPE.filepath)


p<-ggplot(biomass_profile.eigenvector, aes(x=Sweep, y=Value, color=as.factor(SD))) + geom_line() + xlab("Sweep across minor axis") + ylab("Proportion tuber pixels (%)") + scale_color_grey() + theme_bw() + theme(text = element_text(size=15),legend.position = "none") + coord_flip() + scale_x_reverse() + facet_wrap(~PC)

fig.4a_path<-paste(fig.filepath, "/Fig_4a.pdf", sep="")
pdf("Figure_4a.pdf", height=5, width=10)
print(p)
dev.off()


##########################################################################################
## Analysis of color checker standards
##########################################################################################

## First read in data from raw images
scan.cc<-read.csv("color_checker_data_scanner.csv")
blk.cc_NC<-read.csv("color_checker_data_std_not_corrected.csv")
box.cc_NC<-read.csv("color_checker_data_box_not_corrected.csv")

## Lets read in color corrected images
blk.cc<-read.csv("color_checker_data_std.csv")
box.cc<-read.csv("color_checker_data_box.csv")



##########
## Start with scanner images
##########

## Remove image 183_1.jpg (color card not in correct format)
scan.cc<-scan.cc[scan.cc$img_name != '183_1.jpg',]


id.vars<-colnames(scan.cc)[c(1)]
measure.vars<-colnames(scan.cc)[c(2:(ncol(scan.cc)))]

scan.cc_long<-melt(scan.cc,
                # ID variables - all the variables to keep but not split apart on
                id.vars=id.vars,
                # The source columns
                measure.vars=measure.vars,
                # Name of the destination column that will identify the original
                # column that the measurement came from
                variable.name="Trait",
                value.name="measurement"
)


scan.cc_plot<-c()
for(i in 2:25){
  r<-i
  g<-i + 24
  b<-i + 48
  ch<-i-1
  temp<-scan.cc[,c(1,r,g,b)]
  temp$chip<-rep(ch, nrow(temp))
  colnames(temp)[2:4]<-c("Red", "Green", "Blue")
  scan.cc_plot<-rbind(scan.cc_plot, temp)
  
}

## Convert to PC scores
scan.cc_for_pca<-scan.cc_plot[,c(2:4)]

scan.cc_pca_nz<-scan.cc_for_pca[ , which(apply(scan.cc_for_pca, 2, var) != 0)]
scan.cc_pca<-prcomp(scan.cc_pca_nz, scale = TRUE, center = TRUE)

scan_PC<-as.data.frame(scan.cc_pca$x)
scan_PC$img_name<-scan.cc_plot$img_name
scan_PC$chip<-scan.cc_plot$chip

scan.cc_plot<-merge(scan.cc_plot, scan_PC, by=c('img_name', 'chip'), all=T)

scan.cc_plot$Red<-round(scan.cc_plot$Red)
scan.cc_plot$Green<-round(scan.cc_plot$Green)
scan.cc_plot$Blue<-round(scan.cc_plot$Blue)

scan.cc_plot$hex<-rep("NA", nrow(scan.cc_plot))

for(rw in 1:nrow(scan.cc_plot)){
  scan.cc_plot[rw,"hex"]<-rgb(scan.cc_plot[rw, 'Red'],scan.cc_plot[rw, 'Green'],scan.cc_plot[rw, 'Blue'],max=255)
}

cols<-scan.cc_plot$hex
names(cols)<-scan.cc_plot$hex

scan.cc_plot_long<-melt(scan.cc_plot,
                   # ID variables - all the variables to keep but not split apart on
                   id.vars=c("img_name", "chip", "hex"),
                   # The source columns
                   measure.vars=c("Red", "Green", "Blue"),
                   # Name of the destination column that will identify the original
                   # column that the measurement came from
                   variable.name="Color",
                   value.name="Value"
)


## Make a plot of the color checker data
p<-ggplot(scan.cc_plot_long, aes(x=Color, y=Value)) + geom_point(aes(colour = hex)) + theme_bw() + scale_colour_manual(values=cols) + ylim(0,255) + theme(legend.position="none") + facet_wrap(~chip, ncol=6) +ggtitle("Color checker (Scanner)")
#x<-p + stat_summary(fun.y="mean", geom="point", size=6, shape=4, col="red") 
#p +  stat_summary(fun.data=mean_se, geom = "errorbar", width=0.2)
p


fig.S5a_path<-paste(fig.filepath, "/Fig_S5a.pdf", sep="")
pdf(fig.S5a_path, height=6, width=8)
print(p)
dev.off()


## Make a plot derived PC values
p<-ggplot(scan.cc_plot, aes(x=PC1, y=PC2)) + geom_point(aes(colour = hex)) + theme_bw() + scale_colour_manual(values=cols) + theme(legend.position="none") + ggtitle("Color checker PC (Scanner)")

p

fig.S5b_path<-paste(fig.filepath, "/Fig_S5b.pdf", sep="")
pdf(fig.S5b_path, height=6, width=8)
print(p)
dev.off()

## Lets get summary of variation among each PC for each chip
scan_chip.sd<-aggregate(.~ chip, data=scan.cc_plot[,c(2,6:8)], sd)


##########
## Now work on the black background 
##########

##########
## Start with images that haven't been color corrected
##########

## Remove problematic images
blk.cc_NC<-blk.cc_NC[!grepl("117_", blk.cc_NC$img_name),]
blk.cc_NC<-blk.cc_NC[!grepl("177_2", blk.cc_NC$img_name),]


## Lets covert for a more user friendly form
blk.cc_NC_plot<-c()
for(i in 2:25){
  r<-i
  g<-i + 24
  b<-i + 48
  ch<-i-1
  temp<-blk.cc_NC[,c(1,r,g,b)]
  temp$chip<-rep(ch, nrow(temp))
  colnames(temp)[2:4]<-c("Red", "Green", "Blue")
  blk.cc_NC_plot<-rbind(blk.cc_NC_plot, temp)
  
}


## Convert to long form
id.vars<-colnames(blk.cc_NC)[c(1)]
measure.vars<-colnames(blk.cc_NC)[c(2:(ncol(blk.cc_NC)))]

blk.cc_NC_long<-melt(blk.cc_NC,
                  # ID variables - all the variables to keep but not split apart on
                  id.vars=id.vars,
                  # The source columns
                  measure.vars=measure.vars,
                  # Name of the destination column that will identify the original
                  # column that the measurement came from
                  variable.name="Trait",
                  value.name="measurement"
)


## Lets change the chip numbers to reflect the orientation on the scanner
chips<-c(1:24)
replacements<-c(6,12,18,24,5,11,17,23,4,10,16,22,3,9,15,21,2,8,14,20,1,7,13,19)
output<-c()
for(c in 1:length(chips)){
  new_val<-replacements[c]
  temp<-blk.cc_NC_plot[blk.cc_NC_plot$chip == c,]
  temp$chip<-rep(new_val, nrow(temp))
  output<-rbind(output, temp)
}

blk.cc_NC_plot<-output

blk.cc_NC_plot$Red<-round(blk.cc_NC_plot$Red)
blk.cc_NC_plot$Green<-round(blk.cc_NC_plot$Green)
blk.cc_NC_plot$Blue<-round(blk.cc_NC_plot$Blue)

blk.cc_NC_plot$hex<-rep("NA", nrow(blk.cc_NC_plot))

for(rw in 1:nrow(blk.cc_NC_plot)){
  blk.cc_NC_plot[rw,"hex"]<-rgb(blk.cc_NC_plot[rw, 'Red'],blk.cc_NC_plot[rw, 'Green'],blk.cc_NC_plot[rw, 'Blue'],max=255)
}

cols<-blk.cc_NC_plot$hex
names(cols)<-blk.cc_NC_plot$hex

blk.cc_NC_plot_long<-melt(blk.cc_NC_plot,
                       # ID variables - all the variables to keep but not split apart on
                       id.vars=c("img_name", "chip", "hex"),
                       # The source columns
                       measure.vars=c("Red", "Green", "Blue" ),
                       # Name of the destination column that will identify the original
                       # column that the measurement came from
                       variable.name="Color",
                       value.name="Value"
)



p<-ggplot(blk.cc_NC_plot_long, aes(x=Color, y=Value)) + geom_point(aes(colour = hex)) + theme_bw() + scale_colour_manual(values=cols) + ylim(0,255) + theme(legend.position="none") + facet_wrap(~chip, ncol=6) + ggtitle("Color checker (Black background; not color corrected)")
#x<-p + stat_summary(fun.y="mean", geom="point", size=6, shape=4, col="red") 
#p +  stat_summary(fun.data=mean_se, geom = "errorbar", width=0.2)
p

fig.S6a_path<-paste(fig.filepath, "/Fig_S6a.pdf", sep="")
pdf(fig.S6a_path, height=6, width=8)
print(p)
dev.off()




## Lets convert to PC values

## Convert to PC scores
blk.cc_NC_for_pca<-blk.cc_NC_plot[,c(2:4)]

blk.cc_NC_pca_nz<-blk.cc_NC_for_pca[ , which(apply(blk.cc_NC_for_pca, 2, var) != 0)]
blk.cc_NC_pca<-prcomp(blk.cc_NC_pca_nz, scale = TRUE, center = TRUE)


blk_PC<-as.data.frame(blk.cc_NC_pca$x)
blk_PC$img_name<-blk.cc_NC_plot$img_name
blk_PC$chip<-blk.cc_NC_plot$chip


blk.cc_NC_plot.pc<-merge(blk.cc_NC_plot, blk_PC, by=c('img_name', 'chip'), all=T)



blk.cc_NC_plot.pc$Red<-round(blk.cc_NC_plot.pc$Red)
blk.cc_NC_plot.pc$Green<-round(blk.cc_NC_plot.pc$Green)
blk.cc_NC_plot.pc$Blue<-round(blk.cc_NC_plot.pc$Blue)

blk.cc_NC_plot.pc$hex<-rep("NA", nrow(blk.cc_NC_plot.pc))

for(rw in 1:nrow(blk.cc_NC_plot.pc)){
  blk.cc_NC_plot.pc[rw,"hex"]<-rgb(blk.cc_NC_plot.pc[rw, 'Red'],blk.cc_NC_plot.pc[rw, 'Green'],blk.cc_NC_plot.pc[rw, 'Blue'],max=255)
}

cols<-blk.cc_NC_plot.pc$hex
names(cols)<-blk.cc_NC_plot.pc$hex


## Plot derived PC values
p<-ggplot(blk.cc_NC_plot.pc, aes(x=PC1, y=PC2)) + geom_point(aes(colour = hex)) + theme_bw() + scale_colour_manual(values=cols) + theme(legend.position="none") + ggtitle("Color checker PC (Black background; not color corrected)")

fig.S6b_path<-paste(fig.filepath, "/Fig_S6b.pdf", sep="")
pdf(fig.S6b_path, height=6, width=8)
print(p)
dev.off()




## Lets get summary of variation among each PC for each chip
blk_NC_chip.sd<-aggregate(.~ chip, data=blk.cc_NC_plot.pc[,c(2,7:9)], sd)



##########
## Now lets work on images that are color corrected
##########


## Remove problematic images
blk.cc<-blk.cc[!grepl("117_", blk.cc$img_name),]

## Lets covert for a more user friendly form
blk.cc_plot<-c()
for(i in 2:25){
  r<-i
  g<-i + 24
  b<-i + 48
  ch<-i-1
  temp<-blk.cc[,c(1,r,g,b)]
  temp$chip<-rep(ch, nrow(temp))
  colnames(temp)[2:4]<-c("Red", "Green", "Blue")
  blk.cc_plot<-rbind(blk.cc_plot, temp)
  
}


## Convert to long form
id.vars<-colnames(blk.cc)[c(1)]
measure.vars<-colnames(blk.cc)[c(2:(ncol(blk.cc)))]

blk.cc_long<-melt(blk.cc,
                   # ID variables - all the variables to keep but not split apart on
                   id.vars=id.vars,
                   # The source columns
                   measure.vars=measure.vars,
                   # Name of the destination column that will identify the original
                   # column that the measurement came from
                   variable.name="Trait",
                   value.name="measurement"
)


## Lets change the chip numbers to reflect the orientation on the scanner
chips<-c(1:24)
#replacements<-c(9,12,18,24,5,11,17,23,4,10,16,22,3,9,15,21,2,8,14,20,1,7,13,19)
replacements<-c(19,13,8,1,20,14,7,2,21,15,9,3,22,16,10,4,23,17,11,5,24,18,12,6)
output<-c()
for(c in 1:length(chips)){
  new_val<-replacements[c]
  temp<-blk.cc_plot[blk.cc_plot$chip == c,]
  temp$chip<-rep(new_val, nrow(temp))
  output<-rbind(output, temp)
}

blk.cc_plot<-output

blk.cc_plot$Red<-round(blk.cc_plot$Red)
blk.cc_plot$Green<-round(blk.cc_plot$Green)
blk.cc_plot$Blue<-round(blk.cc_plot$Blue)

blk.cc_plot$hex<-rep("NA", nrow(blk.cc_plot))

for(rw in 1:nrow(blk.cc_plot)){
  blk.cc_plot[rw,"hex"]<-rgb(blk.cc_plot[rw, 'Red'],blk.cc_plot[rw, 'Green'],blk.cc_plot[rw, 'Blue'],max=255)
}

cols<-blk.cc_plot$hex
names(cols)<-blk.cc_plot$hex

blk.cc_plot_long<-melt(blk.cc_plot,
                        # ID variables - all the variables to keep but not split apart on
                        id.vars=c("img_name", "chip", "hex"),
                        # The source columns
                        measure.vars=c("Red", "Green", "Blue" ),
                        # Name of the destination column that will identify the original
                        # column that the measurement came from
                        variable.name="Color",
                        value.name="Value"
)



p<-ggplot(blk.cc_plot_long, aes(x=Color, y=Value)) + geom_point(aes(colour = hex)) + theme_bw() + scale_colour_manual(values=cols) + ylim(0,255) + theme(legend.position="none") + facet_wrap(~chip, ncol=6) + ggtitle("Color checker (Black background; color corrected)")
#x<-p + stat_summary(fun.y="mean", geom="point", size=6, shape=4, col="red") 
#p +  stat_summary(fun.data=mean_se, geom = "errorbar", width=0.2)
p

fig.S6c_path<-paste(fig.filepath, "/Fig_S6c.pdf", sep="")
pdf(fig.S6c_path, height=6, width=8)
print(p)
dev.off()



## Remove problematic images
#blk.cc_plot_long[(blk.cc_plot_long$chip == 1) & (blk.cc_plot_long$Value > 100),]
#blk.cc_plot_long[(blk.cc_plot_long$chip == 6) & (blk.cc_plot_long$Value < 100),]
#blk.cc_plot_long[(blk.cc_plot_long$chip == 24) & (blk.cc_plot_long$Value > 100),]
#blk.cc_plot_long[(blk.cc_plot_long$chip == 15) & (blk.cc_plot_long$Value > 100),]

#blk.cc_plot_long<-blk.cc_plot_long[!grepl("117_",blk.cc_plot_long$img_name),]
#blk.cc_plot_long<-blk.cc_plot_long[!grepl("177_",blk.cc_plot_long$img_name),]

## These are now removed above

## Lets convert to PC values

## Convert to PC scores
blk.cc_for_pca<-blk.cc_plot[,c(2:4)]

blk.cc_pca_nz<-blk.cc_for_pca[ , which(apply(blk.cc_for_pca, 2, var) != 0)]
blk.cc_pca<-prcomp(blk.cc_pca_nz, scale = TRUE, center = TRUE)


blk_PC<-as.data.frame(blk.cc_pca$x)
blk_PC$img_name<-blk.cc_plot$img_name
blk_PC$chip<-blk.cc_plot$chip


blk.cc_plot.pc<-merge(blk.cc_plot, blk_PC, by=c('img_name', 'chip'), all=T)



blk.cc_plot.pc$Red<-round(blk.cc_plot.pc$Red)
blk.cc_plot.pc$Green<-round(blk.cc_plot.pc$Green)
blk.cc_plot.pc$Blue<-round(blk.cc_plot.pc$Blue)

blk.cc_plot.pc$hex<-rep("NA", nrow(blk.cc_plot.pc))

for(rw in 1:nrow(blk.cc_plot.pc)){
  blk.cc_plot.pc[rw,"hex"]<-rgb(blk.cc_plot.pc[rw, 'Red'],blk.cc_plot.pc[rw, 'Green'],blk.cc_plot.pc[rw, 'Blue'],max=255)
}

cols<-blk.cc_plot.pc$hex
names(cols)<-blk.cc_plot.pc$hex


## Plot derived PC values
p<-ggplot(blk.cc_plot.pc, aes(x=PC1, y=PC2)) + geom_point(aes(colour = hex)) + theme_bw() + scale_colour_manual(values=cols) + theme(legend.position="none") + ggtitle("Color checker PC (Black background; color corrected)")

fig.S6d_path<-paste(fig.filepath, "/Fig_S6d.pdf", sep="")
pdf(fig.S6d_path, height=6, width=8)
print(p)
dev.off()




## Lets get summary of variation among each PC for each chip
blk_chip.sd<-aggregate(.~ chip, data=blk.cc_plot.pc[,c(2,7:9)], sd)


## Lets find most variable chips on scanner 
scan_chip.sd[order(scan_chip.sd$PC1),]
## Now on black background
blk_chip.sd[order(blk_chip.sd$PC1),]

## Now merge PC1 values to figure out how much
colnames(scan_chip.sd)[2:4]<-paste(colnames(scan_chip.sd)[2:4], ".scan", sep="")
colnames(blk_chip.sd)[2:4]<-paste(colnames(blk_chip.sd)[2:4], ".blk", sep="")
PC_chip.sd<-merge(scan_chip.sd, blk_chip.sd, by=c("chip"))

## Examine relative difference between PC1 on black background vs PC1 on the scanner
PC_chip.sd$PC1.blk/PC_chip.sd$PC1.scan

## Get the average relative difference between PC1 on black background vs PC1 on the scanner
mean(PC_chip.sd$PC1.blk/PC_chip.sd$PC1.scan)


##########
## Now work on the light box color cards 
##########


## Remove problematic images
box.cc_NC<-box.cc_NC[!grepl('145_1_2_box', box.cc_NC$img_name),]

## Lets covert for a more user friendly form
box.cc_NC_plot<-c()
for(i in 2:25){
  r<-i
  g<-i + 24
  b<-i + 48
  ch<-i-1
  temp<-box.cc_NC[,c(1,r,g,b)]
  temp$chip<-rep(ch, nrow(temp))
  colnames(temp)[2:4]<-c("Red", "Green", "Blue")
  box.cc_NC_plot<-rbind(box.cc_NC_plot, temp)
  
}


## Convert to long form
id.vars<-colnames(box.cc_NC)[c(1)]
measure.vars<-colnames(box.cc_NC)[c(2:(ncol(box.cc_NC)))]

box.cc_NC_long<-melt(box.cc_NC,
                  # ID variables - all the variables to keep but not split apart on
                  id.vars=id.vars,
                  # The source columns
                  measure.vars=measure.vars,
                  # Name of the destination column that will identify the original
                  # column that the measurement came from
                  variable.name="Trait",
                  value.name="measurement"
)


## Lets change the chip numbers to reflect the orientation on the scanner
chips<-c(1:24)
replacements<-c(6,12,18,24,5,11,17,23,4,10,16,22,3,9,15,21,2,8,14,20,1,7,13,19)
output<-c()
for(c in 1:length(chips)){
  new_val<-replacements[c]
  temp<-box.cc_NC_plot[box.cc_NC_plot$chip == c,]
  temp$chip<-rep(new_val, nrow(temp))
  output<-rbind(output, temp)
}

box.cc_NC_plot<-output

box.cc_NC_plot$Red<-round(box.cc_NC_plot$Red)
box.cc_NC_plot$Green<-round(box.cc_NC_plot$Green)
box.cc_NC_plot$Blue<-round(box.cc_NC_plot$Blue)

box.cc_NC_plot$hex<-rep("NA", nrow(box.cc_NC_plot))

for(rw in 1:nrow(box.cc_NC_plot)){
  box.cc_NC_plot[rw,"hex"]<-rgb(box.cc_NC_plot[rw, 'Red'],box.cc_NC_plot[rw, 'Green'],box.cc_NC_plot[rw, 'Blue'],max=255)
}

cols<-box.cc_NC_plot$hex
names(cols)<-box.cc_NC_plot$hex

box.cc_NC_plot_long<-melt(box.cc_NC_plot,
                       # ID variables - all the variables to keep but not split apart on
                       id.vars=c("img_name", "chip", "hex"),
                       # The source columns
                       measure.vars=c("Red", "Green", "Blue" ),
                       # Name of the destination column that will identify the original
                       # column that the measurement came from
                       variable.name="Color",
                       value.name="Value"
)



p<-ggplot(box.cc_NC_plot_long, aes(x=Color, y=Value)) + geom_point(aes(colour = hex)) + theme_bw() + scale_colour_manual(values=cols) + ylim(0,255) + theme(legend.position="none") + facet_wrap(~chip, ncol=6) + ggtitle("Color checker (Light box background; not color corrected)")
#x<-p + stat_summary(fun.y="mean", geom="point", size=6, shape=4, col="red") 
#p +  stat_summary(fun.data=mean_se, geom = "errorbar", width=0.2)
p

fig.S7a_path<-paste(fig.filepath, "/Fig_S7a.pdf", sep="")
pdf(fig.S7a_path, height=6, width=8)
print(p)
dev.off()


## Lets convert to PC values

## Convert to PC scores
box.cc_NC_for_pca<-box.cc_NC_plot[,c(2:4)]

box.cc_NC_pca_nz<-box.cc_NC_for_pca[ , which(apply(box.cc_NC_for_pca, 2, var) != 0)]
box.cc_NC_pca<-prcomp(box.cc_NC_pca_nz, scale = TRUE, center = TRUE)


box_PC<-as.data.frame(box.cc_NC_pca$x)
box_PC$img_name<-box.cc_NC_plot$img_name
box_PC$chip<-box.cc_NC_plot$chip


box.cc_NC_plot.pc<-merge(box.cc_NC_plot, box_PC, by=c('img_name', 'chip'), all=T)



box.cc_NC_plot.pc$Red<-round(box.cc_NC_plot.pc$Red)
box.cc_NC_plot.pc$Green<-round(box.cc_NC_plot.pc$Green)
box.cc_NC_plot.pc$Blue<-round(box.cc_NC_plot.pc$Blue)

box.cc_NC_plot.pc$hex<-rep("NA", nrow(box.cc_NC_plot.pc))

for(rw in 1:nrow(box.cc_NC_plot.pc)){
  box.cc_NC_plot.pc[rw,"hex"]<-rgb(box.cc_NC_plot.pc[rw, 'Red'],box.cc_NC_plot.pc[rw, 'Green'],box.cc_NC_plot.pc[rw, 'Blue'],max=255)
}

cols<-box.cc_NC_plot.pc$hex
names(cols)<-box.cc_NC_plot.pc$hex


## Plot derived PC values
p<-ggplot(box.cc_NC_plot.pc, aes(x=PC1, y=PC2)) + geom_point(aes(colour = hex)) + theme_bw() + scale_colour_manual(values=cols) + theme(legend.position="none") + ggtitle("Color checker PC (Light box background; not color corrected)")

fig.S7b_path<-paste(fig.filepath, "/Fig_S7b.pdf", sep="")
pdf(fig.S7b_path, height=6, width=8)
print(p)
dev.off()




## Lets get summary of variation among each PC for each chip
box_NC_chip.sd<-aggregate(.~ chip, data=box.cc_NC_plot.pc[,c(2,7:9)], sd)




## Remove problematic images
box.cc<-box.cc[!grepl('145_1_2_box', box.cc$img_name),]

## Lets covert for a more user friendly form
box.cc_plot<-c()
for(i in 2:25){
  r<-i
  g<-i + 24
  b<-i + 48
  ch<-i-1
  temp<-box.cc[,c(1,r,g,b)]
  temp$chip<-rep(ch, nrow(temp))
  colnames(temp)[2:4]<-c("Red", "Green", "Blue")
  box.cc_plot<-rbind(box.cc_plot, temp)
  
}


## Convert to long form
id.vars<-colnames(box.cc)[c(1)]
measure.vars<-colnames(box.cc)[c(2:(ncol(box.cc)))]

box.cc_long<-melt(box.cc,
                  # ID variables - all the variables to keep but not split apart on
                  id.vars=id.vars,
                  # The source columns
                  measure.vars=measure.vars,
                  # Name of the destination column that will identify the original
                  # column that the measurement came from
                  variable.name="Trait",
                  value.name="measurement"
)


## Lets change the chip numbers to reflect the orientation on the scanner
chips<-c(1:24)
#replacements<-c(9,12,18,24,5,11,17,23,4,10,16,22,3,9,15,21,2,8,14,20,1,7,13,19)
replacements<-c(6,12,18,24,5,11,17,23,4,10,16,22,3,9,15,21,2,8,14,20,1,7,13,19)
output<-c()
for(c in 1:length(chips)){
  new_val<-replacements[c]
  temp<-box.cc_plot[box.cc_plot$chip == c,]
  temp$chip<-rep(new_val, nrow(temp))
  output<-rbind(output, temp)
}

box.cc_plot<-output

box.cc_plot$Red<-round(box.cc_plot$Red)
box.cc_plot$Green<-round(box.cc_plot$Green)
box.cc_plot$Blue<-round(box.cc_plot$Blue)

box.cc_plot$hex<-rep("NA", nrow(box.cc_plot))

for(rw in 1:nrow(box.cc_plot)){
  box.cc_plot[rw,"hex"]<-rgb(box.cc_plot[rw, 'Red'],box.cc_plot[rw, 'Green'],box.cc_plot[rw, 'Blue'],max=255)
}

cols<-box.cc_plot$hex
names(cols)<-box.cc_plot$hex

box.cc_plot_long<-melt(box.cc_plot,
                       # ID variables - all the variables to keep but not split apart on
                       id.vars=c("img_name", "chip", "hex"),
                       # The source columns
                       measure.vars=c("Red", "Green", "Blue" ),
                       # Name of the destination column that will identify the original
                       # column that the measurement came from
                       variable.name="Color",
                       value.name="Value"
)



p<-ggplot(box.cc_plot_long, aes(x=Color, y=Value)) + geom_point(aes(colour = hex)) + theme_bw() + scale_colour_manual(values=cols) + ylim(0,255) + theme(legend.position="none") + facet_wrap(~chip, ncol=6) + ggtitle("Color checker (Light box background; color corrected)")
#x<-p + stat_summary(fun.y="mean", geom="point", size=6, shape=4, col="red") 
#p +  stat_summary(fun.data=mean_se, geom = "errorbar", width=0.2)
p

fig.S7c_path<-paste(fig.filepath, "/Fig_S7c.pdf", sep="")
pdf(fig.S7c_path, height=6, width=8)
print(p)
dev.off()



## Lets convert to PC values

## Convert to PC scores
box.cc_for_pca<-box.cc_plot[,c(2:4)]

box.cc_pca_nz<-box.cc_for_pca[ , which(apply(box.cc_for_pca, 2, var) != 0)]
box.cc_pca<-prcomp(box.cc_pca_nz, scale = TRUE, center = TRUE)


box_PC<-as.data.frame(box.cc_pca$x)
box_PC$img_name<-box.cc_plot$img_name
box_PC$chip<-box.cc_plot$chip


box.cc_plot.pc<-merge(box.cc_plot, box_PC, by=c('img_name', 'chip'), all=T)



box.cc_plot.pc$Red<-round(box.cc_plot.pc$Red)
box.cc_plot.pc$Green<-round(box.cc_plot.pc$Green)
box.cc_plot.pc$Blue<-round(box.cc_plot.pc$Blue)

box.cc_plot.pc$hex<-rep("NA", nrow(box.cc_plot.pc))

for(rw in 1:nrow(box.cc_plot.pc)){
  box.cc_plot.pc[rw,"hex"]<-rgb(box.cc_plot.pc[rw, 'Red'],box.cc_plot.pc[rw, 'Green'],box.cc_plot.pc[rw, 'Blue'],max=255)
}

cols<-box.cc_plot.pc$hex
names(cols)<-box.cc_plot.pc$hex


## Plot derived PC values
p<-ggplot(box.cc_plot.pc, aes(x=PC1, y=PC2)) + geom_point(aes(colour = hex)) + theme_bw() + scale_colour_manual(values=cols) + theme(legend.position="none") + ggtitle("Color checker PC (Light box background; color corrected)")

fig.S7d_path<-paste(fig.filepath, "/Fig_S7d.pdf", sep="")
pdf(fig.S7d_path, height=6, width=8)
print(p)
dev.off()


box_chip.sd<-aggregate(.~ chip, data=box.cc_plot.pc[,c(2,7:9)], sd)


###########
## Lets compare variance (standard deviation between methods)
###########

scan_chip.sd
box_chip.sd
box_NC_chip.sd
blk_chip.sd
blk_NC_chip.sd

## Lets compare non-color corrected images
mean(scan_chip.sd$PC1.scan)
mean(blk_NC_chip.sd$PC1)
mean(box_NC_chip.sd$PC1)

mean(blk_chip.sd$PC1)
mean(box_chip.sd$PC1)

chip<-c(1:24)
scan.order<-rank(scan_chip.sd$PC1.scan)
blk.order<-rank(blk_chip.sd$PC1)
blk_NC.order<-rank(blk_NC_chip.sd$PC1)
box.order<-rank(blk_chip.sd$PC1)
box_NC.order<-rank(box_NC_chip.sd$PC1)

color_order<-cbind(chip, scan.order, blk_NC.order, box_NC.order, blk.order, box.order)
color_order<-as.data.frame(color_order)
color_order$ave_nc<-rowMeans(color_order[,c(2:4)])
color_order[order(color_order$ave_nc, decreasing = T),]
color_order$ave_cc<-rowMeans(color_order[,c(5:6)])
color_order[order(color_order$ave_cc, decreasing = T),]
color_order$ave_all<-rowMeans(color_order[,c(2:6)])
color_order[order(color_order$ave_all, decreasing = T),]

table.dir<-paste(base.dir, "/tables", sep="")
TableS1.name<-paste(table.dir, "/Table_S1.csv", sep="")

write.csv(color_order, file=TableS1.name, quote=F, row.names=T)

color_cor<-cor(color_order[,c(2:6)], use="complete.obs")
corrplot(color_cor, type="upper")

##########################################################################################
## Analysis of tuber flesh
##########################################################################################


#int_color<-scan[,c(2:4,13:ncol(scan))]
int_color<-scan[,c(2:4,13:15)]
int_color<-int_color[int_color$tuber != "marker",]
int_color.cols<-int_color[,c(4:ncol(int_color))]

scan_pca<-int_color.cols[ , which(apply(int_color.cols, 2, var) != 0)]
scan.color.pca<-prcomp(scan_pca, scale = TRUE, center = TRUE)


## Get % variance explained by first 2 PCs
int_pct.explained<-summary(scan.color.pca)$importance[2,1:2] * 100

int_pct_variance<-scan.color.pca$sdev^2/sum(scan.color.pca$sdev^2)
#int_pct_variance<-int_pct_variance[1:6]
int_pct_variance<-int_pct_variance[1:3]
#PCs<-c("PC1", "PC2", "PC3", "PC4","PC5","PC6")
PCs<-c("PC1", "PC2", "PC3")
int_skree<-cbind(PCs, int_pct_variance)
int_skree<-as.data.frame(int_skree)
colnames(int_skree)<-c("PC", "Variance_Explained")
#int_skree$PC<-factor(int_skree$PC, levels=c("PC1", "PC2", "PC3", "PC4","PC5","PC6"))
int_skree$PC<-factor(int_skree$PC, levels=c("PC1", "PC2", "PC3"))
int_skree$Variance_Explained<-as.numeric(as.character(int_skree$Variance_Explained))
int_skree$Variance_Explained<-int_skree$Variance_Explained*100
p<-ggplot(int_skree, aes(x=PC, y=Variance_Explained)) + geom_point(color="red") + ylab("% of Variance explained") + xlab("PC") + theme_bw()

#pdf("skree_plot_flesh_color_PCA.pdf", height=4, width=6)
print(p)
#dev.off()

int_color.PC<-as.data.frame(scan.color.pca$x[,1:3])
int_color.PC$img_name<-int_color$img_name
int_color.PC$clone<-int_color$clone
int_color.PC$rep<-int_color$rep
int_color.PC$tuber<-int_color$tuber


int_color.PC$Red<-int_color$red_ave
int_color.PC$Green<-int_color$green_ave
int_color.PC$Blue<-int_color$blue_ave

int_color.PC$hex<-rep("NA", nrow(int_color.PC))

for(rw in 1:nrow(int_color.PC)){
  int_color.PC[rw,"hex"]<-rgb(int_color.PC[rw, 'Red'],int_color.PC[rw, 'Green'],int_color.PC[rw, 'Blue'],max=255)
}

cols<-int_color.PC$hex
names(cols)<-int_color.PC$hex


## This figure highlights clones with light and dark flesh
int_color_clones<-int_color.PC[int_color.PC$clone %in% c(70, 141, 1, 22, 193, 97, 100,78),]
p<-ggplot() + geom_point(data=int_color_clones, aes(x=PC1, y=PC2, fill=as.factor(clone)), shape=21, size=3) + theme_bw() +  theme(text = element_text(size=15)) +xlim(c(-11, 11)) + ylim(-6,2) + labs(fill = "Clone")
q<-p + geom_point(data = int_color.PC, aes(x=PC1, y=PC2, color = hex), size=0.8) + theme_bw() + scale_color_manual(values=cols, guide="none") #+ theme(legend.position="none", text = element_text(size=15)) 

## This figure highlights clones with light and dark flesh
int_color_clones<-int_color.PC[int_color.PC$clone %in% c(70, 141, 1, 22, 193, 97, 100,78),]
p<-ggplot() + geom_point(data=int_color_clones, aes(x=PC1, y=PC2, fill=as.factor(clone)), shape=21, size=3) + theme_bw() +  theme(text = element_text(size=15)) +xlim(c(-11, 11)) + ylim(-6,2) + labs(fill = "Clone")
q<-p + geom_point(data = int_color.PC, aes(x=PC1, y=PC2, color = hex), size=0.8) + theme_bw() + scale_color_manual(values=cols, guide="none") #+ theme(legend.position="none", text = element_text(size=15)) 


## THis is the figure we want!
int_color_clones2<-int_color.PC[int_color.PC$clone %in% c(70, 141, 1, 193, 97, 100),]
tuber_shape$clone = reorder(tuber_shape$clone, tuber_shape$caliper_ratio, median)
int_color_clones2$clone<-factor(int_color_clones2$clone, levels=c(70, 141, 1, 193, 97, 100))

p<-ggplot(data = int_color.PC, aes(x=PC1, y=PC2, color = hex), size=0.8) + geom_point() + theme_bw() + scale_color_manual(values=cols, guide="none")
q<-p + geom_point(data=int_color_clones2, aes(x=PC1, y=PC2, fill=as.factor(clone)), shape=21, size=3) + scale_fill_manual(values=c("navy", "blue", "deepskyblue", "firebrick4", "red", "pink"),guide="none")
z<-q  + xlab(paste("PC1 (", int_pct.explained[1], "%)", sep="")) + ylab(paste("PC2 (", int_pct.explained[2], "%)", sep=""))
z<-z+ggtitle("Tuber flesh color")



## This figure highlights clones with light and dark flesh
p<-ggplot() + geom_point(data = int_color.PC, aes(x=PC1, y=PC2, color = hex), size=0.8) + theme_bw() + scale_color_manual(values=cols, guide="none") #+ theme(legend.position="none", text = element_text(size=15)) 
q<-p + geom_point(data=int_color_clones, aes(x=PC1, y=PC2, fill=as.factor(clone)), shape=21, size=3) + theme_bw() +  theme(text = element_text(size=15)) + labs(fill = "Clone") # +xlim(c(-11, 11)) + ylim(-6,2) +
z<-q + geom_point(data=int_color_clones, aes(x=PC1, y=PC2, color=hex), size=0.8) + scale_color_manual(values=cols, guide="none")
z<-z  + xlab(paste("PC1 (", int_pct.explained[1], "%)", sep="")) + ylab(paste("PC2 (", int_pct.explained[2], "%)", sep=""))
z<-z+ggtitle("Tuber flesh color")

p<-ggplot(int_color.PC, aes(x=PC1, y=PC2)) + geom_point(aes(colour = hex), size=0.4) + theme_bw() + scale_colour_manual(values=cols) + theme(legend.position="none", text = element_text(size=15))  
#q<- p + ylim(-30, 55) +xlim(-75, 50) + xlab(paste("PC1 (", int_pct.explained[1], "%)", sep="")) + ylab(paste("PC2 (", int_pct.explained[2], "%)", sep=""))
q<- p  + xlab(paste("PC1 (", int_pct.explained[1], "%)", sep="")) + ylab(paste("PC2 (", int_pct.explained[2], "%)", sep=""))
q<-q+ggtitle("Tuber flesh color")

#pdf("PCA_tuber_flesh_color", height=4, width=4)
print(q)
#dev.off()

fig.5b_path<-paste(fig.filepath, "/Fig_5b.pdf", sep="")
pdf(fig.5b_path, height=6, width=8)
#print(q)
print(z)
dev.off()



## Lets calculate mean and variance by clone
int_color.PC<-int_color.PC[,c(4:6,1:3,7:9)]

## Get average of PCs on a per clone basis
int_color.PC.ave<-aggregate(.~ clone, data=int_color.PC[,-c(2,3)], mean)

## Get standard deviation of PCs on a per clone basis
int_color.PC.sd<-aggregate(.~ clone, data=int_color.PC[,-c(2,3,7:9)], sd)
colnames(int_color.PC.sd)[2:ncol(int_color.PC.sd)]<-paste(colnames(int_color.PC.sd)[2:ncol(int_color.PC.sd)], ".sd", sep="")

## Combine both parameters in the same data.frame
int_color.PC.par<-merge(int_color.PC.ave, int_color.PC.sd, by=c("clone"))
int_color.PC.par<-int_color.PC.par[complete.cases(int_color.PC.par),]

## Are they correlated? Yes
cor(int_color.PC.par$PC1, int_color.PC.par$PC1.sd)

## What are the max and min values
int_PC1.max<-max(int_color.PC.par$PC1)
int_PC1.min<-min(int_color.PC.par$PC1)


## Lets make a plot of genotype means
int_color.PC.par$hex<-rep("NA", nrow(int_color.PC.par))
int_color.PC.par$Red<-round(int_color.PC.par$Red)
int_color.PC.par$Green<-round(int_color.PC.par$Green)
int_color.PC.par$Blue<-round(int_color.PC.par$Blue)


for(rw in 1:nrow(int_color.PC.par)){
  int_color.PC.par[rw,"hex"]<-rgb(int_color.PC.par[rw, 'Red'],int_color.PC.par[rw, 'Green'],int_color.PC.par[rw, 'Blue'],max=255)
}

cols<-int_color.PC.par$hex
names(cols)<-int_color.PC.par$hex


#p<-ggplot(int_color.PC.par, aes(x=PC1, y=PC2)) + geom_errorbarh(aes(xmin=PC1-PC1.sd, xmax=PC1+PC1.sd)) + geom_errorbar(aes(ymin=PC2-PC2.sd, ymax=PC2+PC2.sd)) + theme_bw() + geom_point(size=4, shape=19, aes(colour = hex)) + xlab("PC1") + ylab("PC2") + ggtitle("Tuber flesh color (genotype mean)") + xlim(-35, 35) + ylim(-25, 25)  + scale_colour_manual(values=cols) + theme(legend.position="none", text = element_text(size=15))  
p<-ggplot(int_color.PC.par, aes(x=PC1, y=PC2)) + geom_errorbarh(aes(xmin=PC1-PC1.sd, xmax=PC1+PC1.sd)) + geom_errorbar(aes(ymin=PC2-PC2.sd, ymax=PC2+PC2.sd)) + theme_bw() + geom_point(size=4, shape=19, aes(colour = hex)) + xlab("PC1") + ylab("PC2") + ggtitle("Tuber flesh color (genotype mean)")  + scale_colour_manual(values=cols) + theme(legend.position="none", text = element_text(size=15))  


p<-ggplot(int_color.PC.par, aes(x=PC1, y=PC2)) + geom_point(aes(colour = hex), size=4) + theme_bw() + scale_colour_manual(values=cols) + theme(legend.position="none") 

p<-ggplot(int_color.PC.par, aes(x=PC1, y=PC2)) + geom_point() + theme_bw()  + theme(legend.position="none") 

## Calculate H2 for flesh color
int_color.PC_for_h2<-int_color.PC[,c(1:2,4:6)]
h2_int_color<-get_h2(int_color.PC_for_h2)


## Lets do the same for flesh color variance
int_color.sd_for_h2<-int_color.PC[,c(1:2,4:6)]
int_color.sd_for_h2<-aggregate(.~ clone + rep, data=int_color.sd_for_h2, sd)

h2_int_color.sd<-get_h2(int_color.sd_for_h2)



## Lets look at standard error for different sample sizes

## First step is to identify samples with 10 tubers

complete.clones<-names(table(int_color.PC$clone)[table(int_color.PC$clone) > 9])
int_color.PC_complete<-int_color.PC[int_color.PC$clone %in% complete.clones,]

int_color.PC_complete.z<-apply(int_color.PC_complete[,c(4:ncol(int_color.PC_complete))], 2, scale)
int_color.PC_complete.z<-cbind(int_color.PC_complete[,c(1:3)], int_color.PC_complete.z)

std_error_int_color<-c()
for(i in 3:10){
  clones<-unique(int_color.PC_complete.z$clone)
  for (c in clones){
    clone.dat<-int_color.PC_complete.z[int_color.PC_complete.z$clone == c,]
    tubers<-sample(c(1:10), i, replace=F)
    small.dat<-clone.dat[clone.dat$tuber %in% tubers, ]
    PC1.sd.err<-std.error(small.dat$PC1)
    PC2.sd.err<-std.error(small.dat$PC2)
    PC3.sd.err<-std.error(small.dat$PC3)
    temp<-c(i,c,PC1.sd.err,PC2.sd.err,PC3.sd.err)
    std_error_int_color<-rbind(std_error_int_color, temp)
  }
}

std_error_int_color<-as.data.frame(std_error_int_color)
colnames(std_error_int_color)<-c("Replicates","Clone", "PC1.flesh","PC2.flesh","PC3.flesh")

for(i in 1:ncol(std_error_int_color)){
  std_error_int_color[,i]<-as.numeric(as.character(std_error_int_color[,i]))
}

std_error_int_color.ag<-aggregate(.~Replicates, data=std_error_int_color[,-c(2)], mean)

id.vars<-colnames(std_error_int_color.ag)[c(1)]
measure.vars<-colnames(std_error_int_color.ag)[c(2:(ncol(std_error_int_color.ag)))]
std_error_int_color.ag.long<-melt(std_error_int_color.ag,
                              # ID variables - all the variables to keep but not split apart on
                              id.vars=id.vars,
                              # The source columns
                              measure.vars=measure.vars,
                              # Name of the destination column that will identify the original
                              # column that the measurement came from
                              variable.name="Trait",
                              value.name="SE"
)


p<-ggplot(std_error_int_color.ag.long, aes(x=Replicates, y=SE, color=Trait)) + geom_line() + theme_bw()
q<-p+geom_point() + xlab("Number of tubers") + ylab("Standard Error")



#### Co-authors want analysis of how PC1 corresponds with other color space channels
## Table S2

hsv<-convert_colour(int_color.PC[,c(7:9)], from='rgb', to='hsv')
lab<-convert_colour(int_color.PC[,c(7:9)], from='rgb', to='lab')
xyz<-convert_colour(int_color.PC[,c(7:9)], from='rgb', to='xyz')

int_color_cor<-cbind(int_color.PC, hsv)
int_color_cor<-cbind(int_color_cor, lab)
int_color_cor<-cbind(int_color_cor, xyz)

cor(int_color_cor$PC1, int_color_cor$h)
cor(int_color_cor$PC1, int_color_cor$s)
cor(int_color_cor$PC1, int_color_cor$v)
cor(int_color_cor$PC1, int_color_cor$l)
cor(int_color_cor$PC1, int_color_cor$a)
cor(int_color_cor$PC1, int_color_cor$b)
cor(int_color_cor$PC1, int_color_cor$Red)
cor(int_color_cor$PC1, int_color_cor$Green)
cor(int_color_cor$PC1, int_color_cor$Blue)


Hue<-cor(int_color_cor$PC1, int_color_cor$h)
Saturation<-cor(int_color_cor$PC1, int_color_cor$s)
Value<-cor(int_color_cor$PC1, int_color_cor$v)
L<-cor(int_color_cor$PC1, int_color_cor$l)
A<-cor(int_color_cor$PC1, int_color_cor$a)
B<-cor(int_color_cor$PC1, int_color_cor$b)
Red<-cor(int_color_cor$PC1, int_color_cor$Red)
Green<-cor(int_color_cor$PC1, int_color_cor$Green)
Blue<-cor(int_color_cor$PC1, int_color_cor$Blue)

int_PC1.cor<-c(Hue, Saturation, Value, L, A, B)


Hue<-cor(int_color_cor$PC2, int_color_cor$h)
Saturation<-cor(int_color_cor$PC2, int_color_cor$s)
Value<-cor(int_color_cor$PC2, int_color_cor$v)
L<-cor(int_color_cor$PC2, int_color_cor$l)
A<-cor(int_color_cor$PC2, int_color_cor$a)
B<-cor(int_color_cor$PC2, int_color_cor$b)
int_PC2.cor<-c(Hue, Saturation, Value, L, A, B)

int_PC.cor<-rbind(int_PC1.cor,int_PC2.cor)
colnames(int_PC.cor)<-c("Hue", "Saturation", "Value", "L", "A", "B")


##########################################################################################
## Analysis of tuber skin color
##########################################################################################

#ex_color<-std_color[,c(1:6,15:ncol(std_color))]
ex_color<-std_color[,c(1:6,15:17)]
ex_color<-ex_color[ex_color$tuber != "marker",]


## Lets compare sides of each tuber

ex_color.cols_both<-ex_color[,c(7:ncol(ex_color))]

ex_color.cols_nz<-ex_color.cols_both[ , which(apply(ex_color.cols_both, 2, var) != 0)]
ex_color.pca<-prcomp(ex_color.cols_nz, scale = TRUE, center = TRUE)

ex_pct_variance<-ex_color.pca$sdev^2/sum(ex_color.pca$sdev^2)
#ex_pct_variance<-ex_pct_variance[1:6]
#PCs<-c("PC1", "PC2", "PC3", "PC4","PC5","PC6")
ex_pct_variance<-ex_pct_variance[1:3]
PCs<-c("PC1", "PC2", "PC3")
ex_skree<-cbind(PCs, ex_pct_variance)
ex_skree<-as.data.frame(ex_skree)
colnames(ex_skree)<-c("PC", "Variance_Explained")
#ex_skree$PC<-factor(ex_skree$PC, levels=c("PC1", "PC2", "PC3", "PC4","PC5","PC6",))
ex_skree$Variance_Explained<-as.numeric(as.character(ex_skree$Variance_Explained))
ex_skree$Variance_Explained<-ex_skree$Variance_Explained*100
p<-ggplot(ex_skree, aes(x=PC, y=Variance_Explained)) + geom_point(color="red") + ylab("% of Variance explained") + xlab("PC") + theme_bw()

#pdf("skree_plot_external_color_PCA.pdf", height=4, width=6)
print(p)
#dev.off()

#ex_color.PC<-as.data.frame(ex_color.pca$x[,1:5])
ex_color.PC<-as.data.frame(ex_color.pca$x[,1:3])
ex_color.PC$img_name<-ex_color$img_name
ex_color.PC$clone<-ex_color$clone
ex_color.PC$rep<-ex_color$rep
ex_color.PC$tuber<-ex_color$tuber
ex_color.PC$side<-ex_color$side

color_side1<-ex_color.PC[ex_color.PC$side == 1,]
color_side2<-ex_color.PC[ex_color.PC$side == 2,]

colnames(color_side1)[c(1:3)]<-paste(colnames(color_side1)[c(1:3)], "1", sep="_")
colnames(color_side2)[c(1:3)]<-paste(colnames(color_side2)[c(1:3)], "2", sep="_")

color_compareSides<-merge(color_side1[,-c(4,8)], color_side2[,-c(4,8)], by=c('clone', 'rep', 'tuber'), all=T)
color_compareSides<-color_compareSides[complete.cases(color_compareSides),]

ex_col_PC1.cor<-cor(color_compareSides$PC1_1, color_compareSides$PC1_2)
ex_col_PC2.cor<-cor(color_compareSides$PC2_1, color_compareSides$PC2_2)
ex_col_PC3.cor<-cor(color_compareSides$PC3_1, color_compareSides$PC3_2)

ex_col_PC1.cor
ex_col_PC2.cor
ex_col_PC3.cor

p<-ggplot(color_compareSides, aes(x=PC1_1, y=PC1_2)) + geom_point() + theme_bw() + xlab("PC1 Score: Tuber side 1") +  ylab("PC1 Score: Tuber side 2")
p<-ggplot(color_compareSides, aes(x=PC2_1, y=PC2_2)) + geom_point() + theme_bw() + xlab("PC2 Score: Tuber side 1") +  ylab("PC2 Score: Tuber side 2")


## Now lets work on analysis of single side

ex_color.mdl<-ex_color
ex_color<-ex_color[ex_color$side == 1,]

ex_color.cols<-ex_color[,c(7:ncol(ex_color))]

ex_color.cols_nz<-ex_color.cols[ , which(apply(ex_color.cols, 2, var) != 0)]
ex_color.pca<-prcomp(ex_color.cols_nz, scale = TRUE, center = TRUE)

ex_pct_variance<-ex_color.pca$sdev^2/sum(ex_color.pca$sdev^2)
#ex_pct_variance<-ex_pct_variance[1:6]
#PCs<-c("PC1", "PC2", "PC3", "PC4","PC5","PC6")
ex_pct_variance<-ex_pct_variance[1:3]
PCs<-c("PC1", "PC2", "PC3")
ex_skree<-cbind(PCs, ex_pct_variance)
ex_skree<-as.data.frame(ex_skree)
colnames(ex_skree)<-c("PC", "Variance_Explained")
#ex_skree$PC<-factor(ex_skree$PC, levels=c("PC1", "PC2", "PC3", "PC4","PC5","PC6",))
ex_skree$Variance_Explained<-as.numeric(as.character(ex_skree$Variance_Explained))
ex_skree$Variance_Explained<-ex_skree$Variance_Explained*100
p<-ggplot(ex_skree, aes(x=PC, y=Variance_Explained)) + geom_point(color="red") + ylab("% of Variance explained") + xlab("PC") + theme_bw()

#pdf("skree_plot_external_color_PCA.pdf", height=4, width=6)
print(p)
#dev.off()

#ex_color.PC<-as.data.frame(ex_color.pca$x[,1:5])
ex_color.PC<-as.data.frame(ex_color.pca$x[,1:3])
ex_color.PC$img_name<-ex_color$img_name
ex_color.PC$clone<-ex_color$clone
ex_color.PC$rep<-ex_color$rep
ex_color.PC$tuber<-ex_color$tuber
ex_color.PC$side<-ex_color$side


ex_color.PC$Red<-ex_color$red_ave
ex_color.PC$Green<-ex_color$green_ave
ex_color.PC$Blue<-ex_color$blue_ave

ex_color.PC$hex<-rep("NA", nrow(ex_color.PC))

for(rw in 1:nrow(ex_color.PC)){
  ex_color.PC[rw,"hex"]<-rgb(ex_color.PC[rw, 'Red'],ex_color.PC[rw, 'Green'],ex_color.PC[rw, 'Blue'],max=255)
}

cols<-ex_color.PC$hex
names(cols)<-ex_color.PC$hex

PC1.ex.pct<-round(ex_skree[1,2], 1)
PC2.ex.pct<-round(ex_skree[2,2], 1)

## This figure highlights clones with light and dark skin
ex_color_clones<-ex_color.PC[ex_color.PC$clone %in% c(176, 170, 172, 20, 88, 58, 83),]
p<-ggplot() + geom_point(data=ex_color_clones, aes(x=PC1, y=PC2, fill=as.factor(clone)), shape=21, size=3) + theme_bw() +  theme(text = element_text(size=15)) +xlim(c(-6, 6)) + ylim(-3.5,2) + labs(fill = "Clone")
q<-p + geom_point(data = ex_color.PC, aes(x=PC1, y=PC2, color = hex), size=0.8) + theme_bw() + scale_color_manual(values=cols, guide="none") #+ theme(legend.position="none", text = element_text(size=15)) 

ex_color_clones<-ex_color.PC[ex_color.PC$clone %in% c(176, 170, 172, 20, 88, 58, 83),]

p<-ggplot() + geom_point(data = ex_color.PC, aes(x=PC1, y=PC2, color = hex), size=0.8) + theme_bw() + scale_color_manual(values=cols, guide="none") #+ theme(legend.position="none", text = element_text(size=15)) 
q<-p + geom_point(data=ex_color_clones, aes(x=PC1, y=PC2, fill=as.factor(clone)), shape=21, size=3) + theme_bw() +  theme(text = element_text(size=15)) + labs(fill = "Clone") # +xlim(c(-11, 11)) + ylim(-6,2) +
z<-q + geom_point(data=ex_color_clones, aes(x=PC1, y=PC2, color=hex), size=0.8) + scale_color_manual(values=cols, guide="none")
z<- z + xlab(paste("PC1 (", PC1.ex.pct,  " %)",  sep="")) + ylab(paste("PC2 (", PC2.ex.pct," %)", sep=""))
z<-z+ggtitle("Tuber skin color")



## THis is the figure we want!
ex_color_clones2<-ex_color.PC[ex_color.PC$clone %in% c(176,172, 20, 88, 58, 83),]
ex_color_clones2$clone<-factor(ex_color_clones2$clone, levels=c(176,172, 20, 88, 58, 83))

p<-ggplot(data = ex_color.PC, aes(x=PC1, y=PC2, color = hex), size=0.8) + geom_point() + theme_bw() + scale_color_manual(values=cols, guide="none")
q<-p + geom_point(data=ex_color_clones2, aes(x=PC1, y=PC2, fill=as.factor(clone)), shape=21, size=3) + scale_fill_manual(values=c("navy", "blue", "deepskyblue", "firebrick4", "red", "pink"),guide="none")
z<-q  + xlab(paste("PC1 (", PC1.ex.pct, "%)", sep="")) + ylab(paste("PC2 (", PC2.ex.pct, "%)", sep=""))
z<-z+ggtitle("Tuber skin color")





p<-ggplot(ex_color.PC, aes(x=PC1, y=PC2)) + geom_point(aes(colour = hex), size=0.4) + theme_bw() + scale_colour_manual(values=cols) + theme(legend.position="none", text = element_text(size=15))
q<- p + xlab(paste("PC1 (", PC1.ex.pct,  " %)",  sep="")) + ylab(paste("PC2 (", PC2.ex.pct," %)", sep=""))
q<-q+ggtitle("Tuber skin color")

#pdf("PCA_tuber_external_color", height=4, width=4)
print(q)
#dev.off()

fig.5a_path<-paste(fig.filepath, "/Fig_5a.pdf", sep="")
pdf(fig.5a_path, height=6, width=8)
#print(q)
print(z)
dev.off()


## Lets calculate mean and variance by clone
ex_color.PC<-ex_color.PC[,c(5:8,1:3,9:11)]

## Get average of PCs on a per clone basis
ex_color.PC.ave<-aggregate(.~ clone, data=ex_color.PC[,-c(2:4)], mean)

## Get standard deviation of PCs on a per clone basis
ex_color.PC.sd<-aggregate(.~ clone, data=ex_color.PC[,-c(2:4)], sd)
colnames(ex_color.PC.sd)[2:ncol(ex_color.PC.sd)]<-paste(colnames(ex_color.PC.sd)[2:ncol(ex_color.PC.sd)], ".sd", sep="")

## Combine both parameters in the same data.frame
ex_color.PC.par<-merge(ex_color.PC.ave, ex_color.PC.sd, by=c("clone"))
ex_color.PC.par<-ex_color.PC.par[complete.cases(ex_color.PC.par),]

## Are they correlated? Yes
cor(ex_color.PC.par$PC1, ex_color.PC.par$PC1.sd)

## What are the max and min values
ex_color_PC1.max<-max(ex_color.PC.par$PC1)
ex_color_PC1.min<-min(ex_color.PC.par$PC1)


## Lets make a plot of genotype means
ex_color.PC.par$hex<-rep("NA", nrow(ex_color.PC.par))
ex_color.PC.par$Red<-round(ex_color.PC.par$Red)
ex_color.PC.par$Green<-round(ex_color.PC.par$Green)
ex_color.PC.par$Blue<-round(ex_color.PC.par$Blue)


for(rw in 1:nrow(ex_color.PC.par)){
  ex_color.PC.par[rw,"hex"]<-rgb(ex_color.PC.par[rw, 'Red'],ex_color.PC.par[rw, 'Green'],ex_color.PC.par[rw, 'Blue'],max=255)
}

cols<-ex_color.PC.par$hex
names(cols)<-ex_color.PC.par$hex


p<-ggplot(ex_color.PC.par, aes(x=PC1, y=PC2)) + geom_errorbarh(aes(xmin=PC1-PC1.sd, xmax=PC1+PC1.sd)) + geom_errorbar(aes(ymin=PC2-PC2.sd, ymax=PC2+PC2.sd)) + theme_bw() + geom_point(size=4, shape=19, aes(colour = hex)) + xlab("PC1") + ylab("PC2") + ggtitle("Tuber skin color (genotype mean)") + scale_colour_manual(values=cols) + theme(legend.position="none")  


p<-ggplot(ex_color.PC.par, aes(x=PC1, y=PC2)) + geom_point(aes(colour = hex), size=4) + theme_bw() + scale_colour_manual(values=cols) + theme(legend.position="none") 

p<-ggplot(ex_color.PC.par, aes(x=PC1, y=PC2)) + geom_point() + theme_bw()  + theme(legend.position="none") 

## Calculate H2 for flesh color
ex_color.PC_for_h2<-ex_color.PC[,c(1:2,5:10)]
h2_ex_color<-get_h2(ex_color.PC_for_h2)


## Lets do the same for skin color variance
ex_color.sd_for_h2<-ex_color.PC[,c(1:2,5:10)]
ex_color.sd_for_h2<-aggregate(.~ clone + rep, data=ex_color.sd_for_h2, sd)

h2_ex_color.sd<-get_h2(ex_color.sd_for_h2)

##########
## Lets look at standard error for different sample sizes
##########

## First step is to identify samples with 10 tubers

complete.clones<-names(table(ex_color.PC$clone)[table(ex_color.PC$clone) > 9])
ex_color.PC_complete<-ex_color.PC[ex_color.PC$clone %in% complete.clones,]

ex_color.PC_complete.z<-ex_color.PC_complete

std_error_ex_color<-c()
for(i in 3:10){
  clones<-unique(ex_color.PC_complete.z$clone)
  for (c in clones){
    clone.dat<-ex_color.PC_complete.z[ex_color.PC_complete.z$clone == c,]
    tubers<-sample(c(1:10), i, replace=F)
    small.dat<-clone.dat[clone.dat$tuber %in% tubers, ]
    PC1.sd.err<-std.error(small.dat$PC1)
    PC2.sd.err<-std.error(small.dat$PC2)
    PC3.sd.err<-std.error(small.dat$PC3)
    temp<-c(i,c,PC1.sd.err,PC2.sd.err,PC3.sd.err)
    std_error_ex_color<-rbind(std_error_ex_color, temp)
  }
}

std_error_ex_color<-as.data.frame(std_error_ex_color)
colnames(std_error_ex_color)<-c("Replicates","Clone", "PC1.skin","PC2.skin","PC3.skin")

for(i in 1:ncol(std_error_ex_color)){
  std_error_ex_color[,i]<-as.numeric(as.character(std_error_ex_color[,i]))
}

std_error_ex_color.ag<-aggregate(.~Replicates, data=std_error_ex_color[,-c(2)], mean)

id.vars<-colnames(std_error_ex_color.ag)[c(1)]
measure.vars<-colnames(std_error_ex_color.ag)[c(2:(ncol(std_error_ex_color.ag)))]
std_error_ex_color.ag.long<-melt(std_error_ex_color.ag,
                                  # ID variables - all the variables to keep but not split apart on
                                  id.vars=id.vars,
                                  # The source columns
                                  measure.vars=measure.vars,
                                  # Name of the destination column that will identify the original
                                  # column that the measurement came from
                                  variable.name="Trait",
                                  value.name="SE"
)


p<-ggplot(std_error_ex_color.ag.long, aes(x=Replicates, y=SE, color=Trait)) + geom_line() + theme_bw()
q<-p+geom_point() + xlab("Number of tubers") + ylab("Standard Error")



#### Reviewers want analysis of how PC1 corresponds with other color space channels
hsv<-convert_colour(ex_color.PC[,c(8:10)], from='rgb', to='hsv')
lab<-convert_colour(ex_color.PC[,c(8:10)], from='rgb', to='lab')
xyz<-convert_colour(ex_color.PC[,c(8:10)], from='rgb', to='xyz')

ex_color_cor<-cbind(ex_color.PC, hsv)
ex_color_cor<-cbind(ex_color_cor, lab)
ex_color_cor<-cbind(ex_color_cor, xyz)

Hue<-cor(ex_color_cor$PC1, ex_color_cor$h)
Saturation<-cor(ex_color_cor$PC1, ex_color_cor$s)
Value<-cor(ex_color_cor$PC1, ex_color_cor$v)
L<-cor(ex_color_cor$PC1, ex_color_cor$l)
A<-cor(ex_color_cor$PC1, ex_color_cor$a)
B<-cor(ex_color_cor$PC1, ex_color_cor$b)
Red<-cor(ex_color_cor$PC1, ex_color_cor$Red)
Green<-cor(ex_color_cor$PC1, ex_color_cor$Green)
Blue<-cor(ex_color_cor$PC1, ex_color_cor$Blue)

ex_PC1.cor<-c(Hue, Saturation, Value, L, A, B)


Hue<-cor(ex_color_cor$PC2, ex_color_cor$h)
Saturation<-cor(ex_color_cor$PC2, ex_color_cor$s)
Value<-cor(ex_color_cor$PC2, ex_color_cor$v)
L<-cor(ex_color_cor$PC2, ex_color_cor$l)
A<-cor(ex_color_cor$PC2, ex_color_cor$a)
B<-cor(ex_color_cor$PC2, ex_color_cor$b)
ex_PC2.cor<-c(Hue, Saturation, Value, L, A, B)

ex_PC.cor<-rbind(ex_PC1.cor,ex_PC2.cor)
colnames(ex_PC.cor)<-c("Hue", "Saturation", "Value", "L", "A", "B")


## Lets combine skin and flesh 

color_PC.cor<-rbind(ex_PC.cor, int_PC.cor)
PC<-c("PC1", "PC2", "PC1", "PC2")
Tissue<-c("Skin", "Skin", "Flesh", "Flesh")
temp<-cbind(PC, Tissue)
color_PC.cor<-cbind(as.data.frame(temp), color_PC.cor)

table.dir<-paste(base.dir, "/tables", sep="")
TableColorCor.name<-paste(table.dir, "/ColorCor.csv", sep="")

write.csv(color_PC.cor, file=TableColorCor.name, quote=F, row.names=F)

## Lets make a plot of genotype means
ex_color_cor$hex<-rep("NA", nrow(ex_color_cor))
ex_color_cor$Red<-round(ex_color_cor$Red)
ex_color_cor$Green<-round(ex_color_cor$Green)
ex_color_cor$Blue<-round(ex_color_cor$Blue)


for(rw in 1:nrow(ex_color_cor)){
  ex_color_cor[rw,"hex"]<-rgb(ex_color_cor[rw, 'Red'],ex_color_cor[rw, 'Green'],ex_color_cor[rw, 'Blue'],max=255)
}

cols<-ex_color_cor$hex
names(cols)<-ex_color_cor$hex


p<-ggplot(ex_color_cor, aes(x=PC1, y=s)) + geom_point(aes(colour = hex), size=0.4) + theme_bw() + scale_colour_manual(values=cols) + theme(legend.position="none", text = element_text(size=15))
q<- p + xlab(paste("PC1 (", PC1.ex.pct,  " %)",  sep="")) + ylab("Value")
q<-q+ggtitle("Tuber skin color")


##########
## Lets compare with non-color corrected images
##########

std_NC.filepath<-paste(base.dir, "/A08241_potato_measurements_std_shape_not_corrected.csv", sep="")
std_NC<-read.csv(std_NC.filepath)

ex_color_NC<-std_NC[,c(2:7,16:18)]
ex_color_NC<-ex_color_NC[ex_color_NC$tuber != "marker",]

correctedImgs<-unique(ex_color$img_name)

ex_color_NC<-ex_color_NC[ex_color_NC$img_name %in% correctedImgs, ]

ex_color_NC<-ex_color_NC[ex_color_NC$side == 1,]

ex_color_NC.cols<-ex_color_NC[,c(7:ncol(ex_color_NC))]

ex_color_NC.cols_nz<-ex_color_NC.cols[ , which(apply(ex_color_NC.cols, 2, var) != 0)]
ex_color_NC.pca<-prcomp(ex_color_NC.cols_nz, scale = TRUE, center = TRUE)

ex_pct_variance<-ex_color_NC.pca$sdev^2/sum(ex_color_NC.pca$sdev^2)
#ex_pct_variance<-ex_pct_variance[1:6]
#PCs<-c("PC1", "PC2", "PC3", "PC4","PC5","PC6")
ex_pct_variance<-ex_pct_variance[1:3]
PCs<-c("PC1", "PC2", "PC3")
ex_skree<-cbind(PCs, ex_pct_variance)
ex_skree<-as.data.frame(ex_skree)
colnames(ex_skree)<-c("PC", "Variance_Explained")
#ex_skree$PC<-factor(ex_skree$PC, levels=c("PC1", "PC2", "PC3", "PC4","PC5","PC6",))
ex_skree$Variance_Explained<-as.numeric(as.character(ex_skree$Variance_Explained))
ex_skree$Variance_Explained<-ex_skree$Variance_Explained*100
p<-ggplot(ex_skree, aes(x=PC, y=Variance_Explained)) + geom_point(color="red") + ylab("% of Variance explained") + xlab("PC") + theme_bw()

#pdf("skree_plot_external_color_PCA.pdf", height=4, width=6)
print(p)
#dev.off()

#ex_color_NC.PC<-as.data.frame(ex_color_NC.pca$x[,1:5])
ex_color_NC.PC<-as.data.frame(ex_color_NC.pca$x[,1:3])
ex_color_NC.PC$img_name<-ex_color_NC$img_name
ex_color_NC.PC$clone<-ex_color_NC$clone
ex_color_NC.PC$rep<-ex_color_NC$rep
ex_color_NC.PC$tuber<-ex_color_NC$tuber
ex_color_NC.PC$side<-ex_color_NC$side



ex_color_NC.PC$Red<-ex_color_NC$red_ave
ex_color_NC.PC$Green<-ex_color_NC$green_ave
ex_color_NC.PC$Blue<-ex_color_NC$blue_ave

ex_color_NC.PC$hex<-rep("NA", nrow(ex_color_NC.PC))

for(rw in 1:nrow(ex_color_NC.PC)){
  ex_color_NC.PC[rw,"hex"]<-rgb(ex_color_NC.PC[rw, 'Red'],ex_color_NC.PC[rw, 'Green'],ex_color_NC.PC[rw, 'Blue'],max=255)
}

cols<-ex_color_NC.PC$hex
names(cols)<-ex_color_NC.PC$hex

PC1.ex.pct<-round(ex_skree[1,2], 1)
PC2.ex.pct<-round(ex_skree[2,2], 1)

p<-ggplot(ex_color_NC.PC, aes(x=PC1, y=PC2)) + geom_point(aes(colour = hex), size=0.4) + theme_bw() + scale_colour_manual(values=cols) + theme(legend.position="none", text = element_text(size=15))
q<- p + xlab(paste("PC1 (", PC1.ex.pct,  " %)",  sep="")) + ylab(paste("PC2 (", PC2.ex.pct," %)", sep=""))
q<-q+ggtitle("Tuber skin color")

#pdf("PCA_tuber_external_color", height=4, width=4)
print(q)


ex_color_NC.sd<-aggregate(.~ clone, data=ex_color_NC.PC[-c(4,6,7,8,12)], sd)
colnames(ex_color_NC.sd)[2:ncol(ex_color_NC.sd)]<-paste(colnames(ex_color_NC.sd)[2:ncol(ex_color_NC.sd)], ".sd", sep="")
mean(ex_color_NC.sd$PC1.sd)





##########################################################################################
## Make Supplemental Figs 4
##########################################################################################

skree$Trait<-rep("Biomass profile", nrow(skree))
skree<-skree[1:3,]
int_skree$Trait<-rep("Flesh color", nrow(int_skree))
ex_skree$Trait<-rep("Skin color", nrow(ex_skree))

all_skree<-rbind(skree, int_skree)
all_skree<-rbind(all_skree, ex_skree)

p<-ggplot(all_skree, aes(x=PC, y=Variance_Explained, col=Trait)) + geom_point() + facet_wrap(~Trait) + xlab("Principal component (PC)") + ylab("Variance explained (%)") + theme_bw() + theme(text = element_text(size=15),legend.position = "none", axis.text.x = element_text(angle = 90)) + ylim(0,100)


fig.S4_path<-paste(fig.filepath, "/Fig_S4.pdf", sep="")
pdf(fig.S4_path, height=4, width=10)
print(p)
dev.off()






##########################################################################################
## Correlation between traits
##########################################################################################



## Lets combine data.frames
#colnames(tuber_shape)[22:ncol(tuber_shape)]<-paste(colnames(tuber_shape)[22:ncol(tuber_shape)], ".shape", sep="")
colnames(tuber_shape)[25:ncol(tuber_shape)]<-paste(colnames(tuber_shape)[25:ncol(tuber_shape)], ".shape", sep="")


#colnames(ex_color.PC)[5:9]<-paste(colnames(ex_color.PC)[5:9], ".skin", sep="")
#colnames(ex_color.PC)[2]<-c("replicate")
#colnames(int_color.PC)[4:8]<-paste(colnames(int_color.PC)[4:8], ".flesh", sep="")
#colnames(int_color.PC)[2]<-c("replicate")

colnames(ex_color.PC)[5:10]<-paste(colnames(ex_color.PC)[5:10], ".skin", sep="")
colnames(ex_color.PC)[2]<-c("replicate")
colnames(int_color.PC)[4:9]<-paste(colnames(int_color.PC)[4:9], ".flesh", sep="")
colnames(int_color.PC)[2]<-c("replicate")


#tuber_traits<-merge(tuber_shape[,c(1:3,8,9,13:26)], ex_color.PC[,c(1:3,5:9)], by=c("clone", "replicate", "tuber"))
tuber_traits<-merge(tuber_shape[,c(1:3,8,9,16:27)], ex_color.PC[,c(1:3,5:7)], by=c("clone", "replicate", "tuber"))


int_color.PC$clone<-as.factor(int_color.PC$clone)
int_color.PC$tuber<-as.character(int_color.PC$tuber)
tuber_traits<-merge(tuber_traits, int_color.PC[,c(1:6)], by=c("clone", "replicate", "tuber"))

dist_tuber_traits<-tuber_traits[,c(1:3,10,6,12:14,11,7:8,4,5,15:ncol(tuber_traits))]
colnames(dist_tuber_traits)[c(4:13)]<-c("Weight", "Area (MV)", "Length (Caliper)", "Width (Caliper)","L/W ratio (Caliper)", "SVA", "Length (MV)", "Width (MV)", "L/W ratio (MV)", "Eccentricity")

## Convert to long form 


id.vars<-colnames(dist_tuber_traits)[c(1:3)]
measure.vars<-colnames(dist_tuber_traits)[c(4:(ncol(dist_tuber_traits)))]


dist_tuber_traits_long<-melt(dist_tuber_traits,
                      # ID variables - all the variables to keep but not split apart on
                      id.vars=id.vars,
                      # The source columns
                      measure.vars=measure.vars,
                      # Name of the destination column that will identify the original
                      # column that the measurement came from
                      variable.name="Trait",
                      value.name="Value"
)





p<-ggplot(dist_tuber_traits_long, aes(x=Value)) + geom_histogram(color="black", fill="white") + facet_wrap(~Trait, scales="free") + theme_bw() + xlab("") + ylab("")
p


fig.S2_path<-paste(fig.filepath, "/Fig_S2.pdf", sep="")
pdf(fig.S2_path, height=8, width=8)
print(p)
dev.off()




tuber_trait_cor<-cor(tuber_traits[,c(4:ncol(tuber_traits))], use="complete.obs")
corrplot(tuber_trait_cor, type="upper")


## Tuber skin
## Aggregate by mean
#ex_color.PC.ave<-aggregate(.~clone, data=ex_color.PC[c(1:5,7,11:13)], mean)
#ex_color.PC.sd<-aggregate(.~clone, data=ex_color.PC[c(1:5,7,11:13)], sd)
ex_color.PC.ave<-aggregate(.~clone, data=ex_color.PC[c(1,5:7)], mean)
ex_color.PC.sd<-aggregate(.~clone, data=ex_color.PC[c(1,5:7)], sd)

hist(ex_color.PC.ave$PC1, breaks=50)

## What are the lightest clones?
ex_color.PC.ave[order(ex_color.PC.ave$PC1, decreasing=T),]
ex_color.PC.sd[order(ex_color.PC.sd$PC1, decreasing=T),]

## What are the darkest clones?
ex_color.PC.ave[order(ex_color.PC.ave$PC1, decreasing=F),]
ex_color.PC.sd[order(ex_color.PC.sd$PC1, decreasing=F),]

## Aggregate by mean
#int_color.PC.ave<-aggregate(.~clone, data=int_color.PC[c(1:5,7,10:12)], mean)
#int_color.PC.sd<-aggregate(.~clone, data=int_color.PC[c(1:5,7,10:12)], sd)

int_color.PC.ave<-aggregate(.~clone, data=int_color.PC[c(1,4:6)], mean)
int_color.PC.sd<-aggregate(.~clone, data=int_color.PC[c(1,4:6)], sd)

hist(int_color.PC.ave$PC1, breaks=50)

## What are the lightest clones?
int_color.PC.ave[order(int_color.PC.ave$PC1, decreasing=T),]
int_color.PC.sd[order(int_color.PC.sd$PC1, decreasing=T),]

## What are the darkest clones?
int_color.PC.ave[order(int_color.PC.ave$PC1, decreasing=F),]
int_color.PC.sd[order(int_color.PC.sd$PC1, decreasing=F),]


## Lets make a table that summarizes the results

Trait<-c("Weight", "Area (MV)", "Length (Caliper)", "Width (Caliper)","L/W ratio (Caliper)", "SVA", "Length (MV)", "Width (MV)", "L/W ratio (MV)", "Eccentricity", "PC1.shape")
tuber_size.par
tuber_shape.par
ex_color.PC.par
int_color.PC.par



##########################################################################################
## Make a figure of the standard error
##########################################################################################

std_error_size.ag.long
std_error_shape.ag.long

std_error_int_color.ag.long
std_error_ex_color.ag.long

std_error_color<-rbind(std_error_ex_color.ag.long,std_error_int_color.ag.long)



p<-ggplot(std_error_size.ag.long, aes(x=Replicates, y=SE, color=Trait)) + geom_line() + theme_bw()
q<-p+geom_point() + xlab("Number of tubers") + ylab("Standard Error") + scale_color_manual(values=c("firebrick4", "red1","cornflowerblue", "orchid2", "lawngreen","darkblue","orchid4")) + ggtitle("Tuber size measurements")


fig.S3a_path<-paste(fig.filepath, "/Fig_S3a.pdf", sep="")
pdf(fig.S3a_path, height=8, width=8)
print(q)
dev.off()

p<-ggplot(std_error_shape.ag.long, aes(x=Replicates, y=SE, color=Trait)) + geom_line() + theme_bw()
q<-p+geom_point() + xlab("Number of tubers") + ylab("Standard Error") + scale_color_manual(values=c("dodgerblue", "orange","navy", "red", "darkgreen","limegreen","chartreuse"))  + ggtitle("Tuber shape measurements")

fig.S3b_path<-paste(fig.filepath, "/Fig_S3b.pdf", sep="")
pdf(fig.S3b_path, height=8, width=8)
print(q)
dev.off()

## Color
p<-ggplot(std_error_color, aes(x=Replicates, y=SE, color=Trait)) + geom_line() + theme_bw()
q<-p+geom_point() + xlab("Number of tubers") + ylab("Standard Error") + scale_color_manual(values=c("burlywood4", "burlywood3","burlywood1","gray10","gray65","gray95"))  + ggtitle("Tuber color measurements")

fig.S3c_path<-paste(fig.filepath, "/Fig_S3c.pdf", sep="")
pdf(fig.S3c_path, height=8, width=8)
print(q)
dev.off()


##########################################################################################
## Make a table of heritabilities
##########################################################################################

h2_tuber.size
h2_tuber.shape
colnames(h2_tuber.shape)[12:14]<-paste(colnames(h2_tuber.shape)[12:14], 'shape', sep=".")
h2_ex_color
colnames(h2_ex_color)<-paste(colnames(h2_ex_color), "skin", sep=".")
h2_int_color
colnames(h2_int_color)<-paste(colnames(h2_int_color), "flesh", sep=".")

h2_all_traits<-cbind(h2_tuber.size[,c(17,13)], h2_tuber.shape[,c(9,10,11,8,4,5,1,2,12,13,14)], h2_ex_color[,c(1:3)], h2_int_color[,c(1:3)])
colnames(h2_all_traits)[1:10]<-c("Tuber weight", "Tuber area (MV)", "Tuber length (caliper)", "Tuber width (caliper)", "Tuber L/W ratio (caliper)", "SVA", "Tuber length (MV)", "Tuber width (MV)", "Tuber L/W ratio (MV)","Eccentricity")
h2_all_traits<-h2_all_traits[,c(1:3,7,4,8,5,9,6,10,11:ncol(h2_all_traits))]

H2<-t(h2_all_traits)




colnames(H2)<-c("Broad-sense heritability", "Residual error")


#H2_all<-H2_all[c(4,1,5,2,6,3,11,7,10,8,9,12:nrow(H2_all)),]

#rownames(H2_all)[1:11]<-c("Tuber weight (oz)", "Tuber area (mm2) - CV", "Tuber length (caliper)", "Tuber length (CV)", "Tuber width (calipler)", "Tuber width (CV)", "L/W ratio (caliper)", "L/W ratio (CV)", "SVA", "Eccentricity", "Perimeter (mm)")

#H2_all<-round(H2_all, 2)

table.dir<-paste(base.dir, "/tables", sep="")
Table1.name<-paste(table.dir, "/Table_1.csv", sep="")

write.csv(H2, file=Table1.name, quote=F, row.names=T)








### Send data to Jae for mapping
colnames(int_color.PC)[4:9]<-paste(colnames(int_color.PC)[4:9], "flesh", sep=".")
colnames(ex_color.PC)[5:10]<-paste(colnames(ex_color.PC)[5:10], "skin", sep=".")
ex_color.PC<-ex_color.PC[,-c(4)]

color.PC<-merge(ex_color.PC, int_color.PC, by=c("clone", "rep", "tuber"))

tuber_shape2<-tuber_shape
tuber_shape2<-tuber_shape2[,c(1:3,8,9,16:27)]
colnames(tuber_shape2)[c(15:17)]<-paste(colnames(tuber_shape2)[c(15:17)], "shape", sep=".")
colnames(color.PC)[2]<-c("replicate")

tuber_traits<-merge(tuber_shape2, color.PC, by=c("clone", "replicate", "tuber"))

write.csv(tuber_traits, file="TuberTraits_A08241_2022-11-08.csv", quote=F, row.names=F)
