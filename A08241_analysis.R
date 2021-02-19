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
setwd("/Users/max.feldman/Desktop/A08241")
base.dir<-getwd()

## Black background
std.filepath<-paste(base.dir, "/A08241_potato_measurements_std_shape.csv", sep="")
std<-read.csv(std.filepath)
std_color<-std[2:786]

## Extract shape values to seperate data.frame
shape<-std[,c(2:18, 787:ncol(std))]

## Get standard values from black background
std<-std[,c(1:786)]

## Lightbox background
box.filepath<-paste(base.dir, "/A08241_potato_measurements_box.csv", sep="")
box<-read.csv(box.filepath)

## Combine datasets
both<-rbind(std, box)

## Load in scanner measurements
scan.filepath<-paste(base.dir, "/A08241_potato_measurements_scanner.csv", sep="")
scan<-read.csv(scan.filepath)

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

## Lets add a categorical variable named 'light' to the scan data.frame 
## This identifies the collection platform
marker.scan$light<-rep("scanner", nrow(marker.scan))
marker.scan<-marker.scan[names(marker.top_down.s1)]

## Lets make sure the column order are the same between the top-down and scan data.frames
## Then merge them using rbind
marker.all<-rbind(marker.top_down.s1, marker.scan)

## Lets calculate length and width in mm (poker chip size standard has diameter of 37 mm)
marker.all$mm_per_px<-37/marker.all$length
marker.all$length_mm<-marker.all$length * marker.all$mm_per_px
marker.all$width_mm<-marker.all$width * marker.all$mm_per_px
#marker.all$pct_area<-marker.all$area/mean(marker.all$area)
marker.all$rel_area<-marker.all$area * marker.all$mm_per_px
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
q<-x +  stat_summary(fun.data=mean_se, geom = "errorbar", width=0.2) + theme(text = element_text(size=15), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")  + ylab("Distance (mm2)") + xlab("")

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
std$area_mm<-std$area * std$px_per_mm
std$length_mm<-std$length * std$px_per_mm
std$width_mm<-std$width * std$px_per_mm
std$perimeter_mm<-std$perimeter * std$px_per_mm

std_size<-std[,c("img_name", "clone", "rep", "side", "tuber", "area", "perimeter", "length", "width", "ratio", "eccentricity", "red_ave", "green_ave", "blue_ave", "area_mm", "length_mm", "width_mm","perimeter_mm")]
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
tuber_size.ave<-aggregate(.~ clone, data=tuber_size[,c(1,3:16)], mean)
colnames(tuber_size.ave)[2:ncol(tuber_size.ave)]<-paste(colnames(tuber_size.ave[2:ncol(tuber_size.ave)]), "ave", sep=".")

## Standard deviation
tuber_size.sd<-aggregate(.~ clone, data=tuber_size[,c(1,3:16)], sd)
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
tuber_size_for_h2<-tuber_size[,c(1,2,12:14,16,18:19)]

h2_tuber.size<-get_h2(tuber_size_for_h2)

## Lets do the same for tuber size variance
tuber_sd_for_h2<-tuber_size[,c(1:3,12:14,16,18:20)]
tuber_sd_for_h2<-aggregate(.~ clone + replicate, data=tuber_sd_for_h2, sd)
tuber_sd_for_h2<-tuber_sd_for_h2[,-c(2,3)]

h2_tuber.size.sd<-get_h2(tuber_sd_for_h2)

## Make plots of correlation between computer vision and ground truth measurements
#p<-ggplot(tuber_size, aes(x=caliper_length, y=length_mm)) + geom_point(size=0.5, color=brown"") + theme_bw() + xlab("Caliper length (mm)") + ylab("Computer vision length (mm)") + theme(text = element_text(size=15),legend.position = "none") + geom_smooth(method = "lm", se = FALSE, color="red",linetype="dashed") + stat_cor(aes(size = 2, label =  ..rr.label..))
p<-ggplot(tuber_size, aes(x=caliper_length, y=length_mm)) + geom_point(size=0.5, color="grey60") + theme_bw() + xlab("Caliper length (mm)") + ylab("Computer vision length (mm)") + theme(text = element_text(size=15),legend.position = "none") + geom_smooth(method = "lm", se = FALSE, color="black",linetype="dashed") 

fig.1b_path<-paste(fig.filepath, "/Fig_1b.pdf", sep="")
pdf(fig.1b_path, height=4, width=4)
print(p)
dev.off()

fig.1c_path<-paste(fig.filepath, "/Fig_1c.pdf", sep="")

p<-ggplot(tuber_size, aes(x=caliper_width, y=width_mm)) + geom_point(size=0.5, color="grey90") + theme_bw() + xlab("Caliper width (mm)") + ylab("Computer vision width (mm)")  + theme(text = element_text(size=15),legend.position = "none")+ geom_smooth(method = "lm", se = FALSE, color="black",linetype="dashed") 
pdf(fig.1c_path, height=4, width=4)
print(p)
dev.off()

fig.1a_path<-paste(fig.filepath, "/Fig_1a.pdf", sep="")
p<-ggplot(tuber_size, aes(x=weight, y=area_mm)) + geom_point(size=0.5, color=c("grey30")) + theme_bw() + xlab("Weight (oz)") + ylab("Area (mm2)")  + theme(text = element_text(size=15),legend.position = "none") + geom_smooth(method = "lm", se = FALSE, color="black",linetype="dashed") 

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



std_size.model<-std_size[,c(2:4,15:ncol(std_size))]
size.variance_output<-get_fx(std_size)


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
p<-ggplot(tuber_shape, aes(x=caliper_ratio, y=ratio)) + geom_point(size=0.5, color=c("grey30")) + theme_bw() + xlab("L/W ratio (caliper)") + ylab("L/W ratio (computer vision)")  + theme(text = element_text(size=15),legend.position = "none") + geom_smooth(method = "lm", se = FALSE, color="black",linetype="dashed") 

fig.2a_path<-paste(fig.filepath, "/Fig_2a.pdf", sep="")

pdf(fig.2a_path, height=4, width=4)
print(p)
dev.off()

p<-ggplot(tuber_shape, aes(x=sva, y=caliper_ratio)) + geom_point(size=0.5, color=c("grey60")) + theme_bw() + xlab("SVA") + ylab("L/W ratio (caliper)")  + theme(text = element_text(size=15),legend.position = "none") + geom_smooth(method = "lm", se = FALSE, color="black",linetype="dashed") 
fig.2b_path<-paste(fig.filepath, "/Fig_2b.pdf", sep="")

pdf(fig.2b_path, height=4, width=4)
print(p)
dev.off()


p<-ggplot(tuber_shape, aes(x=caliper_ratio, y=eccentricity)) + geom_point(size=0.5, color=c("grey90")) + theme_bw() + xlab("L/W ratio (caliper)") + ylab("Eccentricity")  + theme(text = element_text(size=15),legend.position = "none") + geom_smooth(method = "lm", se = FALSE, color="black",linetype="dashed") 

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
## Keep only one side of the tuber
shape<-shape[shape$side == 1,]

shape.cols<-shape[,c(18:ncol(shape))]
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

p<-ggplot(tuber_shape, aes(x=ratio, y=PC1)) + geom_point(size=0.5, color=c("grey30")) + theme_bw() + xlab("L/W ratio (computer vision)") + ylab("PC1shape")  + theme(text = element_text(size=15),legend.position = "none") + geom_smooth(method = "lm", se = FALSE, color="black",linetype="dashed") 

fig.3a_path<-paste(fig.filepath, "/Fig_3a.pdf", sep="")

pdf(fig.3a_path, height=4, width=4)
print(p)
dev.off()


cor(tuber_shape$ratio, tuber_shape$PC2, use="complete.obs")

p<-ggplot(tuber_shape, aes(x=ratio, y=PC2)) + geom_point(size=0.5, color=c("grey60")) + theme_bw() + xlab("L/W ratio (computer vision)") + ylab("PC2shape")  + theme(text = element_text(size=15),legend.position = "none") + geom_smooth(method = "lm", se = FALSE, color="black",linetype="dashed") 

fig.3b_path<-paste(fig.filepath, "/Fig_3b.pdf", sep="")

pdf(fig.3b_path, height=4, width=4)
print(p)
dev.off()

cor(tuber_shape$ratio, tuber_shape$PC3, use="complete.obs")

p<-ggplot(tuber_shape, aes(x=ratio, y=PC3)) + geom_point(size=0.5, color=c("grey90")) + theme_bw() + xlab("L/W ratio (computer vision)") + ylab("PC3shape")  + theme(text = element_text(size=15),legend.position = "none") + geom_smooth(method = "lm", se = FALSE, color="black",linetype="dashed") 


fig.3c_path<-paste(fig.filepath, "/Fig_3c.pdf", sep="")

pdf(fig.3c_path, height=4, width=4)
print(p)
dev.off()

cor(tuber_shape$ratio, tuber_shape$PC4, use="complete.obs")

cor(tuber_shape$ratio, tuber_shape$PC5, use="complete.obs")

cor(tuber_shape$ratio, tuber_shape$PC6, use="complete.obs")


## Lets plot boxplot of PC2shape
tuber_shape<-tuber_shape[complete.cases(tuber_shape),]
tuber_shape$clone = reorder(tuber_shape$clone, tuber_shape$PC2, median)
p<-ggplot(tuber_shape, aes(x=as.factor(clone), y=PC2)) + geom_boxplot() + theme_bw()  + theme(axis.text.x = element_text(angle = 90)) + ylab("PC2shape") + xlab("Clone") + theme(text = element_text(size=15),legend.position = "none",axis.text.x = element_text(angle = 90, hjust = 1, size=6))
p


fig.3d_path<-paste(fig.filepath, "/Fig_3d.pdf", sep="")
pdf(fig.3d_path, height=3, width=14)
print(p)
dev.off()


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
clone.ratio.min<-min(tuber_shape.ave$ratio)

clone.ratio.sd/clone.ratio

clone.PC2<-mean(tuber_shape.ave$PC2, na.rm=TRUE)
clone.PC2.sd<-mean(tuber_shape.sd$PC2, na.rm=TRUE)

clone.PC2.max<-max(tuber_shape.ave$PC2, na.rm=TRUE)
clone.PC2.min<-min(tuber_shape.ave$PC2, na.rm=TRUE)

clone.PC2.sd/clone.PC2

cor(tuber_shape.sd$PC2, tuber_shape.ave$PC2, use="complete.obs")


tuber_shape_for_h2<-tuber_shape[,c(1:2,8,9,14:ncol(tuber_shape))]

h2_tuber.shape<-get_h2(tuber_shape_for_h2)

## Lets do the same for tuber shape variance
tuber_sd_for_h2<-tuber_shape[,c(1:3,8,9,14:ncol(tuber_shape))]
tuber_sd_for_h2<-aggregate(.~ clone + replicate, data=tuber_sd_for_h2, sd)
tuber_sd_for_h2<-tuber_sd_for_h2[,-c(3)]

h2_tuber.shape.sd<-get_h2(tuber_sd_for_h2)

## Lets reformat for H2 estimates
PC16<-PC16[,c(6,7,8,9,10,11,1:5)]
PC16<-PC16[PC16$side == 1,]


PC16.ave<-tuber_shape.ave[,c(1,20:ncol(tuber_shape.ave))]
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


PC16.sd<-tuber_shape.sd[,c(1,20:ncol(tuber_shape.sd))]

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

p<-ggplot(PC16.ag, aes(x=PC1, y=PC2)) + geom_errorbarh(aes(xmin=PC1-PC1.sd, xmax=PC1+PC1.sd)) + geom_errorbar(aes(ymin=PC2-PC2.sd, ymax=PC2+PC2.sd)) + theme_bw() + geom_point(size=1.5, shape=19, color="tan4") + xlab("PC1 (Variance explained = 22.4%; H2 = 0.56)") + ylab("PC2 (Variance explained = 20.1%; H2 = 0.23)") + ggtitle("Tuber shape characteristics (genotype mean)")  
tempPC1<-PC16.ag[PC16.ag$clone %in% c('19', '37', '122'),]
q<-p + geom_point(data = tempPC1, aes(x=PC1, y=PC2), color="red3", size = 6) 
x<-q + geom_text(data = tempPC1, aes(x=PC1-1, y=PC2-1, label = as.character(clone)), color="red3", size = 4)
tempPC2<-PC16.ag[PC16.ag$clone %in% c('2', '16', '170'),]
p<-x + geom_point(data = tempPC2, aes(x=PC1, y=PC2), color="blue2", size = 6) 
y<-p + geom_text(data = tempPC2, aes(x=PC1-1, y=PC2-1, label = as.character(clone)), color="blue2", size = 4)


#fig.4_path<-paste(fig.filepath, "/Fig_4.pdf", sep="")
#pdf(fig.4_path, height=6, width=6)
#print(p)
#dev.off()
#pdf("PCA_genotype_means.pdf", height=6, width=6)
print(y)
#dev.off()



##### Lets assess fx 
shape.mdl.cols<-shape.mdl[,c(18:ncol(shape.mdl))]
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




##########################################################################################
## Analysis of color checker standards
##########################################################################################

scan.cc<-read.csv("color_checker_data_scanner.csv")
blk.cc<-read.csv("color_checker_data_std.csv")


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


fig.S4a_path<-paste(fig.filepath, "/Fig_S4a.pdf", sep="")
pdf(fig.S4a_path, height=6, width=8)
print(p)
dev.off()


## Make a plot derived PC values
p<-ggplot(scan.cc_plot, aes(x=PC1, y=PC2)) + geom_point(aes(colour = hex)) + theme_bw() + scale_colour_manual(values=cols) + theme(legend.position="none") + ggtitle("Color checker PC (Scanner)")

p

fig.S4c_path<-paste(fig.filepath, "/Fig_S4c.pdf", sep="")
pdf(fig.S4c_path, height=6, width=8)
print(p)
dev.off()

## Lets get summary of variation among each PC for each chip
scan_chip.sd<-aggregate(.~ chip, data=scan.cc_plot[,c(2,6:8)], sd)


##########
## Now work on the black background 
##########

## Remove problematic images
blk.cc<-blk.cc[!grepl("117_", blk.cc$img_name),]
blk.cc<-blk.cc[!grepl("177_", blk.cc$img_name),]

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
replacements<-c(6,12,18,24,5,11,17,23,4,10,16,22,3,9,15,21,2,8,14,20,1,7,13,19)
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



p<-ggplot(blk.cc_plot_long, aes(x=Color, y=Value)) + geom_point(aes(colour = hex)) + theme_bw() + scale_colour_manual(values=cols) + ylim(0,255) + theme(legend.position="none") + facet_wrap(~chip, ncol=6) + ggtitle("Color checker (Black background)")
#x<-p + stat_summary(fun.y="mean", geom="point", size=6, shape=4, col="red") 
#p +  stat_summary(fun.data=mean_se, geom = "errorbar", width=0.2)
p

fig.S4b_path<-paste(fig.filepath, "/Fig_S4b.pdf", sep="")
pdf(fig.S4b_path, height=6, width=8)
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

## Change chip number to reflect values on scanner

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
p<-ggplot(blk.cc_plot.pc, aes(x=PC1, y=PC2)) + geom_point(aes(colour = hex)) + theme_bw() + scale_colour_manual(values=cols) + theme(legend.position="none") + ggtitle("Color checker PC (Black background)")

fig.S4d_path<-paste(fig.filepath, "/Fig_S4d.pdf", sep="")
pdf(fig.S4d_path, height=6, width=8)
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

##########################################################################################
## Analysis of tuber flesh
##########################################################################################


int_color<-scan[,c(2:5,14:ncol(scan))]
int_color<-int_color[int_color$tuber != "marker",]
int_color.cols<-int_color[,c(8:ncol(int_color))]

scan_pca<-int_color.cols[ , which(apply(int_color.cols, 2, var) != 0)]
scan.color.pca<-prcomp(scan_pca, scale = TRUE, center = TRUE)


## Get % variance explained by first 2 PCs
int_pct.explained<-summary(scan.color.pca)$importance[2,1:2] * 100

int_pct_variance<-scan.color.pca$sdev^2/sum(scan.color.pca$sdev^2)
int_pct_variance<-int_pct_variance[1:15]
PCs<-c("PC1", "PC2", "PC3", "PC4","PC5","PC6", "PC7", "PC8", "PC9","PC10","PC11", "PC12", "PC13", "PC14","PC15")
int_skree<-cbind(PCs, int_pct_variance)
int_skree<-as.data.frame(int_skree)
colnames(int_skree)<-c("PC", "Variance_Explained")
int_skree$PC<-factor(int_skree$PC, levels=c("PC1", "PC2", "PC3", "PC4","PC5","PC6", "PC7", "PC8", "PC9","PC10","PC11", "PC12", "PC13", "PC14","PC15"))
int_skree$Variance_Explained<-as.numeric(as.character(int_skree$Variance_Explained))
int_skree$Variance_Explained<-int_skree$Variance_Explained*100
p<-ggplot(int_skree, aes(x=PC, y=Variance_Explained)) + geom_point(color="red") + ylab("% of Variance explained") + xlab("PC") + theme_bw()

#pdf("skree_plot_flesh_color_PCA.pdf", height=4, width=6)
print(p)
#dev.off()

int_color.PC<-as.data.frame(scan.color.pca$x[,1:5])
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

p<-ggplot(int_color.PC, aes(x=PC1, y=PC2)) + geom_point(aes(colour = hex), size=0.4) + theme_bw() + scale_colour_manual(values=cols) + theme(legend.position="none", text = element_text(size=15))  
q<- p + ylim(-30, 55) +xlim(-75, 50) + xlab(paste("PC1 (", int_pct.explained[1], "%)", sep="")) + ylab(paste("PC2 (", int_pct.explained[2], "%)", sep=""))
q<-q+ggtitle("Tuber flesh color")

#pdf("PCA_tuber_flesh_color", height=4, width=4)
print(q)
#dev.off()
fig.5b_path<-paste(fig.filepath, "/Fig_5b.pdf", sep="")
pdf(fig.5b_path, height=6, width=8)
print(q)
dev.off()



## Lets calculate mean and variance by clone
int_color.PC<-int_color.PC[,c(7:9,1:5,10:12)]

## Get average of PCs on a per clone basis
int_color.PC.ave<-aggregate(.~ clone, data=int_color.PC[,-c(2,3)], mean)

## Get standard deviation of PCs on a per clone basis
int_color.PC.sd<-aggregate(.~ clone, data=int_color.PC[,-c(2,3,9:11)], sd)
colnames(int_color.PC.sd)[2:ncol(int_color.PC.sd)]<-paste(colnames(int_color.PC.sd)[2:ncol(int_color.PC.sd)], ".sd", sep="")

## Combine both parameters in the same data.frame
int_color.PC.par<-merge(int_color.PC.ave, int_color.PC.sd, by=c("clone"))
int_color.PC.par<-int_color.PC.par[complete.cases(int_color.PC.par),]

## Are they correlated? Yes
cor(int_color.PC.par$PC1.ave, int_color.PC.par$PC1.sd)

## What are the max and min values
int_PC1.max<-max(int_color.PC.par$PC1.ave)
int_PC1.min<-min(int_color.PC.par$PC1.ave)


## Lets make a plot of genotype means
int_color.PC.par$hex<-rep("NA", nrow(int_color.PC))
int_color.PC.par$Red<-round(int_color.PC.par$Red)
int_color.PC.par$Green<-round(int_color.PC.par$Green)
int_color.PC.par$Blue<-round(int_color.PC.par$Blue)


for(rw in 1:nrow(int_color.PC.par)){
  int_color.PC.par[rw,"hex"]<-rgb(int_color.PC.par[rw, 'Red'],int_color.PC.par[rw, 'Green'],int_color.PC.par[rw, 'Blue'],max=255)
}

cols<-int_color.PC.par$hex
names(cols)<-int_color.PC.par$hex


p<-ggplot(int_color.PC.par, aes(x=PC1, y=PC2)) + geom_errorbarh(aes(xmin=PC1-PC1.sd, xmax=PC1+PC1.sd)) + geom_errorbar(aes(ymin=PC2-PC2.sd, ymax=PC2+PC2.sd)) + theme_bw() + geom_point(size=4, shape=19, aes(colour = hex)) + xlab("PC1") + ylab("PC2") + ggtitle("Tuber flesh color (genotype mean)") + xlim(-35, 35) + ylim(-25, 25)  + scale_colour_manual(values=cols) + theme(legend.position="none", text = element_text(size=15))  


p<-ggplot(int_color.PC.par, aes(x=PC1, y=PC2)) + geom_point(aes(colour = hex), size=4) + theme_bw() + scale_colour_manual(values=cols) + theme(legend.position="none") 

p<-ggplot(int_color.PC.par, aes(x=PC1, y=PC2)) + geom_point() + theme_bw()  + theme(legend.position="none") 

## Calculate H2 for flesh color
int_color.PC_for_h2<-int_color.PC[,c(1:2,4:8)]
h2_int_color<-get_h2(int_color.PC_for_h2)


## Lets do the same for flesh color variance
int_color.sd_for_h2<-int_color.PC[,c(1:2,4:8)]
int_color.sd_for_h2<-aggregate(.~ clone + rep, data=int_color.sd_for_h2, sd)

h2_int_color.sd<-get_h2(int_color.sd_for_h2)


## Lets get variable fx from scans
int_color.variance_output<-get_fx_scan(int_color.PC)


##########################################################################################
## Analysis of tuber skin color
##########################################################################################

ex_color<-std_color[,c(1:6,15:ncol(std_color))]
ex_color<-ex_color[ex_color$tuber != "marker",]
ex_color.mdl<-ex_color
ex_color<-ex_color[ex_color$side == 1,]

ex_color.cols<-ex_color[,c(10:ncol(ex_color))]

ex_color.cols_nz<-ex_color.cols[ , which(apply(ex_color.cols, 2, var) != 0)]
ex_color.pca<-prcomp(ex_color.cols_nz, scale = TRUE, center = TRUE)

ex_pct_variance<-ex_color.pca$sdev^2/sum(ex_color.pca$sdev^2)
ex_pct_variance<-ex_pct_variance[1:15]
PCs<-c("PC1", "PC2", "PC3", "PC4","PC5","PC6", "PC7", "PC8", "PC9","PC10","PC11", "PC12", "PC13", "PC14","PC15")
ex_skree<-cbind(PCs, ex_pct_variance)
ex_skree<-as.data.frame(ex_skree)
colnames(ex_skree)<-c("PC", "Variance_Explained")
ex_skree$PC<-factor(ex_skree$PC, levels=c("PC1", "PC2", "PC3", "PC4","PC5","PC6", "PC7", "PC8", "PC9","PC10","PC11", "PC12", "PC13", "PC14","PC15"))
ex_skree$Variance_Explained<-as.numeric(as.character(ex_skree$Variance_Explained))
ex_skree$Variance_Explained<-ex_skree$Variance_Explained*100
p<-ggplot(ex_skree, aes(x=PC, y=Variance_Explained)) + geom_point(color="red") + ylab("% of Variance explained") + xlab("PC") + theme_bw()

#pdf("skree_plot_external_color_PCA.pdf", height=4, width=6)
print(p)
#dev.off()

ex_color.PC<-as.data.frame(ex_color.pca$x[,1:5])
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

p<-ggplot(ex_color.PC, aes(x=PC1, y=PC2)) + geom_point(aes(colour = hex), size=0.4) + theme_bw() + scale_colour_manual(values=cols) + theme(legend.position="none", text = element_text(size=15))
q<- p + xlab(paste("PC1 (", PC1.ex.pct,  " %)",  sep="")) + ylab(paste("PC2 (", PC2.ex.pct," %)", sep=""))
q<-q+ggtitle("Tuber skin color")

#pdf("PCA_tuber_external_color", height=4, width=4)
print(q)
#dev.off()

fig.5a_path<-paste(fig.filepath, "/Fig_5a.pdf", sep="")
pdf(fig.5a_path, height=6, width=8)
print(q)
dev.off()


## Lets calculate mean and variance by clone
ex_color.PC<-ex_color.PC[,c(7:10,1:5,11:13)]

## Get average of PCs on a per clone basis
ex_color.PC.ave<-aggregate(.~ clone, data=ex_color.PC[,-c(2:4)], mean)

## Get standard deviation of PCs on a per clone basis
ex_color.PC.sd<-aggregate(.~ clone, data=ex_color.PC[,-c(2:4,10:12)], sd)
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


p<-ggplot(ex_color.PC.par, aes(x=PC1, y=PC2)) + geom_errorbarh(aes(xmin=PC1-PC1.sd, xmax=PC1+PC1.sd)) + geom_errorbar(aes(ymin=PC2-PC2.sd, ymax=PC2+PC2.sd)) + theme_bw() + geom_point(size=4, shape=19, aes(colour = hex)) + xlab("PC1") + ylab("PC2") + ggtitle("Tuber skin color (genotype mean)") + xlim(-35, 35) + ylim(-25, 25)  + scale_colour_manual(values=cols) + theme(legend.position="none")  


p<-ggplot(ex_color.PC.par, aes(x=PC1, y=PC2)) + geom_point(aes(colour = hex), size=4) + theme_bw() + scale_colour_manual(values=cols) + theme(legend.position="none") 

p<-ggplot(ex_color.PC.par, aes(x=PC1, y=PC2)) + geom_point() + theme_bw()  + theme(legend.position="none") 

## Calculate H2 for flesh color
ex_color.PC_for_h2<-ex_color.PC[,c(1:2,5:9)]
h2_ex_color<-get_h2(ex_color.PC_for_h2)


## Lets do the same for skin color variance
ex_color.sd_for_h2<-ex_color.PC[,c(1:2,5:9)]
ex_color.sd_for_h2<-aggregate(.~ clone + rep, data=ex_color.sd_for_h2, sd)

h2_ex_color.sd<-get_h2(ex_color.sd_for_h2)


##### Lets assess fx 
ex_color.mdl.cols<-ex_color.mdl[,c(10:ncol(ex_color.mdl))]
ex_color.mdl.cols_nz<-ex_color.mdl.cols[ , which(apply(ex_color.mdl.cols, 2, var) != 0)]
ex_color.mdl.pca<-prcomp(ex_color.mdl.cols_nz, scale = TRUE, center = TRUE)

ex_color.mdl.PC15<-as.data.frame(ex_color.mdl.pca$x[,1:5])
ex_color.mdl.PC15$clone<-ex_color.mdl$clone
ex_color.mdl.PC15$rep<-ex_color.mdl$rep
ex_color.mdl.PC15$side<-ex_color.mdl$side
ex_color.mdl.PC15$tuber<-ex_color.mdl$tuber

## Merge PC values and other measurements of shape

tuber_ex_color.mdl<-merge(ex_color.mdl[,c(2:4,6:9)], ex_color.mdl.PC15, by=c('clone', 'rep', 'side', 'tuber'), all=T)

tuber_ex_color.model<-tuber_ex_color.mdl[,c(1:3,5:ncol(tuber_ex_color.mdl))]
ex_color.variance_output<-get_fx(tuber_ex_color.model)




##########################################################################################
## Make Supplemental Figs 3
##########################################################################################

skree$Trait<-rep("Biomass profile", nrow(skree))
int_skree$Trait<-rep("Flesh color", nrow(int_skree))
ex_skree$Trait<-rep("Skin color", nrow(ex_skree))

all_skree<-rbind(skree, int_skree)
all_skree<-rbind(all_skree, ex_skree)

p<-ggplot(all_skree, aes(x=PC, y=Variance_Explained, col=Trait)) + geom_point() + facet_wrap(~Trait) + xlab("Principal component (PC)") + ylab("Variance explained (%)") + theme_bw() + theme(text = element_text(size=15),legend.position = "none", axis.text.x = element_text(angle = 90))


fig.S3_path<-paste(fig.filepath, "/Fig_S3.pdf", sep="")
pdf(fig.S3_path, height=4, width=10)
print(p)
dev.off()






##########################################################################################
## Correlation between traits
##########################################################################################



## Lets combine data.frames
colnames(tuber_shape)[22:ncol(tuber_shape)]<-paste(colnames(tuber_shape)[22:ncol(tuber_shape)], ".shape", sep="")
colnames(ex_color.PC)[5:9]<-paste(colnames(ex_color.PC)[5:9], ".skin", sep="")
colnames(ex_color.PC)[2]<-c("replicate")
colnames(int_color.PC)[4:8]<-paste(colnames(int_color.PC)[4:8], ".flesh", sep="")
colnames(int_color.PC)[2]<-c("replicate")
tuber_traits<-merge(tuber_shape[,c(1:3,8,9,13:26)], ex_color.PC[,c(1:3,5:9)], by=c("clone", "replicate", "tuber"))
int_color.PC$clone<-as.factor(int_color.PC$clone)
int_color.PC$tuber<-as.character(int_color.PC$tuber)
tuber_traits<-merge(tuber_traits, int_color.PC[,c(1:8)], by=c("clone", "replicate", "tuber"))

dist_tuber_traits<-tuber_traits[,c(1:3,10,6,12:14,11,7:8,4,15:17,20:22,25:27)]
colnames(dist_tuber_traits)[c(4:12)]<-c("Weight", "Area (CV)", "Length (Caliper)", "Width (Caliper)","L/W ratio (Caliper)", "SVA", "Length (CV)", "Width (CV)", "L/W ratio (CV)")

## Convert to long form 


id.vars<-colnames(dist_tuber_traits)[c(1:3)]
measure.vars<-colnames(dist_tuber_traits)[c(4:(ncol(dist_tuber_traits)))]

dist_tuber_traits<-dist_tuber_traits[dist_tuber_traits$PC1.flesh > -45,]

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
ex_color.PC.ave<-aggregate(.~clone, data=ex_color.PC[c(1:5,7,11:13)], mean)
ex_color.PC.sd<-aggregate(.~clone, data=ex_color.PC[c(1:5,7,11:13)], sd)

hist(ex_color.PC.ave$PC1, breaks=50)

## What are the lightest clones?
ex_color.PC.ave[order(ex_color.PC.ave$PC1, decreasing=T),]
ex_color.PC.sd[order(ex_color.PC.sd$PC1, decreasing=T),]

## What are the darkest clones?
ex_color.PC.ave[order(ex_color.PC.ave$PC1, decreasing=F),]
ex_color.PC.sd[order(ex_color.PC.sd$PC1, decreasing=F),]

## Aggregate by mean
int_color.PC.ave<-aggregate(.~clone, data=int_color.PC[c(1:5,7,10:12)], mean)
int_color.PC.sd<-aggregate(.~clone, data=int_color.PC[c(1:5,7,10:12)], sd)

hist(int_color.PC.ave$PC1, breaks=50)

## What are the lightest clones?
int_color.PC.ave[order(int_color.PC.ave$PC1, decreasing=T),]
int_color.PC.sd[order(int_color.PC.sd$PC1, decreasing=T),]

## What are the darkest clones?
int_color.PC.ave[order(int_color.PC.ave$PC1, decreasing=F),]
int_color.PC.sd[order(int_color.PC.sd$PC1, decreasing=F),]

## Plot values of PC1 and PC2
## Get hex values
int_color.PC.ave$Red<-round(int_color.PC.ave$Red)
int_color.PC.ave$Green<-round(int_color.PC.ave$Green)
int_color.PC.ave$Blue<-round(int_color.PC.ave$Blue)



## Now lets combine these two data frames
colnames(shape.PC)[c(1:5,11)]<-c("Shape.PC1","Shape.PC2","Shape.PC3","Shape.PC4","Shape.PC5", "L/W_Ratio")
colnames(ex_color.PC)[c(1:5,11:ncol(ex_color.PC))]<-c("Ex_Col.PC1","Ex_Col.PC2","Ex_Col.PC3","Ex_Col.PC4","Ex_Col.PC5", "Ex_Col.Red", "Ex_Col.Green", "Ex_Col.Blue","Ex_Col.Hex")
shape.base<-std[,c(1:4,6,9:12,14)]
colnames(shape.base)[6:ncol(shape.base)]<-c("Area", "Perimeter", "Length", "Width", "Eccentricity")

external.all<-merge(shape.base, shape.PC, by=c("img_name", "clone", "rep", "tuber","side"),all=T)
external.all<-merge(external.all, ex_color.PC, by=c("img_name", "clone", "rep", "tuber","side"),all=T)

colnames(external.all)[2:5]<-c("Clone", "Rep", "Tuber", "Side")
colnames(external.all)[16]<-c("Ratio")

variance_by_side<-c()
for(i in c(6:24)){
  print(i)
  trait<-colnames(external.all)[i]
  #print(trait)
  #f = paste(trait, "~Clone + Rep + Side", sep="")
  #print(f)
  #mdl<-lm(formula = f, data=external.all)
  
  model<-lmer(external.all[,i]~(1|Clone)+(1|Rep)+(1|Side), data=external.all, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
  re<-as.numeric(VarCorr(model))
  res<-attr(VarCorr(model), "sc")^2
  clone.var<-re[1]
  rep.var<-re[2]
  side.var<-re[3]
  # Total variance is sum of all variances
  tot.var<-sum(re, res)
  # Get proportion of variance for all factors
  clone<-clone.var/tot.var
  rep<-rep.var/tot.var
  side<-side.var/tot.var
  residual<-res/tot.var
  #clones<-c(clones,clone)
  #reps<-c(reps,rep)
  #sides<-c(sides,side)
  #residuals<-c(residuals, res)
  Pct_Variance<-c(clone, rep, side, residual)
  Factor<-c("Clone", "Replicate", "Side", "Residual")
  Trait<-rep(trait, 4)
  output<-cbind(Trait, Factor, Pct_Variance)
  variance_by_side<-rbind(variance_by_side, output)
  
  #mat<-summary(mdl)$coef
  #df<-as.data.frame(mat)
  #df<-df[2:nrow(df),]
  #colnames(df)<-c("Var", "Std.Err", "T-value", "P-value")
  #df$Trait<-rep(trait, nrow(df))
  #df$Source<-rownames(df)
  #rownames(df)<-NULL
  #variance_by_side<-rbind(variance_by_side, df)
}

variance_by_side<-as.data.frame(variance_by_side)
variance_by_side$Pct_Variance<-as.numeric(as.character(variance_by_side$Pct_Variance))

p<-ggplot(variance_by_side, aes(x=Trait, y=Pct_Variance, group=Factor)) + geom_line(aes(colour = Factor)) + theme_bw() 
scale_fill_manual('Treatments', values=c('Genotype', 'Treatment', 'Plot', 'G X Treatment', 'Error')) + ylim(0,1) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

p = ggplot(data=variance_by_side, aes(x=factor(1), y=Pct_Variance, fill = factor(Factor)))
p=p + geom_bar(width = 1, stat="identity")
p=p+facet_grid(facets=. ~ Trait)
p=p+xlab("Trait") + theme(legend.title=element_blank())

########################################
## Now lets get how variable each trait is

test<-aggregate(external.all_long[,8:24], by=list(external.all_long$Clone), mean, na.action=na.pass, na.rm=TRUE)


## Convert external.all to long form                                               
external.all_long<-external.all[,c(1:24)]
id.vars<-colnames(external.all_long)[c(1:5)]
measure.vars<-colnames(external.all_long)[c(7:(ncol(external.all_long)))]

external.all_long<-melt(external.all_long,
                        # ID variables - all the variables to keep but not split apart on
                        id.vars=id.vars,
                        # The source columns
                        measure.vars=measure.vars,
                        # Name of the destination column that will identify the original
                        # column that the measurement came from
                        variable.name="Trait",
                        value.name="Value"
)





##########################################################################################
## Make a table of heritabilities
##########################################################################################

h2_tuber.size
h2_tuber.shape
colnames(h2_tuber.shape)[11:16]<-paste(colnames(h2_tuber.shape)[11:16], 'shape', sep=".")
h2_ex_color
colnames(h2_ex_color)<-paste(colnames(h2_ex_color), "skin", sep=".")
h2_int_color
colnames(h2_int_color)<-paste(colnames(h2_int_color), "flesh", sep=".")

h2_all_traits<-cbind(h2_tuber.size, h2_tuber.shape[,c(1:2,5,7,10:15)], h2_ex_color, h2_int_color)

h2_tuber.size.sd
h2_tuber.shape.sd
colnames(h2_tuber.shape.sd)[11:16]<-paste(colnames(h2_tuber.shape.sd)[11:16], 'shape', sep=".")
h2_ex_color.sd
colnames(h2_ex_color.sd)<-paste(colnames(h2_ex_color.sd), "skin", sep=".")
h2_int_color.sd
colnames(h2_int_color.sd)<-paste(colnames(h2_int_color.sd), "skin", sep=".")

h2_all_traits.sd<-cbind(h2_tuber.size.sd, h2_tuber.shape.sd[,c(1:2,5,7,10:15)], h2_ex_color.sd, h2_int_color.sd)

H2<-t(h2_all_traits)



##########################################################################################
## Make a table of % variance explained by factors
##########################################################################################

colnames(size.variance_output)<-c("Tuber area", "Length", "Width", "Perimeter")
colnames(shape.variance_output)<-c("L/W ratio", "Eccentricity", "PC1.shape", "PC2.shape", "PC3.shape", "PC4.shape", "PC5.shape")
ex_color.variance_output<-ex_color.variance_output[,c(4:8,1:3)]
colnames(ex_color.variance_output)<-c("PC1.skin", "PC2.skin", "PC3.skin", "PC4.skin", "PC5.skin", "Red.skin", "Green.skin","Blue.skin")
colnames(int_color.variance_output)<-paste(colnames(int_color.variance_output), "flesh", sep=".")
Side<-rep("NA", ncol(int_color.variance_output))
int_color.variance_output<-rbind(int_color.variance_output, Side)
int_color.variance_output<-int_color.variance_output[c(1,2,4,3),]

variable_fx<-cbind(size.variance_output, shape.variance_output, ex_color.variance_output, int_color.variance_output)
variable_fx[1,]<-as.numeric(variable_fx[1,])*100
variable_fx[2,]<-as.numeric(variable_fx[2,])*100
variable_fx[3,1:19]<-as.numeric(variable_fx[3,1:19])*100
variable_fx[4,]<-as.numeric(variable_fx[4,])*100

variable_fx[1,]<-round(as.numeric(variable_fx[1,]),1)
variable_fx[2,]<-round(as.numeric(variable_fx[2,]),1)
variable_fx[3,1:19]<-round(as.numeric(variable_fx[3,1:19]),1)
variable_fx[4,]<-round(as.numeric(variable_fx[4,]),1)


write.csv(variable_fx, "Table_S1.csv", quote=F, row.names=T)



H2.sd<-t(h2_all_traits.sd)

H2_all<-cbind(H2[,1], H2.sd[,1])

colnames(H2_all)<-c("Broad-sense heritability", "Broad-sense heritability of trait standard deviation")
H2_all<-H2_all[c(4,1,5,2,6,3,11,7,10,8,9,12:nrow(H2_all)),]

rownames(H2_all)[1:11]<-c("Tuber weight (oz)", "Tuber area (mm2) - CV", "Tuber length (caliper)", "Tuber length (CV)", "Tuber width (calipler)", "Tuber width (CV)", "L/W ratio (caliper)", "L/W ratio (CV)", "SVA", "Eccentricity", "Perimeter (mm)")

H2_all<-round(H2_all, 2)

write.csv(H2_all, file="Table_1.csv", quote=F, row.names=T)




