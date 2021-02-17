## This is a R script to perform analysis presented in Park, Feldman, et al., 2021

## Load libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(lme4)
library(plotrix)

## Load in data
setwd("/Users/max.feldman/Documents/data/2020/A08241_imgs")
base.dir<-getwd()

## Black background
std.filepath<-paste(base.dir, "/std/A08241_potato_measurements_std.csv", sep="")
std<-read.csv(std.filepath)

## Lightbox background
box.filepath<-paste(base.dir, "/box/A08241_potato_measurements_box.csv", sep="")
box<-read.csv(box.filepath)

## Combine datasets
both<-rbind(std, box)

## Load in scanner measurements
scan.filepath<-paste(base.dir, "/scanner/A08241_potato_measurements_scanner.csv", sep="")
scan<-read.csv(scan.filepath)

## Load in ground truth data
gt.filepath<-paste(base.dir, "/ground_truth_data.csv", sep="")
gt<-read.csv(gt.filepath)

## Lets make an directory to output figures for the manuscript
fig.filepath<-paste(base.dir, "/figures", sep="")
dir.create(fig.filepath)


##########################################################################################
## Analysis of size standards
##########################################################################################

## Get only values for the size marker (poker chip) from the top downimaging config.

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

#pdf('variation_of_size_marker.pdf', height=6, width=5)
print(q)
#dev.off()

## Lets make Fig. 1A
rel_area<-marker.all_long[marker.all_long$Trait == 'rel_area',]

p<-ggplot(rel_area, aes(x=Configuration, y=Value, col=Configuration)) + geom_jitter(size = 0.1, width=0.02, height=0) + theme_bw()
x<-p + stat_summary(fun.y="mean", geom="point", size=6, shape=4, col="black") + scale_color_manual(values = c("red", "blue", "black"))
q<-x +  stat_summary(fun.data=mean_se, geom = "errorbar", width=0.2) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Imaging configuration") + ylab("Area (mm2)")


## Fig. 1B
width_mm<-marker.all_long[marker.all_long$Trait == 'width_mm',]

p<-ggplot(width_mm, aes(x=Configuration, y=Value, col=Configuration)) + geom_jitter(size = 0.1, width=0.02, height=0) + theme_bw()
x<-p + stat_summary(fun.y="mean", geom="point", size=6, shape=4, col="black") + scale_color_manual(values = c("red", "blue", "black"))
q<-x +  stat_summary(fun.data=mean_se, geom = "errorbar", width=0.2) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Imaging configuration") + ylab("Distance (mm2)")



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
q<-x +  stat_summary(fun.data=mean_se, geom = "errorbar", width=0.2) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Imaging configuration") + ylab("Difference from average (%)")


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

######################
## Lets evaluate tuber size on black background
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

## Lets get the mean measurement from both sizes
std_size.ag<-aggregate(std_size[,c(2,4:ncol(std_size))], by=list(std_size$clone, std_size$tuber), mean, drop=TRUE)
std_size.ag<-std_size.ag[,c(1,2,6:ncol(std_size.ag))]
colnames(std_size.ag)[1:2]<-c('clone', 'tuber')

tuber_size<-merge(std_size.ag, gt, by=c("clone", "tuber"))

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


tuber_size.ave<-aggregate(tuber_size[,c(1:2,16)], by=list(tuber_size$clone), mean, na.action = na.pass)
tuber_size.ave<-tuber_size.ave[,c(1,4)]
colnames(tuber_size.ave)[1:2]<-c("clone", "weight.ave")
tuber_size.sd<-aggregate(tuber_size[,c(1:2,16)], by=list(tuber_size$clone), sd)
tuber_size.sd<-tuber_size.sd[,c(1,4)]
colnames(tuber_size.sd)[1:2]<-c("clone", "weight.sd")

tuber_size.par<-merge(tuber_size.ave, tuber_size.sd, by=c("clone"))
tuber_size.par<-tuber_size.par[complete.cases(tuber_size.par),]

cor(tuber_size.par$weight.ave, tuber_size.par$weight.sd)

clone.ave<-mean(tuber_size.par$weight.ave)
clone.ave.sd<-mean(tuber_size.par$weight.sd)

clone.ave.max<-max(tuber_size.par$weight.ave)
clone.ave.min<-min(tuber_size.par$weight.ave)

clone.ave.sd/clone.ave


tuber_size_for_h2<-tuber_size[,c(1:3,12:14,16,18:19)]

# Heritability fxn
get_h2<-function(data){

  variance.out<-c()
  H2<-c()
  e2<-c()
  # For each treatment.phenotype calculate variance
  for(i in 4:length(colnames(data))){
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
    colnames(variance)<-colnames(data)[4:length(data)]
    rownames(variance)<-c('Genotype', 'Error')
    return(variance)
}  

h2<-get_h2(tuber_size_for_h2)

## Lets do the same for tuber size variance
tuber_sd_for_h2<-tuber_size[,c(1:4,12:14,16,18:20)]
tuber_sd_for_h2<-aggregate(tuber_sd_for_h2, by=list(tuber_sd_for_h2$clone, tuber_sd_for_h2$replicate), sd)
tuber_sd_for_h2<-tuber_sd_for_h2[,c(1,2,4:10)]
colnames(tuber_sd_for_h2)[1:2]<-c("clone", "replicate")

h2<-get_h2(tuber_sd_for_h2)



model<-lmer(tuber_size$weight~(1|clone)+(1|replicate), data=tuber_size)

re<-as.numeric(VarCorr(model))
res<-attr(VarCorr(model), "sc")^2


tot.var<-sum(re, res)


p<-ggplot(tuber_size, aes(x=caliper_length, y=length_mm)) + geom_point()

p<-ggplot(tuber_size, aes(x=caliper_width, y=width_mm)) + geom_point()

p<-ggplot(tuber_size, aes(x=weight, y=area_mm)) + geom_point()

#######################################################
## Measurement of tuber shape

## Lets start with size standard

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

## Now lets get measurements of tubers on black background
## Get px_per_mm for each img
std$area_mm<-std$area * std$px_per_mm
std$length_mm<-std$length * std$px_per_mm
std$width_mm<-std$width * std$px_per_mm
std$perimeter_mm<-std$perimeter * std$px_per_mm

std_shape<-std[,c("img_name", "clone", "rep", "side", "ratio", "eccentricity", "red_ave", "green_ave", "blue_ave", "area_mm", "length_mm", "width_mm","perimeter_mm")]
std_shape<-std_size[std_size$tuber != 'marker',]

## Lets get the mean measurement from both sizes
std_shape.ag<-aggregate(std_shape[,c(2,5:ncol(std_shape))], by=list(std_shape$clone, std_shape$tuber), mean, drop=TRUE)
std_shape.ag<-std_shape.ag[,c(1,2,5:ncol(std_shape.ag))]
colnames(std_shape.ag)[1:2]<-c('clone', 'tuber')

tuber_shape<-merge(std_shape.ag, gt, by=c("clone", "tuber"))
tuber_shape$caliper_ratio<-tuber_shape$caliper_length/tuber_shape$caliper_width

cor(tuber_shape$sva, tuber_shape$caliper_ratio, use="complete.obs")

cor(tuber_shape$ratio, tuber_shape$caliper_ratio, use="complete.obs")

cor(tuber_shape$ratio, tuber_shape$eccentricity, use="complete.obs")

p<-ggplot(tuber_shape, aes(x=sva, y=caliper_ratio)) + geom_point()

p<-ggplot(tuber_shape, aes(x=ratio, y=caliper_ratio)) + geom_point()

tuber_shape$clone = reorder(tuber_shape$clone, tuber_shape$ratio, median)
p<-ggplot(tuber_shape, aes(x=as.factor(clone), y=ratio)) + geom_boxplot() + theme_bw()  + theme(axis.text.x = element_text(angle = 90)) + ylab("L/W Ratio") + xlab("Clone")
p



## Generate tuber biomass profiles

std.filepath<-paste(base.dir, "/std/A08241_potato_measurements_std_shape.csv", sep="")
std<-read.csv(std.filepath)
std<-std[,2:ncol(std)]

shape<-std[,c(1:6,11:13,786:ncol(std))]
shape<-shape[shape$tuber != "marker",]
shape<-shape[shape$side == 1,]


shape.cols<-shape[,c(10:ncol(shape))]

shape.cols_nz<-shape.cols[ , which(apply(shape.cols, 2, var) != 0)]

shape<-shape[,c(1:6,11:13,786:ncol(shape))]
shape<-shape[shape$tuber != "marker",]


shape.cols<-shape[,c(10:ncol(shape))]

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


tuber_shape<-merge(tuber_shape, PC16, by=c('clone', 'tuber', 'replicate'), all=T)

p<-ggplot(PC16, aes(x=PC1, y=PC2)) + geom_point() + theme_bw() + theme(legend.position="none")  
q<- p + ylim(-30, 55) +xlim(-75, 50) + xlab(paste("PC1 (", pct.explained[1], "%)", sep="")) + ylab(paste("PC2 (", pct.explained[2], "%)", sep=""))
q<-q+ggtitle("PCA of Tuber biomass profile")


cor(tuber_shape$ratio, tuber_shape$PC1, use="complete.obs")

cor(tuber_shape$ratio, tuber_shape$PC2, use="complete.obs")

cor(tuber_shape$ratio, tuber_shape$PC3, use="complete.obs")

cor(tuber_shape$ratio, tuber_shape$PC4, use="complete.obs")

cor(tuber_shape$ratio, tuber_shape$PC5, use="complete.obs")

cor(tuber_shape$ratio, tuber_shape$PC6, use="complete.obs")




tuber_shape


tuber_shape.ave<-aggregate(tuber_shape[,c(1,3:ncol(tuber_shape))], by=list(tuber_shape$clone), mean, na.action = na.pass)
tuber_shape.ave<-tuber_shape.ave[,c(1,4:ncol(tuber_shape.ave))]
colnames(tuber_shape.ave)[1]<-c("clone")
tuber_shape.sd<-aggregate(.~ clone, data=tuber_shape[,c(1,4:ncol(tuber_shape))], sd)

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


h2<-get_h2(tuber_shape[,c(1:25)])



## Lets reformat for H2 estimates
PC16<-PC16[,c(6,7,8,9,10,11,1:5)]
PC16<-PC16[PC16$side == 1,]




#PC16.ave<-aggregate(PC16, by=list(PC16$clone), mean, na.action = na.pass, na.rm=TRUE)
#PC16.ave<-PC16.ave[,-c(1,2,4,6)]

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


PC16.stderr<-aggregate(PC16[,c(2,6:11)], by=list(as.factor(PC16$clone)), std.error)
PC16.stderr<-PC16.stderr[,-c(2)]
colnames(PC16.stderr)[c(1)]<-c('clone')

id.vars<-colnames(PC16.stderr)[c(1:2)]
measure.vars<-colnames(PC16.stderr)[c(3:(ncol(PC16.stderr)))]
PC16.stderr.long<-melt(PC16.stderr,
                       # ID variables - all the variables to keep but not split apart on
                       id.vars=id.vars,
                       # The source columns
                       measure.vars=measure.vars,
                       # Name of the destination column that will identify the original
                       # column that the measurement came from
                       variable.name="Trait",
                       value.name="Value"
)


colnames(PC16.stderr)[2:ncol(PC16.stderr)]<-paste(colnames(PC16.stderr)[2:ncol(PC16.stderr)], "std.err", sep="_")

PC16.ag<-merge(PC16.ave, PC16.stderr, by=c("clone"))

p<-ggplot(PC16.ag, aes(x=PC1, y=PC2)) + geom_errorbarh(aes(xmin=PC1-PC1_std.err, xmax=PC1+PC1_std.err)) + geom_errorbar(aes(ymin=PC2-PC2_std.err, ymax=PC2+PC2_std.err)) + theme_bw() + geom_point(size=1.5, shape=19, color="tan4") + xlab("PC1 (Variance explained = 22.4%; H2 = 0.56)") + ylab("PC2 (Variance explained = 20.1%; H2 = 0.23)") + ggtitle("Tuber shape characteristics (genotype mean)")  
tempPC1<-PC16.ag[PC16.ag$clone %in% c('19', '37', '122'),]
q<-p + geom_point(data = tempPC1, aes(x=PC1, y=PC2), color="red3", size = 6) 
x<-q + geom_text(data = tempPC1, aes(x=PC1-1, y=PC2-1, label = as.character(clone)), color="red3", size = 4)
tempPC2<-PC16.ag[PC16.ag$clone %in% c('2', '16', '170'),]
p<-x + geom_point(data = tempPC2, aes(x=PC1, y=PC2), color="blue2", size = 6) 
y<-p + geom_text(data = tempPC2, aes(x=PC1-1, y=PC2-1, label = as.character(clone)), color="blue2", size = 4)



#pdf("PCA_genotype_means.pdf", height=6, width=6)
print(y)
#dev.off()


## Lets calculate H2
H2<-c()
e2<-c()

# For each treatment.phenotype calculate variance
for(i in 7:length(colnames(PC16))){
  print(colnames(PC16[i]))
  # Use only RILs with all measurements for each treatment.phenotype
  cc.PC16<-PC16[complete.cases(PC16[,i]),c(1:6,i)]
  # Build linear model each cofactor is a random effect
  model<-lmer(cc.PC16[,7]~(1|clone), data=cc.PC16, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
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










#######################################################
## Now lets get color

## Lets start with color checker standards

scan.cc<-read.csv("color_checker_data_scanner_rotated_median.csv")
blk.cc<-read.csv("color_checker_data_std.csv")


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


p<-ggplot(scan.cc_plot_long, aes(x=Color, y=Value)) + geom_point(aes(colour = hex)) + theme_bw() + scale_colour_manual(values=cols) + ylim(0,255) + theme(legend.position="none") + facet_wrap(~chip, ncol=6)
#x<-p + stat_summary(fun.y="mean", geom="point", size=6, shape=4, col="red") 
#p +  stat_summary(fun.data=mean_se, geom = "errorbar", width=0.2)
p



## Plot derived PC values
p<-ggplot(scan.cc_plot, aes(x=PC1, y=PC2)) + geom_point(aes(colour = hex)) + theme_bw() + scale_colour_manual(values=cols) + theme(legend.position="none") 


## Lets get summary of variation among each PC for each chip
scan_chip.sd<-aggregate(.~ chip, data=scan.cc_plot[,c(2,6:8)], sd)



## Now work on the black background 


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


p<-ggplot(blk.cc_plot_long, aes(x=Color, y=Value)) + geom_point(aes(colour = hex)) + theme_bw() + scale_colour_manual(values=cols) + ylim(0,255) + theme(legend.position="none") + facet_wrap(~chip, ncol=6)
#x<-p + stat_summary(fun.y="mean", geom="point", size=6, shape=4, col="red") 
#p +  stat_summary(fun.data=mean_se, geom = "errorbar", width=0.2)
p

blk.cc_plot_long[(blk.cc_plot_long$chip == 1) & (blk.cc_plot_long$Value > 100),]
blk.cc_plot_long[(blk.cc_plot_long$chip == 6) & (blk.cc_plot_long$Value < 100),]
blk.cc_plot_long[(blk.cc_plot_long$chip == 24) & (blk.cc_plot_long$Value > 100),]
blk.cc_plot_long[(blk.cc_plot_long$chip == 15) & (blk.cc_plot_long$Value > 100),]

blk.cc_plot_long<-blk.cc_plot_long[!grepl("117_",blk.cc_plot_long$img_name),]
blk.cc_plot_long<-blk.cc_plot_long[!grepl("177_",blk.cc_plot_long$img_name),]



## Lets convert to PC values

## Convert to PC scores
blk.cc_for_pca<-blk.cc_plot[,c(2:4)]

## Change chip number to reflect values on scanner



blk.cc_pca_nz<-blk.cc_for_pca[ , which(apply(blk.cc_for_pca, 2, var) != 0)]
blk.cc_pca<-prcomp(blk.cc_pca_nz, scale = TRUE, center = TRUE)


blk_PC<-as.data.frame(blk.cc_pca$x)
blk_PC$img_name<-blk.cc_plot$img_name
blk_PC$chip<-blk.cc_plot$chip


blk.cc_plot<-merge(blk.cc_plot, blk_PC, by=c('img_name', 'chip'), all=T)



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
                        measure.vars=c("Red", "Green", "Blue"),
                        # Name of the destination column that will identify the original
                        # column that the measurement came from
                        variable.name="Color",
                        value.name="Value"
)


p<-ggplot(blk.cc_plot_long, aes(x=Color, y=Value)) + geom_point(aes(colour = hex)) + theme_bw() + scale_colour_manual(values=cols) + ylim(0,255) + theme(legend.position="none") + facet_wrap(~chip, ncol=6)
#x<-p + stat_summary(fun.y="mean", geom="point", size=6, shape=4, col="red") 
#p +  stat_summary(fun.data=mean_se, geom = "errorbar", width=0.2)
p



## Plot derived PC values
p<-ggplot(blk.cc_plot, aes(x=PC1, y=PC2)) + geom_point(aes(colour = hex)) + theme_bw() + scale_colour_manual(values=cols) + theme(legend.position="none") 


## Lets get summary of variation among each PC for each chip
blk_chip.sd<-aggregate(.~ chip, data=blk.cc_plot[,c(2,6:8)], sd)






#######################################################
## Now lets work with the tuber values
## Start with internals from scanner
## Tuber flesh color

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

p<-ggplot(int_color.PC, aes(x=PC1, y=PC2)) + geom_point(aes(colour = hex), size=0.4) + theme_bw() + scale_colour_manual(values=cols) + theme(legend.position="none")  
q<- p + ylim(-30, 55) +xlim(-75, 50) + xlab(paste("PC1 (", pct.explained[1], "%)", sep="")) + ylab(paste("PC2 (", pct.explained[2], "%)", sep=""))
q<-q+ggtitle("PCA of tuber flesh color")

#pdf("PCA_tuber_flesh_color", height=4, width=4)
print(q)
#dev.off()

## Lets calculate mean and variance by clone


#####################################################
## Tuber skin color

ex_color<-std[,c(1:6,15:785)]
ex_color<-ex_color[ex_color$tuber != "marker",]


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

p<-ggplot(ex_color.PC, aes(x=PC1, y=PC2)) + geom_point(aes(colour = hex), size=0.4) + theme_bw() + scale_colour_manual(values=cols) + theme(legend.position="none")  
q<- p + xlab(paste("PC1 (", PC1.ex.pct,  " %)",  sep="")) + ylab(paste("PC2 (", PC2.ex.pct," %)", sep=""))
q<-q+ggtitle("PCA of tuber external color")

#pdf("PCA_tuber_external_color", height=4, width=4)
print(q)
#dev.off()


## How about H2 of ex_color 
## reformat for h2

ex_color.PC_h2<-ex_color.PC[ex_color.PC$side == 1,c(7:9,1:5)]


# Heritability fxn
get_h2<-function(data){
  
  variance.out<-c()
  H2<-c()
  e2<-c()
  # For each treatment.phenotype calculate variance
  for(i in 4:length(colnames(data))){
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
  colnames(variance)<-colnames(data)[4:length(data)]
  rownames(variance)<-c('Genotype', 'Error')
  return(variance)
}  

h2<-get_h2(ex_color.PC_h2)



## Lets calculate some summarization statistics

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

## How about H2?




## Get hex values
ex_color.PC.ave$Red<-round(ex_color.PC.ave$Red)
ex_color.PC.ave$Green<-round(ex_color.PC.ave$Green)
ex_color.PC.ave$Blue<-round(ex_color.PC.ave$Blue)

ex_color.PC.ave$hex<-rep("NA", nrow(ex_color.PC.ave))

for(rw in 1:nrow(ex_color.PC.ave)){
  ex_color.PC.ave[rw,"hex"]<-rgb(ex_color.PC.ave[rw, 'Red'],ex_color.PC.ave[rw, 'Green'],ex_color.PC.ave[rw, 'Blue'],max=255)
}

cols<-ex_color.PC.ave$hex
names(cols)<-ex_color.PC.ave$hex

## Plot PC2 as boxplot
ex_color.PC.ave$clone = reorder(ex_color.PC.ave$clone, ex_color.PC.ave$PC1, median)

p<-ggplot(ex_color.PC.ave, aes(x=clone, y=PC1)) + geom_point(aes(colour = hex), size=4) + theme_bw() + scale_colour_manual(values=cols) + theme(legend.position="none")  


## Tuber flesh
## How about H2 of ex_color 
## reformat for h2

int_color.PC_h2<-int_color.PC[,c(7:9,1:5)]


# Heritability fxn
get_h2<-function(data){
  
  variance.out<-c()
  H2<-c()
  e2<-c()
  # For each treatment.phenotype calculate variance
  for(i in 4:length(colnames(data))){
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
  colnames(variance)<-colnames(data)[4:length(data)]
  rownames(variance)<-c('Genotype', 'Error')
  return(variance)
}  

h2<-get_h2(int_color.PC_h2)



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

int_color.PC.ave$hex<-rep("NA", nrow(int_color.PC.ave))

for(rw in 1:nrow(int_color.PC.ave)){
  int_color.PC.ave[rw,"hex"]<-rgb(int_color.PC.ave[rw, 'Red'],int_color.PC.ave[rw, 'Green'],int_color.PC.ave[rw, 'Blue'],max=255)
}

cols<-int_color.PC.ave$hex
names(cols)<-int_color.PC.ave$hex

## Plot PC2 as boxplot
int_color.PC.ave$clone = reorder(int_color.PC.ave$clone, int_color.PC.ave$PC2, median)

p<-ggplot(int_color.PC.ave, aes(x=clone, y=PC2)) + geom_point(aes(colour = hex), size=4) + theme_bw() + scale_colour_manual(values=cols) + theme(legend.position="none")  


























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

































int.filepath<-paste(base.dir, "/A08241_internals.csv", sep = "")
int<-read.csv(int.filepath)
int[] <- lapply(int, as.character)
int[int == '']<-c(0)
int[int == 'X']<-c(1)

int$img_name<-paste(int$img_name, int$tuber, sep="_")
int$binary_image<-paste(int$img_name, "binary.jpg", sep="_")
int$masked_image<-paste(int$img_name, "masked.jpg", sep="_")
int<-int[,c(1,16,15,2:14)]
write.csv(int, file="A08241_internals_machine_learning.csv", quote=F, row.names=F)


table(int$sprouting)
table(int$hollow_heart)
table(int$IBS)
table(int$bruise)
table(int$cracks)
table(int$stemend_browning)
table(int$anthocyanin)
table(int$vascular_discoloration)
table(int$greening)
table(int$unknow_defect)

int[int$greening == '1',]












