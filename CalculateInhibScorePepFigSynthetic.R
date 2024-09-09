##Set working directory to bring in files##
setwd("/Users/Julian/Library/CloudStorage/
      OneDrive-SmithsonianInstitution/Peptides 2020/
      SynthethicPeptide_assay/#Rfiles_Analysis")


setwd("/Users/Carly/OneDrive - Smithsonian Institution/Salamanders2020/Peptides 2020/SynthethicPeptide_assay/#Rfiles_Analysis")
# Load packages
library(dplyr)
library(ggplot2)
library(plyr)

Bd <- read.csv("PlateReadsAll_long.csv")

## Removing outlier wells w/ dplyr filter function
# wells with visible contamination, extremely high readings (>0.1 deviation from well replicates), or read errors
Bd <- Bd %>% filter(!(Plate == 1 & Column == 1 & Row == "G"))  
Bd <- Bd %>% filter(!(Plate == 1 & Column == 1 & Row == "A"))
Bd <- Bd %>% filter(!(Plate == 1 & Column == 12 & Row == "B"))

Bd <- Bd %>% filter(!(Plate == 2 & Column == 9 & Row == "C")) 
Bd <- Bd %>% filter(!(Plate == 3 & Column == 4 & Row == "A"))

Bd <- Bd %>% filter(!(Plate == 4 & Column %in% 3:4 & Row == "H")) 

Bd <- Bd %>% filter(!(Plate == 7 & Column == 9 & Row == "B"))


# Add columns for Transformed OD readings
Bd$T_OD <- NA


# For each plate
# For each row in the isol df 
for(i in 1:nrow(Bd)){
  
  # Take the OD of the row i and transform
  Bd$T_OD[i] <-  log((Bd$OD[i] / (1-Bd$OD[i])+1)) 
  
}



## Look at positive control growth
PC <- Bd[Bd$GroupSampleID == c("NDPC", "PC"),]

(plot12 <-  ggplot(PC, aes(Day, T_OD)) +
    geom_point(aes(color = GroupSampleID)) + 
    geom_smooth(method = loess) +
    facet_grid(~Plate~GroupSampleID, scales = "free")
)


## We want to calculate the slope of growth for each well in our assays

## We want to take an average of those slopes for each peptide concentration 
## or controls

# Then we want to divide the average slope of each sample-Bd combo by the 
# nutrient-depleted positive control for the specific Bd isolate
# and then subtract one and we get our inhibition score!!!!!!
# My old code I wrote for this also subset within Bd isolates (Muletz-Wolz 2017, Frontiers)
# Just keeping that code as is, but we only have 1 isolate

## First let's merge column and row to make a well column
Well <- as.factor(paste(Bd$Row, Bd$Column, sep = ""))
Bd <- cbind(Bd,Well)

## Subset to days based on PC growth
Bd <- subset(Bd, Day %in% c(0,1,4,5,7))


## Make a column to put in calculated slopes
Bd$Inhibit_slope <- NA

## specify unique plates
plate_unique <- unique(Bd$Plate)

## specify our intercept as 0
intercept = 0  


# For each plate
for(i in 1:length(plate_unique)){
  
  Plate_subset <- subset(Bd, Plate == plate_unique[i]) 
  
  iso <- unique(Plate_subset$Isolate)
  
  # For each of unique isolates 
  for(j in 1:length(iso)){
    
    # Subset to only that iso i and assign it to a new dataframe isol
    isol <- subset(Plate_subset, Isolate==iso[j])
    
    well <- unique(isol$Well)
    
    # For each unique well
    for(k in 1:length(well)){
      
      # Subset to that well and assign it to a new dataframe well_subset
      well_subset <- subset(isol, Well == well[k])
      
      # Make a linear model for each bacteria based on the transformed OD
      
      lm_well <- lm(I(well_subset$T_OD - intercept) ~ 0 + well_subset$Day)
      
      ## then take the slope and put into appropriate column in isol df
      ## if we didn't force the intercept slope is in [2]
      ## well_subset$Inhibit_slope <- lm_well$coefficients[2]
      
      well_subset$Inhibit_slope <- lm_well$coefficients[1]
      well_subset$R_squared <- summary(lm_well)$adj.r.squared
      well_subset$R_squared <- ifelse(well_subset$R_squared < 0.2, NA, well_subset$R_squared) 
      well_subset$Inhibit_slope <- ifelse(well_subset$R_squared < 0.2, NA, well_subset$Inhibit_slope) 
      
      # Once that is all completed for that unique isolate, assign it into the original df 
      isol[isol$Well==well[k],"Inhibit_slope"] <- well_subset$Inhibit_slope
      isol[isol$Well==well[k],"R_squared"] <- well_subset$R_squared
      
      
      ## get the average of ND inhibit_slope
      aveND <- mean(isol$Inhibit_slope[isol$GroupSampleID == "NDPC"])
      
      ## for each row in the isol df
      for(l in 1:nrow(isol)){
        
        # Take the inhibit_Slope of the row l and divide by the Avg ND slope and then -1, 
        # Assign to new column called final_inhibition
        isol$final_inhibition[l] <- 1- (isol$Inhibit_slope[l]/aveND)
      }
      
    }
    
    # Once that is all completed for that unique isolate, assign it into the original df 
    Plate_subset[Plate_subset$Isolate==iso[j],"Inhibit_slope"] <- isol$Inhibit_slope
    Plate_subset[Plate_subset$Isolate==iso[j],"final_inhibition"] <- isol$final_inhibition
    Plate_subset[Plate_subset$Isolate==iso[j],"R_squared"] <- isol$R_squared
    
    
  }
  
  Bd[Bd$Plate==plate_unique[i],"Inhibit_slope"] <- Plate_subset$Inhibit_slope
  Bd[Bd$Plate==plate_unique[i],"final_inhibition"] <- Plate_subset$final_inhibition
  Bd[Bd$Plate==plate_unique[i],"R_squared"] <- Plate_subset$R_squared
  
}



## we have replication of the final_inhibition for every day, need to subset only one day

#write.csv(Bd, "InhibitionScoresRaw_Synthetic.csv")

Bd_scores <- Bd[Bd$Day == "7",]

Bd_scores_2 <- na.omit(Bd_scores)

Bd_scores_2$final_inhibition_percent <- Bd_scores_2$final_inhibition * 100

write.csv(Bd_scores_2, "InhibitionScores_Synthetic.csv")


##################


## Double check

Plate_subset <- subset(Bd, Plate == 3) 

well_subsetD5 <- subset(Plate_subset, Well == 'D5')
lm_wellD5 <- lm(I(well_subsetD5$T_OD - intercept) ~ 0 + well_subsetD5$Day)

summary(lm_wellD5)
plot((well_subsetD5$T_OD - intercept) ~ 0 + well_subsetD5$Day)

## Inhibition slope is the coefficient
lm_wellD5$coefficients[1]
Inhibit_slopeD5 <- lm_wellD5$coefficients[1]
summary(lm_wellD5)$adj.r.squared



aveND <- mean(Plate_subset$Inhibit_slope[Plate_subset$GroupSampleID == "NDPC"])



1- (Inhibit_slopeD5/aveND)


## plate 3 well D5 should have inhibit slope as 0.0194 , inhibition score as -0.9187
## and r-squared as 0.998

well_subsetD5


Plate_subset <- subset(Bd, Plate == 4) 

well_subsetA4 <- subset(Plate_subset, Well == 'A4')
lm_wellA4 <- lm(I(well_subsetA4$T_OD - intercept) ~ 0 + well_subsetA4$Day)

summary(lm_wellA4)
plot((well_subsetA4$T_OD - intercept) ~ 0 + well_subsetA4$Day)

## Inhibition slope is the coefficient
lm_wellA4$coefficients[1]
Inhibit_slopeA4 <- lm_wellA4$coefficients[1]
summary(lm_wellA4)$adj.r.squared



aveND <- mean(Plate_subset$Inhibit_slope[Plate_subset$GroupSampleID %in% "NDPC"], na.rm = TRUE)



1- (Inhibit_slopeA4/aveND)


well_subsetA4



#######################        FIGURE        #########################


library(dplyr)
library(ggplot2)
library(plyr)

setwd("/Users/Carly/OneDrive - Smithsonian Institution/Salamanders2020/Peptides 2020/SynthethicPeptide_assay/#Rfiles_Analysis/")




# Some colors I like - don't need to run this!
colors <- c("#7CCD7C", "#FFA500", "#A4D3EE")

colors2 <- c("lightcyan2","palegreen2", "#0072B2", "dodgerblue")

# Inhibition scores w/ added metadata
inhibit_scores <- read.csv("InhibitionScores_Synthetic.csv")
inhibit_scores$PeptideConc <- as.factor(inhibit_scores$PeptideConc)
str(inhibit_scores)

# Subset data - calculate mean and std err,
# there's got to be a more efficient way than this?
meanInhib <- ddply(inhibit_scores, .(GroupSampleID,
                                     PeptideConc),
                   summarize, mean = mean(final_inhibition_percent))

standard_errors <- ddply(inhibit_scores, .(GroupSampleID,
                                           PeptideConc),
                         summarize,
                         stderr = sd(final_inhibition_percent) / 
                           sqrt(length(final_inhibition_percent)))



# Merge standard errors with mean inhibition data, I'll keep the same name
meanInhib <- merge(meanInhib, standard_errors,
                   by = c("GroupSampleID","PeptideConc"))

#write.csv(meanInhib, "SyntheticPeptideInhibition.csv")


### NOTES from Julian

## Tr-Eb4, Tr-Pc5, and Tr-Pc15 were the three peptides that were cloudy at 500 and/or 250 mcg (along with LL-37).
## Tr-Eb1, Tr-Eb4, Tr-Pc15 were recommended to be dissolved in DMSO, 
## however, I didn't notice any cloudiness in Tr-Eb1(according to my notes) but can see it now on the figure.
# Filter out concentrations that were insoluble

## Carly Notes

## PrMb-Eb14 had to have come out of solution and Tr-Eb1 at the lowest concentration oddly. Higher Bd growth than possible
## Remove concentrations where solubility was a problem

meanInhib_soluble <- meanInhib %>% filter(!(GroupSampleID %in% c("LL-37", "Tr-Eb4", "Tr-Pc5", "Tr-Pc15") & PeptideConc %in% c(500, 250)))
meanInhib_soluble <- meanInhib_soluble %>% filter(!(GroupSampleID == "PrMb-Eb14"))
meanInhib_soluble <- meanInhib_soluble %>% filter(!(GroupSampleID == "Tr-Eb1" & PeptideConc == "15.625"))
         
# Set theme for plots
         theme_set(theme_classic() +
                     theme(panel.spacing.x = unit(2, "lines")) +
                     theme(axis.text.x = element_text(vjust = 0, size = 8)) +
                     theme(axis.title.x = element_text(vjust = -1)) +
                     theme(legend.position = "bottom"))
         
         # Scatter for Inhib Scores
         (synth_scatter <- ggplot(meanInhib, aes(x = PeptideConc, y = mean)) +
             geom_point(size = 2.0) +
             geom_errorbar(aes(ymin = mean - stderr,
                               ymax = mean + stderr),
                           width = 0.25) +
             facet_wrap(~ GroupSampleID) +
             labs(y = "Mean Bd inhibiton score", x = "Concentration (ug/mL)") +
             geom_hline(yintercept = 0,
                        linetype = "dashed",
                        linewidth = 0.75, color = "lightgrey") +
             #theme(legend.position = "none") +
             theme(strip.text = element_text(face = "bold",
                                             size = rel(1)),
                   strip.background = element_rect(fill = "lightgrey",
                                                   colour = "black", linewidth = 1)) + ylim(-500,10))
         
#### FINAL FOR soluble
         (synth_scatter2 <- ggplot(meanInhib_soluble, aes(x = PeptideConc, y = mean)) +
             geom_point(size = 2.0) +
             geom_errorbar(aes(ymin = mean - stderr,
                               ymax = mean + stderr),
                           width = 0.25) +
             facet_wrap(~ GroupSampleID) +
             labs(y = "Mean Bd inhibiton score", x = "Concentration (ug/mL)") +
             geom_hline(yintercept = 0,
                        linetype = "dashed",
                        linewidth = 0.5, color = "lightgrey") +
             #theme(legend.position = "none") +
             theme(strip.text = element_text(face = "bold",
                                             size = rel(1)),
                   strip.background = element_rect(fill = "lightgrey",
                                                   colour = "black", linewidth = 1)))+ 
           theme(axis.text.x = element_text(angle=45, vjust=0.6)) 
         





