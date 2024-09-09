

##Set working directory to bring in files##
setwd("/Users/Carly/OneDrive - Smithsonian Institution/Salamanders2020/Peptides 2020/ChallengeAssaysBd-Peptides/#RFilesAnalysis/")


#### LOOK at each plate and make sure no contaminates
#### OR extreme outliers within a sample type
#### remove those here so there are notes of which wells were removed



Bd <- read.csv("PlateReadsAll-SalAMPs.csv")

## Removing outlier wells - some Bd clumping likely occurred. High ODs above 0.5 
## or wells with OD >0.1 deviation from the median

## also had to repeat HOQ-bis-n and used that data from plate 8, instead of plate 5

Bd <- Bd[!(Bd$Plate == 2 & Bd$Column == 2 & Bd$Row == "H"),]
Bd <- Bd[!(Bd$Plate == 2 & Bd$Column == 1 & Bd$Row == "H"),]

Bd <- Bd[!(Bd$Plate == 3 & Bd$Column == 12 & Bd$Row == "D"),]
Bd <- Bd[!(Bd$Plate == 3 & Bd$Column == 1 & Bd$Row == "B"),]
Bd <- Bd[!(Bd$Plate == 3 & Bd$Column == 9 & Bd$Row == "B"),]
Bd <- Bd[!(Bd$Plate == 3 & Bd$Column == 11 & Bd$Row == "C"),]

Bd <- Bd[!(Bd$Plate == 4 & Bd$Column == 3 & Bd$Row == "D"),]
Bd <- Bd[!(Bd$Plate == 4 & Bd$Column == 7 & Bd$Row == "B"),]

## Used rep instead on plate 7
Bd <- Bd[!(Bd$Plate == 5 & Bd$GroupSampleID == "HOQbis-n"),]



## HAVING issues with trying to get this to run with 'skipped' removed here
## Going back to main file and removing skipped
# Bd <- subset(Bd, PeptideConc != "skipped")

#Bd <- read.csv(file = "Plate1-HKcorrected.csv")

Bd$T_OD <- NA


##### NOTE this was not correct in Muletz-Wolz et al. 2017 Frontiers
### Need to put equation again that is correct
### Missing divide by in that paper



# For each plate
# For each row in the isol df 
for(i in 1:nrow(Bd)){
  
  # Take the OD of the row i and divide by the AveND -1 and assign it to Inhibit column in isol df
  Bd$T_OD[i] <-  log((Bd$OD[i] / (1-Bd$OD[i])+1)) 
  
}


## Let's see why we should transform the data this way
hist(Bd$OD)

## This transformation makes the data more linear, better to fit 
## linear models to 
hist(Bd$T_OD)


### Look at positive control growth

PC <- Bd[Bd$GroupSampleID == c("NDPC", "PC"),]

library(ggplot2)


plot12 = ggplot(PC, aes(Day, T_OD))+
                geom_point(aes(color = GroupSampleID)) + geom_smooth(method = loess)

plot12 + facet_grid(~Plate~GroupSampleID, scales = "free")


## Based on these plots, I think we should exclude Days 11 and 12 to fit a LM





## We want to calculate the slope of growth for each well in our assays

## We want to take an average of those slopes for each peptide concentration 
## or controls

# Then we want to divide the average slope of each sample-Bd combo by the 
# nutrient-depleted positive control for the specific Bd isolate
# and then substract one and we get our inhibition score!!!!!!
# My old code I wrote for this also subset within Bd isolates (Muletz-Wolz 2017, Frontiers)
# Just keeping that code as is, but we only have 1 isolate

## First let's merge column and row to make a well column
Well <- as.factor(paste(Bd$Row, Bd$Column, sep = ""))
Bd <- cbind(Bd,Well)

## Subset to days based on PC growth
## Exclude skipped samples
Bd <- subset(Bd, Day %in% c(0,1,4,5,7))

## Had to do some troubleshooting
#Bd <- subset(Bd, Plate %in% c(1,2))
#Bd <- subset(Bd, Plate %in% c(1,3))

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
## This file next isn't necessary with all the reps, just need one day which will now have the slope info
## write.csv(Bd, "InhibitionScoresAllRawSalAMPs.csv")


Bd_scores <- Bd[Bd$Day == "7",]

Bd_scores_2 <- na.omit(Bd_scores)

Bd_scores_2$final_inhibition_percent <- Bd_scores_2$final_inhibition * 100

## write.csv(Bd_scores_2, "InhibitionScoresSalAMPs.csv")

## going in and adding metadata to this output

## Also, Fixed peptide concentration in meta to accurate number based on 1.5 ml dried down from 20 ml quantified

##################


## Double check



Plate_subset <- subset(Bd, Plate == 1) 

well_subsetA4 <- subset(Plate_subset, Well == 'A4')
lm_wellA4 <- lm(I(well_subsetA4$T_OD - intercept) ~ 0 + well_subsetA4$Day)

summary(lm_wellA4)
plot((well_subsetA4$T_OD - intercept) ~ 0 + well_subsetA4$Day)

## Inhibition slope is the coeffient
lm_wellA4$coefficients[1]
Inhibit_slopeA4 <- lm_wellA4$coefficients[1]
summary(lm_wellA4)$adj.r.squared



aveND <- mean(Plate_subset$Inhibit_slope[Plate_subset$GroupSampleID == "NDPC"])



1- (Inhibit_slopeA4/aveND)


## plate 1 well A4 should have inhibit slope as 0.0159, inhibition score as -0.0368 
## and r-squared as 0.9467

well_subsetA4

## We did it

## Let's do one more random well 





Plate_subset <- subset(Bd, Plate == 4) 

well_subsetA4 <- subset(Plate_subset, Well == 'A4')
lm_wellA4 <- lm(I(well_subsetA4$T_OD - intercept) ~ 0 + well_subsetA4$Day)

summary(lm_wellA4)
plot((well_subsetA4$T_OD - intercept) ~ 0 + well_subsetA4$Day)

## Inhibition slope is the coeffient
lm_wellA4$coefficients[1]
Inhibit_slopeA4 <- lm_wellA4$coefficients[1]
summary(lm_wellA4)$adj.r.squared



aveND <- mean(Plate_subset$Inhibit_slope[Plate_subset$GroupSampleID == "NDPC"])



1- (Inhibit_slopeA4/aveND)


## plate 1 well A4 should have inhibit slope as 0.0169, inhibition score as -0.4314 
## and r-squared as 0.9623

well_subsetA4

## Confirmed!


#### FIGURE


# Load packages - need to find a way to use only dplyr (better IMO)
library(ggplot2)
library(plyr)
library(dplyr)
library(gridExtra)

# Some colors I like - don't need to run this!
colors <- c("#7CCD7C", "#FFA500", "#A4D3EE")

colors2 <- c("lightcyan2","palegreen2", "#0072B2", "dodgerblue")

# Inhibition scores w/ added metadata
# Fixed peptide concentration in meta to accurate number based on 1.5 ml dried down from 20 ml quantified
# Also, need to remove one sample's concentration that believe was an error 
## And couple data points that needs to be removed that missed above with outlier OD readings


inhibit_scores <- read.csv("InhibitionScoresSalAMPs_meta.csv")
inhibit_scores <- inhibit_scores[!(inhibit_scores$PeptideConcT == 37.5 & inhibit_scores$GroupSampleID == "DMPcin-n"),]
inhibit_scores <- inhibit_scores[!(inhibit_scores$Plate == 2 & inhibit_scores$Column == 7 &inhibit_scores$Row == "A"),]
inhibit_scores <- inhibit_scores[!(inhibit_scores$Plate == 5 & inhibit_scores$Column == 10 &inhibit_scores$Row == "E"),]




inhibit_scores$PeptideConcT <- as.factor(inhibit_scores$PeptideConcT)
str(inhibit_scores)



# Subset data - calculate mean and std err,
# there's got to be a more efficient way than this?
meanInhib <- ddply(inhibit_scores, .(GroupSampleID3,
                                     PeptideConcT,
                                     Species, 
                                     Injection),
                   summarize, mean = mean(final_inhibition_percent))

standard_errors <- ddply(inhibit_scores, .(GroupSampleID3,
                                           PeptideConcT),
                         summarize,
                         stderr = sd(final_inhibition_percent) / 
                           sqrt(length(final_inhibition_percent)))

# Merge standard errors with mean inhibition data, I'll keep the same name
meanInhib <- merge(meanInhib, standard_errors,
                   by = c("GroupSampleID3","PeptideConcT"))


# Using dplyr function, subset mean inhibition scores by species
bis <- filter(meanInhib, Species == "E. bislineata")
cin <- filter(meanInhib, Species == "P. cinereus")
vir <- filter(meanInhib, Species == "N. viridescens")

# Set theme for plots
theme_set(theme_classic() +
            theme(panel.spacing.x = unit(2, "lines")) +
            theme(axis.text.x = element_text(vjust = 0, size = 10)) +
            theme(axis.title.x = element_text(vjust = -1))) 
#+theme(legend.position = "bottom"))



####### TO DO - make a groupsampleID3 with just Site and No or Yes, add pictures of salamanders
#### Fix the Ach injection at bottom, just remove all together and put in figure legend


## Plots by species ## - Need to find way to run this as a function - efficient

# Bislineata
(bis_scatter <- ggplot(bis, aes(x = PeptideConcT, y = mean, shape = Injection)) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = mean - stderr,
                      ymax = mean + stderr),
                  width = 0.25) +
    facet_wrap(~ GroupSampleID3) +
    labs(y = "Mean Bd inhibiton score", x = "Concentration (ug/mL)") +
    scale_shape_discrete(name = "Ach Injection",
                         breaks = c("n", "y"),
                         labels = c("No", "Yes")) +
    geom_hline(yintercept = 0,
               linetype = "dashed",
               linewidth = 0.5, color = "lightgrey") +
    #theme(legend.position = "none") +
    theme(strip.text = element_text(face = "bold",
                                    size = rel(1)),
          strip.background = element_rect(fill = "lightcyan2",
                                          colour = "black", linewidth = 1))+ ylim(-55,50))
# Cinereus
(cin_scatter <- ggplot(cin, aes(x = PeptideConcT, y = mean, shape = Injection)) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = mean - stderr,
                      ymax = mean + stderr),
                  width = 0.25) +
    facet_wrap(~ GroupSampleID3) +
    labs(y = "Mean Bd inhibiton score", x = "Concentration (ug/mL)") +
    scale_shape_discrete(name = "Ach Injection",
                         breaks = c("n", "y"),
                         labels = c("No", "Yes")) +
    geom_hline(yintercept = 0,
               linetype = "dashed",
               linewidth = 0.5, color = "lightgrey") +
    #theme(legend.position = "none") +
    theme(strip.text = element_text(face = "bold",
                                    size = rel(1)),
          strip.background = element_rect(fill = "palegreen2",
                                          colour = "black", linewidth = 1))+ ylim(-55,50))

# Viridescens
(vir_scatter <- ggplot(vir, aes(x = PeptideConcT, y = mean, shape = Injection)) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = mean - stderr,
                      ymax = mean + stderr),
                  width = 0.25) +
    facet_wrap(~ GroupSampleID3) +
    labs(y = "Mean Bd inhibiton score", x = "Concentration (ug/mL)") +
    scale_shape_discrete(name = "Ach Injection",
                         breaks = c("n", "y"),
                         labels = c("No", "Yes")) +
    geom_hline(yintercept = 0,
               linetype = "dashed",
               linewidth = 0.5, color = "lightgrey") +
    theme(strip.text = element_text(face = "bold",
                                    size = rel(1)),
          strip.background = element_rect(fill = "steelblue1",
                                          colour = "black", linewidth = 1))+ ylim(-55,50))

# In a panel - still a bit cluttered
(panel <-  grid.arrange(bis_scatter + 
                          ggtitle("(A)") + 
                          theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2),
                                                   units = , "cm")),
                        cin_scatter +
                          ggtitle("(B)") + 
                          theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2),
                                                   units = , "cm")),
                        vir_scatter + ggtitle("(C)") +
                          theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2),
                                                   units = , "cm")) +
                          theme(legend.text = element_text(size = 10,
                                                           face = "italic"),
                                legend.position = "bottom")))














