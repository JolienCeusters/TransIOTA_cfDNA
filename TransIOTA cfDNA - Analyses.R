#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
##  TransIOTA ctDNA - final analyses  ##
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

## Load packages
library(data.table)
library(plyr)
library(ggplot2)
library(mice)
library(doBy)
library(ltm)
library(blandr)
library(auRoc)
library(logistf)
library(rms)
library(SDMTools)
library(pROC)
library(doParallel)


#### 1. Load data ####
load("TransIOTA - analyses ctDNA.RData")


#### 2. Create variables ####

## Make ADNEX groups
HistSpecial <- read_excel_allsheets("Histologies TransIOTA - Chiara_CL.xlsx")
HistSpecialLong <- rbindlist(HistSpecial, fill = TRUE)
IDsHis <- HistSpecialLong$`__metadata__.__caseId__`

ctDNAanalyses = ddply(ctDNAanalyses, .(`Patient ID`),
                      function(x){
                        ADNEXgroups =
                          if(!is.na(x$`FOR JOLIEN`)){
                            if(x$`FOR JOLIEN` == "Benign"){
                              "Benign"
                            } # End: benign
                            else{
                              if(x$`FOR JOLIEN` == "Borderline"){
                                "Borderline"
                              } # Borderline
                              else{
                                if(x$`FOR JOLIEN` == "Metastasis"){
                                  "Metastatic"
                                } # Metastatic
                                else{
                                  if(x$`FOR JOLIEN` == "Stage I"){
                                    "Stage I invasive"
                                  } # Stage I
                                  else{
                                    if(x$`FOR JOLIEN` == "Stage II-IV"){
                                      "Stage II-IV invasive"
                                    } # Stage II - IV
                                    else{
                                      "Unknown"
                                    } # End: not Stage II - IV
                                  } # End: not Stage I
                                } # End: not metastatic
                              } # End: not borderline
                            } # End: not benign
                          } # End: if For Jolien is not missing
                        else{
                          if(x$`Patient ID` %in% IDsHis){
                            HistSpecialLong$Group[HistSpecialLong$`__metadata__.__caseId__` == x$`Patient ID`]
                          } # If in list Chiara
                          else{
                            if(x$`Mass Outcome` == "benign"){
                              "Benign"
                            } # Benign
                            else{
                              if(x$`Mass Outcome` == "borderline"){
                                "Borderline"
                              } # Borderline
                              else{
                                if(x$`Mass Outcome` == "infectious_acute_chronic"){
                                  "Benign"
                                } # Infectious acute chronic
                                else{
                                  if(x$`Mass Outcome` == "malignant" | x$`Mass Outcome` == "invasive malignant"){
                                    if(x$`FIGO stage` == 1 | x$`FIGO stage` == "I"){
                                      "Stage I invasive"
                                    } # FIGO stage I
                                    else{
                                      if(x$`FIGO stage` == 2 | x$`FIGO stage` == 3 | x$`FIGO stage` == 4 | x$`FIGO stage` == "II" | x$`FIGO stage` == "III" | x$`FIGO stage` == "IV"){
                                        "Stage II-IV invasive"
                                      } # FIGO stage II or III or IV
                                      else{
                                        "Invasive FIGO stage unknown"
                                      } # Not: FIGO stage II or III or IV
                                    } # Not: FIGO stage I
                                  } # Malignant
                                  else{
                                    if(x$`Mass Outcome` == "metastatic"){
                                      "Metastatic"
                                    } # Metastatic
                                    else{
                                      if(x$`Mass Outcome` == "rare_benign"){
                                        "Benign"
                                      } # Rare benign
                                      else{
                                        "Unknown"
                                      } # Not: rare benign
                                    } # Not: metastatic
                                  } # Not: Malignant
                                } # Not: infectious acute chronic
                              } # Not: borderline
                            } # Not: benign
                          } # Not: in list Chiara
                        } # For Jolien is missing
                        
                        x$ADNEXgroups = ADNEXgroups
                        return(x)
                      } # End function
)
table(ctDNAanalyses$ADNEXgroups)

ctDNAanalyses$ADNEXgroups = factor(ctDNAanalyses$ADNEXgroups, levels = c("Benign", "Borderline", "Stage I invasive", "Stage II-IV invasive", "Metastatic"))

## Binary outcome
ctDNAanalyses$Outcome <- ifelse(ctDNAanalyses$ADNEXgroups == "Benign", "Benign", "Malignant")
table(ctDNAanalyses$Outcome)

ctDNAanalyses$OutcomeBin <- ifelse(ctDNAanalyses$Outcome == "Benign", 0, 1)
table(ctDNAanalyses$OutcomeBin)

## Centercodes
ctDNAanalyses$CenterCode <- ifelse(grepl("Prague", ctDNAanalyses$center), "Prague",
                                   ifelse(grepl("Rome", ctDNAanalyses$center), "Rome", "Leuven"))

## Numeric variable for CA125
ctDNAanalyses$CA125 <- as.numeric(ctDNAanalyses$`CA125 (kU/L)`)


#### 3. Descriptive statistics ####

#### 3.1 Inclusion plot ####
Trans.ds1 <- Trans[, c("__metadata__.__caseId__", "iota7_basic.scan_date")]
colnames(Trans.ds1) <- c("Patient ID", "Exam date")
Trans.ds1$`Exam date` <- as.Date(Trans.ds1$`Exam date`, "%d/%m/%Y")

IOTA5.ds1 <- ds1[, c("Patient ID", "Exam date")]
IOTA5.ds1$`Exam date` <- as.Date(IOTA5.ds1$`Exam date`, "%Y-%m-%d")

Dates <- rbind(IOTA5.ds1, Trans.ds1)

ctDNAinclusion <- merge(ctDNAanalyses, Dates, all.x = TRUE)

ggplot(ctDNAinclusion, aes(x = `Exam date`, y = CenterCode)) +
  geom_point() +
  labs(x = "Exam date", y = "Center")
# 800 x 300


#### 3.2 Descriptive table with patient characteristics and echo variables ####

## Age/Menopausal status
summary(ctDNAanalyses$`Patient age`)
table(ctDNAanalyses$Postmenopausal2)

## Maximal diameter of lesion
summary(ctDNAanalyses$`Lesion largest diameter`)

## Largest diameter of largest solid component (Only for tumors with a solid component)
summary(ctDNAanalyses$soldmax[ctDNAanalyses$solidbin == 1])

## Proportion solid tissue/Presence of solid components
table(ctDNAanalyses$solidbin)

## > 10 cyst locules
table(ctDNAanalyses$loc10)

## Number of papillary projections
table(ctDNAanalyses$papnr)

## Acoustic shadows
table(ctDNAanalyses$shadowsbin)

## Ascites
table(ctDNAanalyses$ascitesbin)

## Metastases
table(ctDNAanalyses$metasbin)

## Bilateral masses
table(ctDNAanalyses$bilatbin)

## Observed CA125
summary(ctDNAanalyses$CA125)

## Ultrasound examiner's subjective impression
table(ctDNAanalyses$Certainty)

## Outcome
table(ctDNAanalyses$ADNEXgroups)
table(ctDNAanalyses$ADNEXgroups, ctDNAanalyses$CenterCode)
table(ctDNAanalyses$Histology, ctDNAanalyses$ADNEXgroups)
ctDNAanalyses[is.na(ctDNAanalyses$Histology),]

## Method 
nrow(ctDNAanalyses[is.na(ctDNAanalyses$GC_Adapter) & !is.na(ctDNAanalyses$gw_z_score),]) # Without GC-adapter (Diether)
table(ctDNAanalyses$CenterCode[is.na(ctDNAanalyses$GC_Adapter) & !is.na(ctDNAanalyses$gw_z_score)])

nrow(ctDNAanalyses[!is.na(ctDNAanalyses$GC_Adapter) & !is.na(ctDNAanalyses$gw_z_score),]) # With GC-adapter (Joris)
table(ctDNAanalyses$CenterCode[!is.na(ctDNAanalyses$GC_Adapter)])

nrow(ctDNAanalyses[is.na(ctDNAanalyses$GC_Adapter) & is.na(ctDNAanalyses$gw_z_score),]) # Samples with no result
table(ctDNAanalyses$CenterCode[is.na(ctDNAanalyses$GC_Adapter) & is.na(ctDNAanalyses$gw_z_score)])

nrow(TransIOTAanalyses[duplicated(TransIOTAanalyses$`Patient ID`, fromLast = T) & !TransIOTAanalyses$`Patient ID` %in% c(15892, 17445),])
table(TransIOTAanalyses$center[duplicated(TransIOTAanalyses$`Patient ID`, fromLast = T) & !TransIOTAanalyses$`Patient ID` %in% c(15892, 17445)]) # Samples analysed twice

## Percentage missing data
summary(ctDNAanalyses$gw_z_score)
summary(ctDNAanalyses$nucleosome_score)
summary(ctDNAanalyses$CA125)


#### 3.3 Summary measures ####

#### 3.3.1 Genome-wide z-score ####

## Overall
summary(ctDNAanalyses$gw_z_score)

box <- ggplot(ctDNAanalyses, aes(y = gw_z_score, x = "")) +
       geom_jitter(width = 0.20, size = 1.5, alpha = 0.3) + # width is the amount of jitter and alpha makes the points more transparant
       geom_boxplot(outlier.shape = NA, colour = "black", fill = NA) # Hiding the outliers
box <- box + labs(y = "Genome-wide z-score", x = "") # Labels for x-axis and y-axis
box <- box + theme_minimal() # White background
box 
# 600 x 800

hist <- ggplot(ctDNAanalyses, aes(x = gw_z_score)) +
        geom_histogram(aes(y = ..count..), alpha=.75) + 
        theme_minimal() +
        labs(x = "Genome wide z-score",
             y = "Frequency")
hist # 800 x 600

## Per outcome
tapply(ctDNAanalyses$gw_z_score, ctDNAanalyses$ADNEXgroups, summary)

library(scales)
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

library(ggallin)
box <- ggplot(ctDNAanalyses, aes(y = gw_z_score, x = ADNEXgroups)) +
       geom_jitter(width = 0.20, size = 1.5, alpha = 0.3) + # width is the amount of jitter and alpha makes the points more transparant
       geom_boxplot(outlier.shape = NA, colour = "black", fill = NA) # Hiding the outliers
box <- box + labs(y = "Genome-wide z-score", x = "") # Labels for x-axis and y-axis
box <- box + theme_minimal() # White background
box <- box + coord_flip()
box
box + scale_y_continuous(trans = pseudolog10_trans, breaks = c(-10, -1, 0, 1, 10, 1000), minor_breaks = NULL, labels = c(-10, -1, 0, 1, 10, 1000))#,
# 800 x 400

box <- ggplot(ctDNAanalyses, aes(y = gw_z_score, x = Outcome)) +
       geom_jitter(width = 0.20, size = 1.5, alpha = 0.3) + # width is the amount of jitter and alpha makes the points more transparant
       geom_boxplot(outlier.shape = NA, colour = "black", fill = NA) # Hiding the outliers
box <- box + labs(y = "Genome-wide z-score", x = "") # Labels for x-axis and y-axis
box <- box + theme_minimal() # White background
box

hist <- ggplot(ctDNAanalyses, aes(x = gw_z_score)) +
        geom_histogram(data=subset(ctDNAanalyses, Outcome == "Benign"),aes(y = ..count.., fill = "Benign"), alpha=.75) + 
        geom_histogram(data=subset(ctDNAanalyses, Outcome == "Malignant"),aes(y = -..count.., fill = "Malignant"), alpha=.75) + 
        theme_minimal() +
        labs(x = "Genome wide z-score",
             y = "Frequency") +
        scale_y_continuous(limits = c(-220, 260), breaks = c(-200, -100, 0, 100, 200), labels = c(200, 100, 0, 100, 200))
hist # 1000 x 600


#### 3.3.2 Nucleosome score ####

## Overall
summary(ctDNAanalyses$nucleosome_score)

box <- ggplot(ctDNAanalyses, aes(y = nucleosome_score, x = "")) +
       geom_jitter(width = 0.20, size = 1.5, alpha = 0.3) + # width is the amount of jitter and alpha makes the points more transparant
       geom_boxplot(outlier.shape = NA, colour = "black", fill = NA) # Hiding the outliers
box <- box + labs(y = "Nucleosome score", x = "") # Labels for x-axis and y-axis
box <- box + theme_minimal() # White background
box # 600 x 800

hist <- ggplot(ctDNAanalyses, aes(x = nucleosome_score)) +
        geom_histogram(aes(y = ..count..), binwidth = 0.01, alpha=.75) +
        theme_minimal() +
        labs(x = "Nucleosome score",
             y = "Frequency") 
hist # 800 x 600

## Per outcome
tapply(ctDNAanalyses$nucleosome_score, ctDNAanalyses$ADNEXgroups, summary)

box <- ggplot(ctDNAanalyses, aes(y = nucleosome_score, x = ADNEXgroups)) +
       geom_jitter(width = 0.20, size = 1.5, alpha = 0.3) + # width is the amount of jitter and alpha makes the points more transparant
       geom_boxplot(outlier.shape = NA, colour = "black", fill = NA) # Hiding the outliers
box <- box + labs(y = "Nucleosome score", x = "") # Labels for x-axis and y-axis
box <- box + theme_minimal() # White background
box <- box + coord_flip()
box # 800 x 400

box <- ggplot(ctDNAanalyses, aes(y = nucleosome_score, x = Outcome)) +
       geom_jitter(width = 0.20, size = 1.5, alpha = 0.3) + # width is the amount of jitter and alpha makes the points more transparant
       geom_boxplot(outlier.shape = NA, colour = "black", fill = NA) # Hiding the outliers
box <- box + labs(y = "Nucleosome score", x = "") # Labels for x-axis and y-axis
box <- box + theme_minimal() # White background
box # 800 x 800

hist <- ggplot(ctDNAanalyses, aes(x = nucleosome_score)) +
        geom_histogram(data=subset(ctDNAanalyses, Outcome == "Benign"),aes(y = ..count.., fill = "Benign"), binwidth = 0.01, alpha=.75) +
        geom_histogram(data=subset(ctDNAanalyses, Outcome == "Malignant"),aes(y = -..count.., fill = "Malignant"), binwidth = 0.01, alpha=.75) +
        theme_minimal() +
        labs(x = "Nucleosome score",
             y = "Frequency") +
        scale_y_continuous(limits = c(-80, 100), breaks = c(-50, 0, 50, 100), labels = c(50, 0, 50, 100))
hist # 1000 x 600


#### 3.4 Correlations ####

#### 3.4.1 Correlation between cfDNA scores ####

cor(ctDNAanalyses[, c("gw_z_score", "nucleosome_score")], method = "spearman", use = "pairwise.complete.obs") # Calculate the correlation


#### 3.4.2 Correlation between ADNEX-variables ####

## Age vs maximum diameter of lesion: Spearman
round(cor(ctDNAanalyses$`Patient age`, ctDNAanalyses$`Lesion largest diameter`, method = "spearman", use="pairwise.complete.obs"), 2)

## Proportion solid tissue vs age: Spearman
round(cor(ctDNAanalyses$propsol, ctDNAanalyses$`Patient age`, method = "spearman", use="pairwise.complete.obs"), 2)

## Proportion solid tissue vs maximum diameter lesion: Spearman
round(cor(ctDNAanalyses$propsol, ctDNAanalyses$`Lesion largest diameter`, method = "spearman", use="pairwise.complete.obs"), 2)

## > 10 locules vs age: Point biserial correlation
round(biserial.cor(ctDNAanalyses$`Patient age`, ctDNAanalyses$loc10, use = "complete.obs"), 2)

## > 10 locules vs maximum diameter lesion: Point biserial correlation
round(biserial.cor(ctDNAanalyses$`Lesion largest diameter`, ctDNAanalyses$loc10, use = "complete.obs"), 2)

## > 10 locules vs proportion solid tissue: Point biserial correlation
round(biserial.cor(ctDNAanalyses$propsol, ctDNAanalyses$loc10, use = "complete.obs"), 2)

## Number papillary projections vs age: Spearman
round(cor(ctDNAanalyses$papnr, ctDNAanalyses$`Patient age`, method = "spearman", use="pairwise.complete.obs"), 2)

## Number papillary projections vs maximum diameter lesion: Spearman
round(cor(ctDNAanalyses$papnr, ctDNAanalyses$`Lesion largest diameter`, method = "spearman", use="pairwise.complete.obs"), 2)

## Number papillary projections vs proportion solid tissue: Spearman
round(cor(ctDNAanalyses$papnr, ctDNAanalyses$propsol, method = "spearman", use="pairwise.complete.obs"), 2)

## Number papillary projections vs > 10 locules: Point biserial correlation
round(biserial.cor(ctDNAanalyses$papnr, ctDNAanalyses$loc10, use = "complete.obs"), 2)

## Acoustic shadows vs age: Point biserial correlation
round(biserial.cor(ctDNAanalyses$`Patient age`, ctDNAanalyses$shadowsbin, use = "complete.obs"), 2)

## Acoustic shadows vs maximum diameter lesion: Point biserial correlation
round(biserial.cor(ctDNAanalyses$`Lesion largest diameter`, ctDNAanalyses$shadowsbin, use = "complete.obs"), 2)

## Acoustic shadows vs proportion solid tissue: Point biserial correlation
round(biserial.cor(ctDNAanalyses$propsol, ctDNAanalyses$shadowsbin, use = "complete.obs"), 2)

## Acoustic shadows vs > 10 locules: phi correlation
round(cor(ctDNAanalyses$shadowsbin, ctDNAanalyses$loc10, method = "pearson", use="pairwise.complete.obs"), 2)

## Acoustic shadows vs number papillary projections: Point biserial correlation
round(biserial.cor(ctDNAanalyses$papnr, ctDNAanalyses$shadowsbin, use = "complete.obs"), 2)

## Ascites vs age: Point biserial correlation
round(biserial.cor(ctDNAanalyses$`Patient age`, ctDNAanalyses$ascitesbin, use = "complete.obs"), 2)

## Ascites vs maximum diameter lesion: Point biserial correlation
round(biserial.cor(ctDNAanalyses$`Lesion largest diameter`, ctDNAanalyses$ascitesbin, use = "complete.obs"), 2)

## Ascites vs proportion solid tissue: Point biserial correlation
round(biserial.cor(ctDNAanalyses$propsol, ctDNAanalyses$ascitesbin, use = "complete.obs"), 2)

## Ascites vs > 10 locules: phi correlation
round(cor(ctDNAanalyses$ascitesbin, ctDNAanalyses$loc10, method = "pearson", use="pairwise.complete.obs"), 2)

## Ascites vs number papillary projections: Point biserial correlation
round(biserial.cor(ctDNAanalyses$papnr, ctDNAanalyses$ascitesbin, use = "complete.obs"), 2)

## Ascites vs acoustic shadows: phi correlation
round(cor(ctDNAanalyses$ascitesbin, ctDNAanalyses$shadowsbin, method = "pearson", use="pairwise.complete.obs"), 2)

## CA125 vs Age: Spearman
round(cor(ctDNAanalyses$CA125, ctDNAanalyses$`Patient age`, method = "spearman", use="pairwise.complete.obs"), 2)

## CA125 vs Maximum diameter of lesion: Spearman
round(cor(ctDNAanalyses$CA125, ctDNAanalyses$`Lesion largest diameter`, method = "spearman", use="pairwise.complete.obs"), 2)

## CA125 vs Proportion of solid tissue: Spearman
round(cor(ctDNAanalyses$CA125, ctDNAanalyses$propsol, method = "spearman", use="pairwise.complete.obs"), 2)

## CA125 vs Presence of more than 10 locules: Point biserial correlation
round(biserial.cor(ctDNAanalyses$CA125, ctDNAanalyses$loc10, use = "complete.obs"), 2)

## CA125 vs Number of papillary projections: Spearman
round(cor(ctDNAanalyses$CA125, ctDNAanalyses$papnr, method = "spearman", use="pairwise.complete.obs"), 2)

## CA125 vs Presence of acoustic shadows: Point biserial correlation
round(biserial.cor(ctDNAanalyses$CA125, ctDNAanalyses$shadowsbin, use = "complete.obs"), 2)

## CA125 vs Presence of ascites: Point biserial correlation
round(biserial.cor(ctDNAanalyses$CA125, ctDNAanalyses$ascitesbin, use = "complete.obs"), 2)


#### 3.4.3 Correlation between ADNEX-variables and proteins ####

## Age: Spearman
CorAge <- cor(ctDNAanalyses[, "Patient age"], ctDNAanalyses[, c("gw_z_score", "nucleosome_score")], method = "spearman", use="pairwise.complete.obs")
round(CorAge, 2)

## Maximum diameter of lesion: Spearman
CorLesdmax <- cor(ctDNAanalyses[, "Lesion largest diameter"], ctDNAanalyses[, c("gw_z_score", "nucleosome_score")], method = "spearman", use="pairwise.complete.obs")
round(CorLesdmax, 2)

## Proportion of solid tissue: Spearman
CorPropsol <- cor(ctDNAanalyses[, "propsol"], ctDNAanalyses[, c("gw_z_score", "nucleosome_score")], method = "spearman", use="pairwise.complete.obs")
round(CorPropsol, 2)

## Presence of more than 10 locules: Point biserial correlation
round(biserial.cor(ctDNAanalyses[, c("gw_z_score")], ctDNAanalyses[, "loc10"], use = "complete.obs"), 2)
round(biserial.cor(ctDNAanalyses[, c("nucleosome_score")], ctDNAanalyses[, "loc10"], use = "complete.obs"), 2)

## Number of papillary projections: Spearman
CorPapnr <- cor(ctDNAanalyses[, "papnr"], ctDNAanalyses[, c("gw_z_score", "nucleosome_score")], method = "spearman", use="pairwise.complete.obs")
round(CorPapnr, 2)

## Presence of acoustic shadows: Point biserial correlation
round(biserial.cor(ctDNAanalyses[, c("gw_z_score")], ctDNAanalyses[, "shadowsbin"], use = "complete.obs"), 2)
round(biserial.cor(ctDNAanalyses[, c("nucleosome_score")], ctDNAanalyses[, "shadowsbin"], use = "complete.obs"), 2)

## Presence of ascites: Point biserial correlation
round(biserial.cor(ctDNAanalyses[, c("gw_z_score")], ctDNAanalyses[, "ascitesbin"], use = "complete.obs"), 2)
round(biserial.cor(ctDNAanalyses[, c("nucleosome_score")], ctDNAanalyses[, "ascitesbin"], use = "complete.obs"), 2)

## CA125: Spearman
CorCA125 <- cor(ctDNAanalyses[, "CA125"], ctDNAanalyses[, c("gw_z_score", "nucleosome_score")], method = "spearman", use="pairwise.complete.obs")
round(CorCA125, 2)

#### 3.5 Bland-Altman analysis ####

blandr.statistics(TransIOTAanalyses$gw_z_score[duplicated(TransIOTAanalyses$`Patient ID`, fromLast = T) & !TransIOTAanalyses$`Patient ID` %in% c(15892, 17445)], TransIOTAanalyses$gw_z_score[duplicated(TransIOTAanalyses$`Patient ID`, fromLast = F) & !TransIOTAanalyses$`Patient ID` %in% c(15892, 17445)])
blandr.draw(TransIOTAanalyses$gw_z_score[duplicated(TransIOTAanalyses$`Patient ID`, fromLast = T) & !TransIOTAanalyses$`Patient ID` %in% c(15892, 17445)], TransIOTAanalyses$gw_z_score[duplicated(TransIOTAanalyses$`Patient ID`, fromLast = F) & !TransIOTAanalyses$`Patient ID` %in% c(15892, 17445)], ciDisplay = FALSE , ciShading = FALSE)
blandr.draw(TransIOTAanalyses$gw_z_score[duplicated(TransIOTAanalyses$`Patient ID`, fromLast = T) & !TransIOTAanalyses$`Patient ID` %in% c(15892, 17445)], TransIOTAanalyses$gw_z_score[duplicated(TransIOTAanalyses$`Patient ID`, fromLast = F) & !TransIOTAanalyses$`Patient ID` %in% c(15892, 17445)])
# 800 x 600

blandr.statistics(TransIOTAanalyses$nucleosome_score[duplicated(TransIOTAanalyses$`Patient ID`, fromLast = T) & !TransIOTAanalyses$`Patient ID` %in% c(15892, 17445)], TransIOTAanalyses$nucleosome_score[duplicated(TransIOTAanalyses$`Patient ID`, fromLast = F) & !TransIOTAanalyses$`Patient ID` %in% c(15892, 17445)])
blandr.draw(TransIOTAanalyses$nucleosome_score[duplicated(TransIOTAanalyses$`Patient ID`, fromLast = T) & !TransIOTAanalyses$`Patient ID` %in% c(15892, 17445)], TransIOTAanalyses$nucleosome_score[duplicated(TransIOTAanalyses$`Patient ID`, fromLast = F) & !TransIOTAanalyses$`Patient ID` %in% c(15892, 17445)], ciDisplay = FALSE , ciShading = FALSE)
blandr.draw(TransIOTAanalyses$nucleosome_score[duplicated(TransIOTAanalyses$`Patient ID`, fromLast = T) & !TransIOTAanalyses$`Patient ID` %in% c(15892, 17445)], TransIOTAanalyses$nucleosome_score[duplicated(TransIOTAanalyses$`Patient ID`, fromLast = F) & !TransIOTAanalyses$`Patient ID` %in% c(15892, 17445)])
# 800 x 600


#### 4. Imputation ####

## Variable for protocol used
ctDNAanalyses$Protocol <- ifelse(is.na(ctDNAanalyses$GC_Adapter), "DL", "JV")
IDs.mis <- c("14409", "14873", "16082", "18155", "18356") # Cases with a missing cfDNA score that should have been analysed in the lab of JV
ctDNAanalyses$Protocol[ctDNAanalyses$Patient.ID %in% IDs.mis] <- "JV"
ctDNAanalyses$Protocol <- as.factor(ctDNAanalyses$Protocol)

## Loglog transformation of CA125
ctDNAanalyses$T_CA125 <- log(log(ctDNAanalyses$CA125 + 1))

## Log transformation of 'maximum diameter of lesion'
ctDNAanalyses$Llesdmax <- log2(ctDNAanalyses$`Lesion largest diameter`) 

## Quadratic term for 'Proportion of solid tissue'
ctDNAanalyses$Qpropsol <- ctDNAanalyses$propsol^2

## Boxcox transformation for genome-wide z-score
Box <- boxcox((ctDNAanalyses$gw_z_score + 1) ~ 1, plotit = FALSE) # lambda = -0.3
Cox = data.frame(Box$x, Box$y)            # Create a data frame with the results
Cox2 = Cox[with(Cox, order(-Cox$Box.y)),] # Order the new data frame by decreasing y
Cox2[1,]                                  # Display the lambda with the greatest log likelihood
lambda_GW = Cox2[1, "Box.x"]              # Extract that lambda
ctDNAanalyses$T_GWZ = ((ctDNAanalyses$gw_z_score + 1) ^ lambda_GW - 1)/lambda_GW   # Transform the original data

ctDNAanalyses <- data.frame(ctDNAanalyses)
VarImp <- c("Patient.age", "Llesdmax", "propsol", "Qpropsol", "loc10", "papnr", "shadowsbin", "ascitesbin", 
            "T_CA125", "T_GWZ", "nucleosome_score", "ADNEXgroups", "Protocol")

hist <- ggplot(ctDNAanalyses, aes(x = T_GWZ)) +
        geom_histogram(aes(y = ..count..), alpha=.75, binwidth = 0.1) + # , binwidth = 0.01
        theme_minimal() +
        labs(x = "Genome wide z-score",
             y = "Frequency")
hist # 800 x 400

# Nucleosome score
hist <- ggplot(ctDNAanalyses, aes(x = nucleosome_score)) +
        geom_histogram(aes(y = ..count..), binwidth = 0.01, alpha=.75) +
        theme_minimal() +
        labs(x = "Nucleosome score",
             y = "Frequency") 
hist # 800 x 400

# CA125
hist <- ggplot(ctDNAanalyses, aes(x = T_CA125)) +
        geom_histogram(aes(y = ..count..), binwidth = 0.1, alpha=.75) +
        theme_minimal() +
        labs(x = "log(CA125)",
             y = "Frequency") 
hist # 800 x 400

## Imputation model
AllVars <- colnames(ctDNAanalyses)

PredMatr = matrix(0, length(AllVars), length(AllVars))
colnames(PredMatr) <- AllVars
rownames(PredMatr) <- AllVars

# PredMatr: A value of 1 indicates that the column variable is used as a predictor to impute the target (row) variable, and a 0 means that it is not used

# Prediction Matrix
PredMatr[c("T_GWZ", "nucleosome_score", "T_CA125"), VarImp] <- 1
# Diagonal of matrix is set to zero
diag(PredMatr) <- 0

Init = mice(ctDNAanalyses, m = 1, maxit = 0, predictorMatrix = PredMatr)
Init$loggedEvents 
Init$method

Method = sapply(colnames(ctDNAanalyses),
                function(x){
                  if(x == "nucleosome_score")
                    "pmm"
                  else if(x == "T_CA125")
                    "pmm"
                  else if(x == "T_GWZ")
                    "pmm"
                  else 
                    ""
                })
Method

imp_ctDNA <- mice(ctDNAanalyses, m = 1, seed = 24, me = Method,
                  predictorMatrix = PredMatr, maxit = 100, vis = "monotone", ridge = 1e-3)
imp_ctDNA$loggedEvents

# Check convergence
plot(imp_ctDNA) # 600 x 400

# Check density plot
densityplot(imp_ctDNA, ~ T_GWZ) # 300 x 300
densityplot(imp_ctDNA, ~ nucleosome_score)
densityplot(imp_ctDNA, ~ T_CA125)

ctDNAimp <-  mice::complete(imp_ctDNA, "long")

ctDNAimp$GWZ <- (lambda_GW * ctDNAimp$T_GWZ + 1)^(1/lambda_GW) - 1
ctDNAimp$CA125 <- exp(exp(ctDNAimp$T_CA125)) - 1
ctDNAimp[, c("gw_z_score", "GWZ")]

## Log transformation of CA125
ctDNAimp$lCA125 <- log2(ctDNAimp$CA125)


#### 5. Analyses ####

source("Scripts/TransIOTA ctDNA - functions.R")


#### 5.1 Univariable analysis ####

#### 5.1.1 Benign vs malignant ####
auc.nonpara.mw(ctDNAimp$GWZ[ctDNAimp$Outcome == "Malignant"], ctDNAimp$GWZ[ctDNAimp$Outcome == "Benign"], method = "pepe")
auc.nonpara.mw(ctDNAimp$nucleosome_score[ctDNAimp$Outcome == "Malignant"], ctDNAimp$nucleosome_score[ctDNAimp$Outcome == "Benign"], method = "pepe")

#### 5.1.2 Pairwise comparisons ####

## Genome-wide z-score
auc.nonpara.mw(ctDNAimp$GWZ[ctDNAimp$ADNEXgroups == "Borderline"], ctDNAimp$GWZ[ctDNAimp$ADNEXgroups == "Benign"], method = "pepe")
auc.nonpara.mw(ctDNAimp$GWZ[ctDNAimp$ADNEXgroups == "Stage I invasive"], ctDNAimp$GWZ[ctDNAimp$ADNEXgroups == "Benign"], method = "pepe")
auc.nonpara.mw(ctDNAimp$GWZ[ctDNAimp$ADNEXgroups == "Stage II-IV invasive"], ctDNAimp$GWZ[ctDNAimp$ADNEXgroups == "Benign"], method = "pepe")
auc.nonpara.mw(ctDNAimp$GWZ[ctDNAimp$ADNEXgroups == "Metastatic"], ctDNAimp$GWZ[ctDNAimp$ADNEXgroups == "Benign"], method = "pepe")
auc.nonpara.mw(ctDNAimp$GWZ[ctDNAimp$ADNEXgroups == "Stage I invasive"], ctDNAimp$GWZ[ctDNAimp$ADNEXgroups == "Borderline"], method = "pepe")
auc.nonpara.mw(ctDNAimp$GWZ[ctDNAimp$ADNEXgroups == "Stage II-IV invasive"], ctDNAimp$GWZ[ctDNAimp$ADNEXgroups == "Borderline"], method = "pepe")
auc.nonpara.mw(ctDNAimp$GWZ[ctDNAimp$ADNEXgroups == "Metastatic"], ctDNAimp$GWZ[ctDNAimp$ADNEXgroups == "Borderline"], method = "pepe")
auc.nonpara.mw(ctDNAimp$GWZ[ctDNAimp$ADNEXgroups == "Stage II-IV invasive"], ctDNAimp$GWZ[ctDNAimp$ADNEXgroups == "Stage I invasive"], method = "pepe")
auc.nonpara.mw(ctDNAimp$GWZ[ctDNAimp$ADNEXgroups == "Metastatic"], ctDNAimp$GWZ[ctDNAimp$ADNEXgroups == "Stage I invasive"], method = "pepe")
auc.nonpara.mw(ctDNAimp$GWZ[ctDNAimp$ADNEXgroups == "Stage II-IV invasive"], ctDNAimp$GWZ[ctDNAimp$ADNEXgroups == "Metastatic"], method = "pepe")

## Nucleosome score
auc.nonpara.mw(ctDNAimp$nucleosome_score[ctDNAimp$ADNEXgroups == "Borderline"], ctDNAimp$nucleosome_score[ctDNAimp$ADNEXgroups == "Benign"], method = "pepe")
auc.nonpara.mw(ctDNAimp$nucleosome_score[ctDNAimp$ADNEXgroups == "Stage I invasive"], ctDNAimp$nucleosome_score[ctDNAimp$ADNEXgroups == "Benign"], method = "pepe")
auc.nonpara.mw(ctDNAimp$nucleosome_score[ctDNAimp$ADNEXgroups == "Stage II-IV invasive"], ctDNAimp$nucleosome_score[ctDNAimp$ADNEXgroups == "Benign"], method = "pepe")
auc.nonpara.mw(ctDNAimp$nucleosome_score[ctDNAimp$ADNEXgroups == "Metastatic"], ctDNAimp$nucleosome_score[ctDNAimp$ADNEXgroups == "Benign"], method = "pepe")
auc.nonpara.mw(ctDNAimp$nucleosome_score[ctDNAimp$ADNEXgroups == "Stage I invasive"], ctDNAimp$nucleosome_score[ctDNAimp$ADNEXgroups == "Borderline"], method = "pepe")
auc.nonpara.mw(ctDNAimp$nucleosome_score[ctDNAimp$ADNEXgroups == "Stage II-IV invasive"], ctDNAimp$nucleosome_score[ctDNAimp$ADNEXgroups == "Borderline"], method = "pepe")
auc.nonpara.mw(ctDNAimp$nucleosome_score[ctDNAimp$ADNEXgroups == "Metastatic"], ctDNAimp$nucleosome_score[ctDNAimp$ADNEXgroups == "Borderline"], method = "pepe")
auc.nonpara.mw(ctDNAimp$nucleosome_score[ctDNAimp$ADNEXgroups == "Stage II-IV invasive"], ctDNAimp$nucleosome_score[ctDNAimp$ADNEXgroups == "Stage I invasive"], method = "pepe")
auc.nonpara.mw(ctDNAimp$nucleosome_score[ctDNAimp$ADNEXgroups == "Metastatic"], ctDNAimp$nucleosome_score[ctDNAimp$ADNEXgroups == "Stage I invasive"], method = "pepe")
auc.nonpara.mw(ctDNAimp$nucleosome_score[ctDNAimp$ADNEXgroups == "Stage II-IV invasive"], ctDNAimp$nucleosome_score[ctDNAimp$ADNEXgroups == "Metastatic"], method = "pepe")


#### 5.2 Multivariable analysis ####

MA.AUC <- matrix(nrow = 8, ncol = 6)
colnames(MA.AUC) <- c("Model", "AUC", "LL_AUC", "UL_AUC", "Optimism", "cAUC")
rownames <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")
MA.AUC <- data.frame(MA.AUC)
MA.AUC$Model <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")

MA.NB5 <- matrix(nrow = 8, ncol = 4)
colnames(MA.NB5) <- c("Model", "NB", "Optimism", "cNB")
rownames <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")
MA.NB5 <- data.frame(MA.NB5)
MA.NB5$Model <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")

MA.NB10 <- matrix(nrow = 8, ncol = 4)
colnames(MA.NB10) <- c("Model", "NB", "Optimism", "cNB")
rownames <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")
MA.NB10 <- data.frame(MA.NB10)
MA.NB10$Model <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")

MA.NB20 <- matrix(nrow = 8, ncol = 4)
colnames(MA.NB20) <- c("Model", "NB", "Optimism", "cNB")
rownames <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")
MA.NB20 <- data.frame(MA.NB20)
MA.NB20$Model <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")

MA.NB30 <- matrix(nrow = 8, ncol = 4)
colnames(MA.NB30) <- c("Model", "NB", "Optimism", "cNB")
rownames <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")
MA.NB30 <- data.frame(MA.NB30)
MA.NB30$Model <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")

MA.NB40 <- matrix(nrow = 8, ncol = 4)
colnames(MA.NB40) <- c("Model", "NB", "Optimism", "cNB")
rownames <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")
MA.NB40 <- data.frame(MA.NB40)
MA.NB40$Model <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")

MA.SS1 <- matrix(nrow = 8, ncol = 11)
colnames(MA.SS1) <- c("Model", "Sens", "LL_Sens", "UL_Sens", "SensOptimism", "cSens", "Spec", "LL_Spec", "UL_Spec", "SpecOptimism", "cSpec")
rownames <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")
MA.SS1 <- data.frame(MA.SS1)
MA.SS1$Model <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")

MA.SS5 <- matrix(nrow = 8, ncol = 11)
colnames(MA.SS5) <- c("Model", "Sens", "LL_Sens", "UL_Sens", "SensOptimism", "cSens", "Spec", "LL_Spec", "UL_Spec", "SpecOptimism", "cSpec")
rownames <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")
MA.SS5 <- data.frame(MA.SS5)
MA.SS5$Model <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")

MA.SS10 <- matrix(nrow = 8, ncol = 11)
colnames(MA.SS10) <- c("Model", "Sens", "LL_Sens", "UL_Sens", "SensOptimism", "cSens", "Spec", "LL_Spec", "UL_Spec", "SpecOptimism", "cSpec")
rownames <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")
MA.SS10 <- data.frame(MA.SS10)
MA.SS10$Model <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")

MA.SS20 <- matrix(nrow = 8, ncol = 11)
colnames(MA.SS20) <- c("Model", "Sens", "LL_Sens", "UL_Sens", "SensOptimism", "cSens", "Spec", "LL_Spec", "UL_Spec", "SpecOptimism", "cSpec")
rownames <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")
MA.SS20 <- data.frame(MA.SS20)
MA.SS20$Model <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")

MA.SS30 <- matrix(nrow = 8, ncol = 11)
colnames(MA.SS30) <- c("Model", "Sens", "LL_Sens", "UL_Sens", "SensOptimism", "cSens", "Spec", "LL_Spec", "UL_Spec", "SpecOptimism", "cSpec")
rownames <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")
MA.SS30 <- data.frame(MA.SS30)
MA.SS30$Model <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")

MA.SS40 <- matrix(nrow = 8, ncol = 11)
colnames(MA.SS40) <- c("Model", "Sens", "LL_Sens", "UL_Sens", "SensOptimism", "cSens", "Spec", "LL_Spec", "UL_Spec", "SpecOptimism", "cSpec")
rownames <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")
MA.SS40 <- data.frame(MA.SS40)
MA.SS40$Model <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")

MA.SS50 <- matrix(nrow = 8, ncol = 11)
colnames(MA.SS50) <- c("Model", "Sens", "LL_Sens", "UL_Sens", "SensOptimism", "cSens", "Spec", "LL_Spec", "UL_Spec", "SpecOptimism", "cSpec")
rownames <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")
MA.SS50 <- data.frame(MA.SS50)
MA.SS50$Model <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")

## Fit the different models

# Model 1
M1 <- logistf(OutcomeBin ~ Patient.age + Llesdmax + propsol + Qpropsol + loc10 + papnr + shadowsbin + ascitesbin, data = ctDNAimp,
              control=logistf.control(maxit = 100))

# Model 2
M2 <- logistf(OutcomeBin ~ Patient.age + Llesdmax + propsol + Qpropsol + loc10 + papnr + shadowsbin + ascitesbin + rcs(GWZ, 3), data = ctDNAimp,
              control=logistf.control(maxit = 100)) 

# Model 3
M3 <- logistf(OutcomeBin ~ Patient.age + Llesdmax + propsol + Qpropsol + loc10 + papnr + shadowsbin + ascitesbin + rcs(nucleosome_score, 3), data = ctDNAimp,
              control=logistf.control(maxit = 100))

# Model 4
M4 <- logistf(OutcomeBin ~ Patient.age + Llesdmax + propsol + Qpropsol + loc10 + papnr + shadowsbin + ascitesbin + rcs(GWZ, 3) + rcs(nucleosome_score, 3), data = ctDNAimp,
              control=logistf.control(maxit = 100))

# Model 5
M5 <- logistf(OutcomeBin ~ Patient.age + Llesdmax + propsol + Qpropsol + loc10 + papnr + shadowsbin + ascitesbin + lCA125, data = ctDNAimp,
              control=logistf.control(maxit = 100))

# Model 6
M6 <- logistf(OutcomeBin ~ Patient.age + Llesdmax + propsol + Qpropsol + loc10 + papnr + shadowsbin + ascitesbin + lCA125 + rcs(GWZ, 3), data = ctDNAimp,
              control=logistf.control(maxit = 100))

# Model 7
M7 <- logistf(OutcomeBin ~ Patient.age + Llesdmax + propsol + Qpropsol + loc10 + papnr + shadowsbin + ascitesbin + lCA125 + rcs(nucleosome_score, 3), data = ctDNAimp,
              control=logistf.control(maxit = 100))

# Model 8
M8 <- logistf(OutcomeBin ~ Patient.age + Llesdmax + propsol + Qpropsol + loc10 + papnr + shadowsbin + ascitesbin + lCA125 + rcs(GWZ, 3) + rcs(nucleosome_score, 3), data = ctDNAimp,
              control=logistf.control(maxit = 100))

## Joint likelihood ratio test
# M1 vs M4
anova(M1, M4)

# M5 vs M8
anova(M5, M8)

#### 5.2.1 Apparent validation ####
Models <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")

for(i in 1:8){
  
  m_M <- (ctDNAimp$OutcomeBin == 1)
  m_B <- (ctDNAimp$OutcomeBin == 0)
  
  ## Calculate the apparent AUC
  f <- get(Models[i])
  auc_mw <- round(auRoc::auc.nonpara.mw(x = f$predict[m_M],
                                        y = f$predict[m_B],
                                        conf.level = 0.95,
                                        method = "pepe"), 3)
  MA.AUC$AUC[i]     <- auc_mw[1]
  MA.AUC$LL_AUC[i]  <- auc_mw[2]
  MA.AUC$UL_AUC[i]  <- auc_mw[3]
  
  predict  <- f$predict
  datapred <- cbind(ctDNAimp, predict)
  
  ## Calculate Net Benefit
  # 5% 
  CM5  <- confusion.matrix(obs = datapred$OutcomeBin, pred = datapred$predict, threshold = 0.05)
  MA.NB5$NB[i]  <- round(CM5[2,2] / nrow(datapred) - CM5[2,1] / nrow(datapred) * (0.05 / (1 - 0.05)), 3) # NB = TP / n - FP / n * (threshold / (1 - threshold))
  # 10%
  CM10 <- confusion.matrix(obs = datapred$OutcomeBin, pred = datapred$predict, threshold = 0.10)
  MA.NB10$NB[i] <- round(CM10[2,2] / nrow(datapred) - CM10[2,1] / nrow(datapred) * (0.10 / (1 - 0.10)), 3)
  # 20%
  CM20 <- confusion.matrix(obs = datapred$OutcomeBin, pred = datapred$predict, threshold = 0.20)
  MA.NB20$NB[i] <- round(CM20[2,2] / nrow(datapred) - CM20[2,1] / nrow(datapred) * (0.20 / (1 - 0.20)), 3)
  # 30%
  CM30 <- confusion.matrix(obs = datapred$OutcomeBin, pred = datapred$predict, threshold = 0.30)
  MA.NB30$NB[i] <- round(CM30[2,2] / nrow(datapred) - CM30[2,1] / nrow(datapred) * (0.30 / (1 - 0.30)), 3)
  # 40%
  CM40 <- confusion.matrix(obs = datapred$OutcomeBin, pred = datapred$predict, threshold = 0.40)
  MA.NB40$NB[i] <- round(CM40[2,2] / nrow(datapred) - CM40[2,1] / nrow(datapred) * (0.40 / (1 - 0.40)), 3)
  
  ## Calculate Sensitivity and specificity
  # 1%
  MA.SS1[i, c("Sens", "LL_Sens", "UL_Sens")] <- Sensitivity(pred = predict, outcome = OutcomeBin, threshold = 0.01, data = datapred)
  MA.SS1[i, c("Spec", "LL_Spec", "UL_Spec")] <- Specificity(pred = predict, outcome = OutcomeBin, threshold = 0.01, data = datapred)
  # 5%
  MA.SS5[i, c("Sens", "LL_Sens", "UL_Sens")] <- Sensitivity(pred = predict, outcome = OutcomeBin, threshold = 0.05, data = datapred)
  MA.SS5[i, c("Spec", "LL_Spec", "UL_Spec")] <- Specificity(pred = predict, outcome = OutcomeBin, threshold = 0.05, data = datapred)
  # 10%
  MA.SS10[i, c("Sens", "LL_Sens", "UL_Sens")] <- Sensitivity(pred = predict, outcome = OutcomeBin, threshold = 0.10, data = datapred)
  MA.SS10[i, c("Spec", "LL_Spec", "UL_Spec")] <- Specificity(pred = predict, outcome = OutcomeBin, threshold = 0.10, data = datapred)
  # 20%
  MA.SS20[i, c("Sens", "LL_Sens", "UL_Sens")] <- Sensitivity(pred = predict, outcome = OutcomeBin, threshold = 0.20, data = datapred)
  MA.SS20[i, c("Spec", "LL_Spec", "UL_Spec")] <- Specificity(pred = predict, outcome = OutcomeBin, threshold = 0.20, data = datapred)
  # 30%
  MA.SS30[i, c("Sens", "LL_Sens", "UL_Sens")] <- Sensitivity(pred = predict, outcome = OutcomeBin, threshold = 0.30, data = datapred)
  MA.SS30[i, c("Spec", "LL_Spec", "UL_Spec")] <- Specificity(pred = predict, outcome = OutcomeBin, threshold = 0.30, data = datapred)
  # 40%
  MA.SS40[i, c("Sens", "LL_Sens", "UL_Sens")] <- Sensitivity(pred = predict, outcome = OutcomeBin, threshold = 0.40, data = datapred)
  MA.SS40[i, c("Spec", "LL_Spec", "UL_Spec")] <- Specificity(pred = predict, outcome = OutcomeBin, threshold = 0.40, data = datapred)
  # 50%
  MA.SS50[i, c("Sens", "LL_Sens", "UL_Sens")] <- Sensitivity(pred = predict, outcome = OutcomeBin, threshold = 0.50, data = datapred)
  MA.SS50[i, c("Spec", "LL_Spec", "UL_Spec")] <- Specificity(pred = predict, outcome = OutcomeBin, threshold = 0.50, data = datapred)
} # loop for each model


## Calculate delta AUC
MA.diffAUC <- matrix(nrow = 6, ncol = 6)
colnames(MA.diffAUC) <- c("Comparison", "dAUC", "LL_dAUC", "UL_dAUC", "Optimism", "cdAUC")
rownames <- c("1vs2", "1vs3", "1vs4", "5vs6", "5vs7", "5vs8")
MA.diffAUC <- data.frame(MA.diffAUC)
MA.diffAUC$Comparison <- c("1vs2", "1vs3", "1vs4", "5vs6", "5vs7", "5vs8")

MA.diffAUC[1, c("dAUC", "LL_dAUC", "UL_dAUC")]    <- dAUC(data = ctDNAimp, predictor1 = M1$predict, predictor2 = M2$predict)[c(1, 3, 4)]
MA.diffAUC[2, c("dAUC", "LL_dAUC", "UL_dAUC")]    <- dAUC(data = ctDNAimp, predictor1 = M1$predict, predictor2 = M3$predict)[c(1, 3, 4)]
MA.diffAUC[3, c("dAUC", "LL_dAUC", "UL_dAUC")]    <- dAUC(data = ctDNAimp, predictor1 = M1$predict, predictor2 = M4$predict)[c(1, 3, 4)]
MA.diffAUC[4, c("dAUC", "LL_dAUC", "UL_dAUC")]    <- dAUC(data = ctDNAimp, predictor1 = M5$predict, predictor2 = M6$predict)[c(1, 3, 4)]
MA.diffAUC[5, c("dAUC", "LL_dAUC", "UL_dAUC")]    <- dAUC(data = ctDNAimp, predictor1 = M5$predict, predictor2 = M7$predict)[c(1, 3, 4)]
MA.diffAUC[6, c("dAUC", "LL_dAUC", "UL_dAUC")]    <- dAUC(data = ctDNAimp, predictor1 = M5$predict, predictor2 = M8$predict)[c(1, 3, 4)]

## Calculate delta NB
MA.diffNB5 <- matrix(nrow = 6, ncol = 6)
colnames(MA.diffNB5) <- c("Comparison", "dNB", "LL_dNB", "UL_dNB", "Optimism", "cdNB")
rownames <- c("1vs2", "1vs3", "1vs4", "5vs6", "5vs7", "5vs8")
MA.diffNB5 <- data.frame(MA.diffNB5)
MA.diffNB5$Comparison <- c("1vs2", "1vs3", "1vs4", "5vs6", "5vs7", "5vs8")

MA.diffNB10 <- matrix(nrow = 6, ncol = 6)
colnames(MA.diffNB10) <- c("Comparison", "dNB", "LL_dNB", "UL_dNB", "Optimism", "cdNB")
rownames <- c("1vs2", "1vs3", "1vs4", "5vs6", "5vs7", "5vs8")
MA.diffNB10 <- data.frame(MA.diffNB10)
MA.diffNB10$Comparison <- c("1vs2", "1vs3", "1vs4", "5vs6", "5vs7", "5vs8")

MA.diffNB20 <- matrix(nrow = 6, ncol = 6)
colnames(MA.diffNB20) <- c("Comparison", "dNB", "LL_dNB", "UL_dNB", "Optimism", "cdNB")
rownames <- c("1vs2", "1vs3", "1vs4", "5vs6", "5vs7", "5vs8")
MA.diffNB20 <- data.frame(MA.diffNB20)
MA.diffNB20$Comparison <- c("1vs2", "1vs3", "1vs4", "5vs6", "5vs7", "5vs8")

MA.diffNB30 <- matrix(nrow = 6, ncol = 6)
colnames(MA.diffNB30) <- c("Comparison", "dNB", "LL_dNB", "UL_dNB", "Optimism", "cdNB")
rownames <- c("1vs2", "1vs3", "1vs4", "5vs6", "5vs7", "5vs8")
MA.diffNB30 <- data.frame(MA.diffNB30)
MA.diffNB30$Comparison <- c("1vs2", "1vs3", "1vs4", "5vs6", "5vs7", "5vs8")

MA.diffNB40 <- matrix(nrow = 6, ncol = 6)
colnames(MA.diffNB40) <- c("Comparison", "dNB", "LL_dNB", "UL_dNB", "Optimism", "cdNB")
rownames <- c("1vs2", "1vs3", "1vs4", "5vs6", "5vs7", "5vs8")
MA.diffNB40 <- data.frame(MA.diffNB40)
MA.diffNB40$Comparison <- c("1vs2", "1vs3", "1vs4", "5vs6", "5vs7", "5vs8")

# 5% 
MA.diffNB5[1, c("dNB", "LL_dNB", "UL_dNB")] <- diffNB(data = ctDNAimp, model1 = M1, model2 = M2, predictor1 = M1$predict, predictor2 = M2$predict, outcome = OutcomeBin, threshold = 0.05) 
MA.diffNB5[2, c("dNB", "LL_dNB", "UL_dNB")] <- diffNB(data = ctDNAimp, model1 = M1, model2 = M3, predictor1 = M1$predict, predictor2 = M3$predict, outcome = OutcomeBin, threshold = 0.05)
MA.diffNB5[3, c("dNB", "LL_dNB", "UL_dNB")] <- diffNB(data = ctDNAimp, model1 = M1, model2 = M4, predictor1 = M1$predict, predictor2 = M4$predict, outcome = OutcomeBin, threshold = 0.05)
MA.diffNB5[4, c("dNB", "LL_dNB", "UL_dNB")] <- diffNB(data = ctDNAimp, model1 = M5, model2 = M6, predictor1 = M5$predict, predictor2 = M6$predict, outcome = OutcomeBin, threshold = 0.05)
MA.diffNB5[5, c("dNB", "LL_dNB", "UL_dNB")] <- diffNB(data = ctDNAimp, model1 = M5, model2 = M7, predictor1 = M5$predict, predictor2 = M7$predict, outcome = OutcomeBin, threshold = 0.05)
MA.diffNB5[6, c("dNB", "LL_dNB", "UL_dNB")] <- diffNB(data = ctDNAimp, model1 = M5, model2 = M8, predictor1 = M5$predict, predictor2 = M8$predict, outcome = OutcomeBin, threshold = 0.05)

# 10%
MA.diffNB10[1, c("dNB", "LL_dNB", "UL_dNB")] <- diffNB(data = ctDNAimp, model1 = M1, model2 = M2, predictor1 = M1$predict, predictor2 = M2$predict, outcome = OutcomeBin, threshold = 0.10)
MA.diffNB10[2, c("dNB", "LL_dNB", "UL_dNB")] <- diffNB(data = ctDNAimp, model1 = M1, model2 = M3, predictor1 = M1$predict, predictor2 = M3$predict, outcome = OutcomeBin, threshold = 0.10)
MA.diffNB10[3, c("dNB", "LL_dNB", "UL_dNB")] <- diffNB(data = ctDNAimp, model1 = M1, model2 = M4, predictor1 = M1$predict, predictor2 = M4$predict, outcome = OutcomeBin, threshold = 0.10)
MA.diffNB10[4, c("dNB", "LL_dNB", "UL_dNB")] <- diffNB(data = ctDNAimp, model1 = M5, model2 = M6, predictor1 = M5$predict, predictor2 = M6$predict, outcome = OutcomeBin, threshold = 0.10)
MA.diffNB10[5, c("dNB", "LL_dNB", "UL_dNB")] <- diffNB(data = ctDNAimp, model1 = M5, model2 = M7, predictor1 = M5$predict, predictor2 = M7$predict, outcome = OutcomeBin, threshold = 0.10)
MA.diffNB10[6, c("dNB", "LL_dNB", "UL_dNB")] <- diffNB(data = ctDNAimp, model1 = M5, model2 = M8, predictor1 = M5$predict, predictor2 = M8$predict, outcome = OutcomeBin, threshold = 0.10)

# 20%
MA.diffNB20[1, c("dNB", "LL_dNB", "UL_dNB")] <- diffNB(data = ctDNAimp, model1 = M1, model2 = M2, predictor1 = M1$predict, predictor2 = M2$predict, outcome = OutcomeBin, threshold = 0.20)
MA.diffNB20[2, c("dNB", "LL_dNB", "UL_dNB")] <- diffNB(data = ctDNAimp, model1 = M1, model2 = M3, predictor1 = M1$predict, predictor2 = M3$predict, outcome = OutcomeBin, threshold = 0.20)
MA.diffNB20[3, c("dNB", "LL_dNB", "UL_dNB")] <- diffNB(data = ctDNAimp, model1 = M1, model2 = M4, predictor1 = M1$predict, predictor2 = M4$predict, outcome = OutcomeBin, threshold = 0.20)
MA.diffNB20[4, c("dNB", "LL_dNB", "UL_dNB")] <- diffNB(data = ctDNAimp, model1 = M5, model2 = M6, predictor1 = M5$predict, predictor2 = M6$predict, outcome = OutcomeBin, threshold = 0.20)
MA.diffNB20[5, c("dNB", "LL_dNB", "UL_dNB")] <- diffNB(data = ctDNAimp, model1 = M5, model2 = M7, predictor1 = M5$predict, predictor2 = M7$predict, outcome = OutcomeBin, threshold = 0.20)
MA.diffNB20[6, c("dNB", "LL_dNB", "UL_dNB")] <- diffNB(data = ctDNAimp, model1 = M5, model2 = M8, predictor1 = M5$predict, predictor2 = M8$predict, outcome = OutcomeBin, threshold = 0.20)

# 30%
MA.diffNB30[1, c("dNB", "LL_dNB", "UL_dNB")] <- diffNB(data = ctDNAimp, model1 = M1, model2 = M2, predictor1 = M1$predict, predictor2 = M2$predict, outcome = OutcomeBin, threshold = 0.30)
MA.diffNB30[2, c("dNB", "LL_dNB", "UL_dNB")] <- diffNB(data = ctDNAimp, model1 = M1, model2 = M3, predictor1 = M1$predict, predictor2 = M3$predict, outcome = OutcomeBin, threshold = 0.30)
MA.diffNB30[3, c("dNB", "LL_dNB", "UL_dNB")] <- diffNB(data = ctDNAimp, model1 = M1, model2 = M4, predictor1 = M1$predict, predictor2 = M4$predict, outcome = OutcomeBin, threshold = 0.30)
MA.diffNB30[4, c("dNB", "LL_dNB", "UL_dNB")] <- diffNB(data = ctDNAimp, model1 = M5, model2 = M6, predictor1 = M5$predict, predictor2 = M6$predict, outcome = OutcomeBin, threshold = 0.30)
MA.diffNB30[5, c("dNB", "LL_dNB", "UL_dNB")] <- diffNB(data = ctDNAimp, model1 = M5, model2 = M7, predictor1 = M5$predict, predictor2 = M7$predict, outcome = OutcomeBin, threshold = 0.30)
MA.diffNB30[6, c("dNB", "LL_dNB", "UL_dNB")] <- diffNB(data = ctDNAimp, model1 = M5, model2 = M8, predictor1 = M5$predict, predictor2 = M8$predict, outcome = OutcomeBin, threshold = 0.30)

# 40%
MA.diffNB40[1, c("dNB", "LL_dNB", "UL_dNB")] <- diffNB(data = ctDNAimp, model1 = M1, model2 = M2, predictor1 = M1$predict, predictor2 = M2$predict, outcome = OutcomeBin, threshold = 0.40)
MA.diffNB40[2, c("dNB", "LL_dNB", "UL_dNB")] <- diffNB(data = ctDNAimp, model1 = M1, model2 = M3, predictor1 = M1$predict, predictor2 = M3$predict, outcome = OutcomeBin, threshold = 0.40)
MA.diffNB40[3, c("dNB", "LL_dNB", "UL_dNB")] <- diffNB(data = ctDNAimp, model1 = M1, model2 = M4, predictor1 = M1$predict, predictor2 = M4$predict, outcome = OutcomeBin, threshold = 0.40)
MA.diffNB40[4, c("dNB", "LL_dNB", "UL_dNB")] <- diffNB(data = ctDNAimp, model1 = M5, model2 = M6, predictor1 = M5$predict, predictor2 = M6$predict, outcome = OutcomeBin, threshold = 0.40)
MA.diffNB40[5, c("dNB", "LL_dNB", "UL_dNB")] <- diffNB(data = ctDNAimp, model1 = M5, model2 = M7, predictor1 = M5$predict, predictor2 = M7$predict, outcome = OutcomeBin, threshold = 0.40)
MA.diffNB40[6, c("dNB", "LL_dNB", "UL_dNB")] <- diffNB(data = ctDNAimp, model1 = M5, model2 = M8, predictor1 = M5$predict, predictor2 = M8$predict, outcome = OutcomeBin, threshold = 0.40)


#### 5.2.2 Enhanced bootstrapping validation ####

# Bootstraps
n_boot <- 1000

# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores)
registerDoParallel(cl)

registerDoSEQ(cl)

# start time
strt <- Sys.time()

# Initial datasets/lists
res_boot <- list()
for(i in 1:n_boot){
  res_boot[[i]] <- list()
}

comp_boot <- list()
for(i in 1:n_boot){
  comp_boot[[i]] <- list()
}

cat("\n\tBootstrapping ...\n", sep = "")
for(i_boot in 1:n_boot){
  cat("\n\tIteration ", i_boot, "/", n_boot, "\n", sep = "")
  set.seed(i_boot)
  # draw a sample with replacement
  m <- sample(x = 1:nrow(ctDNAimp),
              size = nrow(ctDNAimp),
              replace = TRUE)
  x_boot <- ctDNAimp[m, ]
  
  # Fit the different models
  attach(x_boot)
  x_boot$rcs_GWZ <- rcs(GWZ, 3)
  x_boot$rcs_nucleosome <- rcs(nucleosome_score, 3)
  detach(x_boot)
  
  # Model 1
  M1.boot <- logistf(OutcomeBin ~ Patient.age + Llesdmax + propsol + Qpropsol + loc10 + papnr + shadowsbin + ascitesbin, data = x_boot,
                     control=logistf.control(maxit = 2000))
  
  # Model 2
  M2.boot <- logistf(OutcomeBin ~ Patient.age + Llesdmax + propsol + Qpropsol + loc10 + papnr + shadowsbin + ascitesbin + rcs_GWZ, data = x_boot, 
                     control=logistf.control(maxit = 2000)) 
  
  # Model 3
  M3.boot <- logistf(OutcomeBin ~ Patient.age + Llesdmax + propsol + Qpropsol + loc10 + papnr + shadowsbin + ascitesbin + rcs_nucleosome, data = x_boot, 
                     control=logistf.control(maxit = 2000))
  
  # Model 4
  M4.boot <- logistf(OutcomeBin ~ Patient.age + Llesdmax + propsol + Qpropsol + loc10 + papnr + shadowsbin + ascitesbin + rcs_GWZ + rcs_nucleosome, data = x_boot, 
                     control=logistf.control(maxit = 2000))
  
  # Model 5
  M5.boot <- logistf(OutcomeBin ~ Patient.age + Llesdmax + propsol + Qpropsol + loc10 + papnr + shadowsbin + ascitesbin + lCA125, data = x_boot,
                     control=logistf.control(maxit = 2000))
  
  # Model 6
  M6.boot <- logistf(OutcomeBin ~ Patient.age + Llesdmax + propsol + Qpropsol + loc10 + papnr + shadowsbin + ascitesbin + lCA125 + rcs_GWZ, data = x_boot, 
                     control=logistf.control(maxit = 2000))
  
  # Model 7
  M7.boot <- logistf(OutcomeBin ~ Patient.age + Llesdmax + propsol + Qpropsol + loc10 + papnr + shadowsbin + ascitesbin + lCA125 + rcs_nucleosome, data = x_boot, 
                     control=logistf.control(maxit = 2000))
  
  # Model 8
  M8.boot <- logistf(OutcomeBin ~ Patient.age + Llesdmax + propsol + Qpropsol + loc10 + papnr + shadowsbin + ascitesbin + lCA125 + rcs_GWZ + rcs_nucleosome, data = x_boot,
                     control=logistf.control(maxit = 2000))
  
  Models.boot <- c("M1.boot", "M2.boot", "M3.boot", "M4.boot", "M5.boot", "M6.boot", "M7.boot", "M8.boot")
  
  ## Initial dataset for the results
  res.boot <- matrix(nrow = 8, ncol = 41)
  colnames(res.boot) <- c("Model", "AUC.boot", "NB5.boot", "NB10.boot", "NB20.boot", "NB30.boot", "NB40.boot", 
                          "Sens1.boot", "Sens5.boot", "Sens10.boot", "Sens20.boot", "Sens30.boot", "Sens40.boot", "Sens50.boot", 
                          "Spec1.boot", "Spec5.boot", "Spec10.boot", "Spec20.boot", "Spec30.boot", "Spec40.boot", "Spec50.boot", 
                          "AUC.orig", "NB5.orig", "NB10.orig", "NB20.orig", "NB30.orig", "NB40.orig", 
                          "Sens1.orig", "Sens5.orig", "Sens10.orig", "Sens20.orig", "Sens30.orig", "Sens40.orig", "Sens50.orig", 
                          "Spec1.orig", "Spec5.orig", "Spec10.orig", "Spec20.orig", "Spec30.orig", "Spec40.orig", "Spec50.orig")
  rownames(res.boot) <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")
  res.boot <- data.frame(res.boot)
  res.boot$Model <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")
  
  for(i in 1:8){
    ### Bootstrap sample
    m_Mboot <- (x_boot$OutcomeBin == 1)
    m_Bboot <- (x_boot$OutcomeBin == 0)
    
    ## Calculate the apparent AUC
    bootm <- get(Models.boot[i])
    auc_mw <- round(auRoc::auc.nonpara.mw(x = bootm$predict[m_Mboot],
                                          y = bootm$predict[m_Bboot],
                                          conf.level = 0.95,
                                          method = "pepe"), 3)
    res.boot$AUC.boot[i]     <- auc_mw[1]
    
    predict  <- bootm$predict
    boot_pred <- cbind(x_boot, predict)

    ## Calculate Net Benefit
    # 5% 
    CM5.boot  <- confusion.matrix(obs = boot_pred$OutcomeBin, pred = boot_pred$predict, threshold = 0.05)
    res.boot$NB5.boot[i]  <- round(CM5.boot[2,2] / nrow(boot_pred) - CM5.boot[2,1] / nrow(boot_pred) * (0.05 / (1 - 0.05)), 3) # NB = TP / n - FP / n * (threshold / (1 - threshold))
    # 10%
    CM10.boot <- confusion.matrix(obs = boot_pred$OutcomeBin, pred = boot_pred$predict, threshold = 0.10)
    res.boot$NB10.boot[i] <- round(CM10.boot[2,2] / nrow(boot_pred) - CM10.boot[2,1] / nrow(boot_pred) * (0.10 / (1 - 0.10)), 3)
    # 20%
    CM20.boot <- confusion.matrix(obs = boot_pred$OutcomeBin, pred = boot_pred$predict, threshold = 0.20)
    res.boot$NB20.boot[i] <- round(CM20.boot[2,2] / nrow(boot_pred) - CM20.boot[2,1] / nrow(boot_pred) * (0.20 / (1 - 0.20)), 3)
    # 30%
    CM30.boot <- confusion.matrix(obs = boot_pred$OutcomeBin, pred = boot_pred$predict, threshold = 0.30)
    res.boot$NB30.boot[i] <- round(CM30.boot[2,2] / nrow(boot_pred) - CM30.boot[2,1] / nrow(boot_pred) * (0.30 / (1 - 0.30)), 3)
    # 40%
    CM40.boot <- confusion.matrix(obs = boot_pred$OutcomeBin, pred = boot_pred$predict, threshold = 0.40)
    res.boot$NB40.boot[i] <- round(CM40.boot[2,2] / nrow(boot_pred) - CM40.boot[2,1] / nrow(boot_pred) * (0.40 / (1 - 0.40)), 3)
    
    ## Calculate Sensitivity and specificity
    # 1%
    res.boot$Sens1.boot[i] <- Sensitivity(pred = predict, outcome = OutcomeBin, threshold = 0.01, data = boot_pred)[1]
    res.boot$Spec1.boot[i] <- Specificity(pred = predict, outcome = OutcomeBin, threshold = 0.01, data = boot_pred)[1]
    # 5%
    res.boot$Sens5.boot[i] <- Sensitivity(pred = predict, outcome = OutcomeBin, threshold = 0.05, data = boot_pred)[1]
    res.boot$Spec5.boot[i] <- Specificity(pred = predict, outcome = OutcomeBin, threshold = 0.05, data = boot_pred)[1]
    # 10%
    res.boot$Sens10.boot[i] <- Sensitivity(pred = predict, outcome = OutcomeBin, threshold = 0.10, data = boot_pred)[1]
    res.boot$Spec10.boot[i] <- Specificity(pred = predict, outcome = OutcomeBin, threshold = 0.10, data = boot_pred)[1]
    # 20%
    res.boot$Sens20.boot[i] <- Sensitivity(pred = predict, outcome = OutcomeBin, threshold = 0.20, data = boot_pred)[1]
    res.boot$Spec20.boot[i] <- Specificity(pred = predict, outcome = OutcomeBin, threshold = 0.20, data = boot_pred)[1]
    # 30%
    res.boot$Sens30.boot[i] <- Sensitivity(pred = predict, outcome = OutcomeBin, threshold = 0.30, data = boot_pred)[1]
    res.boot$Spec30.boot[i] <- Specificity(pred = predict, outcome = OutcomeBin, threshold = 0.30, data = boot_pred)[1]
    # 40%
    res.boot$Sens40.boot[i] <- Sensitivity(pred = predict, outcome = OutcomeBin, threshold = 0.40, data = boot_pred)[1]
    res.boot$Spec40.boot[i] <- Specificity(pred = predict, outcome = OutcomeBin, threshold = 0.40, data = boot_pred)[1]
    # 50%
    res.boot$Sens50.boot[i] <- Sensitivity(pred = predict, outcome = OutcomeBin, threshold = 0.50, data = boot_pred)[1]
    res.boot$Spec50.boot[i] <- Specificity(pred = predict, outcome = OutcomeBin, threshold = 0.50, data = boot_pred)[1]
    
    ### Original data
    ## Apply model 
    NewGWZ <- cbind(ctDNAimp$GWZ, rcspline.eval(ctDNAimp$GWZ,attr(x_boot$rcs_GWZ,"parms")))
    colnames(NewGWZ) <- c("GWZ", "GWZ'")
    NewNucleo <- cbind(ctDNAimp$nucleosome_score, rcspline.eval(ctDNAimp$nucleosome_score,attr(x_boot$rcs_nucleosome,"parms")))
    colnames(NewNucleo) <- c("nucleosome_score", "nucleosome_score'")
    ctDNAorig <- cbind(ctDNAimp[, c("OutcomeBin", "Patient.age", "Llesdmax", "propsol", "Qpropsol", "loc10", "papnr", "shadowsbin", "ascitesbin", "lCA125")], NewGWZ, NewNucleo)

    M1.boot$pred.orig <- logistf:::predict.logistf(object = M1.boot, newdata = ctDNAorig[, c("Patient.age", "Llesdmax", "propsol", "Qpropsol", "loc10", "papnr", "shadowsbin", "ascitesbin")], type = "response")
    M2.boot$pred.orig <- logistf:::predict.logistf(object = M2.boot, newdata = ctDNAorig[, c("Patient.age", "Llesdmax", "propsol", "Qpropsol", "loc10", "papnr", "shadowsbin", "ascitesbin", "GWZ", "GWZ'")], type = "response")
    M3.boot$pred.orig <- logistf:::predict.logistf(object = M3.boot, newdata = ctDNAorig[, c("Patient.age", "Llesdmax", "propsol", "Qpropsol", "loc10", "papnr", "shadowsbin", "ascitesbin", "nucleosome_score", "nucleosome_score'")], type = "response")
    M4.boot$pred.orig <- logistf:::predict.logistf(object = M4.boot, newdata = ctDNAorig[, c("Patient.age", "Llesdmax", "propsol", "Qpropsol", "loc10", "papnr", "shadowsbin", "ascitesbin", "GWZ", "GWZ'", "nucleosome_score", "nucleosome_score'")], type = "response")
    M5.boot$pred.orig <- logistf:::predict.logistf(object = M5.boot, newdata = ctDNAorig[, c("Patient.age", "Llesdmax", "propsol", "Qpropsol", "loc10", "papnr", "shadowsbin", "ascitesbin", "lCA125")], type = "response")
    M6.boot$pred.orig <- logistf:::predict.logistf(object = M6.boot, newdata = ctDNAorig[, c("Patient.age", "Llesdmax", "propsol", "Qpropsol", "loc10", "papnr", "shadowsbin", "ascitesbin", "lCA125", "GWZ", "GWZ'")], type = "response")
    M7.boot$pred.orig <- logistf:::predict.logistf(object = M7.boot, newdata = ctDNAorig[, c("Patient.age", "Llesdmax", "propsol", "Qpropsol", "loc10", "papnr", "shadowsbin", "ascitesbin", "lCA125", "nucleosome_score", "nucleosome_score'")], type = "response")
    M8.boot$pred.orig <- logistf:::predict.logistf(object = M8.boot, newdata = ctDNAorig[, c("Patient.age", "Llesdmax", "propsol", "Qpropsol", "loc10", "papnr", "shadowsbin", "ascitesbin", "lCA125", "GWZ", "GWZ'", "nucleosome_score", "nucleosome_score'")], type = "response")

    m_Morig <- (ctDNAimp$OutcomeBin == 1)
    m_Borig <- (ctDNAimp$OutcomeBin == 0)
    
    ## Calculate the apparent AUC
    morig <- get(Models.boot[i])
    auc_mw <- round(auRoc::auc.nonpara.mw(x = morig$pred.orig[m_Morig],
                                          y = morig$pred.orig[m_Borig],
                                          conf.level = 0.95,
                                          method = "pepe"), 3)
    res.boot$AUC.orig[i]     <- auc_mw[1]
    
    pred.orig  <- morig$pred.orig
    orig_pred <- cbind(ctDNAimp, pred.orig)
    
    ## Calculate Net Benefit
    # 5% 
    CM5.orig  <- confusion.matrix(obs = orig_pred$OutcomeBin, pred = orig_pred$pred.orig, threshold = 0.05)
    res.boot$NB5.orig[i]  <- round(CM5.orig[2,2] / nrow(orig_pred) - CM5.orig[2,1] / nrow(orig_pred) * (0.05 / (1 - 0.05)), 3) # NB = TP / n - FP / n * (threshold / (1 - threshold))
    # 10%
    CM10.orig <- confusion.matrix(obs = orig_pred$OutcomeBin, pred = orig_pred$pred.orig, threshold = 0.10)
    res.boot$NB10.orig[i] <- round(CM10.orig[2,2] / nrow(orig_pred) - CM10.orig[2,1] / nrow(orig_pred) * (0.10 / (1 - 0.10)), 3)
    # 20%
    CM20.orig <- confusion.matrix(obs = orig_pred$OutcomeBin, pred = orig_pred$pred.orig, threshold = 0.20)
    res.boot$NB20.orig[i] <- round(CM20.orig[2,2] / nrow(orig_pred) - CM20.orig[2,1] / nrow(orig_pred) * (0.20 / (1 - 0.20)), 3)
    # 30%
    CM30.orig <- confusion.matrix(obs = orig_pred$OutcomeBin, pred = orig_pred$pred.orig, threshold = 0.30)
    res.boot$NB30.orig[i] <- round(CM30.orig[2,2] / nrow(orig_pred) - CM30.orig[2,1] / nrow(orig_pred) * (0.30 / (1 - 0.30)), 3)
    # 40%
    CM40.orig <- confusion.matrix(obs = orig_pred$OutcomeBin, pred = orig_pred$pred.orig, threshold = 0.40)
    res.boot$NB40.orig[i] <- round(CM40.orig[2,2] / nrow(orig_pred) - CM40.orig[2,1] / nrow(orig_pred) * (0.40 / (1 - 0.40)), 3)
    
    ## Calculate Sensitivity and specificity
    # 1%
    res.boot$Sens1.orig[i] <- Sensitivity(pred = pred.orig, outcome = OutcomeBin, threshold = 0.01, data = orig_pred)[1]
    res.boot$Spec1.orig[i] <- Specificity(pred = pred.orig, outcome = OutcomeBin, threshold = 0.01, data = orig_pred)[1]
    # 5%
    res.boot$Sens5.orig[i] <- Sensitivity(pred = pred.orig, outcome = OutcomeBin, threshold = 0.05, data = orig_pred)[1]
    res.boot$Spec5.orig[i] <- Specificity(pred = pred.orig, outcome = OutcomeBin, threshold = 0.05, data = orig_pred)[1]
    # 10%
    res.boot$Sens10.orig[i] <- Sensitivity(pred = pred.orig, outcome = OutcomeBin, threshold = 0.10, data = orig_pred)[1]
    res.boot$Spec10.orig[i] <- Specificity(pred = pred.orig, outcome = OutcomeBin, threshold = 0.10, data = orig_pred)[1]
    # 20%
    res.boot$Sens20.orig[i] <- Sensitivity(pred = pred.orig, outcome = OutcomeBin, threshold = 0.20, data = orig_pred)[1]
    res.boot$Spec20.orig[i] <- Specificity(pred = pred.orig, outcome = OutcomeBin, threshold = 0.20, data = orig_pred)[1]
    # 30%
    res.boot$Sens30.orig[i] <- Sensitivity(pred = pred.orig, outcome = OutcomeBin, threshold = 0.30, data = orig_pred)[1]
    res.boot$Spec30.orig[i] <- Specificity(pred = pred.orig, outcome = OutcomeBin, threshold = 0.30, data = orig_pred)[1]
    # 40%
    res.boot$Sens40.orig[i] <- Sensitivity(pred = pred.orig, outcome = OutcomeBin, threshold = 0.40, data = orig_pred)[1]
    res.boot$Spec40.orig[i] <- Specificity(pred = pred.orig, outcome = OutcomeBin, threshold = 0.40, data = orig_pred)[1]
    # 50%
    res.boot$Sens50.orig[i] <- Sensitivity(pred = pred.orig, outcome = OutcomeBin, threshold = 0.50, data = orig_pred)[1]
    res.boot$Spec50.orig[i] <- Specificity(pred = pred.orig, outcome = OutcomeBin, threshold = 0.50, data = orig_pred)[1]
    
  }
  ## Initial dataset for the comparisons
  comp.boot <- matrix(nrow = 6, ncol = 13)
  colnames(comp.boot) <- c("Comparison", "dAUC.boot", "NB5.boot", "NB10.boot", "NB20.boot", "NB30.boot", "NB40.boot", 
                           "dAUC.orig", "NB5.orig", "NB10.orig", "NB20.orig", "NB30.orig", "NB40.orig")
  rownames(comp.boot) <- c("1vs2", "1vs3", "1vs4", "5vs6", "5vs7", "5vs8")
  comp.boot <- data.frame(comp.boot)
  comp.boot$Comparison <- c("1vs2", "1vs3", "1vs4", "5vs6", "5vs7", "5vs8")
  
  ### Bootstrap sample
  ## Calculate delta AUC
  comp.boot$dAUC.boot[1]    <- dAUC(data = x_boot, predictor1 = M1.boot$predict, predictor2 = M2.boot$predict)[1]
  comp.boot$dAUC.boot[2]    <- dAUC(data = x_boot, predictor1 = M1.boot$predict, predictor2 = M3.boot$predict)[1]
  comp.boot$dAUC.boot[3]    <- dAUC(data = x_boot, predictor1 = M1.boot$predict, predictor2 = M4.boot$predict)[1]
  comp.boot$dAUC.boot[4]    <- dAUC(data = x_boot, predictor1 = M5.boot$predict, predictor2 = M6.boot$predict)[1]
  comp.boot$dAUC.boot[5]    <- dAUC(data = x_boot, predictor1 = M5.boot$predict, predictor2 = M7.boot$predict)[1]
  comp.boot$dAUC.boot[6]    <- dAUC(data = x_boot, predictor1 = M5.boot$predict, predictor2 = M8.boot$predict)[1]
  
  ## Calculate delta Net Benefit
  # 5%
  comp.boot$NB5.boot[1] <- diffNB(data = x_boot, model1 = M1.boot, model2 = M2.boot, predictor1 = M1.boot$predict, predictor2 = M2.boot$predict, outcome = OutcomeBin, threshold = 0.05, CI = FALSE)
  comp.boot$NB5.boot[2] <- diffNB(data = x_boot, model1 = M1.boot, model2 = M3.boot, predictor1 = M1.boot$predict, predictor2 = M3.boot$predict, outcome = OutcomeBin, threshold = 0.05, CI = FALSE)
  comp.boot$NB5.boot[3] <- diffNB(data = x_boot, model1 = M1.boot, model2 = M4.boot, predictor1 = M1.boot$predict, predictor2 = M4.boot$predict, outcome = OutcomeBin, threshold = 0.05, CI = FALSE)
  comp.boot$NB5.boot[4] <- diffNB(data = x_boot, model1 = M5.boot, model2 = M6.boot, predictor1 = M5.boot$predict, predictor2 = M6.boot$predict, outcome = OutcomeBin, threshold = 0.05, CI = FALSE)
  comp.boot$NB5.boot[5] <- diffNB(data = x_boot, model1 = M5.boot, model2 = M7.boot, predictor1 = M5.boot$predict, predictor2 = M7.boot$predict, outcome = OutcomeBin, threshold = 0.05, CI = FALSE)
  comp.boot$NB5.boot[6] <- diffNB(data = x_boot, model1 = M5.boot, model2 = M8.boot, predictor1 = M5.boot$predict, predictor2 = M8.boot$predict, outcome = OutcomeBin, threshold = 0.05, CI = FALSE)
  
  # 10%
  comp.boot$NB10.boot[1] <- diffNB(data = x_boot, model1 = M1.boot, model2 = M2.boot, predictor1 = M1.boot$predict, predictor2 = M2.boot$predict, outcome = OutcomeBin, threshold = 0.10, CI = FALSE)
  comp.boot$NB10.boot[2] <- diffNB(data = x_boot, model1 = M1.boot, model2 = M3.boot, predictor1 = M1.boot$predict, predictor2 = M3.boot$predict, outcome = OutcomeBin, threshold = 0.10, CI = FALSE)
  comp.boot$NB10.boot[3] <- diffNB(data = x_boot, model1 = M1.boot, model2 = M4.boot, predictor1 = M1.boot$predict, predictor2 = M4.boot$predict, outcome = OutcomeBin, threshold = 0.10, CI = FALSE)
  comp.boot$NB10.boot[4] <- diffNB(data = x_boot, model1 = M5.boot, model2 = M6.boot, predictor1 = M5.boot$predict, predictor2 = M6.boot$predict, outcome = OutcomeBin, threshold = 0.10, CI = FALSE)
  comp.boot$NB10.boot[5] <- diffNB(data = x_boot, model1 = M5.boot, model2 = M7.boot, predictor1 = M5.boot$predict, predictor2 = M7.boot$predict, outcome = OutcomeBin, threshold = 0.10, CI = FALSE)
  comp.boot$NB10.boot[6] <- diffNB(data = x_boot, model1 = M5.boot, model2 = M8.boot, predictor1 = M5.boot$predict, predictor2 = M8.boot$predict, outcome = OutcomeBin, threshold = 0.10, CI = FALSE)
  
  # 20%
  comp.boot$NB20.boot[1] <- diffNB(data = x_boot, model1 = M1.boot, model2 = M2.boot, predictor1 = M1.boot$predict, predictor2 = M2.boot$predict, outcome = OutcomeBin, threshold = 0.20, CI = FALSE)
  comp.boot$NB20.boot[2] <- diffNB(data = x_boot, model1 = M1.boot, model2 = M3.boot, predictor1 = M1.boot$predict, predictor2 = M3.boot$predict, outcome = OutcomeBin, threshold = 0.20, CI = FALSE)
  comp.boot$NB20.boot[3] <- diffNB(data = x_boot, model1 = M1.boot, model2 = M4.boot, predictor1 = M1.boot$predict, predictor2 = M4.boot$predict, outcome = OutcomeBin, threshold = 0.20, CI = FALSE)
  comp.boot$NB20.boot[4] <- diffNB(data = x_boot, model1 = M5.boot, model2 = M6.boot, predictor1 = M5.boot$predict, predictor2 = M6.boot$predict, outcome = OutcomeBin, threshold = 0.20, CI = FALSE)
  comp.boot$NB20.boot[5] <- diffNB(data = x_boot, model1 = M5.boot, model2 = M7.boot, predictor1 = M5.boot$predict, predictor2 = M7.boot$predict, outcome = OutcomeBin, threshold = 0.20, CI = FALSE)
  comp.boot$NB20.boot[6] <- diffNB(data = x_boot, model1 = M5.boot, model2 = M8.boot, predictor1 = M5.boot$predict, predictor2 = M8.boot$predict, outcome = OutcomeBin, threshold = 0.20, CI = FALSE)
  
  # 30%
  comp.boot$NB30.boot[1] <- diffNB(data = x_boot, model1 = M1.boot, model2 = M2.boot, predictor1 = M1.boot$predict, predictor2 = M2.boot$predict, outcome = OutcomeBin, threshold = 0.30, CI = FALSE)
  comp.boot$NB30.boot[2] <- diffNB(data = x_boot, model1 = M1.boot, model2 = M3.boot, predictor1 = M1.boot$predict, predictor2 = M3.boot$predict, outcome = OutcomeBin, threshold = 0.30, CI = FALSE)
  comp.boot$NB30.boot[3] <- diffNB(data = x_boot, model1 = M1.boot, model2 = M4.boot, predictor1 = M1.boot$predict, predictor2 = M4.boot$predict, outcome = OutcomeBin, threshold = 0.30, CI = FALSE)
  comp.boot$NB30.boot[4] <- diffNB(data = x_boot, model1 = M5.boot, model2 = M6.boot, predictor1 = M5.boot$predict, predictor2 = M6.boot$predict, outcome = OutcomeBin, threshold = 0.30, CI = FALSE)
  comp.boot$NB30.boot[5] <- diffNB(data = x_boot, model1 = M5.boot, model2 = M7.boot, predictor1 = M5.boot$predict, predictor2 = M7.boot$predict, outcome = OutcomeBin, threshold = 0.30, CI = FALSE)
  comp.boot$NB30.boot[6] <- diffNB(data = x_boot, model1 = M5.boot, model2 = M8.boot, predictor1 = M5.boot$predict, predictor2 = M8.boot$predict, outcome = OutcomeBin, threshold = 0.30, CI = FALSE)
  
  # 40%
  comp.boot$NB40.boot[1] <- diffNB(data = x_boot, model1 = M1.boot, model2 = M2.boot, predictor1 = M1.boot$predict, predictor2 = M2.boot$predict, outcome = OutcomeBin, threshold = 0.40, CI = FALSE)
  comp.boot$NB40.boot[2] <- diffNB(data = x_boot, model1 = M1.boot, model2 = M3.boot, predictor1 = M1.boot$predict, predictor2 = M3.boot$predict, outcome = OutcomeBin, threshold = 0.40, CI = FALSE)
  comp.boot$NB40.boot[3] <- diffNB(data = x_boot, model1 = M1.boot, model2 = M4.boot, predictor1 = M1.boot$predict, predictor2 = M4.boot$predict, outcome = OutcomeBin, threshold = 0.40, CI = FALSE)
  comp.boot$NB40.boot[4] <- diffNB(data = x_boot, model1 = M5.boot, model2 = M6.boot, predictor1 = M5.boot$predict, predictor2 = M6.boot$predict, outcome = OutcomeBin, threshold = 0.40, CI = FALSE)
  comp.boot$NB40.boot[5] <- diffNB(data = x_boot, model1 = M5.boot, model2 = M7.boot, predictor1 = M5.boot$predict, predictor2 = M7.boot$predict, outcome = OutcomeBin, threshold = 0.40, CI = FALSE)
  comp.boot$NB40.boot[6] <- diffNB(data = x_boot, model1 = M5.boot, model2 = M8.boot, predictor1 = M5.boot$predict, predictor2 = M8.boot$predict, outcome = OutcomeBin, threshold = 0.40, CI = FALSE)
  
  ### Original data
  ## Calculate delta AUC
  comp.boot$dAUC.orig[1]    <- dAUC(data = ctDNAimp, predictor1 = M1.boot$pred.orig, predictor2 = M2.boot$pred.orig)[1]
  comp.boot$dAUC.orig[2]    <- dAUC(data = ctDNAimp, predictor1 = M1.boot$pred.orig, predictor2 = M3.boot$pred.orig)[1]
  comp.boot$dAUC.orig[3]    <- dAUC(data = ctDNAimp, predictor1 = M1.boot$pred.orig, predictor2 = M4.boot$pred.orig)[1]
  comp.boot$dAUC.orig[4]    <- dAUC(data = ctDNAimp, predictor1 = M5.boot$pred.orig, predictor2 = M6.boot$pred.orig)[1]
  comp.boot$dAUC.orig[5]    <- dAUC(data = ctDNAimp, predictor1 = M5.boot$pred.orig, predictor2 = M7.boot$pred.orig)[1]
  comp.boot$dAUC.orig[6]    <- dAUC(data = ctDNAimp, predictor1 = M5.boot$pred.orig, predictor2 = M8.boot$pred.orig)[1]
  
  ## Calculate delta Net Benefit
  # 5%
  comp.boot$NB5.orig[1] <- diffNB(data = ctDNAimp, model1 = M1.boot, model2 = M2.boot, predictor1 = M1.boot$pred.orig, predictor2 = M2.boot$pred.orig, outcome = OutcomeBin, threshold = 0.05, CI = FALSE)
  comp.boot$NB5.orig[2] <- diffNB(data = ctDNAimp, model1 = M1.boot, model2 = M3.boot, predictor1 = M1.boot$pred.orig, predictor2 = M3.boot$pred.orig, outcome = OutcomeBin, threshold = 0.05, CI = FALSE)
  comp.boot$NB5.orig[3] <- diffNB(data = ctDNAimp, model1 = M1.boot, model2 = M4.boot, predictor1 = M1.boot$pred.orig, predictor2 = M4.boot$pred.orig, outcome = OutcomeBin, threshold = 0.05, CI = FALSE)
  comp.boot$NB5.orig[4] <- diffNB(data = ctDNAimp, model1 = M5.boot, model2 = M6.boot, predictor1 = M5.boot$pred.orig, predictor2 = M6.boot$pred.orig, outcome = OutcomeBin, threshold = 0.05, CI = FALSE)
  comp.boot$NB5.orig[5] <- diffNB(data = ctDNAimp, model1 = M5.boot, model2 = M7.boot, predictor1 = M5.boot$pred.orig, predictor2 = M7.boot$pred.orig, outcome = OutcomeBin, threshold = 0.05, CI = FALSE)
  comp.boot$NB5.orig[6] <- diffNB(data = ctDNAimp, model1 = M5.boot, model2 = M8.boot, predictor1 = M5.boot$pred.orig, predictor2 = M8.boot$pred.orig, outcome = OutcomeBin, threshold = 0.05, CI = FALSE)
  
  # 10%
  comp.boot$NB10.orig[1] <- diffNB(data = ctDNAimp, model1 = M1.boot, model2 = M2.boot, predictor1 = M1.boot$pred.orig, predictor2 = M2.boot$pred.orig, outcome = OutcomeBin, threshold = 0.10, CI = FALSE)
  comp.boot$NB10.orig[2] <- diffNB(data = ctDNAimp, model1 = M1.boot, model2 = M3.boot, predictor1 = M1.boot$pred.orig, predictor2 = M3.boot$pred.orig, outcome = OutcomeBin, threshold = 0.10, CI = FALSE)
  comp.boot$NB10.orig[3] <- diffNB(data = ctDNAimp, model1 = M1.boot, model2 = M4.boot, predictor1 = M1.boot$pred.orig, predictor2 = M4.boot$pred.orig, outcome = OutcomeBin, threshold = 0.10, CI = FALSE)
  comp.boot$NB10.orig[4] <- diffNB(data = ctDNAimp, model1 = M5.boot, model2 = M6.boot, predictor1 = M5.boot$pred.orig, predictor2 = M6.boot$pred.orig, outcome = OutcomeBin, threshold = 0.10, CI = FALSE)
  comp.boot$NB10.orig[5] <- diffNB(data = ctDNAimp, model1 = M5.boot, model2 = M7.boot, predictor1 = M5.boot$pred.orig, predictor2 = M7.boot$pred.orig, outcome = OutcomeBin, threshold = 0.10, CI = FALSE)
  comp.boot$NB10.orig[6] <- diffNB(data = ctDNAimp, model1 = M5.boot, model2 = M8.boot, predictor1 = M5.boot$pred.orig, predictor2 = M8.boot$pred.orig, outcome = OutcomeBin, threshold = 0.10, CI = FALSE)
  
  # 20%
  comp.boot$NB20.orig[1] <- diffNB(data = ctDNAimp, model1 = M1.boot, model2 = M2.boot, predictor1 = M1.boot$pred.orig, predictor2 = M2.boot$pred.orig, outcome = OutcomeBin, threshold = 0.20, CI = FALSE)
  comp.boot$NB20.orig[2] <- diffNB(data = ctDNAimp, model1 = M1.boot, model2 = M3.boot, predictor1 = M1.boot$pred.orig, predictor2 = M3.boot$pred.orig, outcome = OutcomeBin, threshold = 0.20, CI = FALSE)
  comp.boot$NB20.orig[3] <- diffNB(data = ctDNAimp, model1 = M1.boot, model2 = M4.boot, predictor1 = M1.boot$pred.orig, predictor2 = M4.boot$pred.orig, outcome = OutcomeBin, threshold = 0.20, CI = FALSE)
  comp.boot$NB20.orig[4] <- diffNB(data = ctDNAimp, model1 = M5.boot, model2 = M6.boot, predictor1 = M5.boot$pred.orig, predictor2 = M6.boot$pred.orig, outcome = OutcomeBin, threshold = 0.20, CI = FALSE)
  comp.boot$NB20.orig[5] <- diffNB(data = ctDNAimp, model1 = M5.boot, model2 = M7.boot, predictor1 = M5.boot$pred.orig, predictor2 = M7.boot$pred.orig, outcome = OutcomeBin, threshold = 0.20, CI = FALSE)
  comp.boot$NB20.orig[6] <- diffNB(data = ctDNAimp, model1 = M5.boot, model2 = M8.boot, predictor1 = M5.boot$pred.orig, predictor2 = M8.boot$pred.orig, outcome = OutcomeBin, threshold = 0.20, CI = FALSE)
  
  # 30%
  comp.boot$NB30.orig[1] <- diffNB(data = ctDNAimp, model1 = M1.boot, model2 = M2.boot, predictor1 = M1.boot$pred.orig, predictor2 = M2.boot$pred.orig, outcome = OutcomeBin, threshold = 0.30, CI = FALSE)
  comp.boot$NB30.orig[2] <- diffNB(data = ctDNAimp, model1 = M1.boot, model2 = M3.boot, predictor1 = M1.boot$pred.orig, predictor2 = M3.boot$pred.orig, outcome = OutcomeBin, threshold = 0.30, CI = FALSE)
  comp.boot$NB30.orig[3] <- diffNB(data = ctDNAimp, model1 = M1.boot, model2 = M4.boot, predictor1 = M1.boot$pred.orig, predictor2 = M4.boot$pred.orig, outcome = OutcomeBin, threshold = 0.30, CI = FALSE)
  comp.boot$NB30.orig[4] <- diffNB(data = ctDNAimp, model1 = M5.boot, model2 = M6.boot, predictor1 = M5.boot$pred.orig, predictor2 = M6.boot$pred.orig, outcome = OutcomeBin, threshold = 0.30, CI = FALSE)
  comp.boot$NB30.orig[5] <- diffNB(data = ctDNAimp, model1 = M5.boot, model2 = M7.boot, predictor1 = M5.boot$pred.orig, predictor2 = M7.boot$pred.orig, outcome = OutcomeBin, threshold = 0.30, CI = FALSE)
  comp.boot$NB30.orig[6] <- diffNB(data = ctDNAimp, model1 = M5.boot, model2 = M8.boot, predictor1 = M5.boot$pred.orig, predictor2 = M8.boot$pred.orig, outcome = OutcomeBin, threshold = 0.30, CI = FALSE)
  
  # 40%
  comp.boot$NB40.orig[1] <- diffNB(data = ctDNAimp, model1 = M1.boot, model2 = M2.boot, predictor1 = M1.boot$pred.orig, predictor2 = M2.boot$pred.orig, outcome = OutcomeBin, threshold = 0.40, CI = FALSE)
  comp.boot$NB40.orig[2] <- diffNB(data = ctDNAimp, model1 = M1.boot, model2 = M3.boot, predictor1 = M1.boot$pred.orig, predictor2 = M3.boot$pred.orig, outcome = OutcomeBin, threshold = 0.40, CI = FALSE)
  comp.boot$NB40.orig[3] <- diffNB(data = ctDNAimp, model1 = M1.boot, model2 = M4.boot, predictor1 = M1.boot$pred.orig, predictor2 = M4.boot$pred.orig, outcome = OutcomeBin, threshold = 0.40, CI = FALSE)
  comp.boot$NB40.orig[4] <- diffNB(data = ctDNAimp, model1 = M5.boot, model2 = M6.boot, predictor1 = M5.boot$pred.orig, predictor2 = M6.boot$pred.orig, outcome = OutcomeBin, threshold = 0.40, CI = FALSE)
  comp.boot$NB40.orig[5] <- diffNB(data = ctDNAimp, model1 = M5.boot, model2 = M7.boot, predictor1 = M5.boot$pred.orig, predictor2 = M7.boot$pred.orig, outcome = OutcomeBin, threshold = 0.40, CI = FALSE)
  comp.boot$NB40.orig[6] <- diffNB(data = ctDNAimp, model1 = M5.boot, model2 = M8.boot, predictor1 = M5.boot$pred.orig, predictor2 = M8.boot$pred.orig, outcome = OutcomeBin, threshold = 0.40, CI = FALSE)
  
  res_boot[[i_boot]]  <- res.boot
  comp_boot[[i_boot]] <- comp.boot
  

}

stopCluster(cl)

print(Sys.time() - strt)

## Calculate the optimism-corrected results
res_long <- rbindlist(res_boot)
comp_long <- rbindlist(comp_boot)

# AUC
MA.AUC$Optimism <- Optimism(orig.var = res_long$AUC.orig, boot.var = res_long$AUC.boot, Groups = res_long$Model)$Optimism.mean
MA.AUC$cAUC <- paste0(round(MA.AUC$AUC - MA.AUC$Optimism, 2), " (", round(MA.AUC$LL_AUC - MA.AUC$Optimism, 2), " to ", round(MA.AUC$UL_AUC - MA.AUC$Optimism, 2), ")")

# dAUC
MA.diffAUC$Optimism <- Optimism(orig.var = comp_long$dAUC.orig, boot.var = comp_long$dAUC.boot, Groups = comp_long$Comparison)$Optimism.mean
MA.diffAUC$cdAUC <- paste0(round(MA.diffAUC$dAUC - MA.diffAUC$Optimism, 3), " (", round(MA.diffAUC$LL_dAUC - MA.diffAUC$Optimism, 3), " to ", round(MA.diffAUC$UL_dAUC - MA.diffAUC$Optimism, 3), ")")
  
# NB
MA.NB5$Optimism  <- Optimism(orig.var = res_long$NB5.orig, boot.var = res_long$NB5.boot, Groups = res_long$Model)$Optimism.mean
MA.NB5$cNB  <- round(MA.NB5$NB - MA.NB5$Optimism, 3)
MA.NB10$Optimism <- Optimism(orig.var = res_long$NB10.orig, boot.var = res_long$NB10.boot, Groups = res_long$Model)$Optimism.mean
MA.NB10$cNB <- round(MA.NB10$NB - MA.NB10$Optimism, 3)
MA.NB20$Optimism <- Optimism(orig.var = res_long$NB20.orig, boot.var = res_long$NB20.boot, Groups = res_long$Model)$Optimism.mean
MA.NB20$cNB <- round(MA.NB20$NB - MA.NB20$Optimism, 3)
MA.NB30$Optimism <- Optimism(orig.var = res_long$NB30.orig, boot.var = res_long$NB30.boot, Groups = res_long$Model)$Optimism.mean
MA.NB30$cNB <- round(MA.NB30$NB - MA.NB30$Optimism, 3)
MA.NB40$Optimism <- Optimism(orig.var = res_long$NB40.orig, boot.var = res_long$NB40.boot, Groups = res_long$Model)$Optimism.mean
MA.NB40$cNB <- round(MA.NB40$NB - MA.NB40$Optimism, 3)

# Plot for NB
MA.NB5$Threshold  <- 0.05
MA.NB10$Threshold <- 0.10
MA.NB20$Threshold <- 0.20
MA.NB30$Threshold <- 0.30
MA.NB40$Threshold <- 0.40

MA.NB <- rbind(MA.NB5, MA.NB10, MA.NB20, MA.NB30, MA.NB40)

MA.NB$Group <- ifelse(MA.NB$Model %in% c("M1", "M2", "M3", "M4"), "Without CA125", "With CA125")
MA.NB$Group <- factor(MA.NB$Group, levels = c("Without CA125", "With CA125"))
MA.NB$TreatAll <- nrow(ctDNAimp[ctDNAimp$OutcomeBin == 1,]) / nrow(ctDNAimp) - nrow(ctDNAimp[ctDNAimp$OutcomeBin == 0,]) / nrow(ctDNAimp) * (MA.NB$Threshold / (1 - MA.NB$Threshold))

TA1 <- MA.NB[MA.NB$Model == "M1",]
TA2 <- MA.NB[MA.NB$Model == "M5",]
TA <- rbind(TA1, TA2)
TA$Model <- "Treat all"

cols <- c("forestgreen", "firebrick3", "royalblue2", "gold2",
          "springgreen3", "firebrick1", "royalblue1", "goldenrod1")

ggplot(data = MA.NB, aes(x = Threshold, y = cNB, colour = Model)) +
  geom_line(linewidth = 0.75) +
  xlab("Risk threshold") +
  ylab("Net Benefit") +
  theme_minimal() +
  theme(axis.line = element_line(colour = "black")) + 
  theme(panel.spacing = unit(15, "pt"), strip.text = element_text(size = 14),
        axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        legend.title = element_text(size = 12), legend.text = element_text(size = 12)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.5)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.4)) +
  geom_line(data = TA, aes(x = Threshold, y = TreatAll, colour = Model), linewidth = 0.75) +
  scale_color_manual(name = "Model",
                     breaks = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "Treat all"),
                     values = c("M1" = "forestgreen", "M2" = "firebrick3", "M3" = "royalblue2", "M4" = "gold2", "M5" = "springgreen3", "M6" = "firebrick1", "M7" = "royalblue1", "M8" = "goldenrod1", "Treat all" = "grey50")) +
  facet_wrap(~ Group)
# 1000 x 600


# Legend for each facet
library(gridExtra) 
xs <- split(MA.NB,f = MA.NB$Group)
ys <- split(TA,f = TA$Group)

p1 <- ggplot(data = xs$`Without CA125`, aes(x = Threshold, y = cNB, linetype = Model, colour = Model)) +
  geom_line(linewidth = 0.75) +
  xlab("Risk threshold") +
  ylab("Net Benefit") +
  theme_minimal() + 
  theme(axis.line = element_line(colour = "black")) + 
  theme(panel.spacing = unit(15, "pt"), strip.text = element_text(size = 14),
        axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        legend.title = element_text(size = 12), legend.text = element_text(size = 12)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.5)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.4)) +
  geom_line(data = ys$`Without CA125`, aes(x = Threshold, y = TreatAll, linetype = Model, colour = Model), linewidth = 0.75) +
  scale_color_manual(name = "",
                     breaks = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "Treat all"),
                     values = c("gray5", "gray20", "gray10", "gray25", "gray5", "gray20", "gray10", "gray25", "grey50")) +
  scale_linetype_manual(name = "", 
                        breaks = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "Treat all"),
                        values = c("solid", "dotted", "dotdash", "longdash", "solid", "dotted", "dotdash", "longdash", "solid")) +
  facet_wrap(~ Group, nrow = 1)

p2 <- ggplot(data = xs$`With CA125`, aes(x = Threshold, y = cNB, linetype = Model, colour = Model)) +
  geom_line(linewidth = 0.75) +
  xlab("Risk threshold") +
  ylab("Net Benefit") +
  theme_minimal() +
  theme(axis.line = element_line(colour = "black")) + 
  theme(panel.spacing = unit(15, "pt"), strip.text = element_text(size = 14),
        axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        legend.title = element_text(size = 12), legend.text = element_text(size = 12)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.5)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.4)) +
  geom_line(data = ys$`With CA125`, aes(x = Threshold, y = TreatAll, linetype = Model, colour = Model), linewidth = 0.75) +
  scale_color_manual(name = "",
                     breaks = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "Treat all"),
                     values = c("gray5", "gray20", "gray10", "gray25", "gray5", "gray20", "gray10", "gray25", "grey50")) +
  scale_linetype_manual(name = "", 
                        breaks = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "Treat all"),
                        values = c("solid", "dotted", "dotdash", "longdash", "solid", "dotted", "dotdash", "longdash", "solid")) +
  facet_wrap(~ Group, nrow = 1)

grid.arrange(p1,p2, nrow = 1)
# 1000 x 600

# dNB & TT
MA.diffNB5$Optimism  <- Optimism(orig.var = comp_long$NB5.orig, boot.var = comp_long$NB5.boot, Groups = comp_long$Comparison)$Optimism.mean
MA.diffNB5$cdNB  <- MA.diffNB5$dNB - MA.diffNB5$Optimism
MA.diffNB5$cTTO  <- round(1/MA.diffNB5$cdNB, 2)
MA.diffNB5$cdNB  <- paste0(round(MA.diffNB5$dNB - MA.diffNB5$Optimism, 4), " (", round(MA.diffNB5$LL_dNB - MA.diffNB5$Optimism, 4), " to ", round(MA.diffNB5$UL_dNB - MA.diffNB5$Optimism, 4), ")")
MA.diffNB10$Optimism <- Optimism(orig.var = comp_long$NB10.orig, boot.var = comp_long$NB10.boot, Groups = comp_long$Comparison)$Optimism.mean
MA.diffNB10$cdNB <- MA.diffNB10$dNB - MA.diffNB10$Optimism
MA.diffNB10$cTTO <- round(1/MA.diffNB10$cdNB, 2)
MA.diffNB10$cdNB <- paste0(round(MA.diffNB10$dNB - MA.diffNB10$Optimism, 4), " (", round(MA.diffNB10$LL_dNB - MA.diffNB10$Optimism, 4), " to ", round(MA.diffNB10$UL_dNB - MA.diffNB10$Optimism, 4), ")")
MA.diffNB20$Optimism <- Optimism(orig.var = comp_long$NB20.orig, boot.var = comp_long$NB20.boot, Groups = comp_long$Comparison)$Optimism.mean
MA.diffNB20$cdNB <- MA.diffNB20$dNB - MA.diffNB20$Optimism
MA.diffNB20$cTTO <- round(1/MA.diffNB20$cdNB, 2)
MA.diffNB20$cdNB <- paste0(round(MA.diffNB20$dNB - MA.diffNB20$Optimism, 4), " (", round(MA.diffNB20$LL_dNB - MA.diffNB20$Optimism, 4), " to ", round(MA.diffNB20$UL_dNB - MA.diffNB20$Optimism, 4), ")")
MA.diffNB30$Optimism <- Optimism(orig.var = comp_long$NB30.orig, boot.var = comp_long$NB30.boot, Groups = comp_long$Comparison)$Optimism.mean
MA.diffNB30$cdNB <- MA.diffNB30$dNB - MA.diffNB30$Optimism
MA.diffNB30$cTTO <- round(1/MA.diffNB30$cdNB, 2)
MA.diffNB30$cdNB <- paste0(round(MA.diffNB30$dNB - MA.diffNB30$Optimism, 4), " (", round(MA.diffNB30$LL_dNB - MA.diffNB30$Optimism, 4), " to ", round(MA.diffNB30$UL_dNB - MA.diffNB30$Optimism, 4), ")")
MA.diffNB40$Optimism <- Optimism(orig.var = comp_long$NB40.orig, boot.var = comp_long$NB40.boot, Groups = comp_long$Comparison)$Optimism.mean
MA.diffNB40$cdNB <- MA.diffNB40$dNB - MA.diffNB40$Optimism
MA.diffNB40$cTTO <- round(1/MA.diffNB40$cdNB, 2)
MA.diffNB40$cdNB <- paste0(round(MA.diffNB40$dNB - MA.diffNB40$Optimism, 4), " (", round(MA.diffNB40$LL_dNB - MA.diffNB40$Optimism, 4), " to ", round(MA.diffNB40$UL_dNB - MA.diffNB40$Optimism, 4), ")")

# Plot for Test trade-off
MA.diffNB5$Threshold  <- 0.05
MA.diffNB10$Threshold <- 0.10
MA.diffNB20$Threshold <- 0.20
MA.diffNB30$Threshold <- 0.30
MA.diffNB40$Threshold <- 0.40

MA.diffNB <- rbind(MA.diffNB5, MA.diffNB10, MA.diffNB20, MA.diffNB30, MA.diffNB40)

ggplot(data = MA.diffNB, aes(x = Threshold, y = cTTO, colour = Comparison)) +
  geom_line() +
  xlab("Risk threshold") +
  ylab("Test trade-off") +
  theme_minimal() +
  theme(axis.line = element_line(colour = "black")) + 
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0))
# 800 x 400

# Sensitivity and specificity
MA.SS1$SensOptimism  <- Optimism(orig.var = res_long$Sens1.orig, boot.var = res_long$Sens1.boot, Groups = res_long$Model)$Optimism.mean
MA.SS1$cSens <- paste0(round((MA.SS1$Sens - MA.SS1$SensOptimism) * 100, 1), " (", round((MA.SS1$LL_Sens - MA.SS1$SensOptimism) * 100, 1), " to ", round((MA.SS1$UL_Sens - MA.SS1$SensOptimism) * 100, 1), ")")
MA.SS1$SpecOptimism  <- Optimism(orig.var = res_long$Spec1.orig, boot.var = res_long$Spec1.boot, Groups = res_long$Model)$Optimism.mean
MA.SS1$cSpec <- paste0(round((MA.SS1$Spec - MA.SS1$SpecOptimism) * 100, 1), " (", round((MA.SS1$LL_Spec - MA.SS1$SpecOptimism) * 100, 1), " to ", round((MA.SS1$UL_Spec - MA.SS1$SpecOptimism) * 100, 1), ")")

MA.SS5$SensOptimism  <- Optimism(orig.var = res_long$Sens5.orig, boot.var = res_long$Sens5.boot, Groups = res_long$Model)$Optimism.mean
MA.SS5$cSens <- paste0(round((MA.SS5$Sens - MA.SS5$SensOptimism) * 100, 1), " (", round((MA.SS5$LL_Sens - MA.SS5$SensOptimism) * 100, 1), " to ", round((MA.SS5$UL_Sens - MA.SS5$SensOptimism) * 100, 1), ")")
MA.SS5$SpecOptimism  <- Optimism(orig.var = res_long$Spec5.orig, boot.var = res_long$Spec5.boot, Groups = res_long$Model)$Optimism.mean
MA.SS5$cSpec <- paste0(round((MA.SS5$Spec - MA.SS5$SpecOptimism) * 100, 1), " (", round((MA.SS5$LL_Spec - MA.SS5$SpecOptimism) * 100, 1), " to ", round((MA.SS5$UL_Spec - MA.SS5$SpecOptimism) * 100, 1), ")")

MA.SS10$SensOptimism <- Optimism(orig.var = res_long$Sens10.orig, boot.var = res_long$Sens10.boot, Groups = res_long$Model)$Optimism.mean
MA.SS10$cSens <- paste0(round((MA.SS10$Sens - MA.SS10$SensOptimism) * 100, 1), " (", round((MA.SS10$LL_Sens - MA.SS10$SensOptimism) * 100, 1), " to ", round((MA.SS10$UL_Sens - MA.SS10$SensOptimism) * 100, 1), ")")
MA.SS10$SpecOptimism <- Optimism(orig.var = res_long$Spec10.orig, boot.var = res_long$Spec10.boot, Groups = res_long$Model)$Optimism.mean
MA.SS10$cSpec <- paste0(round((MA.SS10$Spec - MA.SS10$SpecOptimism) * 100, 1), " (", round((MA.SS10$LL_Spec - MA.SS10$SpecOptimism) * 100, 1), " to ", round((MA.SS10$UL_Spec - MA.SS10$SpecOptimism) * 100, 1), ")")

MA.SS20$SensOptimism <- Optimism(orig.var = res_long$Sens20.orig, boot.var = res_long$Sens20.boot, Groups = res_long$Model)$Optimism.mean
MA.SS20$cSens <- paste0(round((MA.SS20$Sens - MA.SS20$SensOptimism) * 100, 1), " (", round((MA.SS20$LL_Sens - MA.SS20$SensOptimism) * 100, 1), " to ", round((MA.SS20$UL_Sens - MA.SS20$SensOptimism) * 100, 1), ")")
MA.SS20$SpecOptimism <- Optimism(orig.var = res_long$Spec20.orig, boot.var = res_long$Spec20.boot, Groups = res_long$Model)$Optimism.mean
MA.SS20$cSpec <- paste0(round((MA.SS20$Spec - MA.SS20$SpecOptimism) * 100, 1), " (", round((MA.SS20$LL_Spec - MA.SS20$SpecOptimism) * 100, 1), " to ", round((MA.SS20$UL_Spec - MA.SS20$SpecOptimism) * 100, 1), ")")

MA.SS30$SensOptimism <- Optimism(orig.var = res_long$Sens30.orig, boot.var = res_long$Sens30.boot, Groups = res_long$Model)$Optimism.mean
MA.SS30$cSens <- paste0(round((MA.SS30$Sens - MA.SS30$SensOptimism) * 100, 1), " (", round((MA.SS30$LL_Sens - MA.SS30$SensOptimism) * 100, 1), " to ", round((MA.SS30$UL_Sens - MA.SS30$SensOptimism) * 100, 1), ")")
MA.SS30$SpecOptimism <- Optimism(orig.var = res_long$Spec30.orig, boot.var = res_long$Spec30.boot, Groups = res_long$Model)$Optimism.mean
MA.SS30$cSpec <- paste0(round((MA.SS30$Spec - MA.SS30$SpecOptimism) * 100, 1), " (", round((MA.SS30$LL_Spec - MA.SS30$SpecOptimism) * 100, 1), " to ", round((MA.SS30$UL_Spec - MA.SS30$SpecOptimism) * 100, 1), ")")

MA.SS40$SensOptimism <- Optimism(orig.var = res_long$Sens40.orig, boot.var = res_long$Sens40.boot, Groups = res_long$Model)$Optimism.mean
MA.SS40$cSens <- paste0(round((MA.SS40$Sens - MA.SS40$SensOptimism) * 100, 1), " (", round((MA.SS40$LL_Sens - MA.SS40$SensOptimism) * 100, 1), " to ", round((MA.SS40$UL_Sens - MA.SS40$SensOptimism) * 100, 1), ")")
MA.SS40$SpecOptimism <- Optimism(orig.var = res_long$Spec40.orig, boot.var = res_long$Spec40.boot, Groups = res_long$Model)$Optimism.mean
MA.SS40$cSpec <- paste0(round((MA.SS40$Spec - MA.SS40$SpecOptimism) * 100, 1), " (", round((MA.SS40$LL_Spec - MA.SS40$SpecOptimism) * 100, 1), " to ", round((MA.SS40$UL_Spec - MA.SS40$SpecOptimism) * 100, 1), ")")

MA.SS50$SensOptimism <- Optimism(orig.var = res_long$Sens50.orig, boot.var = res_long$Sens50.boot, Groups = res_long$Model)$Optimism.mean
MA.SS50$cSens <- paste0(round((MA.SS50$Sens - MA.SS50$SensOptimism) * 100, 1), " (", round((MA.SS50$LL_Sens - MA.SS50$SensOptimism) * 100, 1), " to ", round((MA.SS50$UL_Sens - MA.SS50$SensOptimism) * 100, 1), ")")
MA.SS50$SpecOptimism <- Optimism(orig.var = res_long$Spec50.orig, boot.var = res_long$Spec50.boot, Groups = res_long$Model)$Optimism.mean
MA.SS50$cSpec <- paste0(round((MA.SS50$Spec - MA.SS50$SpecOptimism) * 100, 1), " (", round((MA.SS50$LL_Spec - MA.SS50$SpecOptimism) * 100, 1), " to ", round((MA.SS50$UL_Spec - MA.SS50$SpecOptimism) * 100, 1), ")")

# Plot for sensitivity and specificity
MA.SS1$threshold  <- 0.01
MA.SS5$threshold  <- 0.05
MA.SS10$threshold <- 0.10
MA.SS20$threshold <- 0.20
MA.SS30$threshold <- 0.30
MA.SS40$threshold <- 0.40
MA.SS50$threshold <- 0.50

MA.SS <- rbind(MA.SS1, MA.SS5, MA.SS10, MA.SS20, MA.SS30, MA.SS40, MA.SS50)
MA.SS$cSensitivity <- MA.SS$Sens - MA.SS$SensOptimism
MA.SS$CSPecificity <- MA.SS$Spec - MA.SS$SpecOptimism

MA.Sens <- MA.SS[, c("threshold", "cSensitivity", "Model")]
colnames(MA.Sens) <- c("threshold", "Value", "Model")
MA.Sens$Measure <- "Sensitivity"
MA.Spec <- MA.SS[, c("threshold", "CSPecificity", "Model")]
colnames(MA.Spec) <- c("threshold", "Value", "Model")
MA.Spec$Measure <- "Specificity"
SS.plot <- rbind(MA.Sens, MA.Spec)

SS.plot$Group <- ifelse(SS.plot$Model %in% c("M1", "M2", "M3", "M4"), "Without CA125", "With CA125")
SS.plot$Group <- factor(SS.plot$Group, levels = c("Without CA125", "With CA125"))

ggplot(data = SS.plot, aes(x = threshold, y = Value * 100, colour = Model)) +
  geom_line(linewidth = 0.65, aes(linetype = Measure)) +
  xlab("Risk threshold") +
  ylab("Sensitivity/Specificity") +
  theme_minimal() +
  theme(axis.line = element_line(colour = "black")) + 
  theme(panel.spacing = unit(15, "pt"), strip.text = element_text(size = 14),
        axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        legend.title = element_text(size = 12), legend.text = element_text(size = 12)) +
  scale_linetype_manual(values=c("solid", "longdash")) +
  scale_color_manual(values = cols) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  scale_x_continuous(expand = c(0, 0)) +
  facet_wrap(~ Group)
# 1000 x 600
