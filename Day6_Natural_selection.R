
AFR_EUR <- read.delim("AFR_EUR.weir.fst", header=TRUE)

AFR_EAS <- read.delim("AFR_EAS.weir.fst", header=TRUE)

EAS_EUR <- read.delim("EAS_EUR.weir.fst", header=TRUE)
head(AFR_EAS)

AFR_EUR_no_duplictes<- AFR_EUR[!duplicated(AFR_EUR$POS), ]
row_number(AFR_EUR)
row_number(AFR_EUR_no_duplictes)

###There are not duplicates###

is.na(AFR_EUR$WEIR_AND_COCKERHAM_FST)
AFR_EUR_na_counts<- colSums(is.na(AFR_EUR))
na.omit(AFR_EUR) -> AFR_EUR
AFR_EUR_na_counts<- colSums(is.na(AFR_EUR))

### lest do it for all the data sets
na.omit(AFR_EAS) -> AFR_EAS
colSums(is.na(AFR_EAS))

na.omit(EAS_EUR) -> EAS_EUR
colSums(is.na(EAS_EUR))


afr_eur_rmNa_com_cor <- afr_eur_rmNa_com %>% mutate(WEIR_AND_COCKERHAM_FST = ifelse(WEIR_AND_COCKERHAM_FST < 0, 0, WEIR_AND_COCKERHAM_FST))##### filter the files by the overlapping SNP ####

common_ids <- Reduce(intersect, list(AFR_EAS$POS, AFR_EUR$POS, EAS_EUR$POS))

AFR_EAS_com <- AFR_EAS %>% filter(POS %in% (common_ids))
AFR_EUR_com <- AFR_EUR %>% filter(POS %in% (common_ids))
EAS_EUR_com <- EAS_EUR %>% filter(POS %in% (common_ids))

#check that all of the df now have the same number of rows

nrow(AFR_EAS_com) * 2 == (nrow(AFR_EUR_com) + nrow(EAS_EUR_com))

AFR_EUR_com <- AFR_EUR_com %>% mutate(WEIR_AND_COCKERHAM_FST = ifelse(WEIR_AND_COCKERHAM_FST < 0, 0, WEIR_AND_COCKERHAM_FST))

AFR_EAS_com <- AFR_EAS_com %>% mutate(WEIR_AND_COCKERHAM_FST = ifelse(WEIR_AND_COCKERHAM_FST < 0, 0, WEIR_AND_COCKERHAM_FST))

EAS_EUR_com <- EAS_EUR_com %>% mutate(WEIR_AND_COCKERHAM_FST = ifelse(WEIR_AND_COCKERHAM_FST < 0, 0, WEIR_AND_COCKERHAM_FST))

### lest check if it worked ####

#sense check this
min(AFR_EAS_com$WEIR_AND_COCKERHAM_FST) # < 0
min(AFR_EAS_com$WEIR_AND_COCKERHAM_FST) # == 0

#### 3 

plot(AFR_EAS_com$WEIR_AND_COCKERHAM_FST)

AFR_EAS_com_109513601 <- subset(AFR_EAS_com, POS==109513601)
AFR_EUR_com_109513601 <- subset(AFR_EUR_com, POS==109513601)
EAS_EUR_com_109513601 <- subset(EAS_EUR_com, POS==109513601)

### high freqeuncy in AFR_EAS and EAS_EUR ####
### AFr_EAS and EAS_EUR fourth quartile, AFR_EUR

AFR_EAS_com_chro_seg <- AFR_EAS_com$POS[109503601 : 109523601]

### 10000 position around the 109513601 #####
AFR_EAS_com_chro_seg <- subset(AFR_EAS_com,AFR_EAS_com$POS>109503601 & AFR_EAS_com$POS<=109523601)
plot(AFR_EAS_com_chro_seg$WEIR_AND_COCKERHAM_FST, main="AFRICA VS EAST")

ggplot(AFR_EAS_com_chro_seg, aes(x=POS, y=WEIR_AND_COCKERHAM_FST)) + geom_point() + stat_quantile(quantiles =0.95)

AFR_EUR_com_chro_seg <- subset(AFR_EUR_com,AFR_EAS_com$POS>109503601 & AFR_EUR_com$POS<=109523601)
ggplot(AFR_EUR_com_chro_seg, aes(x=POS, y=WEIR_AND_COCKERHAM_FST)) + geom_point() + stat_quantile(quantiles =0.95)

EAS_EUR_com_chro_seg <- subset(EAS_EUR_com,EAS_EUR_com$POS>109503601 & EAS_EUR_com$POS<=109523601)
ggplot(EAS_EUR_com_chro_seg, aes(x=POS, y=WEIR_AND_COCKERHAM_FST, col=cut)) + geom_point() + stat_quantile(quantiles = 0.95)


########## PART II #######

install.packages("rehh")
library("rehh")
EDAR_CHS <- data2haplohh(hap_file = "Chr2_EDAR_CHS_500K.recode.vcf", polarize_vcf = F )

EDAR_LW <- data2haplohh(hap_file = "Chr2_EDAR_LWK_500K.recode.vcf", polarize_vcf = F )


Africa_cal_ehh <- calc_ehh(EDAR_CHS, mrk = "rs3827760")
Asia_cal_ehh <- calc_ehh(EDAR_LW, mrk = "rs3827760")

plot(Africa_cal_ehh, main= "Africa")
plot(Asia_cal_ehh, main = "ASIA")


Africa_furcation <- calc_furcation(EDAR_CHS, mrk = "rs3827760")
plot(Africa_furcation, main= "Africa")

Asian_furcation <- calc_furcation(EDAR_LW, mrk = "rs3827760")
plot(Asian_furcation, main = "ASIA")

afr_all_mkr <- EDAR_CHS@positions[2,1]

### calculation for all the values ###
Africa <- scan_hh(EDAR_CHS)

Asia <- scan_hh(EDAR_LW)

#### check in the selection loci #####

Africa[Africa$POSITION==109513601,]

Asia[Asia$POSITION==109513601,]

### Estimate the iHS in AFR ###

iHS.AFR<-ihh2ihs(Africa, min_maf = 0.02, freqbin = 0.01)

iHS.ASI <-ihh2ihs(Asia, min_maf = 0.02, freqbin = 0.01)

