library(openxlsx)
library("readxl")
library(tidyverse)
library(data.table)
options(scipen = 999)

# Lookup All Species-Antibacterials Combinations For Resistance Cutoffs ------------

# my_data <- read_excel(path = "ANAEuROBE_dataset_matteo_only.xlsx",  sheet="3.Completo")
# 
# names(my_data)
# 
# temp <- my_data %>% select(`Code event`:`Species identification`) %>%
#   bind_cols(
#     my_data %>% select(contains(" MIC"))
#   )
# 
# names(temp)
# 
# temp$`Gram pos bacilli` <- as.numeric(temp$`Gram pos bacilli`)
# temp$`Gram neg bacilli` <- as.numeric(temp$`Gram neg bacilli`)
# temp$`Gram pos cocci` <- as.numeric(temp$`Gram pos cocci`)
# temp$`Gram neg cocci` <- as.numeric(temp$`Gram neg cocci`)
# 
# temp <- temp %>% select(-c(`AST by ANA-ATB microdilution=1`,`AST by MICRONAUT-S MIC test=1`)) %>%
#   gather(Abx, MIC, `Benzilpenicillin MIC`:`Oflox MIC`) 
# 
# temp <- temp %>% select(`Species identification`, Abx, `Gram pos bacilli`:`Gram neg cocci`) %>% distinct()
# 
# length(unique(temp$Group))
# length(unique(temp$Abx))
# length(unique(temp$`Species identification`))
# 
# fwrite(temp, "Lookup_species.csv")

# -----------------------

# Add Resistant/Susceptible Flags using MIC --------------------------------

Lookup_species <- fread("Lookup_species.csv", colClasses = "character")
Lookup_species <- Lookup_species  %>% mutate(`Species identification` =str_replace_all(`Species identification`, "�", " "))  
Lookup_species$`Species identification` <- str_trim(Lookup_species$`Species identification`)
Lookup_species <- Lookup_species %>% distinct()
length(unique(Lookup_species$Abx)) #32
length(unique(Lookup_species$`Species identification`)) #315


Lookup_species$EUCAST_mic_big <- as.numeric(Lookup_species$EUCAST_mic_big)
Lookup_species$CLSI_mic_bigeq <- as.numeric(Lookup_species$CLSI_mic_bigeq)
Lookup_species$CASFM_mic_big <- as.numeric(Lookup_species$CASFM_mic_big)
Lookup_species$EUCAST_disc_lower <- as.numeric(Lookup_species$EUCAST_disc_lower)
Lookup_species$CASFM_disc_lower <- as.numeric(Lookup_species$CASFM_disc_lower)


Lookup_species <- Lookup_species %>% 
  left_join(Lookup_species %>% select(`Species identification`) %>% distinct() %>%
              arrange(`Species identification`) %>% filter(`Species identification`!="") %>%
              drop_na() %>%
              mutate(speciesID=row_number())) 


data.frame(Lookup_species %>% select(`Species identification`, speciesID) %>% distinct()) %>%
  arrange(speciesID)

my_data <- read_excel(path = "ANAEuROBE_dataset_matteo_only.xlsx",  sheet="3.Completo", col_types = "text")
my_data <- my_data  %>% mutate(`Species identification` =str_replace_all(`Species identification`, "�", " "))  
my_data$`Species identification` <- str_trim(my_data$`Species identification`)

length(unique(my_data$`Species identification`)) #315

my_data <- data.frame(my_data %>% select(`Species identification`) %>% distinct() %>%
                        drop_na() %>%
                        arrange(`Species identification`) %>% mutate(speciesID=row_number()) %>%
                        mutate(speciesID=ifelse(speciesID>=44,speciesID+1, speciesID))) %>%
  inner_join(my_data, by=c("Species.identification"="Species identification"))


my_data <- my_data %>% select(-`Ceftriaxone inhib zone diameter...63`)

my_data <- my_data %>% rename("Ceftriaxone inhib zone diameter"="Ceftriaxone inhib zone diameter...69")

names(my_data)

MIC_data <- my_data %>% select(`Code event`, Species.identification, speciesID) %>%
  bind_cols(
    my_data %>% select(contains(" MIC"))
  )

MIC_data <- MIC_data %>% select(-c(`AST by ANA-ATB microdilution=1`,`AST by MICRONAUT-S MIC test=1`)) %>%
  gather(Abx, MIC, `Benzilpenicillin MIC`:`Oflox MIC`) 

unique(MIC_data$MIC)


MIC_data <- MIC_data %>% mutate(R_S=ifelse(grepl("R", MIC), "R",
                                           ifelse(grepl("S", MIC), "S",
                                                  ifelse(grepl("r",MIC), "R",
                                                         ifelse(grepl("s",MIC), "S", 
                                                                ifelse(grepl(">", MIC), "R",
                                                                       ifelse(grepl("<", MIC), "S", NA)))))))


MIC_data <- MIC_data %>% 
  mutate(MIC=ifelse(MIC=="<0.016", "0.008",
                    ifelse(MIC=="<0.019", "0.010",
                           ifelse(MIC=="<=0.016", "0.008",
                                  ifelse(MIC=="1.6E-2", "0.016",
                                         ifelse(MIC=="4.7E-2", "0.047",
                                                ifelse(MIC=="2.3E-2", "0.023",
                                                       ifelse(MIC=="6.4000000000000001E-2", "0.064",
                                                              ifelse(MIC=="9.4E-2", "0.094",
                                                                     ifelse(MIC==">256", "260",
                                                                            ifelse(MIC=="3.2000000000000001E-2", "0.032",
                                                                                   ifelse(MIC=="1.2E-2", "0.012",
                                                                                          ifelse(MIC=="2E-3", "0.002",
                                                                                                 ifelse(MIC=="8.0000000000000002E-3", "0.008",
                                                                                                        ifelse(MIC=="6.0000000000000001E-3", "0.006",
                                                                                                               ifelse(MIC==">32", "48",
                                                                                                                      ifelse(MIC=="3.0000000000000001E-3", "0.003",
                                                                                                                             ifelse(MIC=="4.0000000000000001E-3", "0.004",
                                                                                                                                    ifelse(MIC=="2.4E-2", "0.024",MIC)))))))))))))))))))



MIC_data <- MIC_data %>% mutate(MIC=ifelse(MIC=="≤0.25", "0.12", 
                                           ifelse(MIC==">0.5", "1",
                                                  ifelse(MIC=="≤0.016", "0.008",
                                                         ifelse(MIC=="<0.06", "0.03",
                                                                ifelse(MIC=="3.4000000000000002E-2", "0.03",
                                                                       ifelse(MIC=="<0.002", "0.001",
                                                                              ifelse(MIC=="<0.064", "0.03",
                                                                                     ifelse(MIC=="<0.03", "0.01",
                                                                                            ifelse(MIC=="<0.25", "0.1",
                                                                                                   ifelse(MIC==">2", "3",
                                                                                                          ifelse(MIC=="<256", "128",
                                                                                                                 ifelse(MIC==">8", "12",
                                                                                                                        ifelse(MIC=="00.032", "0.03",
                                                                                                                               ifelse(MIC=="≤4", "2",
                                                                                                                                      ifelse(MIC=="≤0.5", "0.2",
                                                                                                                                             ifelse(MIC=="1-2", "1",
                                                                                                                                                    ifelse(MIC=="TRUE", NA,
                                                                                                                                                           ifelse(MIC=="FALSE", NA,
                                                                                                                                                                  ifelse(MIC=="≤4/2", "2",
                                                                                                                                                                         ifelse(MIC==">8/2", "12",MIC)))))))))))))))))))))

MIC_data <- MIC_data %>% mutate(MIC=ifelse(MIC=="<0,06", "0.03",
                                           ifelse(MIC=="<1", "0.5",
                                                  ifelse(MIC=="<2", "1",
                                                         ifelse(MIC==">0.016", "0.03",
                                                                ifelse(MIC=="<4", "2",
                                                                       ifelse(MIC=="8/2", "8",
                                                                              ifelse(MIC=="1.7000000000000001E-2", "0.017",
                                                                                     ifelse(MIC=="§", "NA",
                                                                                            ifelse(MIC=="<0.16", "0.08",
                                                                                                   ifelse(MIC=="3.7999999999999999E-2", "0.04",
                                                                                                          ifelse(MIC=="≤8/4", "4",
                                                                                                                 ifelse(MIC=="≤1", "0.5",
                                                                                                                        ifelse(MIC=="Not tested", NA,
                                                                                                                               ifelse(MIC==">320", "480",
                                                                                                                                      ifelse(MIC=="<8", "4",MIC))))))))))))))))


MIC_data <- MIC_data %>% mutate(MIC=ifelse(MIC=="45 mm", "45",
                                           ifelse(MIC=="16/4", "16",
                                                  ifelse(MIC==">16/4", "20",
                                                         ifelse(MIC==">128", "200",
                                                                ifelse(MIC==">4", "6",
                                                                       ifelse(MIC=="≤0.06", "0.03",
                                                                              ifelse(MIC=="0,75", "0.7",
                                                                                     ifelse(MIC=="<0.026", "0.01",
                                                                                            ifelse(MIC=="≤0.0625", "0.03",
                                                                                                   ifelse(MIC==".064", "0.064",
                                                                                                          ifelse(MIC==">64", "100",
                                                                                                                 ifelse(MIC=="<0.047", "0.02",
                                                                                                                        ifelse(MIC=="1.2500000000000001E-2", "0.0125",
                                                                                                                               ifelse(MIC=="<0,016", "0.005",
                                                                                                                                      ifelse(MIC=="r", "R",MIC))))))))))))))))


MIC_data <- MIC_data %>% mutate(MIC=ifelse(MIC=="4.8000000000000001E-2", "0.048",
                                           ifelse(MIC=="< 0.016", "0.005",
                                                  ifelse(MIC=="<0.0016", "0.0005",
                                                         ifelse(MIC=="1.6000000000000001E-3", "0.0016",
                                                                ifelse(MIC=="<0.5", "0.2",
                                                                       ifelse(MIC=="<0,002", "0.001",
                                                                              ifelse(MIC=="9.1999999999999998E-2", "0.08",
                                                                                     ifelse(MIC=="> 32", "48",
                                                                                            ifelse(MIC=="2.3E-3", "0.002",
                                                                                                   ifelse(MIC=="<0.02", "0.01",
                                                                                                          ifelse(MIC=="1.9E-2", "0.019",
                                                                                                                 ifelse(MIC=="≤2", "1",MIC)))))))))))))



MIC_data <- MIC_data %>% mutate(MIC=ifelse(MIC=="R", NA,
                                           ifelse(MIC=="S", NA, 
                                                  ifelse(MIC=="r", NA,
                                                         ifelse(MIC=="s",NA,
                                                                ifelse(MIC=="O",NA,MIC))))))


unique(MIC_data$MIC)
MIC_data$MIC <- parse_number(MIC_data$MIC)



MIC_data <- MIC_data %>% left_join(Lookup_species %>% select(-c(`Species identification`, EUCAST_disc_lower, CASFM_disc_lower, ZI)))

unique(MIC_data$Species.identification)
unique(MIC_data$speciesID)


MIC_data <- MIC_data %>% 
  mutate(EUCAST_mic_Resist=ifelse(MIC>EUCAST_mic_big,1,0)) %>%
  mutate(CLIST_mic_Resist=ifelse(MIC>=CLSI_mic_bigeq ,1,0)) %>%
  mutate(CASFM_mic_Resist=ifelse(MIC>CASFM_mic_big,1,0)) 


ignore <- MIC_data %>% select(Species.identification, speciesID, Abx, EUCAST_mic_Resist) %>%
  group_by(Species.identification, speciesID, Abx, EUCAST_mic_Resist) %>% count() %>%
  drop_na() %>% spread(key=EUCAST_mic_Resist, value=n) 

ignore[is.na(ignore)] <- 0

ignore <- ignore %>% mutate(perc=`1`/(`1`+`0`))

ignore <- ignore %>% filter(`0`+`1`>=20) %>% select(-c(`0`, `1`)) %>%
  spread(key=Abx, valu=perc)

# fwrite(ignore, "ignore.csv")

MIC_data %>% group_by(R_S, EUCAST_mic_Resist) %>% count()


summary_MIC_concent <- MIC_data %>% filter(!is.na(MIC)) %>%
  group_by(Species.identification, speciesID, Abx) %>% 
  summarise(mean=mean(MIC, na.rm=T), 
            sd=sd(MIC), 
            mean=mean(MIC), 
            median=median(MIC), 
            Q1=quantile(MIC, 0.25),
            Q3=quantile(MIC, 0.75),
            n=n()) 

R_S_original_counts <- MIC_data %>% filter(!is.na(R_S)) %>%
  group_by(Species.identification, speciesID, Abx, R_S) %>% count() %>%
  spread(key=R_S, value=n) %>%
  mutate(R=ifelse(is.na(R),0,R)) %>%
  mutate(S=ifelse(is.na(S),0,S)) %>%
  mutate(n_n=R+S, perc_r_s=R/(R+S)) %>% select(-c(R,S))


summary_MIC_concent %>% full_join(R_S_original_counts)


EUCAST_resist_counts <- MIC_data %>% filter(!is.na(EUCAST_mic_Resist)) %>%
  group_by(Species.identification, speciesID, Abx, EUCAST_mic_Resist) %>% count() %>%
  spread(key=EUCAST_mic_Resist, value=n) %>%
  mutate(`1`=ifelse(is.na(`1`),0,`1`)) %>%
  mutate(`0`=ifelse(is.na(`0`),0,`0`)) %>%
  mutate(n_eucast=`1`+`0`, perc_r_eucast=`1`/(`1`+`0`)) %>% select(-c(`1`,`0`))

CLSI_resist_counts <- MIC_data %>% filter(!is.na(CLIST_mic_Resist )) %>%
  group_by(Species.identification, speciesID, Abx, CLIST_mic_Resist ) %>% count() %>%
  spread(key=CLIST_mic_Resist , value=n) %>%
  mutate(`1`=ifelse(is.na(`1`),0,`1`)) %>%
  mutate(`0`=ifelse(is.na(`0`),0,`0`)) %>%
  mutate(n_clsi=`1`+`0`, perc_r_clsi=`1`/(`1`+`0`)) %>% select(-c(`1`,`0`))

casfm_resist_counts <- MIC_data %>% filter(!is.na(CASFM_mic_Resist )) %>%
  group_by(Species.identification, speciesID, Abx, CASFM_mic_Resist ) %>% count() %>%
  spread(key=CASFM_mic_Resist , value=n) %>%
  mutate(`1`=ifelse(is.na(`1`),0,`1`)) %>%
  mutate(`0`=ifelse(is.na(`0`),0,`0`)) %>%
  mutate(n_casfm=`1`+`0`, perc_r_casfm=`1`/(`1`+`0`)) %>% select(-c(`1`,`0`))


fwrite(summary_MIC_concent, "summary_MIC_concent.csv")
fwrite(R_S_original_counts, "R_S_original_counts.csv")

fwrite(MIC_data, "MIC_data.csv")

fwrite(EUCAST_resist_counts, "EUCAST_resist_counts.csv")
fwrite(CLSI_resist_counts, "CLSI_resist_counts.csv")
fwrite(casfm_resist_counts, "casfm_resist_counts.csv")


# -----------------------------------

# Add Resistant/Susceptible Flags using Inhib Zone Diam --------------------------------

Lookup_species <- fread("Lookup_species.csv", colClasses = "character")
Lookup_species <- Lookup_species  %>% mutate(`Species identification` =str_replace_all(`Species identification`, "�", " "))  
Lookup_species$`Species identification` <- str_trim(Lookup_species$`Species identification`)
Lookup_species <- Lookup_species %>% distinct()
length(unique(Lookup_species$Abx)) #32
length(unique(Lookup_species$`Species identification`)) #315

Lookup_species <- Lookup_species %>% 
  left_join(Lookup_species %>% select(`Species identification`) %>% distinct() %>%
              arrange(`Species identification`) %>% filter(`Species identification`!="") %>%
              drop_na() %>%
              mutate(speciesID=row_number())) 


Lookup_species$EUCAST_mic_big <- as.numeric(Lookup_species$EUCAST_mic_big)
Lookup_species$CLSI_mic_bigeq <- as.numeric(Lookup_species$CLSI_mic_bigeq)
Lookup_species$CASFM_mic_big <- as.numeric(Lookup_species$CASFM_mic_big)
Lookup_species$EUCAST_disc_lower <- as.numeric(Lookup_species$EUCAST_disc_lower)
Lookup_species$CASFM_disc_lower <- as.numeric(Lookup_species$CASFM_disc_lower)


my_data <- read_excel(path = "ANAEuROBE_dataset_matteo_only.xlsx",  sheet="3.Completo", col_types = "text")
my_data <- my_data  %>% mutate(`Species identification` =str_replace_all(`Species identification`, "�", " "))  
my_data$`Species identification` <- str_trim(my_data$`Species identification`)

length(unique(my_data$`Species identification`)) #315

my_data <- data.frame(my_data %>% select(`Species identification`) %>% distinct() %>%
                        drop_na() %>%
                        arrange(`Species identification`) %>% mutate(speciesID=row_number()) %>%
                        mutate(speciesID=ifelse(speciesID>=44,speciesID+1, speciesID))) %>%
  inner_join(my_data, by=c("Species.identification"="Species identification"))


my_data <- my_data %>% select(-`Ceftriaxone inhib zone diameter...63`)

my_data <- my_data %>% rename("Ceftriaxone inhib zone diameter"="Ceftriaxone inhib zone diameter...69")

names(my_data)


ZONE_data <- my_data %>% select(`Code event`, Species.identification, speciesID) %>%
  bind_cols(
    my_data %>% select(contains(" zone diam"))
  )

ZONE_data <- ZONE_data %>% bind_cols(my_data %>% select(`EUCAST=1; CLSI=0`))


ZONE_data <- ZONE_data %>%
  gather(Abx, ZONE, `Benzilpenicillin inhib zone diameter`:`Oflox zone diam`) 

ZONE_data %>% group_by(ZONE) %>% count() %>% arrange(-n)

ZONE_data <- ZONE_data %>% mutate(ZONE=ifelse(ZONE=="S", NA,
                                              ifelse(ZONE=="R", NA,
                                                     ifelse(ZONE==">35", "50",
                                                            ifelse(ZONE=="na", NA,
                                                                   ifelse(ZONE==">16", "24",
                                                                          ifelse(ZONE==">32","48",
                                                                                 ifelse(ZONE=="8.0000000000000002E-3", "0.008",
                                                                                        ifelse(ZONE=="<0.002", "0.001",
                                                                                               ifelse(ZONE=="4.0000000000000001E-3", "0.004",
                                                                                                      ifelse(ZONE=="2.3E-2", "0.0023",
                                                                                                             ifelse(ZONE=="6.0000000000000001E-3", "0.006",
                                                                                                                    ifelse(ZONE=="9.4E-2", "0.0094", ZONE)))))))))))))



ZONE_data <- ZONE_data %>% mutate(ZONE=ifelse(ZONE==">30", 45,
                                              ifelse(ZONE=="s", NA,
                                                     ifelse(ZONE==">40", "60",
                                                            ifelse(ZONE=="I", NA,
                                                                   ifelse(ZONE=="Not tested", NA,
                                                                          ifelse(ZONE==">25","38",
                                                                                 ifelse(ZONE==">22", "33",
                                                                                        ifelse(ZONE==">20", "30",
                                                                                               ifelse(ZONE=="<40", "20",
                                                                                                      ifelse(ZONE=="6.4000000000000001E-2", "0.064",
                                                                                                             ifelse(ZONE=="1.2E-2", "0.012",
                                                                                                                    ifelse(ZONE=="3.0000000000000001E-3", "0.003", 
                                                                                                                           ifelse(ZONE=="1.6E-2", "0.0016",
                                                                                                                                  ifelse(ZONE=="4.7E-2", "0.047",
                                                                                                                                         ifelse(ZONE=="3.2000000000000001E-2", "0.032",
                                                                                                                                                ifelse(ZONE=="2E-3","0.002",
                                                                                                                                                       ifelse(ZONE==">27","40", ZONE))))))))))))))))))





ZONE_data$ZONE <- as.numeric(ZONE_data$ZONE)

ZONE_data <- ZONE_data %>% drop_na() %>%
  mutate(Abx=str_replace_all(Abx, " inhib zone diameter", "")) %>%
  mutate(Abx=str_replace_all(Abx, " zone diam", "")) %>%
  left_join(Lookup_species %>% select(-`Species identification`) %>% mutate(Abx=str_replace_all(Abx, " MIC", "")) %>% 
              select(-c(EUCAST_mic_big, CLSI_mic_bigeq , CASFM_mic_big)))


data.frame(ZONE_data %>% filter(grepl("TOU",`Code event`)|
                                  grepl("PAR",`Code event`)|
                                  grepl("LIM",`Code event`)) %>% select(Abx) %>% distinct())


# Abx
# 1         Benzilpenicillin
# 2               Ampicillin
# 3              Amoxicillin
# 4  Amoxicillin/clavulanate
# 5  Piperacillin/tazobactam * 
# 6                Cefoxitin
# 7               Cefotaxime
# 8              Clindamycin *
# 9              Ceftriaxone
# 10           Metronidazole *
# 11         Chloramphenicol
# 12               Ertapenem *
# 13                Imipenem *
# 14              Vancomycin
# 15               Linezolid
# 16             Tigecycline
# 17              Rifampicin
# 18            Moxifloxacin
# 19            Tetracycline
# 20                   Genta
# 21                    Levo
# 22                   Oflox


unique(ZONE_data$`EUCAST=1; CLSI=0`)


ZONE_data <- ZONE_data %>% 
  mutate(EUCAST_Diam_Resist=ifelse(ZONE <EUCAST_disc_lower,1,0)) %>%
  mutate(EUCAST_Diam_Resist=ifelse(
    Abx %in% c("Piperacillin/tazobactam", "Clindamycin", "Metronidazole", "Ertapenem", "Imipenem"), EUCAST_Diam_Resist, 
    ifelse(`EUCAST=1; CLSI=0`=="1", EUCAST_Diam_Resist, NA))) %>%
  mutate(CASFM_Diam_Resist=ifelse(ZONE <CASFM_disc_lower,1,0)) %>%
  mutate(CASFM_Diam_Resist=ifelse(
    Abx %in% c("Piperacillin/tazobactam", "Clindamycin", "Metronidazole", "Ertapenem", "Imipenem"), CASFM_Diam_Resist, 
    ifelse(`EUCAST=1; CLSI=0`!="1", CASFM_Diam_Resist, NA)))

ignore <- ZONE_data %>% select(Species.identification, speciesID, Abx, EUCAST_Diam_Resist) %>%
  group_by(Species.identification, speciesID, Abx, EUCAST_Diam_Resist) %>% count() %>%
  drop_na() %>% spread(key=EUCAST_Diam_Resist, value=n) 

ignore[is.na(ignore)] <- 0

ignore <- ignore %>% filter(`1`+`0`>=10) %>% mutate(perc=`1`/(`1`+`0`))

ignore <- ignore %>% filter(`0`+`1`>=20) %>% select(-c(`0`, `1`)) %>%
  spread(key=Abx, valu=perc)

# fwrite(ignore, "ignore_zone.csv")


ZONE_data

summary_ZONE_concent <- ZONE_data %>% filter(!is.na(ZONE)) %>%
  group_by(Species.identification, speciesID, Abx) %>% 
  summarise(mean=mean(ZONE, na.rm=T), 
            sd=sd(ZONE), 
            mean=mean(ZONE), 
            median=median(ZONE), 
            Q1=quantile(ZONE, 0.25),
            Q3=quantile(ZONE, 0.75),
            n=n()) 

EUCAST_resist_counts_zone <- ZONE_data %>% filter(!is.na(EUCAST_Diam_Resist)) %>%
  group_by(Species.identification, speciesID, Abx, EUCAST_Diam_Resist) %>% count() %>%
  spread(key=EUCAST_Diam_Resist, value=n) %>%
  mutate(`1`=ifelse(is.na(`1`),0,`1`)) %>%
  mutate(`0`=ifelse(is.na(`0`),0,`0`)) %>%
  mutate(n_eucast=`1`+`0`, perc_r_eucast=`1`/(`1`+`0`)) %>% select(-c(`1`,`0`))

casfm_resist_counts_zone <- ZONE_data %>% filter(!is.na(CASFM_Diam_Resist  )) %>%
  group_by(Species.identification, speciesID, Abx, CASFM_Diam_Resist ) %>% count() %>%
  spread(key=CASFM_Diam_Resist , value=n) %>%
  mutate(`1`=ifelse(is.na(`1`),0,`1`)) %>%
  mutate(`0`=ifelse(is.na(`0`),0,`0`)) %>%
  mutate(n_casfm=`1`+`0`, perc_r_casfm=`1`/(`1`+`0`)) %>% select(-c(`1`,`0`))


fwrite(summary_ZONE_concent, "summary_ZONE_concent.csv")

fwrite(ZONE_data, "ZONE_data.csv")

fwrite(EUCAST_resist_counts_zone, "EUCAST_resist_counts_zone.csv")
fwrite(casfm_resist_counts_zone, "casfm_resist_counts_zone.csv")


ZI <- ZONE_data %>% filter(grepl(",", ZI)) %>%
  filter(!is.na(EUCAST_Diam_Resist)|!is.na(CASFM_Diam_Resist)) %>%
  separate(ZI, c("Lower", "Upper")) %>%
  select(`Code event`, Species.identification, speciesID, Abx, ZONE, Lower, Upper)


ZI$Lower <- as.numeric(ZI$Lower)
ZI$Upper <- as.numeric(ZI$Upper)


ZI <- ZI %>% mutate(ZI=ifelse(ZONE>Lower&ZONE<Upper,1,0))

ZI %>% group_by(ZI) %>% count()

ZI %>% group_by(Species.identification, speciesID, Abx) %>%
  summarise(mean=mean(ZI),
            n=n()) %>%
  filter(n>=20)

fwrite(ZI, "ZI.csv")




# ---------------------


# Overall figures ---------------------------

my_data <- read_excel(path = "ANAEuROBE_dataset_matteo_only.xlsx",  sheet="3.Completo", col_types = "text")
my_data <- my_data  %>% mutate(`Species identification` =str_replace_all(`Species identification`, "�", " "))  
my_data$`Species identification` <- str_trim(my_data$`Species identification`)

length(unique(my_data$`Species identification`)) #315

my_data <- data.frame(my_data %>% select(`Species identification`) %>% distinct() %>%
                        drop_na() %>%
                        arrange(`Species identification`) %>% mutate(speciesID=row_number()) %>%
                        mutate(speciesID=ifelse(speciesID>=44,speciesID+1, speciesID))) %>%
  inner_join(my_data, by=c("Species.identification"="Species identification"))



names(my_data)
dim(my_data) # 14528   110

my_data <- my_data %>% 
  mutate(`Gram neg bacilli`=ifelse(Species.identification=="Prevotella spp",1,`Gram neg bacilli`)) %>%
  mutate(`Gram neg cocci`=ifelse(Species.identification=="Veillonella spp",1,`Gram neg cocci`))

# Isolates by gram stain

my_data %>% select(`Code event`, Species.identification, speciesID, `Gram pos bacilli`:`Gram neg cocci`) %>%
  filter(!is.na(Species.identification)) %>%
  gather(Group, exp, `Gram pos bacilli`:`Gram neg cocci`) %>%
  filter(exp==1) %>%
  group_by(Group) %>% count()


# isolates per year 

my_data %>% select(`Code event`, Species.identification, speciesID, `2020=1`:`2023=1`) %>%
  filter(!is.na(Species.identification)) %>%
  gather(Group, exp, `2020=1`:`2023=1`) %>%
  filter(exp==1) %>% group_by(`Code event`) %>% filter(Group==max(Group)) %>% ungroup() %>%
  group_by(Group) %>% count()


# Isolates per country

data.frame(my_data %>% select(`Code event`, Species.identification, speciesID, `PT`:`PL`) %>%
             filter(!is.na(Species.identification)) %>%
             gather(Group, exp, `PT`:`PL`) %>%
             filter(exp==1) %>%
             group_by(Group) %>% count() %>% mutate(n=ifelse(Group=="IT", n-1,n)) %>%
             arrange(n))

# Medical ward

data.frame(my_data %>% select(`Code event`, `Species.identification, speciesID,  `Emergency=1`:`ICU=1`) %>%
             filter(!is.na(Species.identification)) %>%
             gather(Group, exp, `Emergency=1`:`ICU=1`) %>%
             filter(exp==1) %>% distinct() %>% 
             group_by(Group) %>%
             count() %>%
             arrange(n))


# Most Common species

TOP50 <- data.frame(my_data %>% select(`Code event`, Species.identification, speciesID, ) %>%
                      filter(!is.na(Species.identification)) %>%
                      group_by(Species.identification) %>% count() %>%
                      arrange(-n) %>% ungroup() %>% mutate(cum=cumsum(n))) %>%
  slice(1:50)

fwrite(TOP50, "TOP50.csv")

# Most Common species by country TOP10 per country

temp <- data.frame(my_data %>% select(`Code event`, Species.identification, speciesID, `PT`:`PL`) %>%
                     filter(!is.na(Species.identification)) %>%
                     mutate(`IT`=ifelse(`GR`==1,NA,`IT`)) %>%
                     gather(Group, exp, `PT`:`PL`) %>%
                     filter(exp==1) %>%
                     group_by(Group, Species.identification, speciesID,) %>% count() %>%
                     rename("sub"="n") %>% ungroup() %>%
                     group_by(Group) %>% mutate(Tot=sum(sub)) %>%
                     mutate(perc=sub/Tot) %>%
                     arrange(Group, -perc ) %>%
                     group_by(Group) %>% slice(1:10) %>% select(-c(sub, Tot)) %>%
                     spread(key=Group, value=perc))


fwrite(temp, "temp.csv")

data.frame(my_data %>% select(`Code event`, Species.identification, speciesID, `Emergency=1`:`ICU=1`) %>%
             filter(!is.na(Species.identification)) %>%
             gather(Group, exp, `Emergency=1`:`ICU=1`) %>%
             filter(exp==1) %>% distinct() %>% 
             group_by(Group,Species.identification, speciesID,) %>%
             count() %>%
             rename("sub"="n") %>% ungroup() %>%
             group_by(Group) %>% mutate(Tot=sum(sub)) %>%
             mutate(perc=sub/Tot) %>%
             arrange(Group, -perc ) %>%
             group_by(Group) %>% slice(1:10) %>% select(-c(sub, Tot)) %>%
             spread(key=Group, value=perc))



# --------------------

# PLot species vs resistance rate MIC ---------------------

Lookup_species <- fread("Lookup_species.csv", colClasses = "character")
Lookup_species <- Lookup_species  %>% mutate(`Species identification` =str_replace_all(`Species identification`, "�", " "))  
Lookup_species$`Species identification` <- str_trim(Lookup_species$`Species identification`)
Lookup_species <- Lookup_species %>% distinct()
length(unique(Lookup_species$Abx)) #32
length(unique(Lookup_species$`Species identification`)) #315


Lookup_species$EUCAST_mic_big <- as.numeric(Lookup_species$EUCAST_mic_big)
Lookup_species$CLSI_mic_bigeq <- as.numeric(Lookup_species$CLSI_mic_bigeq)
Lookup_species$CASFM_mic_big <- as.numeric(Lookup_species$CASFM_mic_big)
Lookup_species$EUCAST_disc_lower <- as.numeric(Lookup_species$EUCAST_disc_lower)
Lookup_species$CASFM_disc_lower <- as.numeric(Lookup_species$CASFM_disc_lower)


Lookup_species <- Lookup_species %>% 
  left_join(Lookup_species %>% select(`Species identification`) %>% distinct() %>%
              arrange(`Species identification`) %>% filter(`Species identification`!="") %>%
              drop_na() %>%
              mutate(speciesID=row_number())) 


data.frame(Lookup_species %>% select(`Species identification`, speciesID) %>% distinct()) %>%
  arrange(speciesID)

my_data <- read_excel(path = "ANAEuROBE_dataset_matteo_only.xlsx",  sheet="3.Completo", col_types = "text")
my_data <- my_data  %>% mutate(`Species identification` =str_replace_all(`Species identification`, "�", " "))  
my_data$`Species identification` <- str_trim(my_data$`Species identification`)

length(unique(my_data$`Species identification`)) #315

my_data <- data.frame(my_data %>% select(`Species identification`) %>% distinct() %>%
                        drop_na() %>%
                        arrange(`Species identification`) %>% mutate(speciesID=row_number()) %>%
                        mutate(speciesID=ifelse(speciesID>=44,speciesID+1, speciesID))) %>%
  inner_join(my_data, by=c("Species.identification"="Species identification"))


my_data <- my_data %>% select(-`Ceftriaxone inhib zone diameter...63`)

my_data <- my_data %>% rename("Ceftriaxone inhib zone diameter"="Ceftriaxone inhib zone diameter...69")

names(my_data)

MIC_data <- my_data %>% select(`Code event`, Species.identification, speciesID) %>%
  bind_cols(
    my_data %>% select(contains(" MIC"))
  )

MIC_data <- MIC_data %>% select(-c(`AST by ANA-ATB microdilution=1`,`AST by MICRONAUT-S MIC test=1`)) %>%
  gather(Abx, MIC, `Benzilpenicillin MIC`:`Oflox MIC`) 

unique(MIC_data$MIC)

MIC_data <- MIC_data %>% mutate(R_S=ifelse(grepl("R", MIC), "R",
                                           ifelse(grepl("S", MIC), "S",
                                                  ifelse(grepl("r",MIC), "R",
                                                         ifelse(grepl("s",MIC), "S", 
                                                                ifelse(grepl(">", MIC), "R",
                                                                       ifelse(grepl("<", MIC), "S", NA)))))))


MIC_data <- MIC_data %>% 
  mutate(MIC=ifelse(MIC=="<0.016", "0.008",
                    ifelse(MIC=="<0.019", "0.010",
                           ifelse(MIC=="<=0.016", "0.008",
                                  ifelse(MIC=="1.6E-2", "0.016",
                                         ifelse(MIC=="4.7E-2", "0.047",
                                                ifelse(MIC=="2.3E-2", "0.023",
                                                       ifelse(MIC=="6.4000000000000001E-2", "0.064",
                                                              ifelse(MIC=="9.4E-2", "0.094",
                                                                     ifelse(MIC==">256", "260",
                                                                            ifelse(MIC=="3.2000000000000001E-2", "0.032",
                                                                                   ifelse(MIC=="1.2E-2", "0.012",
                                                                                          ifelse(MIC=="2E-3", "0.002",
                                                                                                 ifelse(MIC=="8.0000000000000002E-3", "0.008",
                                                                                                        ifelse(MIC=="6.0000000000000001E-3", "0.006",
                                                                                                               ifelse(MIC==">32", "48",
                                                                                                                      ifelse(MIC=="3.0000000000000001E-3", "0.003",
                                                                                                                             ifelse(MIC=="4.0000000000000001E-3", "0.004",
                                                                                                                                    ifelse(MIC=="2.4E-2", "0.024",MIC)))))))))))))))))))



MIC_data <- MIC_data %>% mutate(MIC=ifelse(MIC=="≤0.25", "0.12", 
                                           ifelse(MIC==">0.5", "1",
                                                  ifelse(MIC=="≤0.016", "0.008",
                                                         ifelse(MIC=="<0.06", "0.03",
                                                                ifelse(MIC=="3.4000000000000002E-2", "0.03",
                                                                       ifelse(MIC=="<0.002", "0.001",
                                                                              ifelse(MIC=="<0.064", "0.03",
                                                                                     ifelse(MIC=="<0.03", "0.01",
                                                                                            ifelse(MIC=="<0.25", "0.1",
                                                                                                   ifelse(MIC==">2", "3",
                                                                                                          ifelse(MIC=="<256", "128",
                                                                                                                 ifelse(MIC==">8", "12",
                                                                                                                        ifelse(MIC=="00.032", "0.03",
                                                                                                                               ifelse(MIC=="≤4", "2",
                                                                                                                                      ifelse(MIC=="≤0.5", "0.2",
                                                                                                                                             ifelse(MIC=="1-2", "1",
                                                                                                                                                    ifelse(MIC=="TRUE", NA,
                                                                                                                                                           ifelse(MIC=="FALSE", NA,
                                                                                                                                                                  ifelse(MIC=="≤4/2", "2",
                                                                                                                                                                         ifelse(MIC==">8/2", "12",MIC)))))))))))))))))))))

MIC_data <- MIC_data %>% mutate(MIC=ifelse(MIC=="<0,06", "0.03",
                                           ifelse(MIC=="<1", "0.5",
                                                  ifelse(MIC=="<2", "1",
                                                         ifelse(MIC==">0.016", "0.03",
                                                                ifelse(MIC=="<4", "2",
                                                                       ifelse(MIC=="8/2", "8",
                                                                              ifelse(MIC=="1.7000000000000001E-2", "0.017",
                                                                                     ifelse(MIC=="§", "NA",
                                                                                            ifelse(MIC=="<0.16", "0.08",
                                                                                                   ifelse(MIC=="3.7999999999999999E-2", "0.04",
                                                                                                          ifelse(MIC=="≤8/4", "4",
                                                                                                                 ifelse(MIC=="≤1", "0.5",
                                                                                                                        ifelse(MIC=="Not tested", NA,
                                                                                                                               ifelse(MIC==">320", "480",
                                                                                                                                      ifelse(MIC=="<8", "4",MIC))))))))))))))))


MIC_data <- MIC_data %>% mutate(MIC=ifelse(MIC=="45 mm", "45",
                                           ifelse(MIC=="16/4", "16",
                                                  ifelse(MIC==">16/4", "20",
                                                         ifelse(MIC==">128", "200",
                                                                ifelse(MIC==">4", "6",
                                                                       ifelse(MIC=="≤0.06", "0.03",
                                                                              ifelse(MIC=="0,75", "0.7",
                                                                                     ifelse(MIC=="<0.026", "0.01",
                                                                                            ifelse(MIC=="≤0.0625", "0.03",
                                                                                                   ifelse(MIC==".064", "0.064",
                                                                                                          ifelse(MIC==">64", "100",
                                                                                                                 ifelse(MIC=="<0.047", "0.02",
                                                                                                                        ifelse(MIC=="1.2500000000000001E-2", "0.0125",
                                                                                                                               ifelse(MIC=="<0,016", "0.005",
                                                                                                                                      ifelse(MIC=="r", "R",MIC))))))))))))))))


MIC_data <- MIC_data %>% mutate(MIC=ifelse(MIC=="4.8000000000000001E-2", "0.048",
                                           ifelse(MIC=="< 0.016", "0.005",
                                                  ifelse(MIC=="<0.0016", "0.0005",
                                                         ifelse(MIC=="1.6000000000000001E-3", "0.0016",
                                                                ifelse(MIC=="<0.5", "0.2",
                                                                       ifelse(MIC=="<0,002", "0.001",
                                                                              ifelse(MIC=="9.1999999999999998E-2", "0.08",
                                                                                     ifelse(MIC=="> 32", "48",
                                                                                            ifelse(MIC=="2.3E-3", "0.002",
                                                                                                   ifelse(MIC=="<0.02", "0.01",
                                                                                                          ifelse(MIC=="1.9E-2", "0.019",
                                                                                                                 ifelse(MIC=="≤2", "1",MIC)))))))))))))



MIC_data <- MIC_data %>% mutate(MIC=ifelse(MIC=="R", NA,
                                           ifelse(MIC=="S", NA, 
                                                  ifelse(MIC=="r", NA,
                                                         ifelse(MIC=="s",NA,
                                                                ifelse(MIC=="O",NA,MIC))))))


unique(MIC_data$MIC)
MIC_data$MIC <- parse_number(MIC_data$MIC)


MIC_data <- MIC_data %>% left_join(Lookup_species %>% select(-c(`Species identification`, EUCAST_disc_lower, CASFM_disc_lower, ZI)))

MIC_data <- MIC_data %>% 
  mutate(EUCAST_mic_Resist=ifelse(MIC>EUCAST_mic_big,1,0)) %>%
  mutate(CLIST_mic_Resist=ifelse(MIC>=CLSI_mic_bigeq ,1,0)) %>%
  mutate(CASFM_mic_Resist=ifelse(MIC>CASFM_mic_big,1,0)) 



MIC_data %>% select(Species.identification, speciesID) %>% distinct()
data.frame(Lookup_species %>% select(`Species identification`, speciesID) %>% distinct() %>%
             arrange(speciesID))

temp <- MIC_data %>% select(`Code event`, Species.identification, speciesID, Abx, EUCAST_mic_Resist)

temp <- temp %>% drop_na()


library(pheatmap)
library(dplyr)
library(tidyr)

resistance_summary <- temp %>%
  group_by(Species.identification, speciesID, Abx) %>%
  summarise(
    n_samples = n(),
    resistance_rate = mean(EUCAST_mic_Resist) * 100
  ) %>%
  filter(n_samples >= 20) 

heatmap_data <- resistance_summary %>% ungroup() %>% select(-speciesID) %>%
  select(Species.identification, Abx, resistance_rate) %>%
  pivot_wider(names_from = Abx, values_from = resistance_rate, values_fill = NA)

heatmap_matrix <- as.matrix(heatmap_data[, -1])

rownames(heatmap_matrix) <- heatmap_data$Species.identification

heatmap_matrix[is.na(heatmap_matrix)] <- -99

display_matrix <- round(heatmap_matrix, 0)


plot <- pheatmap(heatmap_matrix,
                 color = colorRampPalette(c("lightgray", "lightblue", "midnightblue"))(50),  
                 cluster_rows = TRUE,  # Cluster species
                 cluster_cols = TRUE,  # Cluster antibiotics
                 na_col = "grey",  # Color for missing values
                 fontsize_row = 10,  # Font size for species
                 fontsize_col = 10,  # Font size for antibiotics
                 display_numbers = display_matrix,  # Show exact resistance rates
                 main = "Antibiotic Resistance Clustering [n>20] \n"
)

ggsave(file="dendo.svg", plot=plot, width=6, height=12)






resistance_summary <- temp %>%
  select(`Code event`, Species.identification, Abx, EUCAST_mic_Resist) %>% drop_na() %>%
  group_by(Species.identification, Abx) %>%
  summarise(
    n_samples = n(),
    resistance_rate = mean(EUCAST_mic_Resist) * 100
  ) %>%
  filter(n_samples >= 20) 

resistance_summary <- resistance_summary %>% mutate(Abx=str_replace_all(Abx, " MIC", ""))

plot <- ggplot(resistance_summary, aes(x = Abx, y = Species.identification)) +
  geom_point(aes(size = n_samples , color = resistance_rate), alpha = 0.9) +
  scale_size(range = c(1, 20)) +  # Adjust the bubble size range
  scale_color_gradient(low = "lightgray", high = "midnightblue") +  # Color scale from susceptible (green) to resistant (red)
  labs(title = "Antibiotic Resistance by Species and Antibiotic \n [n>20]",
       x = "Antibiotic",
       y = "Species",
       size = "Sample Size",
       color = "Proportion Resistant") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  


ggsave(file="bubble.svg", plot=plot, width=6, height=12)




# ----------------
# Plot species vs resistance rate Diameter  --------------

Lookup_species <- fread("Lookup_species.csv", colClasses = "character")
Lookup_species <- Lookup_species  %>% mutate(`Species identification` =str_replace_all(`Species identification`, "�", " "))  
Lookup_species$`Species identification` <- str_trim(Lookup_species$`Species identification`)
Lookup_species <- Lookup_species %>% distinct()
length(unique(Lookup_species$Abx)) #32
length(unique(Lookup_species$`Species identification`)) #315


Lookup_species$EUCAST_mic_big <- as.numeric(Lookup_species$EUCAST_mic_big)
Lookup_species$CLSI_mic_bigeq <- as.numeric(Lookup_species$CLSI_mic_bigeq)
Lookup_species$CASFM_mic_big <- as.numeric(Lookup_species$CASFM_mic_big)
Lookup_species$EUCAST_disc_lower <- as.numeric(Lookup_species$EUCAST_disc_lower)
Lookup_species$CASFM_disc_lower <- as.numeric(Lookup_species$CASFM_disc_lower)


Lookup_species <- Lookup_species %>% 
  left_join(Lookup_species %>% select(`Species identification`) %>% distinct() %>%
              arrange(`Species identification`) %>% filter(`Species identification`!="") %>%
              drop_na() %>%
              mutate(speciesID=row_number())) 



my_data <- read_excel(path = "ANAEuROBE_dataset_matteo_only.xlsx",  sheet="3.Completo", col_types = "text")
my_data <- my_data  %>% mutate(`Species identification` =str_replace_all(`Species identification`, "�", " "))  
my_data$`Species identification` <- str_trim(my_data$`Species identification`)

length(unique(my_data$`Species identification`)) #315

my_data <- data.frame(my_data %>% select(`Species identification`) %>% distinct() %>%
                        drop_na() %>%
                        arrange(`Species identification`) %>% mutate(speciesID=row_number()) %>%
                        mutate(speciesID=ifelse(speciesID>=44,speciesID+1, speciesID))) %>%
  inner_join(my_data, by=c("Species.identification"="Species identification"))


my_data <- my_data %>% select(-`Ceftriaxone inhib zone diameter...63`)

my_data <- my_data %>% rename("Ceftriaxone inhib zone diameter"="Ceftriaxone inhib zone diameter...69")

names(my_data)

ZONE_data <- my_data %>% select(`Code event`, Species.identification, speciesID) %>%
  bind_cols(
    my_data %>% select(contains(" zone diam"))
  )

ZONE_data <- ZONE_data %>% bind_cols(my_data %>% select(`EUCAST=1; CLSI=0`))


ZONE_data <- ZONE_data %>%
  gather(Abx, ZONE, `Benzilpenicillin inhib zone diameter`:`Oflox zone diam`) 

unique(ZONE_data$ZONE)

unique(ZONE_data$Abx)

unique(ZONE_data$ZONE)

ZONE_data %>% group_by(ZONE) %>% count() %>% arrange(-n)

ZONE_data <- ZONE_data %>% mutate(ZONE=ifelse(ZONE=="S", NA,
                                              ifelse(ZONE=="R", NA,
                                                     ifelse(ZONE==">35", "50",
                                                            ifelse(ZONE=="na", NA,
                                                                   ifelse(ZONE==">16", "24",
                                                                          ifelse(ZONE==">32","48",
                                                                                 ifelse(ZONE=="8.0000000000000002E-3", "0.008",
                                                                                        ifelse(ZONE=="<0.002", "0.001",
                                                                                               ifelse(ZONE=="4.0000000000000001E-3", "0.004",
                                                                                                      ifelse(ZONE=="2.3E-2", "0.0023",
                                                                                                             ifelse(ZONE=="6.0000000000000001E-3", "0.006",
                                                                                                                    ifelse(ZONE=="9.4E-2", "0.0094", ZONE)))))))))))))



ZONE_data <- ZONE_data %>% mutate(ZONE=ifelse(ZONE==">30", 45,
                                              ifelse(ZONE=="s", NA,
                                                     ifelse(ZONE==">40", "60",
                                                            ifelse(ZONE=="I", NA,
                                                                   ifelse(ZONE=="Not tested", NA,
                                                                          ifelse(ZONE==">25","38",
                                                                                 ifelse(ZONE==">22", "33",
                                                                                        ifelse(ZONE==">20", "30",
                                                                                               ifelse(ZONE=="<40", "20",
                                                                                                      ifelse(ZONE=="6.4000000000000001E-2", "0.064",
                                                                                                             ifelse(ZONE=="1.2E-2", "0.012",
                                                                                                                    ifelse(ZONE=="3.0000000000000001E-3", "0.003", 
                                                                                                                           ifelse(ZONE=="1.6E-2", "0.0016",
                                                                                                                                  ifelse(ZONE=="4.7E-2", "0.047",
                                                                                                                                         ifelse(ZONE=="3.2000000000000001E-2", "0.032",
                                                                                                                                                ifelse(ZONE=="2E-3","0.002",
                                                                                                                                                       ifelse(ZONE==">27","40", ZONE))))))))))))))))))





ZONE_data$ZONE <- as.numeric(ZONE_data$ZONE)

ZONE_data <- ZONE_data %>% drop_na() %>%
  mutate(Abx=str_replace_all(Abx, " inhib zone diameter", "")) %>%
  mutate(Abx=str_replace_all(Abx, " zone diam", "")) %>%
  left_join(Lookup_species %>% select(-`Species identification`) %>% mutate(Abx=str_replace_all(Abx, " MIC", "")) %>% 
              select(-c(EUCAST_mic_big, CLSI_mic_bigeq , CASFM_mic_big)))


data.frame(ZONE_data %>% filter(grepl("TOU",`Code event`)|
                                  grepl("PAR",`Code event`)|
                                  grepl("LIM",`Code event`)) %>% select(Abx) %>% distinct())


# Abx
# 1         Benzilpenicillin
# 2               Ampicillin
# 3              Amoxicillin
# 4  Amoxicillin/clavulanate
# 5  Piperacillin/tazobactam * 
# 6                Cefoxitin
# 7               Cefotaxime
# 8              Clindamycin *
# 9              Ceftriaxone
# 10           Metronidazole *
# 11         Chloramphenicol
# 12               Ertapenem *
# 13                Imipenem *
# 14              Vancomycin
# 15               Linezolid
# 16             Tigecycline
# 17              Rifampicin
# 18            Moxifloxacin
# 19            Tetracycline
# 20                   Genta
# 21                    Levo
# 22                   Oflox


unique(ZONE_data$`EUCAST=1; CLSI=0`)


ZONE_data <- ZONE_data %>% 
  mutate(EUCAST_Diam_Resist=ifelse(ZONE <EUCAST_disc_lower,1,0)) %>%
  mutate(EUCAST_Diam_Resist=ifelse(
    Abx %in% c("Piperacillin/tazobactam", "Clindamycin", "Metronidazole", "Ertapenem", "Imipenem"), EUCAST_Diam_Resist, 
    ifelse(`EUCAST=1; CLSI=0`=="1", EUCAST_Diam_Resist, NA))) %>%
  mutate(CASFM_Diam_Resist=ifelse(ZONE <CASFM_disc_lower,1,0)) %>%
  mutate(CASFM_Diam_Resist=ifelse(
    Abx %in% c("Piperacillin/tazobactam", "Clindamycin", "Metronidazole", "Ertapenem", "Imipenem"), CASFM_Diam_Resist, 
    ifelse(`EUCAST=1; CLSI=0`!="1", CASFM_Diam_Resist, NA)))


temp <- ZONE_data %>% select(`Code event`, Species.identification, speciesID, Abx, EUCAST_Diam_Resist)
temp <- temp %>% drop_na()



library(pheatmap)
library(dplyr)
library(tidyr)

resistance_summary <- temp %>%
  group_by(Species.identification, Abx) %>%
  summarise(
    n_samples = n(),
    resistance_rate = mean(EUCAST_Diam_Resist) * 100
  ) %>%
  filter(n_samples >= 20) 

heatmap_data <- resistance_summary %>%
  select(Species.identification, Abx, resistance_rate) %>%
  pivot_wider(names_from = Abx, values_from = resistance_rate, values_fill = NA)

heatmap_matrix <- as.matrix(heatmap_data[, -1])

rownames(heatmap_matrix) <- heatmap_data$Species.identification

heatmap_matrix[is.na(heatmap_matrix)] <- -99

display_matrix <- round(heatmap_matrix, 0)

plot <- pheatmap(heatmap_matrix,
                 color = colorRampPalette(c("lightgray", "white", "firebrick"))(50),  
                 cluster_rows = TRUE,  # Cluster species
                 cluster_cols = TRUE,  # Cluster antibiotics
                 na_col = "grey",  # Color for missing values
                 fontsize_row = 10,  # Font size for species
                 fontsize_col = 10,  # Font size for antibiotics
                 display_numbers = display_matrix,  # Show exact resistance rates
                 main = "Antibiotic Resistance Clustering [n>20] \n"
)

ggsave(file="dendo.svg", plot=plot, width=6, height=6)






resistance_summary <- temp %>%
  select(`Code event`, Species.identification Abx, EUCAST_Diam_Resist) %>% drop_na() %>%
  group_by(Species.identification, Abx) %>%
  summarise(
    n_samples = n(),
    resistance_rate = mean(EUCAST_Diam_Resist) * 100
  ) %>%
  filter(n_samples >= 20) 

resistance_summary <- resistance_summary %>% mutate(Abx=str_replace_all(Abx, " MIC", ""))

plot <- ggplot(resistance_summary, aes(x = Abx, y = Species.identification)) +
  geom_point(aes(size = n_samples , color = resistance_rate), alpha = 0.9) +
  scale_size(range = c(1, 20)) +  # Adjust the bubble size range
  scale_color_gradient(low = "lightgray", high = "firebrick") +  # Color scale from susceptible (green) to resistant (red)
  labs(title = "Antibiotic Resistance by Species and Antibiotic \n [n>20]",
       x = "Antibiotic",
       y = "Species",
       size = "Sample Size",
       color = "Proportion Resistant") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  


ggsave(file="bubble.svg", plot=plot, width=6, height=8)

# ---------------------


# Number resistant isolates per country -------

MICs <- read_excel(path = "MIC_Workbook_Oct28.xlsx",  sheet="MIC Clean Data", skip = 0)

my_data <- read_excel(path = "ANAEuROBE_dataset_matteo_only.xlsx",  sheet="3.Completo")

my_data <- data.frame(my_data %>% select(`Code event`,  `PT`:`PL`) %>%
                        gather(Group, exp, `PT`:`PL`) %>%
                        filter(exp==1))

plot <- MICs %>% left_join(my_data, by=c("Code event"="Code.event")) %>%
  select(`Code event`, `Species identification`, Abx, `EUCAST Resistant?`, `CLIST Resistant?`, `CASFM Resistant?`, Group) %>%
  group_by(Group, `EUCAST Resistant?`) %>% count() %>% drop_na() %>%
  spread(key=`EUCAST Resistant?`, value=n) %>%
  mutate(perc=`1`/(`1`+`0`), n=`1`+`0`) %>%
  ggplot(aes(n, perc, colour=Group, fill=Group)) +
  geom_jitter(size=0.1, show.legend = F, colour="black") +
  geom_text(aes(label=paste0(Group,"")), 
            fontface="bold", show.legend = F, colour="black")  +
  theme_minimal() +
  xlab("\n Number of Isolates with \n Known EUCAST resistance MIC cutoff ") +
  ylab("% Resistant Isolates-Antibiotics Tested \n")


ggsave(file="eucas.svg", plot=plot, width=6, height=6)


plot <- MICs %>% left_join(my_data, by=c("Code event"="Code.event")) %>%
  select(`Code event`, `Species identification`, Abx, `EUCAST Resistant?`, `CLIST Resistant?`, `CASFM Resistant?`, Group) %>%
  group_by(Group, `CLIST Resistant?`) %>% count() %>% drop_na() %>%
  spread(key=`CLIST Resistant?`, value=n) %>%
  mutate(perc=`1`/(`1`+`0`), n=`1`+`0`) %>%
  ggplot(aes(n, perc, colour=Group, fill=Group)) +
  geom_jitter(size=0.1,  show.legend = F , colour="black") +
  geom_text(aes(label=paste0(Group,"")), 
            fontface="bold", show.legend = F, colour="black")  +
  theme_minimal() +
  xlab("\n Number of Isolates with \n Known CLSI resistance MIC cutoff ") +
  ylab("% Resistant Isolates-Antibiotics Tested \n")

ggsave(file="clsi.svg", plot=plot, width=6, height=6)


plot <-  MICs %>% left_join(my_data, by=c("Code event"="Code.event")) %>%
  select(`Code event`, `Species identification`, Abx, `EUCAST Resistant?`, `CLIST Resistant?`, `CASFM Resistant?`, Group) %>%
  group_by(Group, `CASFM Resistant?`) %>% count() %>% drop_na() %>%
  mutate(`CASFM Resistant?`=as.numeric(`CASFM Resistant?`)) %>%
  spread(key=`CASFM Resistant?`, value=n) %>%
  mutate(perc=`1`/(`1`+`0`), n=`1`+`0`) %>%
  ggplot(aes(n, perc, colour=Group, fill=Group)) +
  geom_jitter(size=0.1, show.legend = F, colour="black") +
  geom_text(aes(label=paste0(Group,"")), 
            fontface="bold", show.legend = F, colour="black")  +
  theme_minimal() +
  xlab("\n Number of Isolates with \n Known CASFM resistance MIC cutoff ") +
  ylab("% Resistant Isolates-Antibiotics Tested \n") 

ggsave(file="csfm.svg", plot=plot, width=6, height=6)

# -------------

# Summary Table MIC and Diameters ----------------

my_data <- read_excel(path = "MIC_Workbook_Oct28.xlsx",  sheet="Summary MIC Values", skip = 1)
unique(my_data$`Species identification`)
TOP <-  my_data %>% mutate(`Abx MIC tested`=str_replace_all(`Abx MIC tested`, " MIC", ""))

my_data <- read_excel(path = "MIC_Workbook_Oct28.xlsx",  sheet="EUCAST Resist Thresh", skip = 0)
my_data <- my_data%>% mutate(Abx=str_replace_all(Abx, " MIC", ""))
my_data <- my_data %>% rename("N_EUCASAT_thre"="# Nr Isolates")

TOP <- TOP %>% left_join(my_data %>% select(-`Species identification`),  by=c("SpeciesID"="SpeciesID", "Abx MIC tested"="Abx"))

my_data <- read_excel(path = "MIC_Workbook_Oct28.xlsx",  sheet="CLSI Resist Thresh", skip = 0)
my_data <- my_data%>% mutate(Abx=str_replace_all(Abx, " MIC", ""))
my_data <- my_data %>% rename("N_CLSI_thre"="# Nr Isolates")

TOP <- TOP %>% left_join(my_data %>% select(-`Species identification`),  by=c("SpeciesID"="SpeciesID", "Abx MIC tested"="Abx"))

my_data <- read_excel(path = "MIC_Workbook_Oct28.xlsx",  sheet="CASFM Resist Thresh", skip = 0)
my_data <- my_data%>% mutate(Abx=str_replace_all(Abx, " MIC", ""))
my_data <- my_data %>% rename("N_CASFM_thre"="# Nr Isolates")

TOP <- TOP %>% left_join(my_data %>% select(-`Species identification`),  by=c("SpeciesID"="SpeciesID", "Abx MIC tested"="Abx"))

fwrite(TOP, "MIC_Summary_All.csv")



my_data <- read_excel(path = "INHIB_ZONE_Workbook_Oct28.xlsx",  sheet="Summary Zone Diam Values", skip = 1)
TOP <- my_data

my_data  <- read_excel(path = "INHIB_ZONE_Workbook_Oct28.xlsx",  sheet="EUCAST Resist Thresh", skip = 0)
my_data <- my_data %>% rename("N_EUCASAT_thre"="# Nr Isolates")

TOP <- TOP %>% left_join(my_data %>% select(-`Species identification`),  by=c("SpeciesID"="SpeciesID", "Abx MIC tested"="Abx"))

my_data <- read_excel(path = "INHIB_ZONE_Workbook_Oct28.xlsx",  sheet="CASFM Resist Thresh", skip = 0)
my_data <- my_data %>% rename("N_CASFM_thre"="# Nr Isolates")

TOP <- TOP %>% left_join(my_data %>% select(-`Species identification`),  by=c("SpeciesID"="SpeciesID", "Abx MIC tested"="Abx"))

data.frame(TOP)

fwrite(TOP, "Diams_Summary_All.csv")


ZI <- read_excel(path = "INHIB_ZONE_Workbook_Oct28.xlsx",  sheet="ZI", col_types = "text")
ZI$SpeciesID <- as.numeric(ZI$SpeciesID)

Diams_Summary_All <- fread("Diams_Summary_All.csv")

Diams_Summary_All <- Diams_Summary_All %>%
  left_join(ZI %>% group_by(SpeciesID, Abx) %>%
              summarise(mean=mean(as.numeric(ZI)), n=n()) %>%
              rename("# AUT"="n", "Inside AUT"="mean"),
            by=c("SpeciesID"="SpeciesID","Abx MIC tested"="Abx"))


fwrite(Diams_Summary_All, "Diams_Summary_All.csv")
# -------------