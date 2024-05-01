
###############################################################################
#### GLOBAL CASE FATALITY OF BACTERIAL MENINGITIS OVER AN 80-YEAR PERIOD:  ####
#### A SYSTEMATIC REVIEW AND META-ANALYSIS ####################################
###############################################################################

# Ettekoven, Liechti, Brouwer, Bijlsma, van de Beek (2024)

# Prepare ######################################################################
library(tidyverse)
library(gtsummary)
library(meta)
library(metafor)
library(sf)

list_Agegroup <- c(
  "Neonates", "Neonates and children", "Children", "Children and adults",
  "Adults", "Not specified")

list_pathogens <- c(
  "Streppneumoniae", "Neismeningitidis", "HaemophInfluenzae",
  "Listeria", "Eschcoli", "GBS", "Staphaureus", "Enterobacter",
  "Pseudomonas", "Klebsielapneumoniae", "Other", "Unidentified")

list_pathogens_other <- c(
  "Staphaureus", "Enterobacter", "Pseudomonas", "Klebsielapneumoniae", "Other")

pathogen <- "sp" # choose pathogen from: "sp", "nm", "hib"  "ec"  "gbs"  "lm"  ("na" = all)

## Read file ===================================================================
setwd("G:/divd/neu/mening/Fabian/MortBM-MA/5 Manuscript/Submission/github")
data <- read_rds(paste0("data_", pathogen, ".RDS"))
data_merge <- read_rds(paste0("data_merge_", pathogen, ".RDS"))

list_splitdates_levels <- levels(data$Interval)

## Create map ==================================================================
# https://bookdown.org/robinlovelace/geocompr/adv-map.html
# https://thinking-spatial.org/courses/angewandte_geodatenverarbeitung/kurs04/
tab_countries <- data_merge %>%
  group_by(Country) %>%
  summarise(n = n()) %>%
  rename(name = Country) %>%
  arrange(desc(n)) # List of countries with number of studies

World <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
data_sf <- left_join(World, tab_countries) # Join data sets to sf object

### Number of studies ----------------------------------------------------------
tmap::tm_shape(data_sf, projection = "+proj=robin") +
  tmap::tm_fill() +
  tmap::tm_borders(lwd = 1.5) +
  tmap::tm_bubbles(size = "n", col = "blue", 
                   title.size = "Number of studies") +
  tmap::tm_layout(legend.position = c("left", "bottom"))


### Mortality ------------------------------------------------------------------
mortalityrates <- data %>%
  group_by(Country, Interval) %>%
  summarise(mr = mean(CfrCalc)) %>%
  rename(name = Country) %>% # List of countries with mortality rates
  droplevels()

World_all <- World %>%
  mutate(Interval  = list_splitdates_levels[1])
for (w in c(2:length(list_splitdates_levels))) {
  world <- World %>%
    mutate(Interval = list_splitdates_levels[w])
  
  World_all <- World_all %>%
    add_row(world)
}

data_sf_mr <- left_join(World_all, mortalityrates) %>% # Join data sets to sf object
  mutate(Interval = factor(Interval,
                           labels = list_splitdates_levels,
                           levels = list_splitdates_levels))

ggplot() +
  geom_sf(data = data_sf_mr, aes(fill = mr)) +
  scale_fill_gradientn(#limits = c(0.25, 0.75),
    colours=c("white", "darkred"),
    #breaks = c(0.5, 0.5, 0.75),
    name = "Case fatality ratio") +
  #ggtitle(list_splitdates_levels[1]) +
  facet_wrap(. ~ Interval, ncol = 2) +
  theme_classic()  +
  theme(legend.position = c(0.72, 0.1),
        legend.title = element_text( size = rel(.7)), 
        legend.text=element_text(size = rel(.7)),
        legend.key.height = unit(3, 'mm'),
        legend.key.width = unit(4, 'mm'),
        strip.text = element_text(size = rel(.7)))


# Primary analysis #############################################################
## All studies =================================================================
hist(data$Mortalityrate) # explore distribution -> not normally distributed

data_ma <- metafor::escalc(xi = Deaths, ni = Patients, data = data, 
                           measure = "PLO") # logit transformation

ma_all <- metaprop(Deaths, Patients, studlab = Studyauthor, data = data_ma, 
                   sm = "PLO", method.tau = "REML", method.ci = "NAsm", 
                   #control = list(maxiter = 1000),
                   prediction = TRUE, subgroup = Interval)
summary(ma_all)

### Forest plot
metafor::forest(ma_all, 
                sortvar = Periodmean,
                study.results = FALSE, subgroup = k.w > 2,
                print.subgroup.name = FALSE,
                leftcols = c("studlab"), leftlabs = "Subgroup",
                rightlabs = c("CFR", "[95% CI]"),
                prediction = TRUE, print.I2 = TRUE, print.Q = TRUE, 
                print.tau2 = TRUE, common = FALSE, 
                smlab = "",
                col.square = "navy", col.square.lines = "navy",
                col.diamond = "grey", col.diamond.lines = "navy", digits = 2,
                digits.tau2 = 2, digits.TE = 2, digits.se = 2, digits.Q = 1,
                hetlab = "Het.: ", resid.hetlab = "Res.het.: ",
                label.test.subgroup.random = "Subgroup diff.: ",
                prediction.subgroup = k.w > 2,
                fs.hetstat = 8, col.subgroup = "black", xlim = c(0,1))

### Funnel plot
meta::funnel(ma_all, studlab = FALSE, yaxis = "size", common = FALSE, log = "y")
title("All studies")

### Meta-regression
metareg_all <- meta::metareg(ma_all, Periodmean)
summary(metareg_all)
metafor::regplot(metareg_all, pi = TRUE, mod = "Periodmean", transf = transf.ilogit, 
                 legend = FALSE, bg = scales::alpha("#3C5488B2", .3), 
                 xlim = c(1940, 2020), ylim = c(0, 1),
                 ylab = "Case fatality ratio")
title("All studies")

## Moderator analysis (subgroups) ==============================================
list_age <- c(1,3,5)
table_ma <- tribble(
  ~Income, ~Agegroup, ~n, ~re, ~relci, ~reuci, ~ pilci, ~ piuci, 
  ~I2, ~I2uci, ~I2lci) %>%
  add_row(Income = "high and low", Agegroup = "Combined",
          n = ma_all[["k.study"]],
          re = transf.ilogit(ma_all[["TE.random"]]),
          relci = transf.ilogit(ma_all[["lower.random"]]),
          reuci = transf.ilogit(ma_all[["upper.random"]]),
          pilci = transf.ilogit(ma_all[["lower.predict"]]),
          piuci = transf.ilogit(ma_all[["upper.predict"]]),
          I2 = ma_all[["I2"]],
          I2lci = ma_all[["lower.I2"]],
          I2uci = ma_all[["upper.I2"]])

table_mreg <- tribble(~Income, ~Agegroup, ~k, ~beta, ~se, ~ci.lb, ~ci.ub, 
                      ~tau2, ~I2, ~R2, ~QEp, ~QMp) %>%
  add_row(Income = "high and low", Agegroup = "Combined",
          k = metareg_all[["k"]],
          beta = metareg_all$beta[2,1],
          se = metareg_all$se[2],
          ci.lb = metareg_all$ci.lb[2],
          ci.ub = metareg_all$ci.ub[2],
          tau2 = metareg_all[["tau2"]],
          I2 = metareg_all[["I2"]],
          R2 = metareg_all[["R2"]],
          QEp = metareg_all[["QEp"]],
          QMp = metareg_all[["QMp"]])


for (i in levels(data$Income)) {
  for (k in levels(data$Agegroup)[list_age]) {
    data_ma <- subset(data, Income == i & Agegroup == k) 
    
    ### Meta-analysis
    ma_int <- metaprop(Deaths, Patients, studlab = Studyauthor, data = data_ma, 
                       sm = "PLO", method.tau = "REML", method.ci = "NAsm", 
                       control = list(maxiter = 1000),
                       prediction = TRUE, subgroup = Interval, tau.common = FALSE)
    ma_int
    
    table_ma <- table_ma %>%
      add_row(Income = i, Agegroup = k,
              n = ma_int[["k.study"]],
              re = transf.ilogit(ma_int[["TE.random"]]),
              relci = transf.ilogit(ma_int[["lower.random"]]),
              reuci = transf.ilogit(ma_int[["upper.random"]]),
              pilci = transf.ilogit(ma_int[["lower.predict"]]),
              piuci = transf.ilogit(ma_int[["upper.predict"]]),
              I2 = ma_int[["I2"]],
              I2lci = ma_int[["lower.I2"]],
              I2uci = ma_int[["upper.I2"]])
    
    
    ### Forest plot
    tiff(paste0("Rplot_Forest_", i, "_", k, ".tif"), unit = "mm",
         width = 2*107, height = (round((ma_int$k.all)/12, 0)+2)*80, res = 300)
    metafor::forest(ma_int,
                    sortvar = Periodmean,
                    leftcols = c("studlab"), 
                    leftlabs = "Author (publication year) [study period]",
                    rightlabs = c("CFR", "[95% CI]", "Weight"),
                    print.subgroup.name = FALSE,
                    subgroup = k.w > 2,
                    prediction = TRUE, print.I2 = TRUE, print.Q = TRUE, 
                    print.tau2 = TRUE, common = FALSE, 
                    smlab = "",
                    col.square = "navy", col.square.lines = "navy",
                    col.diamond = "grey", col.diamond.lines = "navy",
                    digits.tau2 = 2, digits.TE = 2, digits.se = 2, digits.Q = 1,
                    hetlab = "Het.: ", resid.hetlab = "Res.het.: ",
                    label.test.subgroup.random = "Subgroup diff.: ",
                    prediction.subgroup = k.w > 2,
                    fs.hetstat = 8, col.subgroup = "black")
    dev.off()
    
    ### Funnel plot
    meta::funnel(ma_int, studlab = FALSE, yaxis = "size", common = FALSE)
    title(paste0(k, ", ", i, "-income countries"))
    
    ### Meta-regression
    metareg_int <- meta::metareg(
      ma_int, Periodmean, 
      method = "REML", # If Fisher algorithm does not convert use "DL"
      control = list(maxiter = 1000))
    
    table_mreg <- table_mreg %>%
      add_row(Income = i, Agegroup = k,
              k = metareg_int[["k"]],
              beta = metareg_int$beta[2,1],
              se = metareg_int$se[2],
              ci.lb = metareg_int$ci.lb[2],
              ci.ub = metareg_int$ci.ub[2],
              tau2 = metareg_int[["tau2"]],
              I2 = metareg_int[["I2"]],
              R2 = metareg_int[["R2"]],
              QEp = metareg_int[["QEp"]],
              QMp = metareg_int[["QMp"]])
    
    tiff(paste0("Rplot_metareg_", i, "_", k, ".tif"),
         unit = "mm", width = 2*107, height = 1.5*80, res = 300)
    metafor::regplot(metareg_int, pi = TRUE, mod = "Periodmean", transf = transf.ilogit, 
                     legend = FALSE, bg = scales::alpha("#3C5488B2", .3), 
                     xlim = c(1940, 2020), ylim = c(0, 1),
                     ylab = "Case fatality ratio", xlab = "Year")
    title(paste0(k, ", ", i, "-income countries"))
    dev.off()
  }
}

  table_ma %>%
    mutate_if(is.numeric, round, 3) %>%
    mutate(CFR = paste0(re*100, " (", relci*100, "-", reuci*100, ")")) %>%
    mutate(I2 = round(I2,2)*100) %>%
    select(c("Income", "Agegroup", "Study periods (k)" = "n", 
             "Case fatality ratio (95% CI) [%]" = "CFR", "I^2 [%]" = "I2")) %>%
    filter(Agegroup %in% c("Combined", list_Agegroup[list_age])) 

  table_mreg %>%
    mutate_if(is.numeric, round, 3) %>%
    filter(Agegroup %in% c("Combined", list_Agegroup[list_age])) %>%
    mutate("Change per year (95% CI) [%]" = paste0(
      round(beta, 3)*100, " (", 
      round(ci.lb, 3)*100, " to ", 
      round(ci.ub, 3)*100, ")")) %>%
    mutate(I2 = round(I2, 0)) %>%
    mutate(R2 = round(R2, 0)) %>%
    select(c(1:3, 13, 8:12)) 


## Pathogens ===================================================================
### List of pathogens
list_pathogens_labels <- c(
  "S. pneumoniae", "N. meningitidis", "H. influenzae", "L. monocytogenes", 
  "E. coli", "S. agalactiae", "S. aureus", "Enterobacter spp.", 
  "Pseudomonas spp.", "K. pneumoniae", "Other bacteria", "Unidentified")

list_pathogens_woothers_labels <- c(
  "S. pneumoniae", "N. meningitidis", "H. influenzae",
  "L. monocytogenes", "E. coli", "S. agalactiae", "Other bacteria", "Unidentified")


### Prepare data frame for pathogens -------------------------------------------
table_path2 <- tribble(
  ~Pathogen, ~Income, ~Agegroup, ~Deaths, ~Episodes, ~Interval)

for (p in list_pathogens) {
  for (i in levels(data$Income)) {
    for (k in levels(data$Agegroup)) {
      for (l in levels(data$Interval)) {
        df <- data %>%
          filter(Income == i & Agegroup == k & Interval == l)
        
        table_path2 <- table_path2 %>%
          add_row(Pathogen = p,
                  Income = i,
                  Agegroup = k,
                  Interval = l,
                  Episodes = sum(df[, p], na.rm = TRUE),
                  Deaths = sum(df[, paste0(p, "_A")], na.rm = TRUE)) %>%
          mutate(cfr = Deaths / Episodes *100) %>%
          mutate(cfr_lci = cfr + 1.96*(cfr*((1 - cfr)/Episodes))) %>%
          mutate(cfr_uci = cfr - 1.96*(cfr*((1 - cfr)/Episodes))) %>%
          mutate(Pathogen = factor(Pathogen,
                                   levels = list_pathogens,
                                   labels = list_pathogens)) %>%
          mutate(Income = factor(Income)) %>%
          mutate(Agegroup = factor(Agegroup)) %>%
          mutate(Interval = factor(Interval,
                                   levels = list_splitdates_levels,
                                   labels = list_splitdates_levels))
      }}}}
table_path2


### Tables for pathogen distributions (descriptive) ----------------------------
#### Number of episodes per pathogen
table_path2 %>%
  group_by(Pathogen) %>%
  summarise(n = sum(Episodes)) %>%
  arrange(desc(n)) %>%
  bind_rows(summarise(.,
                      across(where(is.numeric), sum),
                      across(where(is.factor), ~"Total"))) %>%
  mutate(Proportion = round(1/n[Pathogen=="Total"]*n, 3)*100)

#### Number of episodes per interval
table_path2 %>%
  mutate(Pathogen = fct_collapse(
    Pathogen,
    "Identified" = c(list_pathogens[1:11]),
    "Unidentified" = "Unidentified")) %>%
  group_by(Pathogen, Interval) %>%
  summarise(n = sum(Episodes)) %>%
  pivot_wider(names_from = Pathogen, values_from = n) %>%
  bind_rows(summarise(.,
                      across(where(is.numeric), sum),
                      across(where(is.factor), ~ "Total"))) %>%
  mutate(Proportion = round(Identified/(Identified + Unidentified), 3)*100) 

#### Proportion of episodes per pathogen
table_path2 %>%
  droplevels() %>%
  group_by(Pathogen) %>%
  summarise(Deaths = sum(Deaths), Episodes = sum(Episodes)) %>%
  mutate(Proportion = round(Episodes/sum(Episodes)*100, 1)) %>%
  arrange(desc(Episodes)) %>%
  bind_rows(summarise(.,
                      across(where(is.numeric), sum),
                      across(where(is.factor), ~ "Total"))) 

#### Proportion of episodes per pathogen stratified by interval
table_path2 %>%
  group_by(Interval, Pathogen) %>%
  summarise(Deaths = sum(Deaths), Episodes = sum(Episodes)) %>%
  mutate(Proportion = round(Episodes/sum(Episodes)*100, 1)) %>%
  arrange(Pathogen) %>%
  bind_rows(summarise(.,
                        across(where(is.numeric), sum),
                        across(where(is.factor), ~ "Total"))) %>%
  print(n = nrow(.))

#### Proportion of episodes per pathogen stratified by age group
table_path2 %>%
  filter(Agegroup %in% list_Agegroup[c(1,3,5)]) %>%
  droplevels() %>%
  group_by(Agegroup, Pathogen) %>%
  summarise(Deaths = sum(Deaths), Episodes = sum(Episodes)) %>%
  mutate(cfr = Deaths / Episodes *100) %>%
  mutate(cfr_lci = cfr + 1.96*(cfr*((1 - cfr)/Episodes))) %>%
  mutate(cfr_uci = cfr - 1.96*(cfr*((1 - cfr)/Episodes))) %>%
  mutate('CFR (95% CI) [%]' = paste0(
    format(round(cfr, 1), nsmall = 1), " (", 
    format(round(cfr_lci, 1), nsmall = 1), "-", 
    format(round(cfr_uci, 1), nsmall = 1), ")")) %>%
  select(-c(cfr, cfr_lci, cfr_uci)) %>%
  mutate(Proportion = round(Episodes/sum(Episodes)*100, 1)) %>%
  mutate(Pathogen = factor(Pathogen,
                           levels = list_pathogens,
                           labels = list_pathogens_labels)) %>%
  arrange(Agegroup) %>%
  bind_rows(summarise(.,
                      across(where(is.numeric), sum),
                      across(where(is.factor), ~ "Total"))) %>%
  print(n = nrow(.))

#### Proportion of episodes per pathogen stratified by income
table_path2 %>%
  group_by(Income, Pathogen) %>%
  summarise(Deaths = sum(Deaths), Episodes = sum(Episodes)) %>%
  mutate(cfr = Deaths / Episodes *100) %>%
  mutate(cfr_lci = cfr + 1.96*(cfr*((1 - cfr)/Episodes))) %>%
  mutate(cfr_uci = cfr - 1.96*(cfr*((1 - cfr)/Episodes))) %>%
  bind_rows(summarise(.,
                      across(where(is.numeric), sum),
                      across(where(is.factor), ~ "Total"))) %>%
  mutate("Crude case fatality ratio (95% CI) [%]" = paste0(
    format(round(cfr, 1), nsmall = 1), " (", 
    format(round(cfr_lci, 1), nsmall = 1), "-", 
    format(round(cfr_uci, 1), nsmall = 1), ")")) %>%
    select(-c(cfr, cfr_lci, cfr_uci)) %>%
    mutate(Proportion = round(Episodes/sum(Episodes)*100, 1)) %>%
    mutate(Pathogen = factor(
      Pathogen,
      levels = c(list_pathogens, "Total"),
      labels = c(list_pathogens_labels, "Subotal"))) %>%
    arrange(Income) %>%
    print(n = nrow(.)) 


### Plots for pathogen distributions (descriptive) -----------------------------
table_path3 <- table_path2 %>%
  filter(Pathogen != "Unidentified") %>%
  mutate(Pathogen = fct_collapse(
    Pathogen, 
    "Other bact." = c(list_pathogens_other),
    "Streppneumoniae" = "Streppneumoniae",
    "Neismeningitidis" = "Neismeningitidis",
    "HaemophInfluenzae" = "HaemophInfluenzae",
    "Listeria" = "Listeria",
    "Eschcoli"  = "Eschcoli",
    "GBS" = "GBS")) %>%
  droplevels() 

list_col <- c("#FFDB6D", "#C4961A", "#F4EDCA", "#D16103", "#C3D7A4", 
              "#52854C", "#4E84C4", "#7E6148B2", "#DC0000B2", "#293352", 
              "darkgrey")

table_path3 %>%
  droplevels() %>%
  group_by(Pathogen, Interval) %>%
  summarise(Deaths = sum(Deaths), Episodes = sum(Episodes)) %>%
  group_by(Interval) %>%
  mutate(Proportion = Episodes/sum(Episodes)) %>%
  ggplot() +
  aes(x = Interval, y = Proportion, fill = Pathogen) +
  geom_col(color = "black", linewidth = 0.1) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_bw() +
  scale_fill_manual(
    values = list_col[c(1:6, 11)], 
    labels = list_pathogens_woothers_labels) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

table_path3%>%
  droplevels() %>%
  group_by(Pathogen, Interval, Income) %>%
  summarise(Deaths = sum(Deaths), Episodes = sum(Episodes)) %>%
  group_by(Interval, Income) %>%
  mutate(Proportion = Episodes/sum(Episodes)) %>%
  ggplot() +
  aes(x = Interval, y = Proportion, fill = Pathogen) +
  geom_col(color = "black", linewidth = 0.1) +
  theme_bw() +
  scale_fill_manual(
    values = list_col[c(1:6, 11)], 
    labels = list_pathogens_woothers_labels) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ Income)

table_path2 %>%
  filter(Agegroup %in% list_Agegroup[c(1,3,5)]) %>%
  filter(Pathogen != "Unidentified") %>%
  mutate(Agegroup = factor(Agegroup,
                           levels = list_Agegroup[c(1,3,5)])) %>%
  droplevels() %>%
  group_by(Pathogen, Interval, Agegroup) %>%
  summarise(Deaths = sum(Deaths), Episodes = sum(Episodes)) %>%
  group_by(Interval, Agegroup) %>%
  mutate(Proportion = Episodes/sum(Episodes)) %>%
  filter(Episodes > 0) %>% 
  ggplot() +
  aes(x = Interval, y = Proportion, fill = Pathogen) +
  geom_col(color = "black") +
  theme_bw() +
  scale_fill_manual(
    values = list_col, 
    labels = list_pathogens_labels) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ Agegroup)
