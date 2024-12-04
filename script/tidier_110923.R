library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(ciTools)
library(modelr)
library(tidyverse)
library(forcats)
library(ggrepel)

#data from multi-species co-invasions
inv.coms <- read.csv("All_inv_dat_240219.csv", header = T)

### some code just to tidy names and align with the separate dataset on single invasions

strSort <- function(x)
  sapply(lapply(strsplit(x, NULL), sort), paste, collapse="")

inv.coms <- mutate(inv.coms,
                   inv.cfu = count * dil * dil_fact *2)

inv.coms <- mutate(inv.coms,
                   res.cfu = count_other * dil * dil_fact *2) 

inv.coms <- mutate(inv.coms,
                   inv.fit = log(inv.cfu/T0_inv))

inv.coms <- mutate(inv.coms,
                   res.fit = log(res.cfu/T0_other))

inv.coms <- mutate(inv.coms,
                   rel.fit = inv.fit / res.fit) 

inv.coms$res.sp <- gsub('[.]','',inv.coms$res.sp)

inv_dat <- inv.coms[c(1:450), c(2:8, 10, 16:18)]
inv_dat$rel.fit[inv_dat$rel.fit == -Inf] <- 0

inv_dat <- group_by(inv_dat, Sp, res.sp, inv.other, res.div, inv.div) %>%
  summarise_at(vars(rel.fit), list(mean = mean))

summary(inv_dat) 
inv_dat$res.div <- as.numeric(inv_dat$res.div)

inv_dat$res.sp <- gsub('VG',"V",inv_dat$res.sp)
inv_dat$res.sp <- gsub('SR',"S",inv_dat$res.sp)
inv_dat$res.sp <- gsub('PC',"P",inv_dat$res.sp)
inv_dat$res.sp <- gsub('AA',"A",inv_dat$res.sp)
inv_dat$res.sp <- gsub('OD',"O",inv_dat$res.sp)
inv_dat$Sp[inv_dat$Sp == "VG"] <- "V"
inv_dat$Sp[inv_dat$Sp == "SR"] <- "S"
inv_dat$Sp[inv_dat$Sp == "AA"] <- "A"
inv_dat$Sp[inv_dat$Sp == "PC"] <- "P"
inv_dat$Sp[inv_dat$Sp == "OD"] <- "O"

inv_dat2 <- inv_dat

inv_dat2$res.sp <- gsub('[.]','',inv_dat2$res.sp)

inv.coms$res.sp <- gsub('VG',"V",inv.coms$res.sp)
inv.coms$res.sp <- gsub('SR',"S",inv.coms$res.sp)
inv.coms$res.sp <- gsub('PC',"P",inv.coms$res.sp)
inv.coms$res.sp <- gsub('AA',"A",inv.coms$res.sp)
inv.coms$res.sp <- gsub('OD',"O",inv.coms$res.sp)
inv.coms$Sp[inv.coms$Sp == "VG"] <- "V"
inv.coms$Sp[inv.coms$Sp == "SR"] <- "S"
inv.coms$Sp[inv.coms$Sp == "AA"] <- "A"
inv.coms$Sp[inv.coms$Sp == "PC"] <- "P"
inv.coms$Sp[inv.coms$Sp == "OD"] <- "O"

inv.coms2 <- inv.coms

inv.coms2$res.sp <- gsub('[.]','',inv.coms2$res.sp)

fac_labs <- c("(a) Achromobacter sp.", "(b) Ochrobactrum sp.", "(c) Pseudomonas sp.", "(d) Stenotrophomonas sp.", "(e) Variovorax sp.")

names(fac_labs) <- c("AA", "OD", "PC", "SR", "VG")
names(fac_labs) <- c("A", "O", "P", "S", "V")

inv_dat2$res.sp <- strSort(inv_dat2$res.sp)
inv.coms$res.sp <- strSort(inv.coms$res.sp)

inv_dat2$treat <- interaction(inv_dat2$Sp, inv_dat2$res.sp, sep = "_")
inv.coms$treat <- interaction(inv.coms$Sp, inv.coms$res.sp, sep = "_")

inv.coms$rel.fit[inv.coms$rel.fit == -Inf] <- 0

#calculate where relative growth rate measures are significantly different to 1

d_mods <- inv.coms %>%
  nest(-treat) %>%
  mutate(., mod = map(data, ~ t.test(.x$rel.fit, mu = 1)))

d_mods_summary <- d_mods %>%
  mutate(params = map(mod, broom::tidy))

d_mods_summary <- unnest(d_mods_summary, params) %>%
  mutate(padjust = p.adjust(p.value, "fdr"))

d_mods_summary$padjust2 <- round(d_mods_summary$padjust, 5) 

d_mods_summary$sig <- NA
d_mods_summary$sig[d_mods_summary$padjust2 > 0.05] <- 'NS'
d_mods_summary$p.value2 <- round(d_mods_summary$p.value, 5)

d_mods_summary2 <- d_mods_summary[,c(1,14)]

inv_dat2$sig <- NA

inv_dat3 <- merge(d_mods_summary2, inv_dat2, by = 'treat')

inv_dat3$sig.x[inv_dat3$Sp == "A" & inv_dat3$res.sp == "PSV"] <- "NS-s"
inv_dat3$sig.x[inv_dat3$Sp == "A" & inv_dat3$res.sp == "PVO"] <- "NS-s"
inv_dat3$sig.x[inv_dat3$Sp == "A" & inv_dat3$res.sp == "PSVO"] <- "NS-b"
inv_dat3$sig.x[inv_dat3$Sp == "S" & inv_dat3$res.sp == "AVO"] <- "NS-b"

#Work from the single species invasions - Castledine et al. 2024, Microbiology

inv_v <- read.csv("relative_fitness_params.csv", header = T)

#mean invasion success per combination
inv_v2 <- group_by(inv_v, inv_sp, res_sp, res_div) %>%
  summarise_at(vars(rel_fit, inv_m), list(mean = mean, sd = sd)) %>%
  ungroup()

#make resident combinations in same order
strSort <- function(x)
  sapply(lapply(strsplit(x, NULL), sort), paste, collapse="")

inv_v2$res_sp2 <- strSort(inv_v2$res_sp)
inv_dat3$res.sp2 <- strSort(inv_dat3$res.sp)

inv_v2$res_sp2 <- toupper(inv_v2$res_sp2)
inv_v2$Sp <- toupper(inv_v2$inv_sp)

inv_dat3 <- inv_dat3[order(inv_dat3$res.sp2),]
inv_v2 <- inv_v2[order(inv_v2$res_sp2),]

inv_dat3$s_relfit <- inv_v2$rel_fit_mean

#see if relative fitness measures correlated when invading one vs multiple (one duplicate combination is each species into the four resident combination)

is.numeric(inv_dat3$res.div)

m1 <- lm(mean ~ s_relfit * Sp, data = inv_dat3)
summary(m1)
plot(m1) #looks good
m2 <- lm(mean ~ 0 + s_relfit + Sp, data = inv_dat3)
anova(m1,m2,test = "F") #no interaction
m3 <- lm(mean ~ Sp, data = inv_dat3)
anova(m2,m3,test="F") #sig
m4 <- lm(mean ~ s_relfit, data = inv_dat3)
anova(m2,m4,test="F") #sig
summary(m2)

emmeans::emmeans(m2, pairwise ~ Sp)

new <- add_ci(inv_dat3, m2, alpha = 0.5, type = "boot", includeRanef = FALSE, nSims = 2000)

color_scheme <- c(A = "#a06aa0", O = '#db596b', P = '#e79652', S = '#85cdc4', V = '#31646f')

#Figure 2

ggplot() +
  theme_bw() +
  geom_abline(aes(slope = 1, intercept = 0))+
  geom_point(data = inv_dat3, aes(x = s_relfit, y = mean, col = Sp, group = Sp)) + 
  geom_line(data = new, aes(x = s_relfit, y = pred, col = Sp, group = Sp)) +
  geom_ribbon(data = new, aes(x = s_relfit, ymin = LCB0.25, ymax = UCB0.75, group = Sp, fill = Sp), alpha = 0.5) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 12, hjust = 0), axis.text = element_text(size = 11, colour = 'black'), axis.title = element_text(size = 12, colour = 'black'), legend.position = "bottom") +
  ylab("Relative invader growth rate from coinvasions") +
  xlab("Relative invader growth rate from single invasions") +
  labs(colour = "Invading species", fill = "Invading species") +scale_fill_manual(values = color_scheme) +scale_color_manual(values = color_scheme) 

inv_dat3$diff <- inv_dat3$s_relfit - inv_dat3$mean

nrow(inv_dat3[inv_dat3$Sp == 'A' & inv_dat3$diff > 0, ])
nrow(inv_dat3[inv_dat3$Sp == 'O' & inv_dat3$diff > 0, ])
nrow(inv_dat3[inv_dat3$Sp == 'S' & inv_dat3$diff > 0, ])
nrow(inv_dat3[inv_dat3$Sp == 'P' & inv_dat3$diff > 0, ])
nrow(inv_dat3[inv_dat3$Sp == 'V' & inv_dat3$diff > 0, ])

inv.coms$Sp[inv.coms$Sp == "VG"] <- "V"
inv.coms$Sp[inv.coms$Sp == "SR"] <- "S"
inv.coms$Sp[inv.coms$Sp == "AA"] <- "A"
inv.coms$Sp[inv.coms$Sp == "PC"] <- "P"
inv.coms$Sp[inv.coms$Sp == "OD"] <- "O"

inv.coms$res.sp <- gsub('VG',"V",inv.coms$res.sp)
inv.coms$res.sp <- gsub('SR',"S",inv.coms$res.sp)
inv.coms$res.sp <- gsub('PC',"P",inv.coms$res.sp)
inv.coms$res.sp <- gsub('AA',"A",inv.coms$res.sp)
inv.coms$res.sp <- gsub('OD',"O",inv.coms$res.sp)

inv.coms$res.sp <- gsub('[.]','',inv.coms$res.sp)

sing_invs <- read.csv("rel_fitness_sing_inv.csv", header = T)
names(sing_invs)[3] <- "Sp"
names(sing_invs)[4] <- "res.sp"

sing_invs$Sp[sing_invs$Sp == "VG"] <- "V"
sing_invs$Sp[sing_invs$Sp == "SR"] <- "S"
sing_invs$Sp[sing_invs$Sp == "AA"] <- "A"
sing_invs$Sp[sing_invs$Sp == "PC"] <- "P"
sing_invs$Sp[sing_invs$Sp == "OD"] <- "O"

sing_invs$res.sp <- gsub('VG',"V",sing_invs$res.sp)
sing_invs$res.sp <- gsub('SR',"S",sing_invs$res.sp)
sing_invs$res.sp <- gsub('PC',"P",sing_invs$res.sp)
sing_invs$res.sp <- gsub('AA',"A",sing_invs$res.sp)
sing_invs$res.sp <- gsub('OD',"O",sing_invs$res.sp)

sing_invs$res.sp <- gsub('[.]','',sing_invs$res.sp)

sing_invs$res.sp <- strSort(sing_invs$res.sp)
inv.coms$res.sp <- strSort(inv.coms$res.sp)

sing_invs <- sing_invs[,c(3,4,12)]
sing_invs$invasion <- "single"
names(sing_invs)[3] <- "rel.fit"

inv.coms <- inv.coms[,c(3,4,18)]
inv.coms$invasion <- "multi"

all_invs <- merge(sing_invs, inv.coms, all = T)

all_invs$rel.fit[all_invs$rel.fit == -Inf] <- 0

all_invs$treat <- interaction(all_invs$Sp, all_invs$res.sp)

#remove resident diversity 4 as duplicate data between studies
all_invs$res.div <- nchar(all_invs$res.sp)
all_invs <- filter(all_invs, ! res.div == 4)
summary(as.factor(all_invs$treat))

m1 <- lm(rel.fit ~ treat * invasion, data = all_invs)
m2 <- lm(rel.fit ~ treat + invasion, data = all_invs)
anova(m1,m2,test="F") #significant
emmeans::emmeans(m1, pairwise ~ invasion|treat)
#42/70 non-sig

cont <- data.frame(emmeans::emmeans(m1, pairwise ~ invasion|treat)$contrasts)
non_sig <- cont[cont$p.value > 0.05,]
sig <- cont[cont$p.value < 0.05,]

cont <- cont[,c(2,7)]

all_invs2 <- group_by(all_invs, Sp, res.sp, treat) %>%
  summarise_at(vars(rel.fit), list(mean = mean)) %>%
  ungroup()

all_invs3 <- group_by(all_invs, Sp, res.sp, treat, invasion) %>%
  summarise_at(vars(rel.fit), list(mean = mean)) %>%
  ungroup()

all_invs2 <- merge(all_invs2,cont,by = "treat")

all_invs2$sig[all_invs2$p.value < 0.05] <- "*"

all_invs2 <- all_invs2[order(all_invs2$res.sp),]

inv_dat6 <- filter(inv_dat3, ! res.div == 4)

all_invs2$inv <- inv_dat6$sig.x

fac_labs <- c(A = expression((a)~italic(Achromobacter)~sp.), O = expression((b)~italic(Ochrobactrum)~sp.), P = expression((c)~italic(Pseudomonas)~sp.), S = expression((d)~italic(Stenotrophomonas)~sp.), V = expression((e)~italic(Variovorax)~sp.))

fac_labs <- c("(a) Achromobacter sp.", "(b) Ochrobactrum sp.", "(c) Pseudomonas sp.", "(d) Stenotrophomonas sp.", "(e) Variovorax sp.")

names(fac_labs) <- c("A", "O", "P", "S", "V")

all_invs <- mutate(all_invs,
                    Sp2 = forcats::fct_recode(Sp, `(a)~italic("Achromobacter")~sp.` = 'A', `(b)~italic("Ochrobactrum")~sp.` = 'O', `(c)~italic("Pseudomonas")~sp.` = 'P', `(d)~italic("Stenotrophomonas")~sp.` = 'S', `(e)~italic("Variovorax")~sp.` = 'V'))

all_invs2 <- mutate(all_invs2,
                    Sp2 = forcats::fct_recode(Sp, `(a)~italic("Achromobacter")~sp.` = 'A', `(b)~italic("Ochrobactrum")~sp.` = 'O', `(c)~italic("Pseudomonas")~sp.` = 'P', `(d)~italic("Stenotrophomonas")~sp.` = 'S', `(e)~italic("Variovorax")~sp.` = 'V'))

all_invs3 <- mutate(all_invs3,
                   Sp2 = forcats::fct_recode(Sp, `(a)~italic("Achromobacter")~sp.` = 'A', `(b)~italic("Ochrobactrum")~sp.` = 'O', `(c)~italic("Pseudomonas")~sp.` = 'P', `(d)~italic("Stenotrophomonas")~sp.` = 'S', `(e)~italic("Variovorax")~sp.` = 'V'))

all_invs_m <- group_by(all_invs2, Sp) %>%
  summarise(sig_t = sum(p.value < 0.05))

#Figure 3

ggplot() +
  theme_bw() +
  geom_point(data = all_invs, aes(x = rel.fit, y = res.sp, group = interaction(treat,invasion), col = invasion), position = position_dodge(0.75), alpha = 0.5, size = 0.5) +
  geom_point(data = all_invs3, aes(x = mean, y = res.sp, group = interaction(treat,invasion), col = invasion), size = 2, position = position_dodge(0.75)) +
  facet_wrap(~Sp2, scales = "free_y", labeller = label_parsed) + 
  geom_text(data = all_invs2, aes(x = mean, y = res.sp, label = sig), nudge_x = 2, size = 5, nudge_y = -0.25) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 12, hjust = 0), axis.text = element_text(size = 11, colour = 'black'), axis.title = element_text(size = 12, colour = 'black'), legend.position = "bottom", legend.text = element_text(size = 11)) + 
  geom_vline(xintercept=1) +
  palettetown::scale_color_poke(pokemon = "blastoise", spread = 4, labels = c("Co-invasion", "Single invasion")) +
  ylab("Resident species") +
  xlab("Relative invader growth rate (m)") +
  labs(col = 'Invasion type') 

#Figure 3 - panel f

ggplot(all_invs_m, aes(x = Sp, y = sig_t, fill = Sp)) +
  theme_bw()+
  geom_col() +
  scale_fill_manual(values = color_scheme) +
  scale_color_manual(values = color_scheme) +
  theme(axis.text = element_text(size = 11, colour = 'black'), axis.title = element_text(size = 12, colour = 'black'), legend.position = "none")+
  ylab('Number of significant comparisons') +
  xlab('Focal invading species') +
  scale_y_continuous(breaks = c(1:14)) +
  labs(title = '( f )')

d_mods_summary3 <- d_mods_summary[,c(1,4,13,14)]

d_mods_summary3$sig[d_mods_summary3$treat == "A_PSV"] <- "NS-s"
d_mods_summary3$sig[d_mods_summary3$treat == "A_OPV"] <- "NS-s"
d_mods_summary3$sig[d_mods_summary3$treat == "A_OPSV"] <- "NS-b" #exclude this one as duplicate
d_mods_summary3$sig[d_mods_summary3$treat == "S_AOV"] <- "NS-b"

d_mods_summary3 <- d_mods_summary3[-c(1:5),]
d_mods_summary3 <- na.omit(d_mods_summary3)

#add in one-sample t-test results from single invasions

inv.coms$Sp[inv.coms$Sp == "VG"] <- "V"
inv.coms$Sp[inv.coms$Sp == "SR"] <- "S"
inv.coms$Sp[inv.coms$Sp == "AA"] <- "A"
inv.coms$Sp[inv.coms$Sp == "PC"] <- "P"
inv.coms$Sp[inv.coms$Sp == "OD"] <- "O"

inv.coms$res.sp <- gsub('VG',"V",inv.coms$res.sp)
inv.coms$res.sp <- gsub('SR',"S",inv.coms$res.sp)
inv.coms$res.sp <- gsub('PC',"P",inv.coms$res.sp)
inv.coms$res.sp <- gsub('AA',"A",inv.coms$res.sp)
inv.coms$res.sp <- gsub('OD',"O",inv.coms$res.sp)

inv.coms$res.sp <- gsub('[.]','',inv.coms$res.sp)

sing_invs <- read.csv("rel_fitness_sing_inv.csv", header = T)
names(sing_invs)[3] <- "Sp"
names(sing_invs)[4] <- "res.sp"

sing_invs$Sp[sing_invs$Sp == "VG"] <- "V"
sing_invs$Sp[sing_invs$Sp == "SR"] <- "S"
sing_invs$Sp[sing_invs$Sp == "AA"] <- "A"
sing_invs$Sp[sing_invs$Sp == "PC"] <- "P"
sing_invs$Sp[sing_invs$Sp == "OD"] <- "O"

sing_invs$res.sp <- gsub('VG',"V",sing_invs$res.sp)
sing_invs$res.sp <- gsub('SR',"S",sing_invs$res.sp)
sing_invs$res.sp <- gsub('PC',"P",sing_invs$res.sp)
sing_invs$res.sp <- gsub('AA',"A",sing_invs$res.sp)
sing_invs$res.sp <- gsub('OD',"O",sing_invs$res.sp)

sing_invs$res.sp <- gsub('[.]','',sing_invs$res.sp)

sing_invs$res.sp <- strSort(sing_invs$res.sp)
sing_invs$invasion <- "single"
names(sing_invs)[12] <- "rel.fit"

sing_invs$treat <- interaction(sing_invs$Sp, sing_invs$res.sp, sep = "_")

d_mods <- sing_invs %>%
  nest(-treat) %>%
  mutate(., mod = map(data, ~ t.test(.x$rel.fit, mu = 1)))

d_mods_sing <- d_mods %>%
  mutate(params = map(mod, broom::tidy))

d_mods_sing <- unnest(d_mods_sing, params) %>%
  mutate(padjust = p.adjust(p.value, "fdr"))

d_mods_sing$padjust2 <- round(d_mods_sing$padjust, 5) 

d_mods_sing <- d_mods_sing[,c(1,4,13)]
d_mods_sing$invasion <- "single"
d_mods_summary3$invasion <- "multi"

d_all <- merge(d_mods_sing, d_mods_summary3, by = "treat")
#first set of columns are single invasion estimates
#second set are co-invasion estimates

d_all <- d_all[,c(-4,-7,-8)]

d_all[c('Sp', 'Residents')] <- str_split_fixed(d_all$treat, '_', 2)
d_all <- d_all[,-1]
d_all <- d_all[,c(5,6,1:4)]

d_all$estimate.x <- round(d_all$estimate.x, 3)
d_all$estimate.y <- round(d_all$estimate.y, 3)
d_all$padjust2.x <- round(d_all$padjust2.x, 3)
d_all$padjust2.y <- round(d_all$padjust2.y, 3)

d_all$padjust2.x[d_all$padjust2.x == 0] <- '<0.001'
d_all$padjust2.y[d_all$padjust2.y == 0] <- '<0.001'

library(flextable)

#Table 1

men_flex <- flextable(d_all) %>%
  set_header_labels(estimate.x = "Single invasion\nMean", padjust2.x = "\np-value", estimate.y = "Co-invasion\nMean", padjust2.y = "\np-value") %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 16, part = 'all') %>%
  vline(j = c(2,4,6), border = officer::fp_border(color="black")) %>%
  autofit() 
