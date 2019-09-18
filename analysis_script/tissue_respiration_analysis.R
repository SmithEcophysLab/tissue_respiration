#############################
#####################
# load libraries
#####################
#############################

library(lme4)
library(car)
library(emmeans)
library(tidyverse)
library(dplyr)
library(gtable)
library(grid)
library(simr)

#############################
#####################
# load respiration data
#####################
#############################

# load up fits
stem_fits = read.csv('../tresp_curve_fits/stem_fits_v3.csv')
root_fits = read.csv('../tresp_curve_fits/root_fits_v3.csv')
leaf_fits = read.csv('../tresp_curve_fits/leaf_fits_v3.csv')

# combine fits
all_analysis1 = full_join(stem_fits, root_fits, by = c('genus', 'species', 'GrowthT', 'Rep'))
all_analysis = full_join(all_analysis1, leaf_fits, by = c('genus', 'species', 'GrowthT', 'Rep'))

#############################
#####################
# clean respiration data
#####################
#############################

stem_mass = all_analysis[ , c(2:5, 11:15)]
root_mass = all_analysis[ , c(2:5, 22:26)]
leaf_mass = all_analysis[ , c(2:5, 33:37)]

stem_mass$tissue = 'stem'
root_mass$tissue = 'root'
leaf_mass$tissue = 'leaf'

colnames(stem_mass) = c("genus", "species", "GrowthT", "Rep", "a_mass", "b_mass", "c_mass", "RMSE_mass", "r2_mass", "tissue")
colnames(leaf_mass) = c("genus", "species", "GrowthT", "Rep", "a_mass", "b_mass", "c_mass", "RMSE_mass", "r2_mass", "tissue")
colnames(root_mass) = c("genus", "species", "GrowthT", "Rep", "a_mass", "b_mass", "c_mass", "RMSE_mass", "r2_mass", "tissue")

all_mass_sub = rbind(stem_mass, root_mass, leaf_mass)

all_mass_sub$individual = paste(all_mass$genus, all_mass$species, all_mass$GrowthT, all_mass$Rep, sep = '_')

# add in plant type characteristics
all_mass_sub$lifespan = 'annual'
all_mass_sub$lifespan[all_mass_sub$genus == 'Pinus' | all_mass_sub$genus == 'Betula'] = 'perennial'

# fit statistics
# mean(subset(all_mass_sub, tissue == 'leaf')$r2_mass, na.rm = T)
# mean(subset(all_mass_sub, tissue == 'stem')$r2_mass, na.rm = T)
# mean(subset(all_mass_sub, tissue == 'root')$r2_mass, na.rm = T)
# mean(subset(all_mass_sub, tissue == 'leaf')$RMSE_mass, na.rm = T)
# mean(subset(all_mass_sub, tissue == 'stem')$RMSE_mass, na.rm = T)
# mean(subset(all_mass_sub, tissue == 'root')$RMSE_mass, na.rm = T)

# count individuals per species per GrowthT
# all_mass_sub$GrowthTfac = as.factor(all_mass_sub$GrowthT)
# all_mass_sub_group = group_by(all_mass_sub, GrowthTfac, genus, species, tissue)
# all_mass_sub_count = summarise(all_mass_sub_group, n = n())
# write.csv(subset(all_mass_sub_count, tissue == 'leaf'), '/Users/nicksmith/Documents/Research/Thesis/Growth chamber/Stem_Root/Analysis/all_mass_sub_count.csv') # only leaf subset because leaf, stem, root are all the same

# separate by photosynthetic and non-photosynthetic tissue
all_mass_sub$ps_tis = 'no'
all_mass_sub$ps_tis[all_mass_sub$species == 'max' & all_mass_sub$tissue == 'stem'] = 'yes'
all_mass_sub$ps_tis[all_mass_sub$species == 'sativa' & all_mass_sub$tissue == 'stem'] = 'yes'
all_mass_sub$ps_tis[all_mass_sub$tissue == 'leaf'] = 'yes'

all_mass_sub$tissue_ps = paste(all_mass_sub$tissue, all_mass_sub$ps_tis, sep = '_')

#############################
#####################
# add in photosynthetic data
#####################
#############################

gc_data=read.csv('../raw_data/gc_data_merged.csv')
gc_data$individual = paste(gc_data$Genus, gc_data$Species, gc_data$GrowthT, gc_data$Rep, sep = '_')
gc_data_rm = gc_data[-c(189:199, 334:362), ] #remove all non-measured Zea mays and Pinus sylvestris

# but how do we make sure this gets populated for all leaf, stem, and root?
all_mass_sub_ps_leaf = left_join(subset(all_mass_sub, tissue == 'leaf'), gc_data_rm, by = c('individual'))
all_mass_sub_ps_stem = left_join(subset(all_mass_sub, tissue == 'stem'), gc_data_rm, by = c('individual'))
all_mass_sub_ps_root = left_join(subset(all_mass_sub, tissue == 'root'), gc_data_rm, by = c('individual'))

all_mass_sub_ps = rbind(all_mass_sub_ps_leaf, all_mass_sub_ps_stem, all_mass_sub_ps_root)

all_mass_sub_ps$vcmax_acc_area = exp(all_mass_sub_ps$a_Vcmax + all_mass_sub_ps$b_Vcmax * all_mass_sub_ps$GrowthT.y + all_mass_sub_ps$c_Vcmax * (all_mass_sub_ps$GrowthT.y ^2))
all_mass_sub_ps$vcmax_acc_mass = all_mass_sub_ps$vcmax_acc_area * ((all_mass_sub_ps$Larea_P * 0.0001) / all_mass_sub_ps$Lmass_P)

all_mass_sub_ps$rd_acc_mass = exp(all_mass_sub_ps$a_mass + all_mass_sub_ps$b_mass * all_mass_sub_ps$GrowthT.x + all_mass_sub_ps$c_mass * (all_mass_sub_ps$GrowthT.x ^2))

all_mass_sub_ps$rdvcmax_acc_mass = all_mass_sub_ps$rd_acc_mass / all_mass_sub_ps$vcmax_acc_mass

#############################
#####################
# run mixed models
#####################
#############################
#MM analyses
a_lm = lmer(a_mass ~ GrowthT * tissue_ps + species + (1|individual) , data = all_mass_sub)
summary(a_lm)
Anova(a_lm) # 
cld(emmeans(a_lm, ~species))
a_lm_lsm = emtrends(a_lm, ~ tissue_ps, var = 'GrowthT')
test(a_lm_lsm) # 
cld(a_lm_lsm) # 
contrast(a_lm_lsm, list(c1 = c(-1, 1, 1, -1))) # test non-photosynthetic versus photosynthetic slopes
plot(residuals(a_lm)~fitted(a_lm))

b_lm = lmer(b_mass ~ GrowthT * tissue_ps + species + (1|individual) , data = all_mass_sub)
summary(b_lm)
Anova(b_lm) # T x tissue
b_lm_lsm = emtrends(b_lm, ~ tissue_ps, var = 'GrowthT')
test(b_lm_lsm) # 
cld(b_lm_lsm) # 
contrast(b_lm_lsm, list(c1 = c(-1, 1, 1, -1))) # test non-photosynthetic versus photosynthetic slopes
plot(residuals(b_lm)~fitted(b_lm))

c_lm = lmer(c_mass ~ GrowthT * tissue_ps + species + (1|individual) , data = all_mass_sub)
summary(c_lm)
Anova(c_lm) # tissue and T x tissue
c_lm_lsm = emtrends(c_lm, ~ tissue_ps, var = 'GrowthT')
test(c_lm_lsm) # 
cld(c_lm_lsm) # 
contrast(c_lm_lsm, list(c1 = c(-1, 1, 1, -1))) # test non-photosynthetic versus photosynthetic slopes
plot(residuals(c_lm)~fitted(c_lm))

rdvcmax_lm = lmer(log(rdvcmax_acc_mass) ~ GrowthT.x * tissue_ps + species + 
                    (1|individual), 
                  data = subset(all_mass_sub_ps, rdvcmax_acc_mass < 1 & rdvcmax_acc_mass > 0))
Anova(rdvcmax_lm)
rdvcmax_lm_lsm = emtrends(rdvcmax_lm, ~ tissue_ps, var = 'GrowthT.x')
test(rdvcmax_lm_lsm) # 
cld(rdvcmax_lm_lsm) # 
contrast(rdvcmax_lm_lsm, list(c1 = c(-1, 1, 1, -1))) # test non-photosynthetic versus photosynthetic slopes
plot(residuals(rdvcmax_lm)~fitted(rdvcmax_lm))
cld(emmeans(rdvcmax_lm, ~ species))

#############################
#####################
# calculate Acclimhomeo values
#####################
#############################
# leaf
leaf_tatl_15 = exp(summary(lsmeans(a_lm, ~tissue_ps, at = list(GrowthT = 15)))[1, 2] + summary(lsmeans(b_lm, ~tissue_ps, at = list(GrowthT = 15)))[1, 2] * 15 + summary(lsmeans(c_lm, ~tissue_ps, at = list(GrowthT = 15)))[1, 2] * 15 * 15)
leaf_tatl_20 = exp(summary(lsmeans(a_lm, ~tissue_ps, at = list(GrowthT = 20)))[1, 2] + summary(lsmeans(b_lm, ~tissue_ps, at = list(GrowthT = 20)))[1, 2] * 20 + summary(lsmeans(c_lm, ~tissue_ps, at = list(GrowthT = 20)))[1, 2] * 20 * 20)
leaf_tatl_25 = exp(summary(lsmeans(a_lm, ~tissue_ps, at = list(GrowthT = 25)))[1, 2] + summary(lsmeans(b_lm, ~tissue_ps, at = list(GrowthT = 25)))[1, 2] * 25 + summary(lsmeans(c_lm, ~tissue_ps, at = list(GrowthT = 25)))[1, 2] * 25 * 25)
leaf_tatl_30 = exp(summary(lsmeans(a_lm, ~tissue_ps, at = list(GrowthT = 30)))[1, 2] + summary(lsmeans(b_lm, ~tissue_ps, at = list(GrowthT = 30)))[1, 2] * 30 + summary(lsmeans(c_lm, ~tissue_ps, at = list(GrowthT = 30)))[1, 2] * 30 * 30)
leaf_tatl_35 = exp(summary(lsmeans(a_lm, ~tissue_ps, at = list(GrowthT = 35)))[1, 2] + summary(lsmeans(b_lm, ~tissue_ps, at = list(GrowthT = 35)))[1, 2] * 35 + summary(lsmeans(c_lm, ~tissue_ps, at = list(GrowthT = 35)))[1, 2] * 35 * 35)

leaf_homeo = c(leaf_tatl_15 / leaf_tatl_25, leaf_tatl_20 / leaf_tatl_25, leaf_tatl_25 / leaf_tatl_30, leaf_tatl_25 / leaf_tatl_35)

root_tatl_15 = exp(summary(lsmeans(a_lm, ~tissue_ps, at = list(GrowthT = 15)))[2, 2] + summary(lsmeans(b_lm, ~tissue_ps, at = list(GrowthT = 15)))[2, 2] * 15 + summary(lsmeans(c_lm, ~tissue_ps, at = list(GrowthT = 15)))[2, 2] * 15 * 15)
root_tatl_20 = exp(summary(lsmeans(a_lm, ~tissue_ps, at = list(GrowthT = 20)))[2, 2] + summary(lsmeans(b_lm, ~tissue_ps, at = list(GrowthT = 20)))[2, 2] * 20 + summary(lsmeans(c_lm, ~tissue_ps, at = list(GrowthT = 20)))[2, 2] * 20 * 20)
root_tatl_25 = exp(summary(lsmeans(a_lm, ~tissue_ps, at = list(GrowthT = 25)))[2, 2] + summary(lsmeans(b_lm, ~tissue_ps, at = list(GrowthT = 25)))[2, 2] * 25 + summary(lsmeans(c_lm, ~tissue_ps, at = list(GrowthT = 25)))[2, 2] * 25 * 25)
root_tatl_30 = exp(summary(lsmeans(a_lm, ~tissue_ps, at = list(GrowthT = 30)))[2, 2] + summary(lsmeans(b_lm, ~tissue_ps, at = list(GrowthT = 30)))[2, 2] * 30 + summary(lsmeans(c_lm, ~tissue_ps, at = list(GrowthT = 30)))[2, 2] * 30 * 30)
root_tatl_35 = exp(summary(lsmeans(a_lm, ~tissue_ps, at = list(GrowthT = 35)))[2, 2] + summary(lsmeans(b_lm, ~tissue_ps, at = list(GrowthT = 35)))[2, 2] * 35 + summary(lsmeans(c_lm, ~tissue_ps, at = list(GrowthT = 35)))[2, 2] * 35 * 35)

root_homeo = c(root_tatl_15 / root_tatl_25, root_tatl_20 / root_tatl_25, root_tatl_25 / root_tatl_30, root_tatl_25 / root_tatl_35)

stem_no_tatl_15 = exp(summary(lsmeans(a_lm, ~tissue_ps, at = list(GrowthT = 15)))[3, 2] + summary(lsmeans(b_lm, ~tissue_ps, at = list(GrowthT = 15)))[3, 2] * 15 + summary(lsmeans(c_lm, ~tissue_ps, at = list(GrowthT = 15)))[3, 2] * 15 * 15)
stem_no_tatl_20 = exp(summary(lsmeans(a_lm, ~tissue_ps, at = list(GrowthT = 20)))[3, 2] + summary(lsmeans(b_lm, ~tissue_ps, at = list(GrowthT = 20)))[3, 2] * 20 + summary(lsmeans(c_lm, ~tissue_ps, at = list(GrowthT = 20)))[3, 2] * 20 * 20)
stem_no_tatl_25 = exp(summary(lsmeans(a_lm, ~tissue_ps, at = list(GrowthT = 25)))[3, 2] + summary(lsmeans(b_lm, ~tissue_ps, at = list(GrowthT = 25)))[3, 2] * 25 + summary(lsmeans(c_lm, ~tissue_ps, at = list(GrowthT = 25)))[3, 2] * 25 * 25)
stem_no_tatl_30 = exp(summary(lsmeans(a_lm, ~tissue_ps, at = list(GrowthT = 30)))[3, 2] + summary(lsmeans(b_lm, ~tissue_ps, at = list(GrowthT = 30)))[3, 2] * 30 + summary(lsmeans(c_lm, ~tissue_ps, at = list(GrowthT = 30)))[3, 2] * 30 * 30)
stem_no_tatl_35 = exp(summary(lsmeans(a_lm, ~tissue_ps, at = list(GrowthT = 35)))[3, 2] + summary(lsmeans(b_lm, ~tissue_ps, at = list(GrowthT = 35)))[3, 2] * 35 + summary(lsmeans(c_lm, ~tissue_ps, at = list(GrowthT = 35)))[3, 2] * 35 * 35)

stem_no_homeo = c(stem_no_tatl_15 / stem_no_tatl_25, stem_no_tatl_20 / stem_no_tatl_25, stem_no_tatl_25 / stem_no_tatl_30, stem_no_tatl_25 / stem_no_tatl_35)

stem_yes_tatl_15 = exp(summary(lsmeans(a_lm, ~tissue_ps, at = list(GrowthT = 15)))[4, 2] + summary(lsmeans(b_lm, ~tissue_ps, at = list(GrowthT = 15)))[4, 2] * 15 + summary(lsmeans(c_lm, ~tissue_ps, at = list(GrowthT = 15)))[4, 2] * 15 * 15)
stem_yes_tatl_20 = exp(summary(lsmeans(a_lm, ~tissue_ps, at = list(GrowthT = 20)))[4, 2] + summary(lsmeans(b_lm, ~tissue_ps, at = list(GrowthT = 20)))[4, 2] * 20 + summary(lsmeans(c_lm, ~tissue_ps, at = list(GrowthT = 20)))[4, 2] * 20 * 20)
stem_yes_tatl_25 = exp(summary(lsmeans(a_lm, ~tissue_ps, at = list(GrowthT = 25)))[4, 2] + summary(lsmeans(b_lm, ~tissue_ps, at = list(GrowthT = 25)))[4, 2] * 25 + summary(lsmeans(c_lm, ~tissue_ps, at = list(GrowthT = 25)))[4, 2] * 25 * 25)
stem_yes_tatl_30 = exp(summary(lsmeans(a_lm, ~tissue_ps, at = list(GrowthT = 30)))[4, 2] + summary(lsmeans(b_lm, ~tissue_ps, at = list(GrowthT = 30)))[4, 2] * 30 + summary(lsmeans(c_lm, ~tissue_ps, at = list(GrowthT = 30)))[4, 2] * 30 * 30)
stem_yes_tatl_35 = exp(summary(lsmeans(a_lm, ~tissue_ps, at = list(GrowthT = 35)))[4, 2] + summary(lsmeans(b_lm, ~tissue_ps, at = list(GrowthT = 35)))[4, 2] * 35 + summary(lsmeans(c_lm, ~tissue_ps, at = list(GrowthT = 35)))[4, 2] * 35 * 35)

stem_yes_homeo = c(stem_yes_tatl_15 / stem_yes_tatl_25, stem_yes_tatl_20 / stem_yes_tatl_25, stem_yes_tatl_25 / stem_yes_tatl_30, stem_yes_tatl_25 / stem_yes_tatl_35)

homeo = cbind(c(15, 20, 30, 35), leaf_homeo, stem_yes_homeo, stem_no_homeo, root_homeo)
homeo_mean = c(0, mean(leaf_homeo), mean(stem_yes_homeo), mean(stem_no_homeo), mean(root_homeo))
homeo = rbind(homeo, homeo_mean)
colnames(homeo) = c('Ta', 'leaf_homeo', 'stem_yes_homeo', 'stem_no_homeo', 'root_homeo')

