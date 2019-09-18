# read in raw flux data
stem = read.csv('../raw_data/stem_raw.csv')
root = read.csv('../raw_data/root_raw.csv')
leaf = read.csv('../raw_data/leaf_raw.csv')

stem_sub = subset(stem, Stem_Rarea <= 0 & Stem_Rmass <= 0)
root_sub = subset(root, Root_Rarea <= 0 & Root_Rmass <= 0)
individuals_stem = levels(stem_sub$individual)
individuals_root = levels(root_sub$individual)
leaf_sub = subset(leaf, Rarea >= 0 & Rmass >= 0)
individuals_leaf = levels(leaf_sub$individual)

# fit each curve and put in dataframe
stem_fits = c()
for(i in 1:length(individuals_stem)){
	
	fit = lm(log(-Stem_Rarea) ~ Stem_Tmeas + Stem_Tmeas2, data = subset(stem, individual == individuals_stem[i]))
	summary = summary(fit)
	
	fit_mass = lm(log(-Stem_Rmass) ~ Stem_Tmeas + Stem_Tmeas2, data = subset(stem, individual == individuals_stem[i]))
	summary_mass = summary(fit_mass)
	
	data = c(as.character(subset(stem, individual == individuals_stem[i])$genus[1]), as.character(subset(stem, individual == individuals_stem[i])$species[1]), subset(stem, individual == individuals_stem[i])$GrowthT[1], subset(stem, individual == individuals_stem[i])$Rep[1], summary$coefficients[1], summary$coefficients[2], summary$coefficients[3], sqrt(mean(summary$residuals^2)), summary$adj.r.squared, summary_mass$coefficients[1], summary_mass$coefficients[2], summary_mass$coefficients[3], sqrt(mean(summary_mass$residuals^2)), summary_mass$adj.r.squared)
	stem_fits = rbind(stem_fits, data)
	
}

colnames(stem_fits) = c('genus', 'species', 'GrowthT', 'Rep', 'a_stem_area', 'b_stem_area', 'c_stem_area', 'RMSE_stem_area', 'r2_stem_area', 'a_stem_mass', 'b_stem_mass', 'c_stem_mass', 'RMSE_stem_mass', 'r2_stem_mass')

root_fits = c()
for(i in 1:length(individuals_root)){
	
	fit = lm(log(-Root_Rarea) ~ Root_Tmeas + Root_Tmeas2, data = subset(root, individual == individuals_root[i]))
	summary = summary(fit)
	
	fit_mass = lm(log(-Root_Rmass) ~ Root_Tmeas + Root_Tmeas2, data = subset(root, individual == individuals_root[i]))
	summary_mass = summary(fit_mass)
	
	data = c(as.character(subset(root, individual == individuals_root[i])$genus[1]), as.character(subset(root, individual == individuals_root[i])$species[1]), subset(root, individual == individuals_root[i])$GrowthT[1], subset(root, individual == individuals_root[i])$Rep[1], summary$coefficients[1], summary$coefficients[2], summary$coefficients[3], sqrt(mean(summary$residuals^2)), summary$adj.r.squared, summary_mass$coefficients[1], summary_mass$coefficients[2], summary_mass$coefficients[3], sqrt(mean(summary_mass$residuals^2)), summary_mass$adj.r.squared)
	root_fits = rbind(root_fits, data)
	
}

colnames(root_fits) = c('genus', 'species', 'GrowthT', 'Rep', 'a_root_area', 'b_root_area', 'c_root_area', 'RMSE_root_area', 'r2_root_area', 'a_root_mass', 'b_root_mass', 'c_root_mass', 'RMSE_root_mass', 'r2_root_mass')

leaf_fits = c()
for(i in 1:length(individuals_leaf)){
	
	fit = lm(log(Rarea) ~ LeafT + LeafT2, data = subset(leaf, individual == individuals_leaf[i]))
	summary = summary(fit)
	
	fit_mass = lm(log(Rmass) ~ LeafT + LeafT2, data = subset(leaf, individual == individuals_leaf[i]))
	summary_mass = summary(fit_mass)
	
	data = c(as.character(subset(leaf, individual == individuals_leaf[i])$Genus[1]), as.character(subset(leaf, individual == individuals_leaf[i])$Species[1]), subset(leaf, individual == individuals_leaf[i])$GrowthT[1], subset(leaf, individual == individuals_leaf[i])$Rep[1], summary$coefficients[1], summary$coefficients[2], summary$coefficients[3], sqrt(mean(summary$residuals^2)), summary$adj.r.squared, summary_mass$coefficients[1], summary_mass$coefficients[2], summary_mass$coefficients[3], sqrt(mean(summary_mass$residuals^2)), summary_mass$adj.r.squared)
	leaf_fits = rbind(leaf_fits, data)
	
}

colnames(leaf_fits) = c('genus', 'species', 'GrowthT', 'Rep', 'a_leaf_area', 'b_leaf_area', 'c_leaf_area', 'RMSE_leaf_area', 'r2_leaf_area', 'a_leaf_mass', 'b_leaf_mass', 'c_leaf_mass', 'RMSE_leaf_mass', 'r2_leaf_mass')

# write new data frames

#write.csv(stem_fits, '../tresp_curve_fits/stem_fits_v3.csv')
#write.csv(root_fits, '../tresp_curve_fits/root_fits_v3.csv')
#write.csv(leaf_fits, '../tresp_curve_fits/leaf_fits_v3.csv')




