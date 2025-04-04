# AP. Woodward, 2025, University of Georgia.
# Implements some descriptive analyses of various outcomes from canine bladder cancer cases.
# To quantify some sampling errors for proportions the Wilson confidence interval is used, following the recommendation of (https://doi.org/10.1214/ss/1009213286).
# A simple hierarchical model is used to probe for differences in metastasis risk by variant, using weak prior information, via package 'brms' (https://doi.org/10.18637/jss.v080.i01).
# A hierarchical survival model to quantify the effect of variant on survival time was computationally unstable and not trustworthy.
# Survival is instead described visually using the Kaplan-Meier method, across variants. 

# Load the required packages.
library(brms)
library(ggplot2)
library(matrixStats)
library(survival)
library(ggpubr)
library(reshape2)
library(binom)

# Conduct some data cleaning (can't have survival times of zero, so we will set those to an arbitrary small value) 
bladder_data_cleaned$SURV[bladder_data_cleaned$SURV == 0] <- 0.01

# Conduct a hierarchical logistic regression model for the presence of metastasis (METASTASIS = 1), with risk varying by VARIANT, using weak prior information.
bladder_data_metastasis <- bladder_data_cleaned[!(is.na(bladder_data_cleaned$METASTASIS)),]
bladder_data_metastasis$VARIANT <- factor(bladder_data_metastasis$VARIANT)
bladder_mod_metastasis <- brm(METASTASIS ~ 1 + (1|VARIANT), bladder_data_metastasis, family = bernoulli(link = 'logit'), prior = set_prior('normal(0,1)', class = 'sd') + set_prior('normal(0,2)', class = 'Intercept'))
bladder_mod_metastasis_pred <- cbind(data.frame(VARIANT = levels(bladder_data_metastasis$VARIANT)), as.data.frame(fitted(bladder_mod_metastasis, newdata = data.frame(VARIANT = levels(bladder_data_metastasis$VARIANT)), robust = TRUE, probs = c(0.05,0.95))))
ggplot(data = bladder_mod_metastasis_pred, aes(x = VARIANT, y = Estimate)) + geom_point() +  geom_errorbar(aes(ymin = Q5, ymax = Q95)) + theme(axis.text.x = element_text(angle = 90)) + ylim(c(0,1))

# Generate some postprocessing visualization from the metastasis model. The barplot is presented on the scale of the observed counts (to clearly demonstrate the sample size).
metastasis_freq <- data.frame(VARIANT = row.names(table(bladder_data_metastasis$VARIANT, bladder_data_metastasis$METASTASIS)), negative = (table(bladder_data_metastasis$VARIANT, bladder_data_metastasis$METASTASIS))[,1], positive = (table(bladder_data_metastasis$VARIANT, bladder_data_metastasis$METASTASIS))[,2], total = rowSums(table(bladder_data_metastasis$VARIANT, bladder_data_metastasis$METASTASIS)))
metastasis_freq <- melt(metastasis_freq, id.vars = 'VARIANT')
bladder_mod_metastasis_predcount <- data.frame(VARIANT = levels(bladder_data_metastasis$VARIANT))
bladder_mod_metastasis_predcount <- cbind(bladder_mod_metastasis_predcount, as.data.frame(bladder_mod_metastasis_pred[,2:5])*matrix((metastasis_freq$value[((metastasis_freq$variable == 'total') & !(metastasis_freq$VARIANT == 'Large nested UC'))]), ncol = 4, nrow = 5))
bladder_mod_metastasis_predprop  <- data.frame(VARIANT = levels(bladder_data_metastasis$VARIANT))
bladder_mod_metastasis_predprop <- cbind(bladder_mod_metastasis_predprop, as.data.frame(bladder_mod_metastasis_pred[,2:5]))
metastasis_prop_plot  <- ggplot(data = metastasis_freq[(!(metastasis_freq$variable == 'total')),], aes (x = VARIANT, y = value, fill = variable)) + geom_bar(stat = 'identity', position = 'fill', alpha = 0.6) + theme_bw() + coord_flip() + theme(legend.title = element_blank(), panel.grid.minor = element_blank()) + xlab('') + ylab('Metastasis (proportion)') + geom_errorbar(data = bladder_mod_metastasis_predprop, aes(x = VARIANT, ymin = Q5, ymax = Q95), inherit.aes = FALSE, alpha = 0.7, width = 0.5) + geom_point(data = bladder_mod_metastasis_predprop, aes(x = VARIANT, y = Estimate), inherit.aes = FALSE, alpha = 0.7) + scale_fill_grey(labels = c('Metastasis-negative','Metastasis-positive')) + scale_y_continuous(breaks = seq(0,1,0.5))
metastasis_count_plot <- ggplot(data = metastasis_freq[(!(metastasis_freq$variable == 'total')),], aes (x = VARIANT, y = value, fill = variable)) + geom_bar(stat = 'identity', alpha = 0.6) + theme_bw() + coord_flip() + theme(legend.title = element_blank(), panel.grid.minor = element_blank()) + xlab('') + ylab('Metastasis (count)') + geom_errorbar(data = bladder_mod_metastasis_predcount, aes(x = VARIANT, ymin = Q5, ymax = Q95), inherit.aes = FALSE, alpha = 0.7, width = 0.5) + geom_point(data = bladder_mod_metastasis_predcount, aes(x = VARIANT, y = Estimate), inherit.aes = FALSE, alpha = 0.7) + scale_fill_grey(labels = c('Metastasis-negative','Metastasis-positive')) + scale_y_continuous(breaks = seq(0,10,2))
metastasis_plots <- ggarrange(metastasis_count_plot, metastasis_prop_plot, ncol = 1, common.legend = TRUE, legend = 'bottom')
ggsave(metastasis_plots, file = 'metastasis_plots.svg', units = 'mm', height = 200, width = 120)

# Implement and visualize a Kaplan-Meier analysis of the survival time.
bladder_KM_mod <- survfit(Surv(SURV, DEATH) ~ 1, data = bladder_data_cleaned)
bladder_KM_plot <- ggsurvplot(bladder_KM_mod, data = bladder_data_cleaned, ggtheme = theme_bw(), palette = 'black', legend = 'none', alpha = 0.5, xlab = 'Time (months)')
bladder_surv_combined <- ggarrange(bladder_KM_plot$plot, bladder_surv_points, ncol = 1)
ggsave(bladder_surv_combined, filename = 'bladder_surv_combined.svg', units = 'mm', height = 150, width = 150)

# Implement a hierarchical model for survival time, using the lognormal survival model.
#     To have any chance of achieving a reasonable fit this model must allow variability of 'sigma' by 'VARIANT', but this is too unstable with the small dataset and should not be used.
bladder_mod_survival <- brm(bf(SURV|cens(CENS) ~ 1 + (1|VARIANT), sigma ~ (1|VARIANT)), bladder_data_cleaned, family = lognormal(), prior = set_prior('normal(0,1)', class = 'sd') + set_prior('normal(1.8,0.5)', class = 'Intercept', dpar = 'mu') + set_prior('normal(0,1)', class = 'Intercept', dpar = 'sigma'), control = list(adapt_delta = 0.99))
bladder_surv_median <- as.data.frame(exp(((coef(bladder_mod_survival, robust = TRUE, probs = c(0.05,0.25,0.75,0.95)))$VARIANT)[,,1]))
bladder_surv_median$VARIANT <- row.names(bladder_surv_median)
ggplot() + geom_errorbar(data = bladder_surv_median, aes(x = VARIANT, ymin = Q5, ymax = Q95), inherit.aes = FALSE, alpha = 0.6) + geom_point(data = bladder_data_cleaned, aes(x = VARIANT, y = SURV, shape = factor(CENS)), fill = 'white', alpha = 0.7) + scale_shape_manual(values = c(19,21)) + scale_y_log10() + theme(legend.position = 'none')  + theme_bw() + coord_cartesian(ylim = c(0.01,100)) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), legend.position = 'none') + theme(axis.text.x = element_text(angle = 90)) + coord_flip() + ylab('Survival time (months)') + xlab('') 

# Generate a 95% confidence interval for the proportion 4/31 using method = 'Wilson'.
binom.wilson(4,31)

# Generate a 95% confidence interval for the proportion 8/31 using method = 'Wilson'.
binom.wilson(8,31)

# Generate a 95% confidence interval for the proportion 2/31 using method = 'Wilson'.
binom.wilson(2,31)

# Generate a 95% confidence interval for the proportion 2/31 using method = 'Wilson'.
binom.wilson(1,31)

# Export the working data.
write.csv(bladder_data_cleaned, file = 'bladder_data_cleaned.csv', row.names = FALSE)
