# This code is designed to be ran line-by-line in an R interpreter rather than as a script

library(reshape2)

setwd(<insert working directory here>)

dat <- read.csv('normalized_distances_fr10.csv')

dat_melted <- reshape2::melt(dat, id.vars = "Gene")

head(dat_melted)

results <- lm(value~relevel(factor(Gene),ref="ApoA-II"), data=dat_melted)

summary(results)

# Call:
# lm(formula = value ~ relevel(factor(Gene), ref = "ApoA-II"),
#     data = dat_melted)
# 
# Residuals:
#      Min       1Q   Median       3Q      Max
# -0.43305 -0.07292  0.01587  0.05663  0.49697
# 
# Coefficients:
#                                                Estimate Std. Error t value
# (Intercept)                                     1.44387    0.03103  46.535
# relevel(factor(Gene), ref = "ApoA-II")ApoA-I   -0.39374    0.04180  -9.420
# relevel(factor(Gene), ref = "ApoA-II")ApoA-IV  -0.41485    0.04296  -9.658
# relevel(factor(Gene), ref = "ApoA-II")ApoA-V   -0.39384    0.04496  -8.759
# relevel(factor(Gene), ref = "ApoA-II")ApoC-I   -0.77368    0.04700 -16.463
# relevel(factor(Gene), ref = "ApoA-II")ApoC-II  -0.30668    0.04496  -6.821
# relevel(factor(Gene), ref = "ApoA-II")ApoC-III -0.30668    0.04496  -6.821
# relevel(factor(Gene), ref = "ApoA-II")АроЕ     -0.51810    0.04254 -12.178
#                                                Pr(>|t|)
# (Intercept)                                     < 2e-16 ***
# relevel(factor(Gene), ref = "ApoA-II")ApoA-I    < 2e-16 ***
# relevel(factor(Gene), ref = "ApoA-II")ApoA-IV   < 2e-16 ***
# relevel(factor(Gene), ref = "ApoA-II")ApoA-V   2.12e-15 ***
# relevel(factor(Gene), ref = "ApoA-II")ApoC-I    < 2e-16 ***
# relevel(factor(Gene), ref = "ApoA-II")ApoC-II  1.59e-10 ***
# relevel(factor(Gene), ref = "ApoA-II")ApoC-III 1.59e-10 ***
# relevel(factor(Gene), ref = "ApoA-II")АроЕ      < 2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1455 on 167 degrees of freedom
#   (57 observations deleted due to missingness)
# Multiple R-squared:  0.649,	Adjusted R-squared:  0.6343
# F-statistic: 44.12 on 7 and 167 DF,  p-value: < 2.2e-16

