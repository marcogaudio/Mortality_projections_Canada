# data and models

library(demography)
library(StMoMo)
library(MTS)
library(devtools)
library(writexl)
library(readxl)
source("Functions.r") # uploading functions from an other script in the directory
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(DT)
library(fanplot)
library(lattice)


# starting age
a.min <- 0
# maximum age
a.max <- 100
# extreme age
omega <-110
# age interval
A.fit <-c(a.min:a.max)
# fitting time interval
y.fit.min <- 1968
y.fit.max <- 2018
Y.fit <- c(y.fit.min:y.fit.max)
# projection interval
y.pred <- 100
Y.pred <- c((y.fit.max+1):(y.fit.max+y.pred))

n.sim <- 2500
n.boot <- 250

# vector of interest rates (flat)
int_vect <- c(rep(0.02,omega-a.min))
# vector of discount factors
v_vect <- cumprod((1+int_vect)^-1)


# Demographic Data about Denmark, loaded from HMD.
Denmark_data <- hmd.mx(country = "DNK",
                       username = paste0("marcogaudio1997@outlook.it"),
                       password = 1623581009, label = "Denmark")



# distinguish between male and female
DNKmStMoMoC <- StMoMoData(Denmark_data, series = "male")
DNKfStMoMoC <- StMoMoData(Denmark_data, series = "female")

# from central exposed to risk to initial exposed to risk
DNKmStMoMoI <- central2initial(DNKmStMoMoC)
DNKfStMoMoI <- central2initial(DNKfStMoMoC)

# specify the weights we give to the data: we need to pass age range and time
# range.
wxt <- genWeightMat(ages = A.fit, years = Y.fit, clip = 4)


# due to the fact that we have 0 deaths in 2018 for male pop., the log of the quantity below returns -inf. So, this code is to fix this problem, basically it is assumed 2 deaths at age 3 and 8 for last fitted year
DNKmStMoMoC$Dxt[4,184] <- 2
DNKmStMoMoC$Dxt[8,184] <- 2
DNKfStMoMoC$Dxt[5,184] <- 2
DNKfStMoMoC$Dxt[11,184] <- 2

DNKmStMoMoC$Dxt[which(DNKmStMoMoC$Dxt == 0)] = 2

DNKmStMoMoI <- central2initial(DNKmStMoMoC)
DNKfStMoMoI <- central2initial(DNKfStMoMoC)



# we DNK get death rates for male and female pop.
DNKmRates <- DNKmStMoMoC$Dxt/DNKmStMoMoC$Ext

#matrix of the males death rates
DNKmRates <- DNKmRates[A.fit+1,tail(DNKmStMoMoC$years+1,length(Y.fit))-DNKmStMoMoC$years[1]]

DNKfRates <- DNKfStMoMoC$Dxt/DNKfStMoMoC$Ext

#matrix of the males death rates
DNKfRates <- DNKfRates[A.fit+1,tail(DNKfStMoMoC$years+1,length(Y.fit))-DNKfStMoMoC$years[1]]



# fitting models

# log link
LCfit_DNKm <- fit(lc(link = "log"),
                  data = DNKmStMoMoC,
                  ages.fit = a.min:a.max,
                  years.fit = y.fit.min:y.fit.max,
                  wxt = wxt)

LCfit_DNKf <- fit(lc(link = "log"),
                  data = DNKfStMoMoC,
                  ages.fit = a.min:a.max,
                  years.fit = y.fit.min:y.fit.max,
                  wxt = wxt)

# in RH model we need to specify an initial value for alpha and beta.
# one of the problem of this model is that it is slow to be fitted.
# we use the starting value of LC to initialize the algo.
RHfit_DNKm <- fit(rh(link = "logit", cohortAgeFun="1"),
                  data = DNKmStMoMoI,
                  ages.fit = a.min:a.max,
                  years.fit=y.fit.min:y.fit.max,
                  wxt = wxt,
                  start.ax = LCfit_DNKm$ax,
                  start.bx = LCfit_DNKm$bx,
                  start.kt = LCfit_DNKm$kt)

RHfit_DNKf <- fit(rh(link = "logit", cohortAgeFun="1"),
                  data = DNKfStMoMoI,
                  ages.fit = a.min:a.max,years.fit=y.fit.min:y.fit.max,
                  wxt = wxt,
                  start.ax = LCfit_DNKf$ax,
                  start.bx = LCfit_DNKf$bx,
                  start.kt = LCfit_DNKf$kt)

APCfit_DNKm <- fit(apc(link = "logit"),
                   data = DNKmStMoMoI,
                   ages.fit = a.min:a.max,
                   years.fit=y.fit.min:y.fit.max,
                   wxt = wxt)

APCfit_DNKf <- fit(apc(link = "logit"),
                   data = DNKfStMoMoI,
                   ages.fit = a.min:a.max,
                   years.fit = y.fit.min:y.fit.max,
                   wxt = wxt)

CBDfit_DNKm <- fit(cbd(),
                   data = DNKmStMoMoI,
                   ages.fit = a.min:a.max,
                   years.fit = y.fit.min:y.fit.max,
                   wxt = wxt)

CBDfit_DNKf <- fit(cbd(),
                   data = DNKfStMoMoI,
                   ages.fit = a.min:a.max,
                   years.fit = y.fit.min:y.fit.max,
                   wxt = wxt)

M7fit_DNKm <- fit(m7(link = "logit"),
                  data = DNKmStMoMoI,
                  ages.fit = a.min:a.max,
                  years.fit = y.fit.min:y.fit.max,
                  wxt = wxt)

M7fit_DNKf <- fit(m7(link = "logit"),
                  data = DNKfStMoMoI,
                  ages.fit = a.min:a.max,
                  years.fit = y.fit.min:y.fit.max,
                  wxt = wxt)

# the function PLAT is defined in the scirpt functions.R
PLATfit_DNKm <- fit(PLAT, data = DNKmStMoMoI,
                    ages.fit = a.min:a.max,
                    years.fit = y.fit.min:y.fit.max,
                    wxt = wxt)

PLATfit_DNKf <- fit(PLAT, data = DNKfStMoMoI,
                    ages.fit = a.min:a.max,
                    years.fit = y.fit.min:y.fit.max,
                    wxt = wxt)


# models comparison
LL <- tibble(Gender = c("Male", "Female"),
             Lee_Carter = c(LCfit_DNKm$loglik, LCfit_DNKf$loglik),
             Renshaw_Haberman = c(RHfit_DNKm$loglik, RHfit_DNKm$loglik),
             APC = c(APCfit_DNKm$loglik, APCfit_DNKm$loglik),
             CBD = c(CBDfit_DNKm$loglik, CBDfit_DNKf$loglik),
             M7 = c(M7fit_DNKm$loglik, M7fit_DNKf$loglik),
             Plat = c(PLATfit_DNKm$loglik, PLATfit_DNKm$loglik))
LL


info_criteria_male <- tibble(
  Criteria = c("BIC", "AIC", "n.parameters"),
  Lee_Carter = c(BIC(LCfit_DNKm), AIC(LCfit_DNKm),
                 LCfit_DNKm$npar),
  Renshaw_Haberman = c(BIC(RHfit_DNKm), AIC(RHfit_DNKm),
                       RHfit_DNKm$npar),
  APC = c(BIC(APCfit_DNKm), AIC(APCfit_DNKm), APCfit_DNKm$npar),
  CBD = c(BIC(CBDfit_DNKm), AIC(CBDfit_DNKm), CBDfit_DNKm$npar),
  M7 = c(BIC(M7fit_DNKm), AIC(M7fit_DNKm), M7fit_DNKm$npar),
  Plat = c(BIC(PLATfit_DNKm), AIC(PLATfit_DNKm),
           PLATfit_DNKm$npar))

info_criteria_female <- tibble(
  Criteria = c("BIC", "AIC", "n.parameters"),
  Lee_Carter = c(BIC(LCfit_DNKf), AIC(LCfit_DNKf),
                 LCfit_DNKf$npar),
  Renshaw_Haberman = c(BIC(RHfit_DNKf), AIC(RHfit_DNKf),
                       RHfit_DNKf$npar),
  APC = c(BIC(APCfit_DNKf), AIC(APCfit_DNKf), APCfit_DNKf$npar),
  CBD = c(BIC(CBDfit_DNKf), AIC(CBDfit_DNKf), CBDfit_DNKf$npar),
  M7 = c(BIC(M7fit_DNKf), AIC(M7fit_DNKf), M7fit_DNKf$npar),
  Plat = c(BIC(PLATfit_DNKf), AIC(PLATfit_DNKf),
           PLATfit_DNKf$npar))
info_criteria_male
info_criteria_female



## Residuals

LCres_DNKm <- residuals(LCfit_DNKm)
LCres_DNKf <- residuals(LCfit_DNKf)


# this function DNK be found in the R script Function.R
reshape_heatmap <- function(residuals, ages){

  data_heatmap <- tibble(as.data.frame(residuals)) %>%
    pivot_longer(cols = `1968`:`2018`) %>%
    rename(years = "name",
           residuals = "value") %>%
    arrange(years) %>%
    mutate(ages = rep(ages,51))
  data_heatmap$years <- as.numeric(data_heatmap$years)

  return(data_heatmap)
}

data_LC_male <- reshape_heatmap(LCres_DNKm$residuals,
                                LCres_DNKm$ages)
data_LC_female <- reshape_heatmap(LCres_DNKf$residuals,
                                  LCres_DNKf$ages)

male_plot <- ggplot(data = data_LC_male, aes(x = years,
                                             y = ages,
                                             fill = residuals)) +
  geom_tile() +
  scale_fill_distiller(palette = "Blues") +
  labs(subtitle = "Lee-Carter residuals, Male") +
  theme_economist() +
  scale_x_continuous(breaks=seq(1970, 2010, 10)) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

female_plot <- ggplot(data = data_LC_female, aes(x = years,
                                                 y = ages,
                                                 fill = residuals)) +
  geom_tile() +
  scale_fill_distiller(palette = "Blues") +
  labs(subtitle = "Lee-Carter residuals, Female" )+
  theme_economist() +
  scale_x_continuous(breaks=seq(1970, 2010, 10)) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggarrange(female_plot, male_plot)


APCres_DNKm <- residuals(APCfit_DNKm)
APCres_DNKf <- residuals(APCfit_DNKf)


# this function DNK be found in the R script Function.R
reshape_heatmap <- function(residuals, ages){

  data_heatmap <- tibble(as.data.frame(residuals)) %>%
    pivot_longer(cols = `1968`:`2018`) %>%
    rename(years = "name",
           residuals = "value") %>%
    arrange(years) %>%
    mutate(ages = rep(ages,51))
  data_heatmap$years <- as.numeric(data_heatmap$years)

  return(data_heatmap)
}

data_APC_male <- reshape_heatmap(APCres_DNKm$residuals,
                                 APCres_DNKm$ages)
data_APC_female <- reshape_heatmap(APCres_DNKf$residuals,
                                   APCres_DNKf$ages)

male_plot <- ggplot(data = data_APC_male, aes(x = years,
                                              y = ages,
                                              fill = residuals)) +
  geom_tile() +
  scale_fill_distiller(palette = "Blues") +
  labs(subtitle = "Lee-Carter residuals, Male") +
  theme_economist() +
  scale_x_continuous(breaks=seq(1970, 2010, 10)) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

female_plot <- ggplot(data = data_APC_female, aes(x = years,
                                                  y = ages,
                                                  fill = residuals)) +
  geom_tile() +
  scale_fill_distiller(palette = "Blues") +
  labs(subtitle = "Lee-Carter residuals, Female" )+
  theme_economist() +
  scale_x_continuous(breaks=seq(1970, 2010, 10)) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggarrange(female_plot, male_plot)



PLATres_DNKm <- residuals(PLATfit_DNKm)
PLATres_DNKf <- residuals(PLATfit_DNKf)

data_PLAT_male <- reshape_heatmap(PLATres_DNKm$residuals,
                                  PLATres_DNKm$ages)
data_PLAT_female <- reshape_heatmap(PLATres_DNKf$residuals,
                                    PLATres_DNKf$ages)

male_plot <- ggplot(data = data_PLAT_male, aes(x = years,
                                               y = ages,
                                               fill = residuals)) +
  geom_tile() +
  scale_fill_distiller(palette = "Blues") +
  labs(subtitle = "PLAT residuals, Male") +
  theme_economist() +
  scale_x_continuous(breaks=seq(1970, 2010, 10)) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

female_plot <- ggplot(data = data_PLAT_female, aes(x = years,
                                                   y = ages,
                                                   fill = residuals)) +
  geom_tile() +
  scale_fill_distiller(palette = "Blues") +
  labs(subtitle = "PLAT residuals, Female") +
  theme_economist() +
  scale_x_continuous(breaks=seq(1970, 2010, 10)) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggarrange(female_plot, male_plot)



res_cy_m <- ggplot(data = data_LC_male, aes(x = years, y = residuals)) +
  geom_point() +
  theme_economist() +
  labs(subtitle = "Lee-Carter residuals for calendar year, Male" )+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

res_cy_f <- ggplot(data = data_LC_female, aes(x = years, y = residuals)) +
  geom_point() +
  theme_economist() +
  labs(subtitle = "Lee-Carter residuals for calendar year, Female" )+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




# residuals lee carter by calendar year, male and female

res_age_m <- ggplot(data = data_LC_male, aes(x = ages, y = residuals)) +
  geom_point() +
  theme_economist() +
  labs(subtitle = "Lee-Carter residuals for ages, Male" )+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

res_age_f <- ggplot(data = data_LC_female, aes(x = ages, y = residuals)) +
  geom_point() +
  theme_economist() +
  labs(subtitle = "Lee-Carter residuals for ages, Female" )+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



ggarrange(res_age_f,res_age_m, res_cy_m, res_cy_f)


res_cy_m <- ggplot(data = data_APC_male, aes(x = years, y = residuals)) +
  geom_point() +
  theme_economist() +
  labs(subtitle = "APC residuals for calendar year, Male" )+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

res_cy_f <- ggplot(data = data_APC_female, aes(x = years, y = residuals)) +
  geom_point() +
  theme_economist() +
  labs(subtitle = "APC residuals for calendar year, Female" )+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




# residuals lee carter by calendar year, male and female

res_age_m <- ggplot(data = data_APC_male, aes(x = ages, y = residuals)) +
  geom_point() +
  theme_economist() +
  labs(subtitle = "APC residuals for ages, Male" )+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

res_age_f <- ggplot(data = data_APC_female, aes(x = ages, y = residuals)) +
  geom_point() +
  theme_economist() +
  labs(subtitle = "APC residuals for ages, Female" )+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



ggarrange(res_age_f,res_age_m, res_cy_m, res_cy_f)


res_cy_m <- ggplot(data = data_PLAT_male, aes(x = years, y = residuals)) +
  geom_point() +
  theme_economist() +
  labs(subtitle = "PLAT residuals for calendar year, Male" )+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

res_cy_f <- ggplot(data = data_PLAT_female, aes(x = years, y = residuals)) +
  geom_point() +
  theme_economist() +
  labs(subtitle = "PLAT residuals for calendar year, Female" )+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




# residuals lee carter by calendar year, male and female

res_age_m <- ggplot(data = data_PLAT_male, aes(x = ages, y = residuals)) +
  geom_point() +
  theme_economist() +
  labs(subtitle = "PLAT residuals for ages, Male" )+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

res_age_f <- ggplot(data = data_PLAT_female, aes(x = ages, y = residuals)) +
  geom_point() +
  theme_economist() +
  labs(subtitle = "PLAT residuals for ages, Female" )+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



ggarrange(res_age_f,res_age_m, res_cy_m, res_cy_f)


# Parameters representation

# even in this case, we need to reshape the data
lee_carter_params <- tibble(ages =  as.numeric(LCfit_DNKm$ages),
                            ax_male = LCfit_DNKm$ax,
                            bx_male = as.numeric(LCfit_DNKm$bx),
                            ax_female = LCfit_DNKf$ax,
                            bx_female = as.numeric(LCfit_DNKf$bx))

kt_lee_carter <- tibble(kt_male = as.numeric(LCfit_DNKm$kt),
                        kt_female = as.numeric(LCfit_DNKf$kt),
                        years = LCfit_DNKm$years)



ax_m <- ggplot(data = lee_carter_params, aes(x = ages, y = ax_male)) +
  geom_line() +
  xlab("AGE") +
  ylab("ax") +
  labs(subtitle = "alpha plot, Male")  +
  theme_economist() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

bx_m <- ggplot(data = lee_carter_params, aes(x = ages, y = bx_male)) +
  geom_line() +
  xlab("AGE") +
  ylab("bx") +
  labs(subtitle = "beta plot, Male") +
  theme_economist() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


kt_m <- ggplot(data = kt_lee_carter, aes(x = years , y = kt_male)) +
  geom_line() +
  xlab("AGE") +
  ylab("bx") +
  labs(subtitle = "kt plot, Male") +
  theme_economist() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ax_f <- ggplot(data = lee_carter_params, aes(x = ages, y = ax_female)) +
  geom_line() +
  xlab("AGE") +
  ylab("ax") +
  labs(subtitle = "alpha plot, Female") +
  theme_economist()

bx_f <- ggplot(data = lee_carter_params, aes(x = ages, y = bx_female)) +
  geom_line() +
  xlab("AGE") +
  ylab("bx") +
  labs(subtitle = "beta plot, Female") +
  theme_economist() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


kt_f <- ggplot(data = kt_lee_carter, aes(x = years , y = kt_female)) +
  geom_line() +
  xlab("AGE") +
  ylab("bx") +
  labs(subtitle = "kt plot, Female") +
  theme_economist() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggarrange(ax_m, bx_m, kt_m, ax_f, bx_f, kt_f)

### APC
APC_params <- tibble(ages = as.numeric(APCfit_DNKm$ages),
                     ax_male = APCfit_DNKm$ax,
                     ax_female = APCfit_DNKf$ax)

kt_APC <- tibble(kt1_male = as.numeric(APCfit_DNKm$kt),
                 kt1_female = as.numeric(APCfit_DNKf$kt),
                 years = APCfit_DNKm$years)

gc_APC <- tibble(gc_male = APCfit_DNKm$gc,
                 gc_female = APCfit_DNKf$gc,
                 cohorts = APCfit_DNKm$cohorts)

ax_m <- ggplot(data = APC_params, aes(x = ages, y = ax_male)) +
  geom_line() +
  labs(subtitle = "alpha plot, Male") +
  theme_economist() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ax_f <- ggplot(data = APC_params, aes(x = ages, y = ax_female)) +
  geom_line() +
  labs(subtitle = "alpha plot, Female") +
  theme_economist() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



kt1_m <- ggplot(data = kt_APC, aes(x = years, y = kt1_male)) +
  geom_line() +
  labs(subtitle = "kt plot, Male", x = "YEARS") +
  theme_economist() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


kt1_f <- ggplot(data = kt_APC, aes(x = years, y = kt1_female)) +
  geom_line() +
  labs(subtitle = "kt plot, Female", x = "YEARS") +
  theme_economist() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


gc_m <- ggplot(data = gc_APC, aes(x = cohorts, y = gc_male)) +
  geom_line() +
  labs(subtitle = "gamma plot, Male") +
  theme_economist() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



gc_f <- ggplot(data = gc_APC, aes(x = cohorts, y = gc_female)) +
  geom_line() +
  labs(subtitle = "gamma plot, Female") +
  theme_economist() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggarrange(ax_m, ax_f, kt1_m, kt1_f, gc_m, gc_f)



### Plat
PLAT_params <- tibble(ages = as.numeric(PLATfit_DNKm$ages),
                      ax_male = PLATfit_DNKm$ax,
                      ax_female = PLATfit_DNKf$ax)

kt_PLAT <- tibble(kt1_male = as.numeric(PLATfit_DNKm$kt[1,]),
                  kt2_male = as.numeric(PLATfit_DNKm$kt[2,]),
                  kt1_female = as.numeric(PLATfit_DNKf$kt[1,]),
                  kt2_female = as.numeric(PLATfit_DNKf$kt[2,]),
                  years = PLATfit_DNKm$years)

gc_PLAT <- tibble(gc_male = PLATfit_DNKm$gc,
                  gc_female = PLATfit_DNKf$gc,
                  cohorts = PLATfit_DNKm$cohorts)


ax_m <- ggplot(data = PLAT_params, aes(x = ages, y = ax_male)) +
  geom_line() +
  labs(subtitle = "alpha plot, Male") +
  theme_economist() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ax_f <- ggplot(data = PLAT_params, aes(x = ages, y = ax_female)) +
  geom_line() +
  labs(subtitle = "alpha plot, Female") +
  theme_economist() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


kt1_m <- ggplot(data = kt_PLAT, aes(x = years, y = kt1_male)) +
  geom_line() +
  labs(subtitle = "Kt1 plot, Male") +
  theme_economist() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


kt1_f <- ggplot(data = kt_PLAT, aes(x = years, y = kt1_female)) +
  geom_line() +
  labs(subtitle = "Kt1 plot, Female") +
  theme_economist() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


kt2_m <- ggplot(data = kt_PLAT, aes(x = years, y = kt2_male)) +
  geom_line() +
  labs(subtitle = "Kt2 plot, Male") +
  theme_economist() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


kt2_f <- ggplot(data = kt_PLAT, aes(x = years, y = kt2_female)) +
  geom_line() +
  labs(subtitle = "Kt2 plot, Female") +
  theme_economist() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



gc_m <- ggplot(data = gc_PLAT, aes(x = cohorts, y = gc_male)) +
  geom_line() +
  labs(subtitle = "gamma plot, Male") +
  theme_economist() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



gc_f <- ggplot(data = gc_PLAT, aes(x = cohorts, y = gc_female)) +
  geom_line() +
  labs(subtitle = "gamma plot, Female") +
  theme_economist() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggarrange(ax_m, ax_f, kt1_m, kt1_f, kt2_m, kt2_f, gc_m, gc_f)






