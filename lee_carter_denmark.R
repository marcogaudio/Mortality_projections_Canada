# Projections LC model, Denmark

source("data_and_models.R")

# male forecasts
LCfor_DNKm <- forecast(LCfit_DNKm, h=y.pred, kt.method = "iarima",
                       kt.order = NULL)
#female forecasts
LCfor_DNKf <- forecast(LCfit_DNKf, h=y.pred, kt.method = "iarima",
                       kt.order = NULL)

plot(LCfor_DNKm, only.kt = TRUE)
plot(LCfor_DNKf, only.kt = TRUE)


data_death_rates <- tibble(Ages = c(a.min:a.max),
                           last_pred_m=log(LCfor_DNKm$rates[,y.pred]),
                           last_fit_m=log((LCfit_DNKm$Dxt/LCfit_DNKm$Ext)[,length(Y.fit)]),
                           last_pred_f=log(LCfor_DNKf$rates[,y.pred]),
                           last_fit_f=log((LCfit_DNKf$Dxt/LCfit_DNKf$Ext)[,length(Y.fit)]))

male_death_rates <- ggplot(data = data_death_rates, aes(x = Ages,
                                                        y = last_pred_m), color = "blue") +
  geom_line() +

  xlab("Ages") +
  ylab("Death rates") +
  theme_economist() +
  geom_line(data = data_death_rates, aes(x =Ages,
                                         y=last_fit_m),
            color = "darkcyan", linetype = "dotted") +
  labs(subtitle = "evolution of death rates, Female",
       color = "Legend")



female_death_rates <- ggplot(data = data_death_rates, aes(x = Ages,
                                                          y = last_pred_f), color = "blue") +
  geom_line() +

  xlab("Ages") +
  ylab("Death rates") +
  theme_economist() +
  geom_line(data = data_death_rates, aes(x =Ages,
                                         y=last_fit_f),
            color = "darkcyan", linetype = "dotted") +
  labs(subtitle = "evolution of death rates, Male",
       color = "Legend")


ggarrange(male_death_rates, female_death_rates)



# male death rates
rates_LCfit_DNKm <- fitted(LCfit_DNKm, type = "rates")
rates_LCfor_DNKm <- LCfor_DNKm$rates
rates_LC_DNKm <- cbind(DNKmRates,rates_LCfor_DNKm) # observed vs forcasted rates

# starting from rates, we DNK obtain the death probs
q_LC_DNKm <- 1- exp(-rates_LC_DNKm) # for both past and future

# in order to consider also this pop, we use an extrapolation
q_LC_DNKm.ext  <- extrapolation.fit(q_LC_DNKm)

# female death rates
rates_LCfit_DNKf <- fitted(LCfit_DNKf, type = "rates")
rates_LCfor_DNKf <- LCfor_DNKf$rates

rates_LC_DNKf <- cbind(DNKfRates,rates_LCfor_DNKf)
q_LC_DNKf <- 1- exp(-rates_LC_DNKf)
q_LC_DNKf.ext  <- extrapolation.fit(q_LC_DNKf)

# create the dataframes

q_LC.DNKm<-tibble(as.data.frame(q_LC_DNKm.ext))
q_LC.DNKf<-tibble(as.data.frame(q_LC_DNKf.ext))



# ssurvival probs and Life expectancy

# one-year survival probability
p_LC_DNKm.ext <- 1-q_LC_DNKm.ext
# n year survival probability
p0n_LC_DNKm.ext <- apply(p_LC_DNKm.ext, 2, cumprod)
# from death prob we DNK also obtain life expectancy
ex_LC_DNKm.ext <- life.exp(q_LC_DNKm.ext)

# the same we DNK do for females
p_LC_DNKf.ext <- 1-q_LC_DNKf.ext
p0n_LC_DNKf.ext <- apply(p_LC_DNKf.ext, 2, cumprod)
ex_LC_DNKf.ext <- life.exp(q_LC_DNKf.ext)


E_birth_m <- ggplot(data=NULL, aes(x = c(1968:2118) , y = ex_LC_DNKm.ext[1,])) +
  geom_line() +
  labs(subtitle = "Life Expectancy at birth, Male") +
  xlab("Years") +
  ylab("Life expectancy") +
  theme_economist()

E_65_m <- ggplot(data=NULL, aes(x = c(1968:2118) , y = ex_LC_DNKm.ext[66,])) +
  geom_line() +
  labs(subtitle = "Life Expectancy at 65, Male") +
  xlab("Years") +
  ylab("Life expectancy") +
  theme_economist()


E_birth_f <- ggplot(data=NULL, aes(x = c(1968:2118) , y = ex_LC_DNKf.ext[1,])) +
  geom_line() +
  labs(subtitle = "Life Expectancy at birth, Female") +
  xlab("Years") +
  ylab("Life expectancy") +
  theme_economist()

E_65_f <- ggplot(data=NULL, aes(x = c(1968:2118) , y = ex_LC_DNKf.ext[66,])) +
  geom_line() +
  labs(subtitle = "Life Expectancy at 65, Female") +
  xlab("Years") +
  ylab("Life expectancy") +
  theme_economist()

ggarrange(E_birth_m, E_65_m, E_birth_f, E_65_f)


LCsim_DNKm.mrwd <- simulate(LCfit_DNKm, nsim = n.sim,
                            h=y.pred
                            # ,
                            # kt.method = "iarima",
                            # kt.order = NULL
)
rates_LC_DNKm.st <- LCsim_DNKm.mrwd$rates
q_LC_DNKm.st <- 1- exp(-rates_LC_DNKm.st)
q_LC_DNKm.st.ext <-  extrapolation.sim(q_LC_DNKm.st)

# female
LCsim_DNKf.mrwd <- simulate(LCfit_DNKf, nsim = n.sim,
                            h=y.pred, kt.method = "iarima",
                            kt.order = NULL)
rates_LC_DNKf.st <- LCsim_DNKf.mrwd$rates
q_LC_DNKf.st <- 1- exp(-rates_LC_DNKf.st)
q_LC_DNKf.st.ext <-  extrapolation.sim(q_LC_DNKf.st)

# fun charts
probs = c(2.5, 10, 25, 50, 75, 90, 97.5)
par(mfrow=c(1, 2))
matplot(LCfit_DNKm$years,
        t(q_LC_DNKm[c("65", "75", "85"),
                    c(1:length(Y.fit)) ]),
        xlim = c(y.fit.min, (y.fit.min+y.pred)),
        ylim = c(0.001, 0.2),
        pch = 20, col = "black",
        log = "y",
        xlab = "year",
        ylab = "male mortality rate (log scale)")
fan(t(LCsim_DNKm.mrwd$rates["65", , ]),
    start = y.fit.max+1,
    probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("black", "white")),
    ln = NULL)
fan(t(LCsim_DNKm.mrwd$rates["75", , ]),
    start = y.fit.max+1,
    probs = probs,
    n.fan = 4,
    fan.col = colorRampPalette(c("red", "white")),
    ln = NULL)
fan(t(LCsim_DNKm.mrwd$rates["85", , ]),
    start = y.fit.max+1,
    probs = probs,
    n.fan = 4,
    fan.col = colorRampPalette(c("blue", "white")),
    ln = NULL)
text(y.fit.min+4,
     q_LC_DNKm[c("65", "75", "85"), 30],
     labels = c("x = 65", "x = 75", "x = 85"))
matplot(LCfit_DNKf$years,
        t(q_LC_DNKf[c("65", "75", "85"),
                    c(1:length(Y.fit)) ]),
        xlim = c(y.fit.min, (y.fit.min+y.pred)),
        ylim = c(0.001, 0.2),
        pch = 20, col = "black",
        log = "y", xlab = "year", ylab = "female mortality rate (log scale)")
fan(t(LCsim_DNKf.mrwd$rates["65", , ]),
    start = y.fit.max+1,
    probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("black", "white")),
    ln = NULL)
fan(t(LCsim_DNKf.mrwd$rates["75", , ]),
    start = y.fit.max+1,
    probs = probs,
    n.fan = 4,
    fan.col = colorRampPalette(c("red", "white")),
    ln = NULL)
fan(t(LCsim_DNKf.mrwd$rates["85", , ]),
    start = y.fit.max+1,
    probs = probs,
    n.fan = 4,
    fan.col = colorRampPalette(c("blue", "white")), ln = NULL)
text(y.fit.min+4,
     q_LC_DNKf[c("65", "75", "85"), 30],
     labels = c("x = 65", "x = 75", "x = 85"))


### Bootstrap
LCboot_DNKm <- bootstrap(LCfit_DNKm,
                         nBoot = n.boot,
                         type = "semiparametric")
LCsim_DNKm.boot <- simulate(LCboot_DNKm, nsim = n.sim/n.boot, h = y.pred)
rates_LC_DNKm.boot.st <- LCsim_DNKm.boot$rates
q_LC_DNKm.boot.st <- 1- exp(-rates_LC_DNKm.boot.st)
q_LC_DNKm.boot.st.ext <-  extrapolation.sim(q_LC_DNKm.boot.st)

LCboot_DNKf <- bootstrap(LCfit_DNKf,
                         nBoot = n.boot,
                         type = "semiparametric")
LCsim_DNKf.boot <- simulate(LCboot_DNKf, nsim = n.sim/n.boot, h = y.pred)
rates_LC_DNKf.boot.st <- LCsim_DNKf.boot$rates
q_LC_DNKf.boot.st <- 1- exp(-rates_LC_DNKf.boot.st)
q_LC_DNKf.boot.st.ext <-  extrapolation.sim(q_LC_DNKf.boot.st)


# Confidence intervals at 90%

# we DNK represent the confidence intervals for the death probs and store them in a df.
# 1 year death probability
conf.lev <- 0.9
q_LC_DNKm.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_LC_DNKm.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_LC_DNKm.boot.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_LC_DNKm.boot.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))

for (j in 1:(omega-a.min)){
  for (k in 1:y.pred){
    q_LC_DNKm.q95[j,k] <- quantile(q_LC_DNKm.st.ext[j,k,],
                                   probs=0.5+conf.lev/2)
    q_LC_DNKm.q05[j,k] <- quantile(q_LC_DNKm.st.ext[j,k,],
                                   probs=(1-conf.lev)/2)
    q_LC_DNKm.boot.q95[j,k] <- quantile(q_LC_DNKm.boot.st.ext[j,k,],
                                        probs=0.5+conf.lev/2)
    q_LC_DNKm.boot.q05[j,k] <- quantile(q_LC_DNKm.boot.st.ext[j,k,],
                                        probs=(1-conf.lev)/2)
  }
}

q_LC_DNKf.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_LC_DNKf.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_LC_DNKf.boot.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_LC_DNKf.boot.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
for (j in 1:(omega-a.min)){
  for (k in 1:y.pred){
    q_LC_DNKf.q95[j,k] <- quantile(q_LC_DNKf.st.ext[j,k,], probs=0.5+conf.lev/2)
    q_LC_DNKf.q05[j,k] <- quantile(q_LC_DNKf.st.ext[j,k,], probs=(1-conf.lev)/2)
    q_LC_DNKf.boot.q95[j,k] <- quantile(q_LC_DNKf.boot.st.ext[j,k,], probs=0.5+conf.lev/2)
    q_LC_DNKf.boot.q05[j,k] <- quantile(q_LC_DNKf.boot.st.ext[j,k,], probs=(1-conf.lev)/2)
  }
}



### Actuarial quantities

#### death probabilities at 65

# data preparation
death_probs_ci_m <- tibble(q95 = q_LC_DNKm.q95["65",],
                           q05 = q_LC_DNKm.q05["65",],
                           q95boot = q_LC_DNKm.boot.q95["65",],
                           q05boot = q_LC_DNKm.boot.q05["65",],
                           years =seq(2019,2118),
)


death_probs_ci_m <- death_probs_ci_m %>% gather("value", "rates", 1:4)

death_prob_m <- tibble(years = seq(1968,2118),
                       value = "qx",
                       rates = q_LC_DNKm.ext["65",]
)

death_probs_ci_f <- tibble(q95 = q_LC_DNKf.q95["65",],
                           q05 = q_LC_DNKf.q05["65",],
                           q95boot = q_LC_DNKf.boot.q95["65",],
                           q05boot = q_LC_DNKf.boot.q05["65",],
                           years =seq(2019,2118),
)


death_probs_ci_f <- death_probs_ci_f %>%  gather("value", "rates", 1:4)

death_prob_f <- tibble(years = seq(1968,2118),
                       value = "qx",
                       rates = q_LC_DNKf.ext["65",]
)

data_deaths_m <- rbind(death_prob_m, death_probs_ci_m)
data_deaths_f <- rbind(death_prob_f, death_probs_ci_f)

d_p_m <- ggplot(data = data_deaths_m) +
  geom_line( aes(x = years, y = rates, colour = value)) +
  theme_economist() +
  labs(subtitle = "death probabilities at 65 + C.I., Male") +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


d_p_f <- ggplot(data = data_deaths_f) +
  geom_line( aes(x = years, y = rates, colour = value)) +
  theme_economist() +
  labs(subtitle = "death probabilities at 65 + C.I., Female") +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



ggarrange(d_p_m, d_p_f)


#### death probabilities at age 20


# data preparation
death_probs_ci_m <- tibble(q95 = q_LC_DNKm.q95["20",],
                           q05 = q_LC_DNKm.q05["20",],
                           q95boot = q_LC_DNKm.boot.q95["20",],
                           q05boot = q_LC_DNKm.boot.q05["20",],
                           years =seq(2019,2118),
)


death_probs_ci_m <- death_probs_ci_m %>%  gather("value", "rates", 1:4)

death_prob_m <- tibble(years = seq(1968,2118),
                       value = "qx",
                       rates = q_LC_DNKm.ext["20",]
)

death_probs_ci_f <- tibble(q95 = q_LC_DNKf.q95["20",],
                           q05 = q_LC_DNKf.q05["20",],
                           q95boot = q_LC_DNKf.boot.q95["20",],
                           q05boot = q_LC_DNKf.boot.q05["20",],
                           years =seq(2019,2118),
)


death_probs_ci_f <- death_probs_ci_f %>%  gather("value", "rates", 1:4)

death_prob_f <- tibble(years = seq(1968,2118),
                       value = "qx",
                       rates = q_LC_DNKf.ext["20",]
)

data_deaths_m <- rbind(death_prob_m, death_probs_ci_m)
data_deaths_f <- rbind(death_prob_f, death_probs_ci_f)

d_p_m <- ggplot(data = data_deaths_m) +
  geom_line( aes(x = years, y = rates, colour = value)) +
  theme_economist() +
  labs(subtitle = "death probabilities at 20 + C.I., Male") +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


d_p_f <- ggplot(data = data_deaths_f) +
  geom_line( aes(x = years, y = rates, colour = value)) +
  theme_economist() +
  labs(subtitle = "death probabilities at 20 + C.I., Female") +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



ggarrange(d_p_m, d_p_f)




#### Present values of stochastic annuities


p_LC_DNKm.st.ext <- 1-q_LC_DNKm.st.ext # 1 year survival probability
ann_LC_DNKm.st <- annuity.st(q_LC_DNKm.st.ext,v_vect,0.9)
ann_LC_DNKm.mean <- ann_LC_DNKm.st[[1]]
ann_LC_DNKm.q95 <- ann_LC_DNKm.st[[2]]
ann_LC_DNKm.q05 <- ann_LC_DNKm.st[[3]]

p_LC_DNKf.st.ext <- 1-q_LC_DNKf.st.ext # 1 year survival probability
ann_LC_DNKf.st <- annuity.st(q_LC_DNKf.st.ext,v_vect,0.9)
ann_LC_DNKf.mean <- ann_LC_DNKf.st[[1]]
ann_LC_DNKf.q95 <- ann_LC_DNKf.st[[2]]
ann_LC_DNKf.q05 <- ann_LC_DNKf.st[[3]]

# data preparation

annuity_male <- tibble(annuity_mean = ann_LC_DNKm.mean["65",],
                       q05 = ann_LC_DNKm.q05["65",],
                       q95 = ann_LC_DNKm.q95["65",],
                       years = seq(2019,2118))

annuity_male <- annuity_male %>% gather("lab", "value", 1:3)

annuity_female <- tibble(annuity_mean = ann_LC_DNKf.mean["65",],
                         q05 = ann_LC_DNKf.q05["65",],
                         q95 = ann_LC_DNKf.q95["65",],
                         years = seq(2019,2118))

annuity_female <- annuity_female %>% gather("lab", "value", 1:3)


ann_2050_m <- tibble(annuity_mean = ann_LC_DNKm.mean[,"2050"],
                     q05 = ann_LC_DNKm.q05[,"2050"],
                     q95 = ann_LC_DNKm.q95[,"2050"],
                     ages = seq(0,109))

annuity_2050_m <- ann_2050_m %>% gather("lab", "value", 1:3)

ann_2050_f <- tibble(annuity_mean = ann_LC_DNKf.mean[,"2050"],
                     q05 = ann_LC_DNKf.q05[,"2050"],
                     q95 = ann_LC_DNKf.q95[,"2050"],
                     ages = seq(0,109))

annuity_2050_f <- ann_2050_f %>% gather("lab", "value", 1:3)



an_male <- ggplot(data = annuity_male) +
  geom_line(aes(x = years, y = value, color = lab)) +
  theme_economist() +
  labs(subtitle = "stochastic annuity at 65 + C.I., Male") +
  theme(legend.position = "right")

an_fem <- ggplot(data = annuity_female) +
  geom_line(aes(x = years, y = value, color = lab)) +
  theme_economist() +
  labs(subtitle = "stochastic annuity at 65 + C.I., Female") +
  theme(legend.position = "right")


an_male_2050 <- ggplot(data = annuity_2050_m) +
  geom_line(aes(x = ages, y = value, color = lab)) +
  theme_economist() +
  labs(subtitle = "stochastic annuity 2050 + C.I., Male") +
  theme(legend.position = "right")

an_female_2050 <- ggplot(data = annuity_2050_f) +
  geom_line(aes(x = ages, y = value, color = lab)) +
  theme_economist() +
  labs(subtitle = "stochastic annuity 2050 + C.I., Female") +
  theme(legend.position = "right")

ggarrange(an_male, an_fem, an_male_2050, an_female_2050)



