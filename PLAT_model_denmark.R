## Plat model.

source("data_and_models.R")

# male forecasts
PLATfor_DNKm <- forecast(PLATfit_DNKm, h=y.pred, gc.order = c(2, 0, 0))
#female forecasts
PLATfor_DNKf <- forecast(PLATfit_DNKf, h=y.pred, gc.order = c(2, 0, 0))

# to be changed!!
plot(PLATfor_DNKm, only.kt = TRUE)
plot(PLATfor_DNKf, only.kt = TRUE)
plot(PLATfor_DNKm, only.gc = TRUE)
plot(PLATfor_DNKf, only.gc = TRUE)


### Death rates: last observed vs last projected


data_death_rates <- tibble(Ages = c(a.min:a.max),
                           last_pred_m=log(PLATfor_DNKm$rates[,y.pred]),
                           last_fit_m=log((PLATfit_DNKm$Dxt/PLATfit_DNKm$Ext)[,length(Y.fit)]),
                           last_pred_f=log(PLATfor_DNKf$rates[,y.pred]),
                           last_fit_f=log((PLATfit_DNKf$Dxt/PLATfit_DNKf$Ext)[,length(Y.fit)]))

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
rates_PLATfit_DNKm <- fitted(PLATfit_DNKm, type = "rates")
rates_PLATfor_DNKm <- PLATfor_DNKm$rates

# starting from rates, we DNK obtain the death probs
rates_PLAT_DNKm <- cbind(DNKmRates,rates_PLATfor_DNKm)
q_PLAT_DNKm <- 1- exp(-rates_PLAT_DNKm) # for both past and future

# in order to consider also this pop, we use an extrapolation
q_PLAT_DNKm.ext  <- extrapolation.fit(q_PLAT_DNKm)

# female death rates
rates_PLATfit_DNKf <- fitted(PLATfit_DNKf, type = "rates")
rates_PLATfor_DNKf <- PLATfor_DNKf$rates
rates_PLAT_DNKf <- cbind(DNKfRates,rates_PLATfor_DNKf)
q_PLAT_DNKf <- 1- exp(-rates_PLAT_DNKf)
q_PLAT_DNKf.ext  <- extrapolation.fit(q_PLAT_DNKf)

# create the dataframes

q_PLAT.DNKm<-tibble(as.data.frame(q_PLAT_DNKm.ext))
q_PLAT.DNKf<-tibble(as.data.frame(q_PLAT_DNKf.ext))



### Survival Probabilities and Life expectancies.




# one-year survival probability
p_PLAT_DNKm.ext <- 1-q_PLAT_DNKm.ext
# n year survival probability
p0n_PLAT_DNKm.ext <- apply(p_PLAT_DNKm.ext, 2, cumprod)
# from death prob we DNK also obtain life expectancy
ex_PLAT_DNKm.ext <- life.exp(q_PLAT_DNKm.ext)

# the same we DNK do for females
p_PLAT_DNKf.ext <- 1-q_PLAT_DNKf.ext
p0n_PLAT_DNKf.ext <- apply(p_PLAT_DNKf.ext, 2, cumprod)
ex_PLAT_DNKf.ext <- life.exp(q_PLAT_DNKf.ext)



E_birth_m <- ggplot(data=NULL, aes(x = c(1968:2118) , y = ex_PLAT_DNKm.ext[1,])) +
  geom_line() +
  labs(subtitle = "Life Expectancy at birth, Male") +
  xlab("Years") +
  ylab("Life expectancy") +
  theme_economist()

E_65_m <- ggplot(data=NULL, aes(x = c(1968:2118) , y = ex_PLAT_DNKm.ext[66,])) +
  geom_line() +
  labs(subtitle = "Life Expectancy at 65, Male") +
  xlab("Years") +
  ylab("Life expectancy") +
  theme_economist()


E_birth_f <- ggplot(data=NULL, aes(x = c(1968:2118) , y = ex_PLAT_DNKf.ext[1,])) +
  geom_line() +
  labs(subtitle = "Life Expectancy at birth, Female") +
  xlab("Years") +
  ylab("Life expectancy") +
  theme_economist()

E_65_f <- ggplot(data=NULL, aes(x = c(1968:2118) , y = ex_PLAT_DNKf.ext[66,])) +
  geom_line() +
  labs(subtitle = "Life Expectancy at 65, Female") +
  xlab("Years") +
  ylab("Life expectancy") +
  theme_economist()

ggarrange(E_birth_m, E_65_m, E_birth_f, E_65_f)




### Evolution of death rates

#male
PLATsim_DNKm.mrwd <- simulate(PLATfit_DNKm,
                              nsim = n.sim,
                              h=y.pred,
                              gc.order = c(2, 0, 0))

rates_PLAT_DNKm.st <- PLATsim_DNKm.mrwd$rates
q_PLAT_DNKm.st <- 1- exp(-rates_PLAT_DNKm.st)
q_PLAT_DNKm.st.ext <-  extrapolation.sim(q_PLAT_DNKm.st)

# female
PLATsim_DNKf.mrwd <- simulate(PLATfit_DNKf,
                              nsim = n.sim,
                              h=y.pred,
                              gc.order = c(2, 0, 0))

rates_PLAT_DNKf.st <- PLATsim_DNKf.mrwd$rates
q_PLAT_DNKf.st <- 1- exp(-rates_PLAT_DNKf.st)
q_PLAT_DNKf.st.ext <-  extrapolation.sim(q_PLAT_DNKf.st)

# fun charts
probs = c(2.5, 10, 25, 50, 75, 90, 97.5)
par(mfrow=c(1, 2))
matplot(PLATfit_DNKm$years,
        t(q_PLAT_DNKm[c("65", "75", "85"),
                      c(1:length(Y.fit)) ]),
        xlim = c(y.fit.min, (y.fit.min+y.pred)),
        ylim = c(0.001, 0.2),
        pch = 20, col = "black",
        log = "y",
        xlab = "year",
        ylab = "male mortality rate (log scale)")
fan(t(PLATsim_DNKm.mrwd$rates["65", , ]),
    start = y.fit.max+1,
    probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("black", "white")),
    ln = NULL)
fan(t(PLATsim_DNKm.mrwd$rates["75", , ]),
    start = y.fit.max+1,
    probs = probs,
    n.fan = 4,
    fan.col = colorRampPalette(c("red", "white")),
    ln = NULL)
fan(t(PLATsim_DNKm.mrwd$rates["85", , ]),
    start = y.fit.max+1,
    probs = probs,
    n.fan = 4,
    fan.col = colorRampPalette(c("blue", "white")),
    ln = NULL)
text(y.fit.min+4,
     q_PLAT_DNKm[c("65", "75", "85"), 30],
     labels = c("x = 65", "x = 75", "x = 85"))
matplot(PLATfit_DNKf$years,
        t(q_PLAT_DNKf[c("65", "75", "85"),
                      c(1:length(Y.fit)) ]),
        xlim = c(y.fit.min, (y.fit.min+y.pred)),
        ylim = c(0.001, 0.2),
        pch = 20, col = "black",
        log = "y", xlab = "year", ylab = "female mortality rate (log scale)")
fan(t(PLATsim_DNKf.mrwd$rates["65", , ]),
    start = y.fit.max+1,
    probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("black", "white")),
    ln = NULL)
fan(t(PLATsim_DNKf.mrwd$rates["75", , ]),
    start = y.fit.max+1,
    probs = probs,
    n.fan = 4,
    fan.col = colorRampPalette(c("red", "white")),
    ln = NULL)
fan(t(PLATsim_DNKf.mrwd$rates["85", , ]),
    start = y.fit.max+1,
    probs = probs,
    n.fan = 4,
    fan.col = colorRampPalette(c("blue", "white")), ln = NULL)
text(y.fit.min+4,
     q_PLAT_DNKf[c("65", "75", "85"), 30],
     labels = c("x = 65", "x = 75", "x = 85"))



### Bootstrap



PLATboot_DNKm <- bootstrap(PLATfit_DNKm,
                           nBoot = n.boot,
                           type = "semiparametric")
PLATsim_DNKm.boot <- simulate(PLATboot_DNKm, nsim = n.sim/n.boot, h = y.pred)
rates_PLAT_DNKm.boot.st <- PLATsim_DNKm.boot$rates
q_PLAT_DNKm.boot.st <- 1- exp(-rates_PLAT_DNKm.boot.st)
q_PLAT_DNKm.boot.st.ext <-  extrapolation.sim(q_PLAT_DNKm.boot.st)

PLATboot_DNKf <- bootstrap(PLATfit_DNKf,
                           nBoot = n.boot,
                           type = "semiparametric")
PLATsim_DNKf.boot <- simulate(PLATboot_DNKf, nsim = n.sim/n.boot, h = y.pred)
rates_PLAT_DNKf.boot.st <- PLATsim_DNKf.boot$rates
q_PLAT_DNKf.boot.st <- 1- exp(-rates_PLAT_DNKf.boot.st)
q_PLAT_DNKf.boot.st.ext <-  extrapolation.sim(q_PLAT_DNKf.boot.st)


# Confidence intervals at 90%

# we DNK represent the confidence intervals for the death probs and store them in a df.
# 1 year death probability
conf.lev <- 0.9
q_PLAT_DNKm.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_PLAT_DNKm.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_PLAT_DNKm.boot.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_PLAT_DNKm.boot.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))

for (j in 1:(omega-a.min)){
  for (k in 1:y.pred){
    q_PLAT_DNKm.q95[j,k] <- quantile(q_PLAT_DNKm.st.ext[j,k,],
                                     probs=0.5+conf.lev/2)
    q_PLAT_DNKm.q05[j,k] <- quantile(q_PLAT_DNKm.st.ext[j,k,],
                                     probs=(1-conf.lev)/2)
    q_PLAT_DNKm.boot.q95[j,k] <- quantile(q_PLAT_DNKm.boot.st.ext[j,k,],
                                          probs=0.5+conf.lev/2)
    q_PLAT_DNKm.boot.q05[j,k] <- quantile(q_PLAT_DNKm.boot.st.ext[j,k,],
                                          probs=(1-conf.lev)/2)
  }
}

q_PLAT_DNKf.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_PLAT_DNKf.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_PLAT_DNKf.boot.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_PLAT_DNKf.boot.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
for (j in 1:(omega-a.min)){
  for (k in 1:y.pred){
    q_PLAT_DNKf.q95[j,k] <- quantile(q_PLAT_DNKf.st.ext[j,k,], probs=0.5+conf.lev/2)
    q_PLAT_DNKf.q05[j,k] <- quantile(q_PLAT_DNKf.st.ext[j,k,], probs=(1-conf.lev)/2)
    q_PLAT_DNKf.boot.q95[j,k] <- quantile(q_PLAT_DNKf.boot.st.ext[j,k,], probs=0.5+conf.lev/2)
    q_PLAT_DNKf.boot.q05[j,k] <- quantile(q_PLAT_DNKf.boot.st.ext[j,k,], probs=(1-conf.lev)/2)
  }
}



### Actuarial quantities {.tabset}


#### death probabilities at 65


# data preparation
death_probs_ci_m <- tibble(q95 = q_PLAT_DNKm.q95["65",],
                           q05 = q_PLAT_DNKm.q05["65",],
                           q95boot = q_PLAT_DNKm.boot.q95["65",],
                           q05boot = q_PLAT_DNKm.boot.q05["65",],
                           years =seq(2019,2118),
)


death_probs_ci_m <- death_probs_ci_m %>%  gather("value", "rates", 1:4)

death_prob_m <- tibble(years = seq(1968,2118),
                       value = "qx",
                       rates = q_PLAT_DNKm.ext["65",]
)

death_probs_ci_f <- tibble(q95 = q_PLAT_DNKf.q95["65",],
                           q05 = q_PLAT_DNKf.q05["65",],
                           q95boot = q_PLAT_DNKf.boot.q95["65",],
                           q05boot = q_PLAT_DNKf.boot.q05["65",],
                           years =seq(2019,2118),
)


death_probs_ci_f <- death_probs_ci_f %>%  gather("value", "rates", 1:4)

death_prob_f <- tibble(years = seq(1968,2118),
                       value = "qx",
                       rates = q_PLAT_DNKf.ext["65",]
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
death_probs_ci_m <- tibble(q95 = q_PLAT_DNKm.q95["20",],
                           q05 = q_PLAT_DNKm.q05["20",],
                           q95boot = q_PLAT_DNKm.boot.q95["20",],
                           q05boot = q_PLAT_DNKm.boot.q05["20",],
                           years =seq(2019,2118),
)


death_probs_ci_m <- death_probs_ci_m %>%  gather("value", "rates", 1:4)

death_prob_m <- tibble(years = seq(1968,2118),
                       value = "qx",
                       rates = q_PLAT_DNKm.ext["20",]
)

death_probs_ci_f <- tibble(q95 = q_PLAT_DNKf.q95["20",],
                           q05 = q_PLAT_DNKf.q05["20",],
                           q95boot = q_PLAT_DNKf.boot.q95["20",],
                           q05boot = q_PLAT_DNKf.boot.q05["20",],
                           years =seq(2019,2118),
)


death_probs_ci_f <- death_probs_ci_f %>%  gather("value", "rates", 1:4)

death_prob_f <- tibble(years = seq(1968,2118),
                       value = "qx",
                       rates = q_PLAT_DNKf.ext["20",]
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



p_PLAT_DNKm.st.ext <- 1-q_PLAT_DNKm.st.ext # 1 year survival probability
ann_PLAT_DNKm.st <- annuity.st(q_PLAT_DNKm.st.ext,v_vect,0.9)
ann_PLAT_DNKm.mean <- ann_PLAT_DNKm.st[[1]]
ann_PLAT_DNKm.q95 <- ann_PLAT_DNKm.st[[2]]
ann_PLAT_DNKm.q05 <- ann_PLAT_DNKm.st[[3]]

p_PLAT_DNKf.st.ext <- 1-q_PLAT_DNKf.st.ext # 1 year survival probability
ann_PLAT_DNKf.st <- annuity.st(q_PLAT_DNKf.st.ext,v_vect,0.9)
ann_PLAT_DNKf.mean <- ann_PLAT_DNKf.st[[1]]
ann_PLAT_DNKf.q95 <- ann_PLAT_DNKf.st[[2]]
ann_PLAT_DNKf.q05 <- ann_PLAT_DNKf.st[[3]]

# data preparation

annuity_male <- tibble(annuity_mean = ann_PLAT_DNKm.mean["65",],
                       q05 = ann_PLAT_DNKm.q05["65",],
                       q95 = ann_PLAT_DNKm.q95["65",],
                       years = seq(2019,2118))

annuity_male <- annuity_male %>% gather("lab", "value", 1:3)

annuity_female <- tibble(annuity_mean = ann_PLAT_DNKf.mean["65",],
                         q05 = ann_PLAT_DNKf.q05["65",],
                         q95 = ann_PLAT_DNKf.q95["65",],
                         years = seq(2019,2118))

annuity_female <- annuity_female %>% gather("lab", "value", 1:3)


ann_2050_m <- tibble(annuity_mean = ann_PLAT_DNKm.mean[,"2050"],
                     q05 = ann_PLAT_DNKm.q05[,"2050"],
                     q95 = ann_PLAT_DNKm.q95[,"2050"],
                     ages = seq(0,109))

annuity_2050_m <- ann_2050_m %>% gather("lab", "value", 1:3)

ann_2050_f <- tibble(annuity_mean = ann_PLAT_DNKf.mean[,"2050"],
                     q05 = ann_PLAT_DNKf.q05[,"2050"],
                     q95 = ann_PLAT_DNKf.q95[,"2050"],
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

