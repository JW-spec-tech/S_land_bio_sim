library(ggformula)
library(mgcv)
library(tibble)
library(patch)
library(dvmisc)
library(dplyr)

dat = tibble( x= runif(500), y= runif(500),
              biomass = rTweedie(mu = 10*dnorm(x-0.5,sd = 0.5)*dnorm(y-0.5,sd = 0.5),p = 1.5,phi = 1.1))
gam_mod <- gam(biomass~s(x,y,bs = "gp"), family = tw, data=dat)

dat_grid <- as.data.frame(expand_grid(x=seq(0,1,length=100), y= seq(0,1,length=100)))%>%
            dplyr::mutate(fit  = predict.gam(gam_mod,type = "response",newdata = .),
                          true = 10*dnorm(x-0.5,sd = 0.5)*dnorm(y-0.5,sd = 0.5))

dat = tibble( x= runif(500), y= runif(500),
              biomass = rTweedie(mu = 10*dnorm(x-0.5,sd = 0.5)*dnorm(y-0.5,sd = 0.5),p = 1.5,phi = 1.1))
gam_mod <- gam(biomass~s(x,y,bs = "gp"), family = tw, data=dat)



sum(dat_grid$fit)
sum(dat_grid$true)



sum(dat_grid$fit)
sum(dat_grid$true)

library(patchwork)
fit_plot <- gf_raster(lat~long, data= dat_grid, fill = ~fit_simple_gam)+ scale_fill_viridis_c()
true_plot <- gf_raster(lat~long, data= dat_grid, fill = ~fit_spatialonly_gam)+ scale_fill_viridis_c()


true_plot + fit_plot

fit_coef <- coef(simple_gam)
fit_VCV <- simple_gam$Vp
fit_lp_mat <- predict(simple_gam, newdata = dat_grid, type = "lpmatrix")

fit_coef_sims <- rmvn(500,mu =fit_coef, V =  fit_VCV)

fit_function_sims <- fit_lp_mat %*% t(fit_coef_sims)
sim_biomasses <- colSums(exp(fit_function_sims))
biomass_ci <- quantile(sim_biomasses,probs =  c(0.025, 0.975))
