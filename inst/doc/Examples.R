## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
set.seed(17)  # For reproducibility's sake

## ------------------------------------------------------------------------
library(mstate)
data(ebmt3)
head(ebmt3)

## ------------------------------------------------------------------------
tmat <- trans.illdeath()
tmat

## ------------------------------------------------------------------------
long <- msprep(time=c(NA, 'prtime', 'rfstime'), 
               status=c(NA, 'prstat', 'rfsstat'), 
               data=ebmt3, 
               trans=tmat, 
               keep=c('age', 'dissub'))
head(long)

## ------------------------------------------------------------------------
library(flexsurv)
models <- lapply(1:3, function(i) {
    flexsurvreg(Surv(time, status) ~ age + dissub, data=long, dist='weibull')
})

## ------------------------------------------------------------------------
newdata <- data.frame(age="20-40", dissub="AML")

## ----example1------------------------------------------------------------
library(multistateutils)
predict_transitions(models, newdata, tmat, times=365)

## ----example2------------------------------------------------------------
pmatrix.simfs(models, tmat, newdata=newdata, t=365)

## ----example3------------------------------------------------------------
predict_transitions(models, newdata, tmat, times=365, ci=TRUE, M=10)

## ----example4------------------------------------------------------------
pmatrix.simfs(models, tmat, newdata=newdata, t=365, ci=TRUE, M=9)

## ------------------------------------------------------------------------
library(microbenchmark)
microbenchmark("multistateutils"=predict_transitions(models, newdata, tmat, times=365),
               "flexsurv"=pmatrix.simfs(models, tmat, newdata=newdata, t=365), times=10)

## ----example5------------------------------------------------------------
predict_transitions(models, newdata, tmat, times=seq(9)*365)

## ----example6------------------------------------------------------------
do.call('rbind', lapply(seq(9)*365, function(t) {
    pmatrix.simfs(models, tmat, newdata=newdata, t=t)
}))

## ----benchmarkmultipletimes----------------------------------------------
microbenchmark("multistateutils"=predict_transitions(models, newdata, tmat, times=seq(9)*365),
               "flexsurv"={do.call('rbind', lapply(seq(9)*365, function(t) {
                            pmatrix.simfs(models, tmat, newdata=newdata, t=t)}))
               }, times=10)

## ------------------------------------------------------------------------
predict_transitions(models, newdata, tmat, times=365, start_times = 365/2)

## ------------------------------------------------------------------------
predict_transitions(models, newdata, tmat, times=365, 
                    start_times = c(0.25, 0.5, 0.75) * 365)

## ------------------------------------------------------------------------
predict_transitions(models, newdata, tmat, times=seq(2)*365, 
                    start_times = c(0.25, 0.5, 0.75) * 365)

## ------------------------------------------------------------------------
microbenchmark("time"=predict_transitions(models, newdata, tmat, 
                                          times=seq(2)*365, 
                                          start_times = c(0.25, 0.5, 0.75)*365),
               times=10)

## ------------------------------------------------------------------------
newdata_multi <- data.frame(age=c("20-40", ">40"), dissub=c("AML", "CML"))

## ----exampleinds1--------------------------------------------------------
predict_transitions(models, newdata_multi, tmat, times=365)

## ----exampleinds2, error=T, eval=T---------------------------------------
pmatrix.simfs(models, tmat, newdata=newdata_multi, t=365)

## ----exampleinds3--------------------------------------------------------
do.call('rbind', lapply(seq(nrow(newdata_multi)), function(i) {
    pmatrix.simfs(models, tmat, newdata=newdata_multi[i, ], t=365)
}))

## ------------------------------------------------------------------------
models_arrival <- lapply(1:3, function(i) {
    if (i == 3) {
        flexsurvreg(Surv(time, status) ~ age + dissub + Tstart, data=long, dist='weibull')
    } else {
        
        flexsurvreg(Surv(time, status) ~ age + dissub, data=long, dist='weibull')
    }
})

## ------------------------------------------------------------------------
models_arrival[[3]]

## ------------------------------------------------------------------------
newdata_arrival <- data.frame(age="20-40", dissub="AML", Tstart=0)

## ------------------------------------------------------------------------
predict_transitions(models_arrival, newdata_arrival, tmat, times=365, tcovs='Tstart')

## ------------------------------------------------------------------------
predict_transitions(models_arrival, newdata_arrival, tmat, times=365)

## ----tcovs1--------------------------------------------------------------
pmatrix.simfs(models_arrival, tmat, newdata=newdata_arrival, t=365, tcovs='Tstart')

## ----tcovs2--------------------------------------------------------------
pmatrix.simfs(models_arrival, tmat, newdata=newdata_arrival, t=365)

## ------------------------------------------------------------------------
models_mix <- lapply(1:3, function(i) {
    if (i == 1) {
        flexsurvreg(Surv(time, status) ~ age + dissub, data=long, dist='weibull')
    } else if (i == 2) {
        flexsurvreg(Surv(time, status) ~ age + dissub, data=long, dist='exp')
    } else {
        flexsurvreg(Surv(time, status) ~ age + dissub, data=long, dist='lnorm')
    }
})

## ------------------------------------------------------------------------
predict_transitions(models_mix, newdata, tmat, times=365)

## ------------------------------------------------------------------------
pmatrix.simfs(models_mix, tmat, newdata=newdata, t=365)

## ------------------------------------------------------------------------
length_of_stay(models, 
               newdata=newdata,
               tmat, times=365.25*3,
               start=1)

## ------------------------------------------------------------------------
totlos.simfs(models, tmat, t=365.25*3, start=1, newdata=newdata)

## ------------------------------------------------------------------------
length_of_stay(models, 
               newdata=data.frame(age=c(">40", ">40"),
                                  dissub=c('CML', 'AML')),
               tmat, times=c(1, 3, 5)*365.25,
               start=c('healthy', 'illness'))

## ------------------------------------------------------------------------
time_points <- seq(0, 10, by=2) * 365.25
plot_predicted_pathway(models, tmat, newdata, time_points, 1)

