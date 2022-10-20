## Load survival package
library(survival)
## List datasets in survival package
data(package = "survival")

## Load lung data
data(lung)

## Show first 6 rows
head(lung)

## Add survival object
lung$SurvObj <- with(lung, Surv(time, status == 2))

## Check data
head(lung)

## Kaplan-Meier estimator. The "log-log" confidence interval is preferred.
km.as.one <- survfit(SurvObj ~ 1, data = lung, conf.type = "log-log")
km.by.sex <- survfit(SurvObj ~ sex, data = lung, conf.type = "log-log")

## Show object
km.as.one

plot(km.as.one)

plot(km.as.one, conf = F, mark.time = F)
