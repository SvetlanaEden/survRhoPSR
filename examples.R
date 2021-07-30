#source("sourcesFunctions.R")
source("/Users/svetlanaeden/stuff/RStuff/survRhoPSR/sourcesFunctions.R")

library(survival)
data(diabetes, package = "SurvCorr")

### We call SurvCorr package to access the diabetic retinopathy data. The first five rows of the dataset are printed below.
diabetes[1:5,]

### Estimator $\rho_{PSR}$ provides a way to approximate the overall Spearman's correlation. The function unadjusted.CorPSRs() computes $\widehat{\rho}_{PSR}$ and its $95\%$ confidence interval:
res = unadjusted.CorPSRs(X = diabetes$TIME1, Y = diabetes$TIME2,
                         deltaX = diabetes$STATUS1, deltaY = diabetes$STATUS2)
round(res[ c("est", "lower.CI", "upper.CI") ], 3)
   
### The correlation can be affected by different factors, e.g., age, what eye was treated, and treatment type. We can adjust for these other variables by estimating the partial correlation, $\widehat{\rho}_{PSR\cdot Z}$. To compute the partial correlation, we first need to fit separate models for both of the time to event outcomes on the covariates. Here we fit Cox models:
survObjX = Surv(diabetes$TIME1, diabetes$STATUS1)
survObjY = Surv(diabetes$TIME2, diabetes$STATUS2)
modX = coxph(survObjX ~ TRT_EYE + AGE_DX + LASER,     
               data=diabetes, method = "breslow", timefix = FALSE)
modY = coxph(survObjY ~ TRT_EYE + AGE_DX + LASER,
              data=diabetes, method = "breslow", timefix = FALSE)

### The function partial.corPSRs() then takes these model objects as input and computes the partial correlation estimate, $\widehat{\rho}_{PSR\cdot \pmb{Z}}$. In this example, the correlation slightly increased after adjusting for age, treatment eye, and type of treatment:
round(partial.corPSRs(modX, modY)[ c("est", "lower.CI", "upper.CI") ], 3)

### When the survival probabilities are modeled using Cox proportional hazards model, the confidence interval of $\widehat{\rho}_{PSR\cdot \pmb{Z}}$ is computed using score equations from the Cox regression partial likelihood. Using partial likelihood results in slightly underestimated variability of $\rho_{PSR}$ because it does not take into account the variability of the baseline hazard. The user can choose the full likelihood option instead:
round(partial.corPSRs(modX, modY, likelihood = "full")[c("est", "lower.CI", "upper.CI") ], 3)

### The confidence intervals obtained from the full and partial likelihoods are almost the same. Note that Cox regression should be fit with timefix = FALSE because of issues related to the floating point round-off error of the time to event stored in the Cox regression object.

### Instead of Cox proportional hazards, parametric survival models can be used, for example, the log-logistic model:
modX = survreg(survObjX ~ AGE_DX, data=diabetes, dist = "loglogistic")
modY = survreg(survObjY ~ AGE_DX, data=diabetes, dist = "loglogistic")
round(partial.corPSRs(modX, modY)[ c("est", "lower.CI", "upper.CI") ], 3)

### The correlation computed using log-logistic model is very similar.

### In addition to the partial correlation, there is an interest in computing correlation conditional on other variables, $\rho_{PSR|Z}$. Suppose we are interested in estimating the rank correlation conditional on age at diagnosis. The function conditional.corPSRs() can be used to estimate $\rho_{PSR|Z}$. The function requires inputting models of the times-to-event conditional on age; Cox models are a natural choice and are implemented below. Since age is a continuous variable, conditional.corPSRs() allows the user to specify how to model the rank correlation. Specifically, $\rho_{PSR|Z}$ can be modeled linearly (numKnots = 0) or using restricted cubic splines with mumKnots set to a whole number between $3$ and $7$. The location of the spline knots is defined in terms of quantiles suggested by Harrell (see for example, package rms). The code below estimates and plots $\rho_{PSR|Z}$ with restricted cubic splines with three knots.
modXC = coxph(survObjX ~ AGE_DX, data=diabetes, method = "breslow", timefix = FALSE)
modYC = coxph(survObjY ~ AGE_DX, data=diabetes, method = "breslow", timefix = FALSE)
z = diabetes[["AGE_DX"]]
newZ = diabetes[["AGE_DX"]]
par(mfrow = c(1, 2))
plotCol = "#22222222"
for(n_knots in c(0, 3)){
  resultXY = conditional.corPSRs(modXC, modYC, z, newZ, numKnots = n_knots)
  plotData = data.frame(x = newZ, y = resultXY$est, 
   yLower = resultXY$lower.CI, yUpper = resultXY$upper.CI)
  plotData = plotData[order(plotData$x),]
  plot(plotData$x, plotData$y, type = "n", ylim = c(0, 1))
  points(plotData$x, plotData$y, pch = 19, col = plotCol, cex = .3)  
  lines(plotData$x, plotData$yLower, col = plotCol)
  lines(plotData$x, plotData$yUpper, col = plotCol)
  abline(h=0, lty = 3, col = plotCol)
}

### Note that the current version of function conditional.corPSRs() computes correlation conditional only on one variable.
### Note also that argument z contains variable $Z$, and newZ contains only those values of $Z$, for which $\rho_{PSR | Z}$ is computed. The plot below shows $\widehat{\rho}_{PSR | Z}$ computed as a function of age, where age is modeled as a linear (left panel) and quadratic (right panel) variable. The panels generally show the rank correlation increasing with age, suggesting that the association between times to retinopathy for different eyes in the same patient is likely stronger for older patients. Larger sample sizes are needed to determine the functional form of the correlation with greater precision.\\

### The same function conditional.corPSRs() can also compute the partial-conditional correlation. We may want to estimate the rank correlation conditional on age after adjusting for which eye was treated TRT_EYE and the type of treatment LASER. This can be estimated using the same code except inputting objects from models that include TRT_EYE and LASER covariates:
plotCol = "#44444444"
modXC = coxph(survObjX ~ TRT_EYE + AGE_DX + LASER, data=diabetes,
                  method = "breslow", timefix = FALSE)
modYC = coxph(survObjY ~ TRT_EYE + AGE_DX + LASER, data=diabetes,
                  method = "breslow", timefix = FALSE)
z = diabetes[["AGE_DX"]]
newZ = diabetes[["AGE_DX"]]
par(mfrow = c(1, 2))
for(n_knots in c(0, 3)){
      resultXY = conditional.corPSRs(modXC, modYC, z, newZ, numKnots = n_knots)
      plotData = data.frame(x = newZ, y = resultXY$est, 
        yLower = resultXY$lower.CI, yUpper = resultXY$upper.CI)
      plotData = plotData[order(plotData$x),]
      plot(plotData$x, plotData$y, type = "n", ylim = c(0, 1))
      points(plotData$x, plotData$y, pch = 19, col = plotCol, cex = .3)  
      lines(plotData$x, plotData$yLower, col = plotCol)
      lines(plotData$x, plotData$yUpper, col = plotCol)
      abline(h=0, lty = 3, col = plotCol)
}

### There is not much difference between the correlation conditional on age computed earlier and the partial correlation conditional on age.

### To check if the correlation is the same, whether the right eye or left eye is treated or for the two types of treatment, we compute the partial rank correlation adjusted for age and conditional on the treated eye and treatment type.
par(mfrow = c(1, 2))
nameList = list(TRT_EYE = c("Right", "Left"), LASER = c("Xenon", "Argon"))
xList = list(TRT_EYE = 2:1, LASER = 1:2)
newZ = c(1, 2)
for(var_i in c("TRT_EYE", "LASER")){
  z = diabetes[[var_i]]
  resultXY = conditional.corPSRs(modXC, modYC, z, newZ, numKnots = 0)
  resultXY[[var_i]] = newZ
  resultXY[["NAMES"]] = nameList[[var_i]][newZ]
  plot(0, 0, type = "n", xlim = c(0, 3), ylim = range(-0.2, 1),
      axes = FALSE, xlab = "", ylab = "")
  abline(h = 0, col = "gray", lty = 3)
  box()
  points(xList[[var_i]], resultXY[, "est"], pch = 18)
  segments(x0 = xList[[var_i]], y0 = resultXY$lower.CI,
          x1 = xList[[var_i]], y1 = resultXY$upper.CI)
  axis(side = c(2), at = seq(-0.2, 1, .2))
  axis(side = c(4), at = seq(-0.2, 1, .2))
  mtext(resultXY$NAMES, side = 1, line = 1, at = xList[[var_i]])
}

### The plot above shows estimates of the partial rank correlation conditional on eye (left panel) and treatment (right panel). Interestingly, it appears that the correlation between the times to retinopathy for different eyes within the same patient was higher when the left eye was treated than when the right eye was treated. In contrast, there appears to be no difference in correlation between the two treatment groups.
###
