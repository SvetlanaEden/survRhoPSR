
###################################################################################
##### prodTS.new.SKE() computes the standard error of conditional estimate
#####    using large sample approximation, author: Qi Liu and Svetlana Eden
#####    Qi originally wrote this function. Svetlana slightly modified it 
#####    to use splines with z and to get predicted values and CIs right away
###################################################################################

prodTS.new.SKE = function(x.resid, y.resid, z, numKnots = 3, newZ = z,
                       xres2.method=c("emp", "constant", "model"),
                       yres2.method=c("emp", "constant", "model"),
                       xz.dl.dtheta, yz.dl.dtheta,
                       xz.d2l.dtheta.dtheta, yz.d2l.dtheta.dtheta,
                       dxresid.dtheta, dyresid.dtheta, debugMode = TRUE, quadraticTerms = FALSE){

  newZ = as.numeric(newZ)
  ### Define the design matrix for z - variable of interest
  ###                         and newZ - values of z at wich we want conditional corr.
  if (!any(numKnots == c(-1, 0, 3, 4, 5, 6, 7))) stop("The number of knots should be either 0, 3, 4, or 5")
  if(numKnots == -1){
    designMatrOfZ = matrixForZ = matrix(1, ncol = 1, nrow = length(z))
    newDesignMatrOfZ = matrixForNewZ = matrix(1, ncol = 1, nrow = length(newZ))
  }else{
    if(numKnots == 0){
      ### perhaps you want to change it to numKnots <= 2 b/c you cannot have 1 or 2
      matrixForZ = matrix(z, ncol = 1)
      matrixForNewZ = matrix(newZ, ncol = 1)
    }else{
      splinePercentileList = list("3" = c(.1, .5, .9), "4" = c(.05, .35, .65, .95), "5" = c(.05, .275, .5, .725, .95), "6" = c(.05, .23, .41, .59, .77, .95), "7" = c(.05, .1833, .3417, .5, .6583, .8167, .975)) ### see Regression Modeling Strategies by F. Harrell.
      knots = quantile(z, splinePercentileList[[as.character(numKnots)]], na.rm = TRUE)
      matrixForZ = componentsOfSplines(z, knots)
      matrixForNewZ = componentsOfSplines(newZ, knots)
    }
    designMatrOfZ <- cbind(1, matrixForZ)
    newDesignMatrOfZ <- cbind(1, matrixForNewZ)
    
    ### Bryan wanted a quadratic term:
    ### This is a weird patch, but that will for now
    if(quadraticTerms & numKnots == 3){
      matrixForZ = cbind(z, z^2)
      designMatrOfZ <- matrix(c(rep(1, length(z)), matrixForZ), ncol = 3)
      newDesignMatrOfZ <- matrix(c(rep(1, length(newZ)), newZ, newZ^2), ncol = 3)
    }
  }
  
  npar.xz <- dim(xz.dl.dtheta)[2]
  npar.yz <- dim(yz.dl.dtheta)[2]

  score.prod <- lm.scores(x.resid*y.resid, matrixForZ)

  npar.prod <- dim(score.prod$dl.dtheta)[2]
  prod.dl.dtheta <- score.prod$dl.dtheta
  prod.d2l.dtheta.dtheta <- score.prod$d2l.dtheta.dtheta
  N <- dim(score.prod$dl.dtheta)[1]
  
  if(xres2.method[1]=="model"){
    score.xres <- lm.scores(x.resid, matrixForZ)
    score.xres2 <- lm.scores(x.resid^2, matrixForZ)
    npar.xres2 <- dim(score.xres2$dl.dtheta)[2]
    xres2.dl.dtheta <- score.xres2$dl.dtheta
    xres2.d2l.dtheta.dtheta <- score.xres2$d2l.dtheta.dtheta
    
  } else if (xres2.method[1]=="emp"){
    xres2.dl.dtheta <- as.matrix(x.resid^2-mean(x.resid^2) )
    npar.xres2 <- 1
    xres2.d2l.dtheta.dtheta <- -N
    
  } else if (xres2.method[1]=="constant"){
    npar.xres2 <- 0
    xres2.dl.dtheta <- NULL
    xres2.d2l.dtheta.dtheta <- NULL
  }
  
  if(yres2.method[1]=="model"){
    score.yres <- lm.scores(y.resid, matrixForZ)
    score.yres2 <- lm.scores(y.resid^2, matrixForZ)
    npar.yres2 <- dim(score.yres2$dl.dtheta)[2]
    yres2.dl.dtheta <- score.yres2$dl.dtheta
    yres2.d2l.dtheta.dtheta <- score.yres2$d2l.dtheta.dtheta
    
  } else if (yres2.method[1]=="emp"){
    yres2.dl.dtheta <- as.matrix( y.resid^2-mean(y.resid^2)  )
    npar.yres2 <- 1
    yres2.d2l.dtheta.dtheta <- -N
    
  } else if (yres2.method[1]=="constant"){
    npar.yres2 <- 0
    yres2.dl.dtheta <- NULL
    yres2.d2l.dtheta.dtheta <- NULL
  }
  
  Ntheta <- npar.xz + npar.yz + npar.prod + npar.xres2 + npar.yres2
  bigphi <- cbind(xz.dl.dtheta,
                  yz.dl.dtheta,
                  prod.dl.dtheta,
                  xres2.dl.dtheta,
                  yres2.dl.dtheta)
  
  A <- matrix(0, Ntheta, Ntheta)
  ### initiate the diagonal of A 
  A[(1:npar.xz), (1:npar.xz)] <- as.matrix(xz.d2l.dtheta.dtheta)
  A[npar.xz + (1:npar.yz), npar.xz + (1:npar.yz)] <- as.matrix(yz.d2l.dtheta.dtheta)
  A[npar.xz + npar.yz + (1: npar.prod), npar.xz + npar.yz + (1: npar.prod)] <- prod.d2l.dtheta.dtheta
  if (npar.xres2>0){
    A[npar.xz + npar.yz + npar.prod+ (1: npar.xres2), npar.xz + npar.yz + npar.prod+ (1: npar.xres2)] <- xres2.d2l.dtheta.dtheta
  }
  if(npar.yres2>0){
    A[npar.xz + npar.yz + npar.prod+ npar.xres2+(1: npar.yres2), npar.xz + npar.yz + npar.prod+ npar.xres2+(1: npar.yres2)] <- yres2.d2l.dtheta.dtheta
  }
  par.1 <- t(y.resid * designMatrOfZ) %*% t(dxresid.dtheta)
  par.2 <- t(x.resid * designMatrOfZ) %*% t(dyresid.dtheta)
  A[npar.xz + npar.yz + (1: npar.prod), 1:npar.xz] <- par.1
  A[npar.xz + npar.yz + (1: npar.prod), npar.xz + (1:npar.yz)] <- par.2
  
  if (xres2.method[1]=="model"){
    par.3 <- t(2*x.resid * designMatrOfZ) %*% t(dxresid.dtheta)
  } else if (xres2.method[1]=="emp"){
    par.3 <- t(2*x.resid ) %*% t(dxresid.dtheta)
  } else if(xres2.method[1]=="constant"){
    par.3 <- NULL
  }
  
  if(npar.xres2>0){
    A[npar.xz + npar.yz + npar.prod+ (1: npar.xres2), 1:npar.xz] <- par.3
  }
  
  if (yres2.method[1]=="model"){
    par.4 <- t(2*y.resid * designMatrOfZ) %*% t(dyresid.dtheta)
  } else if(yres2.method[1]=="emp"){
    par.4 <- t(2*y.resid ) %*% t(dyresid.dtheta)
  } else if(yres2.method[1]=="constant"){
    par.4 <- NULL
  }
  
  if(npar.yres2>0){
    A[npar.xz + npar.yz + npar.prod+ npar.xres2+(1:npar.yres2), npar.xz + (1:npar.yz)] <- par.4
  }
  
  if(debugMode){
    # A = diag(1, nrow(A))
    require(matlib)
    matrB = t(bigphi) %*% bigphi
    SS = Ginv(A) %*% matrB %*% t(Ginv(A))
    warning("Used generalized inverse")
  }else{
    SS <- solve(A, t(bigphi))
  }
  var.theta <- tcrossprod(SS, SS)

  prod.coef <- score.prod$mod$coef
  names(prod.coef) <- paste("prod:", names(prod.coef))
  if(xres2.method[1]=="model"){
    xres2.coef <- score.xres2$mod$coef
  } else if (xres2.method[1]=="emp"){
    xres2.coef <- mean(x.resid^2)
  } else if (xres2.method[1]=="constant"){
    xres2.coef <- 1/3
  }
  names(xres2.coef) <- paste("xres2:", names(xres2.coef))
  
  if(yres2.method[1]=="model"){
    yres2.coef <- score.yres2$mod$coef
  } else if (yres2.method[1]=="emp"){
    yres2.coef <- mean(y.resid^2)
  } else if (yres2.method[1]=="constant"){
    yres2.coef <- 1/3
  }
  names(yres2.coef) <- paste("yres2:", names(yres2.coef))
  coef <- rbind(prod.coef, xres2.coef, yres2.coef)
  
  var.coef <- var.theta[(npar.xz+npar.yz+1): Ntheta , (npar.xz+npar.yz+1): Ntheta]
  
  ####### compute predicted values
  # regResult = ols(x.resid*y.resid ~ rcs(z, numKnots))
  # predProd = model.matrix(regResult) %*% regResult$coefficients
  #score.prod$mod
  if(FALSE){
    regYX = ols(scoreRes1$presid*scoreRes2$presid ~ rcs(z, numKnots))
    regX2 = ols((scoreRes1$presid)^2 ~ rcs(z, numKnots))
    regY2 = ols((scoreRes2$presid)^2 ~ rcs(z, numKnots))
    pXY = model.matrix(regXY) %*% regXY$coefficients
    pX2 = model.matrix(regX2) %*% regX2$coefficients
    pY2 = model.matrix(regY2) %*% regY2$coefficients
    pRho =  pXY/sqrt(pX2*pY2)
  }
  predYX = matrix(prod.coef, nrow = 1) %*% t(newDesignMatrOfZ)
  predX2 = matrix(xres2.coef, nrow = 1) %*% t(newDesignMatrOfZ)
  predY2 = matrix(yres2.coef, nrow = 1) %*% t(newDesignMatrOfZ)
  predRho = predYX/sqrt(predX2*predY2)
  
  ####### compute derivatives for delta method:
  ### npar.prod + npar.xres2 + npar.yres2
  derivForDeltaMethod = cbind(newDesignMatrOfZ, newDesignMatrOfZ, newDesignMatrOfZ)
  ### divide by the denominator
  derivForDeltaMethod = derivForDeltaMethod/matrix(rep(sqrt(predX2 * predY2), npar.prod + npar.xres2 + npar.yres2), nrow = nrow(derivForDeltaMethod), ncol = npar.prod + npar.xres2 + npar.yres2)

  ### multiply by -(1/2) and by predYX
  derivForDeltaMethod[, npar.prod + (1 : (npar.xres2 + npar.yres2))] = -(1/2)*derivForDeltaMethod[, npar.prod + (1 : (npar.xres2 + npar.yres2))] * matrix(rep(predYX, npar.xres2 + npar.yres2), nrow = nrow(derivForDeltaMethod), ncol = npar.xres2 + npar.yres2)
  
  ### divide by predX2
  derivForDeltaMethod[, npar.prod + (1 : npar.xres2)] = derivForDeltaMethod[, npar.prod + (1 : npar.xres2)] / matrix(rep(predX2, npar.xres2), nrow = nrow(derivForDeltaMethod), ncol = npar.xres2) 
  
  ### divide by predY2
  derivForDeltaMethod[, npar.prod + npar.xres2 + (1 : npar.yres2)] = derivForDeltaMethod[, npar.prod + npar.xres2 + (1 : npar.yres2)] / matrix(rep(predY2, npar.yres2), nrow = nrow(derivForDeltaMethod), ncol = npar.yres2) 
  
  resultingVar = diag(derivForDeltaMethod %*% var.coef %*% t(derivForDeltaMethod))
  
  ci = 0.95
  pointEstAndCIs = data.frame(pointEst = as.vector(predRho), lower = as.vector(predRho - abs(qnorm(0.5*(1-ci)))*sqrt(resultingVar)), upper = as.vector(predRho + abs(qnorm(0.5*(1-ci)))*sqrt(resultingVar)))
    
  result <- list(prod=prod.coef,
                 xres2=xres2.coef,
                 yres2=yres2.coef,
                 var.coef=var.coef,
                 pointEstAndCIs = pointEstAndCIs)
                 
  return(result)
}

##############################################################################
##### similar to prodTS.new.SKE(), but unlike it, this function uses the full formula 
##### for the correlation
##############################################################################
prodTS.new.SKE.full = function(x.resid, y.resid, z, numKnots = 3, newZ = z,
                       xres2.method=c("emp", "constant", "model"),
                       yres2.method=c("emp", "constant", "model"),
                       xz.dl.dtheta, yz.dl.dtheta,
                       xz.d2l.dtheta.dtheta, yz.d2l.dtheta.dtheta,
                       dxresid.dtheta, dyresid.dtheta, debugMode = FALSE){
                         
                         ### argument debugMode is not in use

  if (!any(numKnots == c(-1, 0, 3, 4, 5, 6, 7))) stop("The number of knots should be either 0, 3, 4, or 5")
  if(numKnots == -1){
    designMatrOfZ = matrixForZ = matrix(1, ncol = 1, nrow = length(z))
    newDesignMatrOfZ = matrixForNewZ = matrix(1, ncol = 1, nrow = length(newZ))
  }else{
    if(numKnots == 0){
      matrixForZ = matrix(z, ncol = 1)
      matrixForNewZ = matrix(newZ, ncol = 1)
    }else{
      splinePercentileList = list("3" = c(.1, .5, .9), "4" = c(.05, .35, .65, .95), "5" = c(.05, .275, .5, .725, .95), "6" = c(.05, .23, .41, .59, .77, .95), "7" = c(.05, .1833, .3417, .5, .6583, .8167, .975)) ### see Regression Modeling Strategies by F. Harrell.
      knots = quantile(z, splinePercentileList[[as.character(numKnots)]], na.rm = TRUE)
      matrixForZ = componentsOfSplines(z, knots)
      matrixForNewZ = componentsOfSplines(newZ, knots)
    }
    designMatrOfZ <- cbind(1, matrixForZ)
    newDesignMatrOfZ <- cbind(1, matrixForNewZ)
  }
  xres.method = "model"
  yres.method = "model"

  npar.xz <- dim(xz.dl.dtheta)[2]
  npar.yz <- dim(yz.dl.dtheta)[2]
  
  score.prod <- lm.scores(x.resid*y.resid, matrixForZ)

  npar.prod <- dim(score.prod$dl.dtheta)[2]
  prod.dl.dtheta <- score.prod$dl.dtheta
  prod.d2l.dtheta.dtheta <- score.prod$d2l.dtheta.dtheta
  N <- dim(score.prod$dl.dtheta)[1]
  
  if(xres2.method[1]=="model"){
    score.xres <- lm.scores(x.resid, matrixForZ)
    score.xres2 <- lm.scores(x.resid^2, matrixForZ)

    npar.xres <- dim(score.xres$dl.dtheta)[2]
    xres.dl.dtheta <- score.xres$dl.dtheta
    xres.d2l.dtheta.dtheta <- score.xres$d2l.dtheta.dtheta
    npar.xres2 <- dim(score.xres2$dl.dtheta)[2]
    xres2.dl.dtheta <- score.xres2$dl.dtheta
    xres2.d2l.dtheta.dtheta <- score.xres2$d2l.dtheta.dtheta
    
  } else if (xres2.method[1]=="emp"){
    xres2.dl.dtheta <- as.matrix(x.resid^2-mean(x.resid^2) )
    npar.xres2 <- 1
    xres2.d2l.dtheta.dtheta <- -N
    
  } else if (xres2.method[1]=="constant"){
    npar.xres2 <- 0
    xres2.dl.dtheta <- NULL
    xres2.d2l.dtheta.dtheta <- NULL
  }
  
  if(yres2.method[1]=="model"){
    score.yres <- lm.scores(y.resid, matrixForZ)
    score.yres2 <- lm.scores(y.resid^2, matrixForZ)
    npar.yres <- dim(score.yres$dl.dtheta)[2]
    yres.dl.dtheta <- score.yres$dl.dtheta
    yres.d2l.dtheta.dtheta <- score.yres$d2l.dtheta.dtheta
    npar.yres2 <- dim(score.yres2$dl.dtheta)[2]
    yres2.dl.dtheta <- score.yres2$dl.dtheta
    yres2.d2l.dtheta.dtheta <- score.yres2$d2l.dtheta.dtheta
    
  } else if (yres2.method[1]=="emp"){
    yres2.dl.dtheta <- as.matrix( y.resid^2-mean(y.resid^2)  )
    npar.yres2 <- 1
    yres2.d2l.dtheta.dtheta <- -N
    
  } else if (yres2.method[1]=="constant"){
    npar.yres2 <- 0
    yres2.dl.dtheta <- NULL
    yres2.d2l.dtheta.dtheta <- NULL
  }
  
  Ntheta <- npar.xz + npar.yz + npar.prod + npar.xres2 + npar.yres2 + npar.xres + npar.yres
  bigphi <- cbind(xz.dl.dtheta,
                  yz.dl.dtheta,
                  prod.dl.dtheta,
                  xres2.dl.dtheta,
                  yres2.dl.dtheta,
                  xres.dl.dtheta,
                  yres.dl.dtheta)
  
  A <- matrix(0, Ntheta, Ntheta)
  ### initiate the diagonal of A 
  A[(1:npar.xz), (1:npar.xz)] <- as.matrix(xz.d2l.dtheta.dtheta)
  A[npar.xz + (1:npar.yz), npar.xz + (1:npar.yz)] <- as.matrix(yz.d2l.dtheta.dtheta)
  A[npar.xz + npar.yz + (1: npar.prod), npar.xz + npar.yz + (1: npar.prod)] <- prod.d2l.dtheta.dtheta
  if (npar.xres2>0){
    A[npar.xz + npar.yz + npar.prod+ (1: npar.xres2), npar.xz + npar.yz + npar.prod+ (1: npar.xres2)] <- xres2.d2l.dtheta.dtheta
  }
  if(npar.yres2>0){
    A[npar.xz + npar.yz + npar.prod+ npar.xres2+(1:npar.yres2), npar.xz + npar.yz + npar.prod+ npar.xres2+(1: npar.yres2)] <- yres2.d2l.dtheta.dtheta
  }
  if (npar.xres>0){
    A[npar.xz + npar.yz + npar.prod+ npar.xres2 + npar.yres2 + (1:npar.xres), npar.xz + npar.yz + npar.prod+ npar.xres2 + npar.yres2 + (1:npar.xres)] <- xres.d2l.dtheta.dtheta
  }
  if(npar.yres>0){
    A[npar.xz + npar.yz + npar.prod+ npar.xres2 + npar.yres2 + npar.xres + (1:npar.yres), npar.xz + npar.yz + npar.prod+ npar.xres2 + npar.yres2 + npar.xres + (1:npar.yres)] <- yres.d2l.dtheta.dtheta
  }

  par.1 <- t(y.resid * designMatrOfZ) %*% t(dxresid.dtheta)
  par.2 <- t(x.resid * designMatrOfZ) %*% t(dyresid.dtheta)
  A[npar.xz + npar.yz + (1: npar.prod), 1:npar.xz] <- par.1
  A[npar.xz + npar.yz + (1: npar.prod), npar.xz + (1:npar.yz)] <- par.2
  
  if (xres2.method[1]=="model"){
    par.3 <- t(2*x.resid * designMatrOfZ) %*% t(dxresid.dtheta)
  } else if (xres2.method[1]=="emp"){
    par.3 <- t(2*x.resid ) %*% t(dxresid.dtheta)
  } else if(xres2.method[1]=="constant"){
    par.3 <- NULL
  }
  
  if(npar.xres2>0){
    A[npar.xz + npar.yz + npar.prod+ (1: npar.xres2), 1:npar.xz] <- par.3
  }
  
  if (yres2.method[1]=="model"){
    par.4 <- t(2*y.resid * designMatrOfZ) %*% t(dyresid.dtheta)
  } else if(yres2.method[1]=="emp"){
    par.4 <- t(2*y.resid ) %*% t(dyresid.dtheta)
  } else if(yres2.method[1]=="constant"){
    par.4 <- NULL
  }
  
  if(npar.yres2>0){
    A[npar.xz + npar.yz + npar.prod+ npar.xres2 + (1: npar.yres2), npar.xz + (1:npar.yz)] <- par.4
  }
  
  if (xres.method[1]=="model"){
    par.5 <- t(designMatrOfZ) %*% t(dxresid.dtheta)
  } ############ other methods are not implemented
   # else if (xres2.method[1]=="emp"){
#     par.5 <- t(2*x.resid ) %*% t(dxresid.dtheta)
#   } else if(xres2.method[1]=="constant"){
#     par.5 <- NULL
#   }
  
  if(npar.xres>0){
    A[npar.xz + npar.yz + npar.prod + npar.xres2 + npar.yres2 + (1: npar.xres), 1:npar.xz] <- par.5
  }
  
  if (yres2.method[1]=="model"){
    par.6 <- t(designMatrOfZ) %*% t(dyresid.dtheta)
  } ############ other methods are not implemented
  # else if(yres2.method[1]=="emp"){
  #   par.6 <- t(2*y.resid ) %*% t(dyresid.dtheta)
  # } else if(yres2.method[1]=="constant"){
  #   par.6 <- NULL
  # }
  
  if(npar.yres2>0){
    A[npar.xz + npar.yz + npar.prod + npar.xres2 + npar.yres2 + npar.xres + (1:npar.yres), npar.xz + (1:npar.yz)] <- par.6
  }
  
  SS <- solve(A, t(bigphi))
  var.theta <- tcrossprod(SS, SS)
  
  prod.coef <- score.prod$mod$coef
  xres.coef <- score.xres$mod$coef
  yres.coef <- score.yres$mod$coef
  
  names(prod.coef) <- paste("prod:", names(prod.coef))
  if(xres2.method[1]=="model"){
    xres2.coef <- score.xres2$mod$coef
  } else if (xres2.method[1]=="emp"){
    xres2.coef <- mean(x.resid^2)
  } else if (xres2.method[1]=="constant"){
    xres2.coef <- 1/3
  }
  names(xres2.coef) <- paste("xres2:", names(xres2.coef))
  
  if(yres2.method[1]=="model"){
    yres2.coef <- score.yres2$mod$coef
  } else if (yres2.method[1]=="emp"){
    yres2.coef <- mean(y.resid^2)
  } else if (yres2.method[1]=="constant"){
    yres2.coef <- 1/3
  }
  names(yres2.coef) <- paste("yres2:", names(yres2.coef))
  coef <- rbind(prod.coef, xres2.coef, yres2.coef, xres.coef, yres.coef)
  
  var.coef <- var.theta[(npar.xz+npar.yz+1): Ntheta , (npar.xz+npar.yz+1): Ntheta]
  
  ####### compute predicted values
  # regResult = ols(x.resid*y.resid ~ rcs(z, numKnots))
  # predProd = model.matrix(regResult) %*% regResult$coefficients
  #score.prod$mod
  if(FALSE){
    regYX = ols(scoreRes1$presid*scoreRes2$presid ~ rcs(z, numKnots))
    regX2 = ols((scoreRes1$presid)^2 ~ rcs(z, numKnots))
    regY2 = ols((scoreRes2$presid)^2 ~ rcs(z, numKnots))
    regX = ols((scoreRes1$presid) ~ rcs(z, numKnots))
    regY = ols((scoreRes2$presid) ~ rcs(z, numKnots))
    pYX = model.matrix(regYX) %*% regYX$coefficients
    pX2 = model.matrix(regX2) %*% regX2$coefficients
    pY2 = model.matrix(regY2) %*% regY2$coefficients
    pX = model.matrix(regX) %*% regX$coefficients
    pY = model.matrix(regY) %*% regY$coefficients
    pRho =  pYX/sqrt(pX2*pY2)
    # pRho =  pYX/sqrt(pX2*pY2)
  }
  predYX = matrix(prod.coef, nrow = 1) %*% t(newDesignMatrOfZ)
  predX2 = matrix(xres2.coef, nrow = 1) %*% t(newDesignMatrOfZ)
  predY2 = matrix(yres2.coef, nrow = 1) %*% t(newDesignMatrOfZ)
  predX = matrix(xres.coef, nrow = 1) %*% t(newDesignMatrOfZ)
  predY = matrix(yres.coef, nrow = 1) %*% t(newDesignMatrOfZ)
  C = (predYX - predX*predY)
  D = sqrt((predX2 - predX^2)*(predY2 - predY^2))
  predRho = C/D
  
  ####### compute derivatives for delta method:
  ### npar.prod + npar.xres2 + npar.yres2
  derivForDeltaMethod = cbind(newDesignMatrOfZ, newDesignMatrOfZ, newDesignMatrOfZ, newDesignMatrOfZ, newDesignMatrOfZ)
  
  ### divide by the denominator
  derivForDeltaMethod = derivForDeltaMethod/matrix(rep(D, npar.prod + npar.xres2 + npar.yres2 + npar.xres + npar.yres), nrow = nrow(derivForDeltaMethod), ncol = npar.prod + npar.xres2 + npar.yres2 + npar.xres + npar.yres)

  ### multiply by -(1/2) and by C
  derivForDeltaMethod[, npar.prod + (1 : (npar.xres2 + npar.yres2))] = -(1/2)*derivForDeltaMethod[, npar.prod + (1 : (npar.xres2 + npar.yres2))] * matrix(rep(C, npar.xres2 + npar.yres2), nrow = nrow(derivForDeltaMethod), ncol = npar.xres2 + npar.yres2)
  
  ### divide by predX2
  derivForDeltaMethod[, npar.prod + (1 : npar.xres2)] = derivForDeltaMethod[, npar.prod + (1 : npar.xres2)] / matrix(rep(predX2, npar.xres2), nrow = nrow(derivForDeltaMethod), ncol = npar.xres2) 
  
  ### divide by predY2
  derivForDeltaMethod[, npar.prod + npar.xres2 + (1 : npar.yres2)] = derivForDeltaMethod[, npar.prod + npar.xres2 + (1 : npar.yres2)] / matrix(rep(predY2, npar.yres2), nrow = nrow(derivForDeltaMethod), ncol = npar.yres2) 
  
  ### multiply by -predY + C/(predX2 - predX^2)
  derivForDeltaMethod[, npar.prod + npar.xres2 + npar.yres2 + (1:npar.xres)] = derivForDeltaMethod[, npar.prod + npar.xres2 + npar.yres2 + (1:npar.xres)] / matrix(rep(-predY - C/(predX2 - predX^2), npar.xres), nrow = nrow(derivForDeltaMethod), ncol = npar.xres) 
  
  ### multiply by -predX + C/(predY2 - predY^2)
  derivForDeltaMethod[, npar.prod + npar.xres2 + npar.yres2 + npar.xres + (1:npar.yres)] = derivForDeltaMethod[, npar.prod + npar.xres2 + npar.yres2 + npar.xres + (1:npar.yres)] / matrix(rep(-predX - C/(predY2 - predY^2), npar.yres), nrow = nrow(derivForDeltaMethod), ncol = npar.yres) 
  
  resultingVar = diag(derivForDeltaMethod %*% var.coef %*% t(derivForDeltaMethod))
  
  ci = 0.95
  pointEstAndCIs = data.frame(pointEst = as.vector(predRho), lower = as.vector(predRho - abs(qnorm(0.5*(1-ci)))*sqrt(resultingVar)), upper = as.vector(predRho + abs(qnorm(0.5*(1-ci)))*sqrt(resultingVar)))
    
  result <- list(prod=prod.coef,
                 xres2=xres2.coef,
                 yres2=yres2.coef,
                 xres=xres.coef,
                 yres=yres.coef,
                 var.coef=var.coef,
                 pointEstAndCIs = pointEstAndCIs)
                 
  return(result)
}

##################################### conditional cor for simulated data (used by
###       conditionalCorSimulatedData)
##################################### conditional cor for simulated data
##################################### conditional cor for simulated data
### this function is used by conditionalCorSimulatedData
computeCondCor = function(dataForCondCor, variableName = "calendarTime", whichFun = prodTS.new.SKE, 
                          numKnots = 3, method = "exp", newZ = seq(0, 3, .1), 
                          plotStuff = FALSE, quadraticTerms = FALSE){
  ### Returns conditional correlation for given variable "calendarTime"
  ###    It uses function prodTS.new.SKE() (see parameter whichFun)
  ###    that returns conditional correlation
  
  debugMode = FALSE
  formulaForVars = c(calendarTime = "calendarTime")
  formula1 = as.formula(paste0("Surv(timeToRegChangeEvent, regChangeEvent) ~ ", paste0(formulaForVars, collapse = " + ")))
  formula2 = as.formula(paste0("Surv(timeToVirFailEvent, virFailEvent) ~ ", paste0(formulaForVars, collapse = " + ")))
  
  # print(formula1)
  if(method == "exp"){
    modC1 = survreg(formula1, dist="exponential", data = dataForCondCor, score=TRUE)
    modC2 = survreg(formula2, dist="exponential", data = dataForCondCor, score=TRUE)
    partialCor = mestimation.survreg(object1 = modC1, object2 = modC2)
    scoreRes1 = mestimation.survregPSR(object = modC1, inverseA = FALSE)
    scoreRes2 = mestimation.survregPSR(object = modC2, inverseA = FALSE)
  }else{
    modC1 = coxph(formula1, data=dataForCondCor, timefix = FALSE, x = TRUE, y = TRUE)
    modC2 = coxph(formula2, data=dataForCondCor, timefix = FALSE)
    partialCor = mestimation.coxph(object1 = modC1, object2 = modC2, method = "full")
  
    if(method == "old"){
      scoreRes1 = mestimation.coxph_Old(object = modC1)
      scoreRes2 = mestimation.coxph_Old(object = modC2)
    }else{
      if(method == "full"){
        scoreRes1 = mestimation.coxph.full.likelihood(object = modC1, continuous = FALSE)
        scoreRes2 = mestimation.coxph.full.likelihood(object = modC2, continuous = FALSE)
      }else{
        scoreRes1 = mestimation.coxph.how.did.i.not.think.of.it(object = modC1)
        scoreRes2 = mestimation.coxph.how.did.i.not.think.of.it(object = modC1)
      }
    }
  }### end of if(method == "exp")
  if(det(scoreRes1$d2l.dtheta.dtheta) == 0 & !debugMode){
    scoreRes1 = mestimation.coxph_Old(object = modC1)
    warning("Reverted to M-estimation assuming known hazard for object1")
  }
  if(det(scoreRes2$d2l.dtheta.dtheta) == 0 & !debugMode){
    scoreRes2 = mestimation.coxph_Old(object = modC2)
    warning("Reverted to M-estimation assuming known hazard for object2")
  }
  z = dataForCondCor[[variableName]]
  #newZ = dataForCondCor[[variableName]]

  if(numKnots >= 0){
    result <- whichFun(x.resid = scoreRes1$presid, y.resid = scoreRes2$presid, z = z, newZ = newZ, numKnots = numKnots, xres2.method="model", yres2.method="model", xz.dl.dtheta = t(scoreRes1$dl.dtheta), yz.dl.dtheta = t(scoreRes2$dl.dtheta), xz.d2l.dtheta.dtheta = scoreRes1$d2l.dtheta.dtheta, yz.d2l.dtheta.dtheta = scoreRes2$d2l.dtheta.dtheta, dxresid.dtheta = scoreRes1$dpresid.dtheta, dyresid.dtheta = scoreRes2$dpresid.dtheta, debugMode = debugMode, quadraticTerms = quadraticTerms)
    # x.resid = scoreRes1$presid; y.resid = scoreRes2$presid; z = z; newZ = newZ; numKnots = numKnots; xres2.method="model"; yres2.method="model"; xz.dl.dtheta = t(scoreRes1$dl.dtheta); yz.dl.dtheta = t(scoreRes2$dl.dtheta); xz.d2l.dtheta.dtheta = scoreRes1$d2l.dtheta.dtheta; yz.d2l.dtheta.dtheta = scoreRes2$d2l.dtheta.dtheta; dxresid.dtheta = scoreRes1$dpresid.dtheta; dyresid.dtheta = scoreRes2$dpresid.dtheta; debugMode = debugMode
    plotData = data.frame(x = newZ, y = result$pointEstAndCIs$pointEst, yLower = result$pointEstAndCIs$lower, yUpper = result$pointEstAndCIs$upper, partialY = rep(partialCor$TS, length(newZ)), parialYLower = rep(partialCor$CI_TS[1], length(newZ)), parialYUpper = rep(partialCor$CI_TS[2], length(newZ))  )
  }else{
    # plotData = data.frame(x = newZ, y = rep(partialCor$TS, length(newZ)), yLower = rep(partialCor$CI_TS[1], length(newZ)), yUpper = rep(partialCor$CI_TS[2], length(newZ)) )
    plotData = data.frame(x = newZ, partialY = rep(partialCor$TS, length(newZ)), parialYLower = rep(partialCor$CI_TS[1], length(newZ)), parialYUpper = rep(partialCor$CI_TS[2], length(newZ)) )
  }
  plotData = plotData[order(plotData$x),]

  if(plotStuff){
    #plot(plotData$x, plotData$y, ylim = range(plotData[, c("y", "yLower", "yUpper")]))
    plotCol = "#22222222"
    points(plotData$x, plotData$y, pch = 19, col = plotCol, cex = .3)  
    lines(plotData$x, plotData$yLower, col = plotCol)
    lines(plotData$x, plotData$yUpper, col = plotCol)
    abline(h=0, lty = 3, col = plotCol)
  }
  plotData
}

##################################### conditional corr: making the data that Bryan wanted
##################################### conditional corr: making the data that Bryan wanted
##################################### conditional corr: making the data that Bryan wanted
conditionalCorSimulatedData = function(subjNum = 500, simNum = 20, newZ = seq(from = 0.01, to = 3, by = .5), 
                                       methods = c("full", "exp"), seeds = NULL, whichFun = prodTS.new.SKE, 
                                       fileNameStarter = "conditionalCorrelation", quadraticTerms = FALSE,
                                       censVec = c(.99, .7), family){
  corOfRes = rep(NA, 9)
  simdata = data.frame(calendarTime = rep(NA, subjNum))
  randomCensoring = TRUE
  knotsVec = c(3, 0, -1)
  #censVec = c(.99, .7) ### if .99 - no censoring; if .7 - censoring 30 percent
  dependCens = TRUE
  # correlationPatterns = list(constant = function(x){rep(.2, length(x))}, linearBryan = function(x){4*x*qexp(.95)/90}, quadraticBryanTwo = function(x){ (.4*x ) *  ( 1.2 - .4*x)}, quadraticBryanTwo = function(x){ (.5*x ) *  ( 1.5 - .5*x)}, runifQexp = function(x){(x*qexp(.95))/18}, sinFunction = function(x){(sin(4*x-pi/2)+1)/4} )
  correlationPatterns = list(constant = function(x){rep(.2, length(x))}, linearBryan = function(x){4*x*qexp(.95)/90}, quadraticBryanTwo = function(x){ (.4*x ) *  ( 1.2 - .4*x) + .001}, runifQexp = function(x){(x*qexp(.95))/18}, sinFunction = function(x){(sin(4*x-pi/2)+1)/4} )
  functionsToPlot = c("constant", "linearBryan", "quadraticBryanTwo")
  fitText = c("Quadratic Fit", "Linear Fit", "Constant Fit")
  
  familyLow = tolower(family)
  familyCap = gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", familyLow, perl = TRUE)
  
  awesomeArray = array(0, c(4, length(newZ), simNum, length(methods), length(functionsToPlot), length(knotsVec), length(censVec)), dimnames = list(cor = c("true", "est", "lower", "upper"), newZPoints = newZ, simNum = 1:simNum, method = methods, func = functionsToPlot, fit = knotsVec, cens = censVec))
  if(is.null(seeds)){
    set.seed(1)
    seeds = sample(1:(simNum*10), simNum)
  }
  
  ###### save approx percent censoring:
  apprPercent = matrix(NA, nrow = simNum, ncol = 2)
  
  goodSimulation = TRUE
  sim_i = 1
  seed_i = 1
  while(sim_i <= simNum & seed_i <= length(seeds)){
    # set.seed(seeds[seed_i])
    cat("sim:", sim_i, "   seed_i:", seed_i, "   seed:",  seeds[seed_i], "\n")
    ########## begin try
    for(mi in methods){
      for(ck in censVec){
        for (ni in functionsToPlot){
          
          if(goodSimulation){
            simdata$calendarTime = runif(subjNum)*qexp(.95)
            corrFunctionGrid = correlationPatterns[[ni]](as.numeric(newZ))
            # ### temporary fix
            # corrFunctionGrid[corrFunctionGrid < 0] = 0
            # ### temporary fix
            corrValuesDataGrid = mapCorValueToParamInCopula(corValues = corrFunctionGrid, family = familyCap)
            corrFunction = correlationPatterns[[ni]](simdata$calendarTime)
            ### temporary fix
            corrFunction[corrFunction < 0] = 0
            ### temporary fix
            corrValuesData = mapCorValueToParamInCopula(corValues = corrFunction, family = familyCap)
            ### Let's simulate some conditional correlation
            res12 = rCopulaConditional(nSim = subjNum, thetaPars = corrValuesData$predTheta, lambda1 = 1, lambda2 = 1, percentCensoring = 0, family = familyLow)
            simdata$timeToRegChangeEvent = res12[, "t1"]
            simdata$timeToVirFailEvent = res12[, "t2"]
            if(!randomCensoring){
              censoringTime = rep(quantile(simdata$calendarTime + simdata$timeToRegChangeEvent, ck), nrow(res12))
            }else{
              if (ck == .99){
                censoringTime = rep(Inf, nrow(res12))
              }else{
                if(FALSE){
                  beta = (1/(1-ck) - 1)
                  censoringTime = simdata$calendarTime + rexp(nrow(res12), rate = 1/beta)
                }else{
                  censMean1 = beta = (1/(1-ck) - 1)
                  censTime1 = rexp(subjNum, 1/rep(censMean1, nrow(simdata)))
                  censTime2 = rexp(subjNum, 1/rep(censMean1, nrow(simdata)))
                  if(dependCens){
                    censoringTime = censTime2 = censTime1
                  }
                }
              }
            }
            if(FALSE){
              simdata$regChangeEvent = as.numeric(simdata$calendarTime + simdata$timeToRegChangeEvent <= censoringTime)
              simdata$virFailEvent = as.numeric(simdata$calendarTime + simdata$timeToVirFailEvent <= censoringTime)
              simdata$timeToRegChangeEvent[simdata$regChangeEvent == 0] = censoringTime[simdata$regChangeEvent == 0]
              simdata$timeToVirFailEvent[simdata$virFailEvent == 0] = censoringTime[simdata$virFailEvent == 0]         
            }else{
              simdata$regChangeEvent = as.numeric(simdata$timeToRegChangeEvent <= censoringTime)
              simdata$virFailEvent = as.numeric(simdata$timeToVirFailEvent <= censoringTime)
              simdata$timeToRegChangeEvent[simdata$regChangeEvent == 0] = censoringTime[simdata$regChangeEvent == 0]
              simdata$timeToVirFailEvent[simdata$virFailEvent == 0] = censoringTime[simdata$virFailEvent == 0]         
            }
            apprPercent[sim_i, ] = c(1 - mean(simdata$regChangeEvent), 1 - mean(simdata$virFailEvent))
            # cat("Approx. percent cens. =", apprPercent[sim_i, 1], apprPercent[sim_i, 2], "\n")
            for (ki in knotsVec){
              if(goodSimulation){
                tryOutput <- try(tmp <- computeCondCor(dataForCondCor = simdata, variableName = "calendarTime", numKnots = ki, whichFun = whichFun, method = mi,  newZ = newZ, plotStuff = FALSE, quadraticTerms = quadraticTerms), silent = TRUE)
                if(is.null(dim(tryOutput))){
                  goodSimulation = FALSE
                  cat("Bad simulation i =", sim_i,  ",  seed_i =", seed_i, "seed = ", seeds[seed_i], "\n")
                }
              }### end of 1st if (goodSimulation) inside ki
              if(goodSimulation){
                if(ki >= 0){
                  awesomeArray[c("est", "lower", "upper"), , as.character(sim_i), method = mi, func = ni, fit = as.character(ki), cens = as.character(ck)] = t(tmp[, c("y", "yLower", "yUpper")])
                }else{
                  awesomeArray[c("est", "lower", "upper"), , as.character(sim_i), method = mi, func = ni, fit = as.character(ki), cens = as.character(ck)] = t(tmp[, c("partialY", "parialYLower", "parialYUpper")])
                }
                awesomeArray[c("true"), , as.character(sim_i), mi, ni, as.character(ki), as.character(ck)] = corrValuesDataGrid$corValues
              } ### end of 2nd if (goodSimulation) inside ki
            } ### end of for (ki)
          
          } ### end of if (goodSimulation)  outside ki, inside ni
          
        } ### end of for (ni)
      } ### end of for (ck)
    } ### end of for (mi)
    
    if(goodSimulation){
      sim_i = sim_i + 1
    }
    goodSimulation = TRUE
    seed_i = seed_i + 1
    
  } ### end of while sim_i
  fileName = paste0(fileNameStarter, "_subj_", subjNum, "_sim_", simNum, ".rda")
  save(awesomeArray, file = fileName, compress = TRUE)
  save(apprPercent, file = "apprPercent.rda", compress = TRUE)
  fileName
}

##################################### conditional corr: making the data that Bryan wanted
##################################### conditional corr: making the data that Bryan wanted
##################################### conditional corr: making the data that Bryan wanted
conditionalCorSimulatedDataTakeTwo = function(subjNum = 500, simNum = 20, newZ = seq(from = 0.01, to = 3, by = .5), 
                                              methods = c("full", "exp"), seeds = NULL, whichFun = prodTS.new.SKE, 
                                              fileNameStarter = "conditionalCorrelation", quadraticTerms = FALSE){
  corOfRes = rep(NA, 9)
  simdata = data.frame(calendarTime = rep(NA, subjNum))
  randomCensoring = TRUE
  knotsVec = c(3, 0, -1)
  censVec = c(.99, .7) ### if .99 - no censoring; if .7 - censoring 30 percent
  dependCens = TRUE
  correlationPatterns = list(constant = function(x){rep(.2, length(x))}, linearBryan = function(x){4*x*qexp(.95)/90}, quadraticBryanTwo = function(x){ (.4*x ) *  ( 1.2 - .4*x) + .001}, runifQexp = function(x){(x*qexp(.95))/18}, sinFunction = function(x){(sin(4*x-pi/2)+1)/4} )
  functionsToPlot = c("constant", "linearBryan", "quadraticBryanTwo")
  fitText = c("Quadratic Fit", "Linear Fit", "Constant Fit")

  awesomeArray = array(0, c(4, length(newZ), simNum, length(methods), length(functionsToPlot), length(knotsVec), length(censVec)), dimnames = list(cor = c("true", "est", "lower", "upper"), newZPoints = newZ, simNum = 1:simNum, method = methods, func = functionsToPlot, fit = knotsVec, cens = censVec))
  if(is.null(seeds)){
    set.seed(1)
    seeds = sample(1:(simNum*10), simNum)
  }
  
  ###### save approx percent censoring:
  apprPercent = matrix(NA, nrow = simNum, ncol = 2)
  
  goodSimulation = TRUE
  sim_i = 1
  seed_i = 1
  while(sim_i <= simNum & seed_i <= length(seeds)){
    # set.seed(seeds[seed_i])
    # cat("sim:", sim_i, "   seed_i:", seed_i, "   seed:",  seeds[seed_i], "\n")
    ########## begin try
    for(mi in methods){
      for(ck in censVec){
        for (ni in functionsToPlot){
          
          if(goodSimulation){
            simdata$calendarTime = runif(subjNum)*qexp(.95)
            corrFunctionGrid = correlationPatterns[[ni]](as.numeric(newZ))
            # ### temporary fix
            corrValuesDataGrid = mapCorValueToParamInCopula(corValues = corrFunctionGrid, family = "Clayton")
            corrFunction = correlationPatterns[[ni]](simdata$calendarTime)
            ### temporary fix
            corrFunction[corrFunction < 0] = 0
            ### temporary fix
            corrValuesData = mapCorValueToParamInCopula(corValues = corrFunction, family = "Clayton")
            ### Let's simulate some conditional correlation
            res12 = rCopulaConditional(nSim = subjNum, thetaPars = corrValuesData$predTheta, lambda1 = 1, lambda2 = 1, percentCensoring = 0, family = "clayton")
            simdata$timeToRegChangeEvent = res12[, "t1"]
            simdata$timeToVirFailEvent = res12[, "t2"]
            if(!randomCensoring){
              censoringTime = rep(quantile(simdata$calendarTime + simdata$timeToRegChangeEvent, ck), nrow(res12))
            }else{
              if (ck == .99){
                censoringTime = rep(Inf, nrow(res12))
              }else{                
                censMean1 = beta = (1/(1-ck) - 1)
                censTime1 = rexp(subjNum, 1/rep(censMean1, nrow(simdata)))
                censTime2 = rexp(subjNum, 1/rep(censMean1, nrow(simdata)))
                if(dependCens){
                  censoringTime = censTime2 = censTime1
                }        
              }
            }
            simdata$regChangeEvent = as.numeric(simdata$timeToRegChangeEvent <= censoringTime)
            simdata$virFailEvent = as.numeric(simdata$timeToVirFailEvent <= censoringTime)
            simdata$timeToRegChangeEvent[simdata$regChangeEvent == 0] = censoringTime[simdata$regChangeEvent == 0]
            simdata$timeToVirFailEvent[simdata$virFailEvent == 0] = censoringTime[simdata$virFailEvent == 0]         
            apprPercent[sim_i, ] = c(1 - mean(simdata$regChangeEvent), 1 - mean(simdata$virFailEvent))
            for (ki in knotsVec){
              cat("ki", ki, "\n")
              if(goodSimulation){
                ### dataForCondCor = simdata; variableName = "calendarTime"; numKnots = ki; whichFun = whichFun; method = mi;  newZ = newZ; plotStuff = FALSE; quadraticTerms = quadraticTerms
                tryOutput <- try(tmp <- computeCondCor(dataForCondCor = simdata, variableName = "calendarTime", numKnots = ki, whichFun = whichFun, method = mi,  newZ = newZ, plotStuff = FALSE, quadraticTerms = quadraticTerms), silent = TRUE)
                if(is.null(dim(tryOutput))){
                  goodSimulation = FALSE
                  cat("Bad simulation i =", sim_i,  ",  seed_i =", seed_i, "seed = ", seeds[seed_i], "\n")
                }
              }### end of 1st if (goodSimulation) inside ki
              if(goodSimulation){
                if(ki >= 0){
                  awesomeArray[c("est", "lower", "upper"), , as.character(sim_i), method = mi, func = ni, fit = as.character(ki), cens = as.character(ck)] = t(tmp[, c("y", "yLower", "yUpper")])
                }else{
                  awesomeArray[c("est", "lower", "upper"), , as.character(sim_i), method = mi, func = ni, fit = as.character(ki), cens = as.character(ck)] = t(tmp[, c("partialY", "parialYLower", "parialYUpper")])
                }
                awesomeArray[c("true"), , as.character(sim_i), mi, ni, as.character(ki), as.character(ck)] = corrValuesDataGrid$corValues
              } ### end of 2nd if (goodSimulation) inside ki
            } ### end of for (ki)
          
          } ### end of if (goodSimulation)  outside ki, inside ni
          
        } ### end of for (ni)
      } ### end of for (ck)
    } ### end of for (mi)
    
    if(goodSimulation){
      sim_i = sim_i + 1
    }
    goodSimulation = TRUE
    seed_i = seed_i + 1
    
  } ### end of while sim_i
  fileName = paste0(fileNameStarter, "_subj_", subjNum, "_sim_", simNum, ".rda")
  save(awesomeArray, file = fileName, compress = TRUE)
  save(apprPercent, file = "apprPercent.rda", compress = TRUE)
  fileName
}

###################################################
### This function matches corValues to copula parameters needed to simulate conditional cor.
###################################################
mapCorValueToParamInCopula = function(corValues, family = "Frank"){
  ### This function figures out the value of the copula parameter depending on what correlation 
  ### we need.
  ### Currently, it only takes absolute correlation values from 0 to .5
  library(cubature)        # load the package "cubature"
  claytonsFamily <- function(x, theta = .5) {
    (max(x[1]^(-theta) + x[2]^(-theta) - 1,   0))^(-1/theta)
  }
  franksFamily <- function(x, theta = .5) {
    (-1/theta)*log(  1  + (exp(-theta*x[1]) - 1) * (exp(-theta*x[2]) - 1)/(exp(-theta) - 1)  )
  }
  
  cat(min(corValues), "\n")
  
  if(min(corValues)< -.6 | max(corValues) > .6) stop("Argument _corValue_ should be between -.6 and .6")
  #if(min(corValues) <= 0 & family == "Clayton") stop("Clayton family can generate only positive correlation")
  if(!all(family %in%  c("Clayton", "Frank"))) stop("Only Clayton or Frank family is implemented")
  
  if(family == "Clayton"){
    copFun = claytonsFamily
    y = seq(0.001, 3, by=.1)
    # y = seq(0.0, 3, by=.1)
  }else{
    copFun = franksFamily
    y = seq(-.001, 3, by=.1)
  }
  x = sapply(y, function(x){12*adaptIntegrate(copFun, lowerLimit = c(0, 0), upperLimit = c(1, 1), theta = x)$integral - 3})
  fitThetaModel = lm(y ~ poly(x, 4, raw=TRUE))
  predTheta = predict(fitThetaModel, newdata = data.frame(x=abs(corValues)), type= "response")
  
  if(family == "Clayton"){
    ### Clayton's family is not allowed to have 0's or less
    predTheta[corValues == 0] = 0
    predTheta[corValues  < 0] = NA
  }else{
    ### We use the fact that for Franks family cor(-theta) = -cor(theta),
    predTheta[corValues < 0] = -predTheta[corValues < 0]
  }
  res = data.frame(predTheta = predTheta, corValues = corValues)
  res
}

###################################################
### generates conditionaly dependent variables
###################################################
rCopulaConditional = function(nSim, thetaPars, lambda1=1, lambda2=1, percentCensoring = 0.1, family = "clayton"){
  ### Generates a given family of copulas with different parametres in the same sample
  ### and with marginal exponential distributions
  ### Compared to function rclayton(), this is an added flexibility that allows to simulate 
  ### conditional correlation. Otherwise function rclayton() does the same thing
  ### Note, that we could just use this function instead of rclayton(),
  ### but this function will be slower because it calls sapply, creates lists and stuff
  ### thetaPars = is a vector of parameters that has a length of nSim
  
  if(length(thetaPars) != nSim) stop("Length of thetaPars has to be equal to nSim")
  if(family == "clayton"){
    if(any(thetaPars == 0)){  ###### trying to decide which one is better
      t1t2List = t(sapply(thetaPars, function(p){rarchi(n = 1, "clayton", d = 2, theta = p)}))
    }else{
      cop12List = sapply(thetaPars, function(p) {archmCopula(family = "clayton", param = p)})
      t1t2List = t(sapply(cop12List, function(cop12){rCopula(1, cop12)}))
    }
  }else{
    if(family != "frank") stop("Only Clayton's and Frank's families were implemented\n")
    cop12List = sapply(thetaPars, function(p) {archmCopula(family = "frank", param = p)})
    t1t2List = t(sapply(cop12List, function(cop12){rCopula(1, cop12)}))
  }
  t1 = -log(1-t1t2List[,1])/lambda1
  t2 = -log(1-t1t2List[,2])/lambda2
  
  delta1 = delta2 = rep(1, nSim)
  if(percentCensoring > 0){
    if(FALSE){
      ######### the old way according to Shih and Louis
      delta1 = sample(x = c(0, 1), size = nSim, replace = TRUE, prob = c(percentCensoring, 1-percentCensoring))
      delta2 = sample(x = c(0, 1), size = nSim, replace = TRUE, prob = c(percentCensoring, 1-percentCensoring))
      t1[delta1==0] = runif(sum(delta1==0))*2.3
      t2[delta2==0] = runif(sum(delta2==0))*2.3
    }else{
      ######### the new way according to Fan, Hsu, Prentice
      censoringTime1 = rexp(n = nSim, rate = lambda1/2)
      censoringTime2 = rexp(n = nSim, rate = lambda2/2)
      delta1 = as.numeric(t1 <= censoringTime1) 
      delta2 = as.numeric(t2 <= censoringTime2) 
      t1[delta1==0] = censoringTime1[delta1==0]
      t2[delta2==0] = censoringTime2[delta2==0]
    }
  }
  
  res = cbind(t1, delta1, t2, delta2)
  colnames(res) = c("t1", "delta1", "t2", "delta2")
  res
}

###############################################
############################################### M-estimation for survival regression
###############################################
mestimation.coxph <- function(object1, object2, method = "full"){
  if(method == "old"){
    scoreRes1 = mestimation.coxph_Old(object = object1)
    scoreRes2 = mestimation.coxph_Old(object = object2)
  }
  if(method == "full"){
    scoreRes1 = mestimation.coxph.full.likelihood(object = object1, continuos = FALSE)
    scoreRes2 = mestimation.coxph.full.likelihood(object = object2, continuos = FALSE)
    if(det(scoreRes1$d2l.dtheta.dtheta) == 0){
      scoreRes1 = mestimation.coxph_Old(object = object1)
      cat("FUNCTION mestimation.coxph: Reverted to M-estimation assuming known hazard for object1\n")
      warning("FUNCTION mestimation.coxph: Reverted to M-estimation assuming known hazard for object1")
    }
    if(det(scoreRes2$d2l.dtheta.dtheta) == 0){
      scoreRes2 = mestimation.coxph_Old(object = object2)
      cat("FUNCTION mestimation.coxph: Reverted to M-estimation assuming known hazard for object2\n")
      warning("FUNCTION mestimation.coxph: Reverted to M-estimation assuming known hazard for object2")
    }
  }
  if(method == "stute"){
    scoreRes1 = mestimation.coxph.how.did.i.not.think.of.it(object = object1)
    scoreRes2 = mestimation.coxph.how.did.i.not.think.of.it(object = object2)
  }
  cor1 = corTS.survreg.modified(xresid = scoreRes1$presid, yresid = scoreRes2$presid,
                                xz.dl.dtheta = t(scoreRes1$dl.dtheta),
                                yz.dl.dtheta = t(scoreRes2$dl.dtheta),
                                xz.d2l.dtheta.dtheta = scoreRes1$d2l.dtheta.dtheta,
                                yz.d2l.dtheta.dtheta = scoreRes2$d2l.dtheta.dtheta,
                                dxresid.dthetax = scoreRes1$dpresid.dtheta,
                                dyresid.dthetay = scoreRes2$dpresid.dtheta, inverseA = FALSE)
  # cor1 = corTS.survreg(xresid = scoreRes1$presid, yresid = scoreRes2$presid,
  #              xz.dl.dtheta = t(scoreRes1$dl.dtheta),
  #              yz.dl.dtheta = t(scoreRes2$dl.dtheta),
  #              xz.d2l.dtheta.dtheta = scoreRes1$d2l.dtheta.dtheta,
  #              yz.d2l.dtheta.dtheta = scoreRes2$d2l.dtheta.dtheta,
  #              dxresid.dthetax = scoreRes1$dpresid.dtheta,
  #              dyresid.dthetay = scoreRes2$dpresid.dtheta, inverseA = FALSE)
  cor1
}

###############################################
############################################### M-estimation for Cox using full likelihood
###############################################
mestimation.coxph.full.likelihood = function(object, time, continuous = TRUE, ...){
  ### in general, argument time can be extracted from object: time = object$y[, 1]
  ### but there were numerical issue with it 
  ################################# M-estimation for Cox  
  Y = object$y[,1]
  nSubj = length(Y)
  delta = object$y[,2]
  resid <- residuals(object, type="martingale")
  orderedY = Y[order(Y)]
  orderedDelta = delta[order(Y)]
  orderedUniqueFailureTimes = unique(orderedY[orderedDelta == 1])
  orderedFailureTimesWithZero = c(0, orderedUniqueFailureTimes)
  nHazardParam = length(orderedUniqueFailureTimes)
  orderedMartRes = resid[order(Y)]
  
  orderMatrix = matrix(c(1:nSubj, (1:nSubj)[order(Y)]), ncol=2)
  colnames(orderMatrix) = c("new", "old")
  originalOrder = orderMatrix[order(orderMatrix[,2]), 1]
  
  X1 = model.matrix(object)
  if(dim(X1)[2] == 1){
    orderedX1 = as.matrix(X1[order(Y),], ncol = 1)
  }else{
    orderedX1 = X1[order(Y),]
  }
  Xbeta =  X1 %*% matrix(object$coefficients, nrow = length(object$coefficients))
  expXbeta = exp(Xbeta)
  orderedExpXbeta = expXbeta[order(Y)]
  
  nParam = dim(object$var)[1]
  if(is.null(nParam)){
    nParam = 1
  }
  
  H0 = basehaz(object, centered=FALSE)
  # uniqueBlHaz = H0[, ]
  # timeDiff = diff(c(0, uniqueBlHaz[, "time"]))
  # lambdas = diff(c(0, uniqueBlHaz[, "hazard"]))/timeDiff
  ### the following code is for "discrete" time,  but it turns out we need it for everyting
  # #if(!continuous){
  #   H0$delta = 1
  #   H0_minus = hazardStretcher(H0)
  #   H0_minus$delta = orderedDelta
  #   H0 = H0_minus[, c("hazard", "time")]
  # #}
  
  #H0_minus = hazardStretcher_take_two(H0, orderedY)
  H0_minus = hazardStretcher_take_two(H0, Y)
  H0_minus$delta = orderedDelta
  uniqueBlHaz = H0_minus[H0_minus$delta == 1,]
  uniqueBlHaz = unique(uniqueBlHaz[, c("hazard", "time")])
  timeDiff = diff(c(0, uniqueBlHaz[, "time"]))
  lambdas = diff(c(0, uniqueBlHaz[, "hazard"]))/timeDiff
  
  ######################################## dl.dtheta
  ######################################## dl.dtheta
  ######################################## dl.dtheta
  dl.dtheta = matrix(NA, nrow = nParam + nHazardParam, ncol = nSubj)
  ######################################## dl.dtheta   for beta:
  for (j in 1:nParam){
    dl.dtheta[j, ] = (orderedDelta - H0_minus[,"hazard"]*orderedExpXbeta)*orderedX1[,j]
  }  
  ######################################## dl.dtheta   for lambda:
  dl.dlambda = matrix(0, nrow=nHazardParam, ncol = nSubj)
  for (k in 1:nHazardParam){
    hadEventAtT = (orderedDelta == 1 & orderedY == orderedFailureTimesWithZero[k+1])
    hadEventAfterT = (orderedY > orderedFailureTimesWithZero[k+1] | (orderedDelta == 0 & orderedY == orderedFailureTimesWithZero[k+1]))
    dl.dlambda[k, hadEventAtT] = orderedDelta[hadEventAtT]/lambdas[k] - timeDiff[k] * orderedExpXbeta[hadEventAtT]
    dl.dlambda[k, hadEventAfterT] = -timeDiff[k] * orderedExpXbeta[hadEventAfterT]
  }
  ######################################## merge dl.dtheta for betand and for lambda:
  dl.dtheta[(nParam+1):(nParam+nHazardParam), ] = dl.dlambda
  
  ######################################## d2l.dtheta.dtheta (matrix A)
  ######################################## d2l.dtheta.dtheta (matrix A)
  ######################################## d2l.dtheta.dtheta (matrix A)
  d2l.dtheta.dtheta = diag(nParam + nHazardParam)
  ######################################## d2l.dtheta.dtheta  for beta/beta:
  for (i in 1:nParam){
    for (j in 1:i){
      element = H0_minus[,"hazard"] * orderedExpXbeta * orderedX1[,i] * orderedX1[,j]
      d2l.dtheta.dtheta[i, j] = d2l.dtheta.dtheta[j, i] = -sum(element)
    }
  }  
  ######################################## d2l.dtheta.dtheta  for lambda/lambda:
  for (k in (1:nHazardParam)){
    hadEventAtT = (orderedY == orderedFailureTimesWithZero[k+1])
    ############################################
    element = orderedDelta/(lambdas[k]^2)
    d2l.dtheta.dtheta[nParam+k, nParam+k] = -sum(element[hadEventAtT])
  }
  
  ######################################## d2l.dtheta.dtheta  for beta/lambda:
  for (j in 1:nParam){
    for (k in 1:nHazardParam){
      greaterOrEqual = (orderedY >= orderedFailureTimesWithZero[k+1])
      element = timeDiff[k] * orderedExpXbeta * orderedX1[, j]
      d2l.dtheta.dtheta[j, nParam + k] = d2l.dtheta.dtheta[nParam + k, j] = -sum(element[greaterOrEqual])
    }
  }  
  
  ######################################## prob. scale residuals (PSR)
  presid = presid.coxph(object, continuous = continuous)
  orderedPresid = presid[order(Y)]         #### ORDERING BY Y
  
  ######################################## PSR derivatives by beta and lambda_k
  ######################################## PSR derivatives by beta and lambda_k
  dpresid.dtheta = matrix(0, nrow = nParam + nHazardParam, ncol = nSubj)
  if(continuous){
    ######################################## PSR derivatives by beta  
    expression1 = (1 + orderedDelta)*exp(orderedMartRes - orderedDelta)*(orderedDelta-orderedMartRes)
    dpresid.dtheta[1:nParam, ] = matrix(expression1, nrow = nParam, ncol = nSubj, byrow = TRUE) * t(orderedX1)
    
    ######################################## PSR derivatives by lambda_k
    expression2 = (1 + orderedDelta)*exp(orderedMartRes - orderedDelta) * orderedExpXbeta
    for(k in 1:nHazardParam){
      condGreaterOrEqual = (orderedY >= orderedFailureTimesWithZero[k+1])
      dpresid.dtheta[nParam + k, condGreaterOrEqual]  = timeDiff[k] * expression2[condGreaterOrEqual]
    }
  }else{
    element1 = exp(-H0_minus[, "hazard"] * orderedExpXbeta)
    # element2 = exp(-c(0, H0[1:(nrow(H0)-1),1]) * orderedExpXbeta)
    element2 = exp(-H0_minus[, "Hminus"] * orderedExpXbeta)
    # element1[H0[,1] == 0] = 1
    # element2[c(0, H0[1:(nrow(H0)-1),1]) == 0] = 1
    element3 = element1 * orderedExpXbeta
    element4 = element2 * orderedExpXbeta
    
    ######################################## PSR derivatives by beta  
    if(nParam>=1){
      for (j in 1:nParam){
        dpresid.dtheta[j, ] = (element1*H0_minus[, "hazard"] + orderedDelta * element2 * H0_minus[, "Hminus"]) * orderedExpXbeta * orderedX1[, j]
      }
    }
    
    ######################################## PSR derivatives by lambda_k
    for(k in 1:nHazardParam){
      # condGreater = (orderedY > orderedFailureTimesWithZero[k+1])
      # condEqual = (orderedY == orderedFailureTimesWithZero[k+1])
      # dpresid.dtheta[nParam + k, condGreater]  = (timeDiff[k] * (element3 + orderedDelta * element4) )[condGreater]
      # dpresid.dtheta[nParam + k, condEqual]  = (timeDiff[k] * element3)[condEqual]
      condGreaterOrEqual = (orderedY >= orderedFailureTimesWithZero[k+1])
      dpresid.dtheta[nParam + k, condGreaterOrEqual]  = (timeDiff[k] * (element3 + orderedDelta * element4) )[condGreaterOrEqual]
    }
  }
  
  res = list(dl.dtheta = dl.dtheta[, originalOrder],
             d2l.dtheta.dtheta = d2l.dtheta.dtheta,
             resid = NULL,
             dresid.dtheta = NULL,
             presid = orderedPresid[originalOrder],
             presid.k= NULL,
             dpresid.dtheta = dpresid.dtheta[, originalOrder],
             dpresid.dtheta.k = NULL)
  res
}

###############################################
############################################### Hazard Stretcher:
### it turns out that coxph does not stretch hazard for multiple events per time point
###############################################
hazardStretcher_take_two = function(data, time){
  ### data - data with hazard, time  what is supplied by cox model
  ### time - is the original time variable (with duplicates and everything...)
  
  ### data = data.frame(hazard = c(0, 0, .11, .11, .11, .22, .33, .33), time = c(1.1, 1.2, 2.2, 2.22, 2.3, 3, 4, 5.5)); time = c(1.1, 1.2, 2.2, 2.22, 2.22, 2.3, 3, 4, 4, 4, 5.5)
  ### data = data.frame(hazard = c(0, 0, .11, .11, .11, .22, .33, .44), delta = c(0, 0, 1, 0, 0, 1, 1, 1), time = c(1.1, 1.2, 2, 2.2, 2.3, 3, 4, 5.5))
  ### data = data.frame(hazard = c(.11, .11, .22, .33, .44), delta = c(1, 0, 1, 1, 1), time = c(2, 2.2, 3, 4, 5.5))
  ### data = data.frame(hazard = c(.11, .11, .22, .33, .44), delta = c(1, 0, 1, 1, 1), time = c(2, 2.2, 3, 4, 5.5))
  
  ######################### important note:
  ### when this funciton is applied to coxph baseline hazard (model1), 
  ### it is important that time comes from object$y[,1] because of issue with
  ### number representation (see argument timefix in coxph.control and Terri Therneau's email)
  
  if(!all(c("hazard", "time") %in% names(data))){
    stop("The data has to have the following names: hazard", "time")
  }
  
  orderedTime = time[order(time)]
  data = data[order(data$time),]
  uniqueHazard = unique(data[, c("hazard", "time")])
  #uniqueHazard = rbind(c(0, 1), uniqueHazard)
  uniqueHazard = rbind(c(0, 0), uniqueHazard)
  data$Hminus = NA
  
  ##################### stretch the data:
  ##################### define hazard for every point of the original time:
  orig_i = 1
  stretchedData = data.frame(time = orderedTime, hazard = NA)
  
  for(i in 1:nrow(uniqueHazard)){
    while((uniqueHazard$time[i] == orderedTime[orig_i]) & (orig_i <= length(orderedTime))){
      stretchedData$hazard[orig_i] = uniqueHazard$hazard[i]
      orig_i = orig_i + 1
    }
  }
  
  i_minus = 1
  i_unique = 2
  while(stretchedData$time[i_minus] <= uniqueHazard$time[nrow(uniqueHazard)] & i_minus <= nrow(stretchedData) & i_unique < nrow(uniqueHazard) ){
    ### fill in censored observations 
    if(stretchedData$time[i_minus] < uniqueHazard$time[i_unique+1]){
      stretchedData$Hminus[i_minus] = uniqueHazard$hazard[i_unique-1]
      i_minus = i_minus + 1
    }else{
      stretchedData$Hminus[i_minus] = uniqueHazard$hazard[i_unique]
      i_minus = i_minus + 1
      i_unique = i_unique + 1
    }
    #cat(i_minus, i_unique, "\n")
  }
  
  ### if there are any observations that are censored after the last failure point.
  while(i_minus <= nrow(stretchedData)){
    stretchedData$Hminus[i_minus] = uniqueHazard$hazard[nrow(uniqueHazard) - 1]
    i_minus = i_minus + 1
  }
  
  stretchedData
}

############################################## the code of Qi: just in case
###coxph()
#' @export
presid.coxph <- function(object, continuous = FALSE) {
  if(continuous){ ### this code does not take into account the discontinuity of the surv.fun.
    # time <- object$y[,1]
    delta <- object$y[,2]
    resid <- residuals(object, type="martingale")
    1 - exp(resid - delta) - delta*exp(resid - delta)
  }else{
    Y = object$y[,1]
    nSubj = length(Y)
    delta = object$y[,2]
    orderedY = Y[order(Y)]
    orderedDelta = delta[order(Y)]
    orderedUniqueFailureTimes = unique(orderedY[orderedDelta == 1])
    
    X1 = model.matrix(object)
    if(dim(X1)[2] == 1){
      orderedX1 = as.matrix(X1[order(Y),], ncol = 1)
    }else{
      orderedX1 = X1[order(Y),]
    }
    Xbeta =  X1 %*% matrix(object$coefficients, nrow = length(object$coefficients))
    expXbeta = exp(Xbeta)
    orderedExpXbeta = expXbeta[order(Y)]
    
    H0 = basehaz(object, centered=FALSE)
    H0 = hazardStretcher_take_two(H0, orderedY)
    uniqueBlHaz = H0[orderedDelta ==  1, ]
    timeDiff = diff(c(0, uniqueBlHaz[, "time"]))
    lambdas = diff(c(0, uniqueBlHaz[, "hazard"]))/timeDiff
    
    element1 = exp(-H0[, "hazard"] * orderedExpXbeta)
    #element2 = exp(-c(0, H0[1:(nrow(H0)-1), "hazard"]) * orderedExpXbeta)
    element2 = exp(-H0[, "Hminus"] * orderedExpXbeta)
    
    element1[H0[, "hazard"] == 0] = 1
    #element2[c(0, H0[1:(nrow(H0)-1), "hazard"]) == 0] = 1
    element2[H0[, "Hminus"] == 0] = 1
    
    ######################################## restore original order
    orderMatrix = matrix(c(1:nSubj, (1:nSubj)[order(Y)]), ncol=2)
    colnames(orderMatrix) = c("new", "old")
    originalOrder = orderMatrix[order(orderMatrix[,2]), 1]
    
    res = 1 - element1 - orderedDelta * element2
    res[originalOrder]
  }
}

###################################################################
######### Partial correlation usingM-estimation the modification of 
######### Qi's code that have correct sign for matrix A and that
######### that is built in a way that one does not have to invert
######### the part of matrix A that is obtained from the models
######### because it may sometime not invert and the models 
######### provide its inverse anyway ...
###################################################################
corTS.survreg.modified <- function(xresid, yresid,
                                   xz.dl.dtheta, yz.dl.dtheta,
                                   xz.d2l.dtheta.dtheta, yz.d2l.dtheta.dtheta,
                                   dxresid.dthetax, dyresid.dthetay, inverseA = TRUE,
                                   fisher=FALSE, confLevel = 0.95){
  
  TS = cor(xresid, yresid)
  
  xresid2 = xresid^2
  yresid2 = yresid^2
  xbyyresid = xresid * yresid
  mean.xresid = mean(xresid)
  mean.yresid = mean(yresid)
  mean.xbyyresid = mean(xbyyresid)
  
  bigphi = cbind(xz.dl.dtheta,
                 yz.dl.dtheta,
                 mean.xresid - xresid,
                 mean.yresid - yresid,
                 mean.xbyyresid - xbyyresid,
                 mean(xresid2)-xresid2,
                 mean(yresid2)-yresid2,
                 0)
  
  npar.xz = dim(xz.dl.dtheta)[2]
  npar.yz = dim(yz.dl.dtheta)[2]
  Ntheta = npar.xz + npar.yz + 6
  N = dim(xz.dl.dtheta)[1]
  
  A = matrix(0,Ntheta,Ntheta)
  A[1:npar.xz, 1:npar.xz] = xz.d2l.dtheta.dtheta
  A[npar.xz+(1:npar.yz), npar.xz+(1:npar.yz)] = yz.d2l.dtheta.dtheta
  
  A[Ntheta-6+(1:6), Ntheta-6+(1:6)] = diag(N, 6)
  
  bigpartial = rbind(c(dxresid.dthetax %*% rep(1, N), rep(0, npar.yz)),
                     c(rep(0, npar.xz), dyresid.dthetay %*% rep(1, N)),
                     c(dxresid.dthetax %*% yresid, dyresid.dthetay %*% xresid),
                     c(dxresid.dthetax %*% (2*xresid), rep(0, npar.yz)),
                     c(rep(0, npar.xz), dyresid.dthetay %*% (2*yresid)))
  
  A[Ntheta-6+(1:5), 1:(npar.xz+npar.yz)] = bigpartial
  
  ## TS also equals numTS / sqrt(varprod) = numTS * revsvp
  numTS = mean.xbyyresid - mean.xresid * mean.yresid
  var.xresid = mean(xresid2) - mean.xresid^2
  var.yresid = mean(yresid2) - mean.yresid^2
  varprod = var.xresid * var.yresid
  revsvp = 1/sqrt(varprod)
  dTS.dvarprod = numTS * (-0.5) * revsvp^3
  
  smallpartial = N *
    c(-mean.yresid * revsvp + dTS.dvarprod * (-2*mean.xresid*var.yresid),
      -mean.xresid * revsvp + dTS.dvarprod * (-2*mean.yresid*var.xresid),
      revsvp,
      dTS.dvarprod * var.yresid,
      dTS.dvarprod * var.xresid)
  
  A[Ntheta, Ntheta-6+(1:5)] = -smallpartial
  
  if (inverseA){ ### this option means that the part of matrix A (A11) is supplied in inverse form
    A11_inv = A[1:(npar.xz+npar.yz), 1:(npar.xz+npar.yz)]
    A21 = A[Ntheta-6+(1:6), 1:(npar.xz+npar.yz)]
    B_inv = solve(A[Ntheta-6+(1:6), Ntheta-6+(1:6)])
    A_inv = A
    A_inv[Ntheta-6+(1:6), 1:(npar.xz+npar.yz)] = - B_inv %*% A21 %*% A11_inv
    A_inv[Ntheta-6+(1:6), Ntheta-6+(1:6)] = B_inv
    var.theta = A_inv %*% t(bigphi) %*% bigphi %*% t(A_inv)
  }else{
    SS = solve(A, t(bigphi))
    var.theta = tcrossprod (SS, SS)
  }
  varTS = var.theta[Ntheta, Ntheta]  ### this is actually a squared standard error
  pvalTS = 2 * pnorm( -abs(TS)/sqrt(varTS))
  addToValue = abs(qnorm((1-confLevel)/2))*sqrt(varTS)
  CI_TS = TS + c(-1, 1)*addToValue
  
  
  ####Fisher's transformation
  TS_f <- log( (1+TS)/(1-TS) )
  varTS_f <- varTS*(2/(1-TS^2))^2
  pvalTS_f <- 2 * pnorm( -abs(TS_f)/sqrt(varTS_f))
  addToValue_f = abs(qnorm((1-confLevel)/2))*sqrt(varTS_f)
  CI_TS_f = TS_f + c(-1, 1)*addToValue_f
  
  list(TS=TS, varTS=varTS, pvalTS=pvalTS, CI_TS = CI_TS, pvalTS=pvalTS_f, CI_TS_f = CI_TS_f, var.thetaPSR = diag(var.theta)[dim(var.theta)[1] + c(-5 : 0)])
}

###################################################################
### Spline components
###################################################################
componentsOfSplines = function(x, knots){
  ### gives sline components
  positive = function(x){
    ### used by spline functions
    x[x <= 0]=0
    x
  }
  k = length(knots)
  res = matrix(NA, nrow=length(x), ncol=k-1)
  res[, 1] = x
  kd = (knots[k] - knots[1])^(2/3)
  for (j in 1:(k-2)){
    res[, j+1] = positive(((x-knots[j])/kd)^3) - positive(((x-knots[k-1])/kd)^3) * (knots[k]-knots[j])/(knots[k]-knots[k-1]) + positive(((x-knots[k])/kd)^3) * (knots[k-1]-knots[j])/(knots[k]-knots[k-1])
  }
  res
}

###################################################################
### Qi's code that is necessary for conditinal correlation
###################################################################
lm.scores <- function(y, X){
  N = length(y)  
  mod = lm(y~X)
  smod = summary(mod)
  resid = smod$residuals
  
  d2l.dtheta.dtheta = -crossprod(cbind(1, X))
  
  dl.dtheta <- resid*cbind(1, X)
  presid = 2*pnorm((y - mod$fitted.values)/smod$sigma) -1
  dresid.dtheta = t(cbind(-1, -X))
  dpresid.dtheta = t(cbind(-2*dnorm((y - mod$fitted.values)/smod$sigma)/smod$sigma,
                           -2*dnorm((y - mod$fitted.values)/smod$sigma)/smod$sigma *
                             X))
  
  f.y<-density(resid)
  fy.ry <- NULL
  presid.k <- NULL
  for (i in 1:length(resid)){
    fy.ry[i] <- f.y$y[which(abs(f.y$x-resid[i])==min(abs(f.y$x-resid[i])))]
    presid.k[i] <- sum(resid<resid[i])/length(resid) - sum(resid>resid[i])/length(resid)
  }
  dpresid.dtheta.k <- t(cbind(-2*fy.ry,
                              -2*fy.ry*X))
  list(mod = mod, 
       dl.dtheta = dl.dtheta,
       d2l.dtheta.dtheta = d2l.dtheta.dtheta,
       resid = resid,
       dresid.dtheta = dresid.dtheta,
       presid = presid,
       presid.k= presid.k,
       dpresid.dtheta = dpresid.dtheta,
       dpresid.dtheta.k = dpresid.dtheta.k)
}


##################################### conditional corr: preparing table data
#####################################conditional corr: preparing table data
##################################### conditional corr: preparing table data
preparePowerForCondCorr = function(simNumForRhoStar = 100, fileNameWithAwesomeArray, fileNameStarter = "powerForConditional",
                                   family){
  
  familyLow = tolower(family)
  familyCap = gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", familyLow, perl = TRUE)
  
  cat(familyCap, "\n")
  
  load(fileNameWithAwesomeArray)
  newZ = dimnames(awesomeArray)$newZPoints
  functionsToPlot = dimnames(awesomeArray)$func
  knotsVec = dimnames(awesomeArray)$fit
  censVec = as.numeric(dimnames(awesomeArray)$cens)
  methods = dimnames(awesomeArray)$method
  smallerArray = array(0, c(4, length(newZ), 2, length(functionsToPlot), length(knotsVec), length(censVec)), dimnames = list(cor = c("bias", "coverage", "biasStar", "coverageStar"), newZPoints = newZ, method = c("full", "exp"), func = functionsToPlot, fit = knotsVec, cens = censVec))
  for (ni in functionsToPlot){
    for(ck in censVec){
      trueCorValues = awesomeArray[c("true"), , 1, dimnames(awesomeArray)$method[1], ni, 1, as.character(ck)]
      ### temporary fix
      trueCorValues[trueCorValues < 0] = 0
      ### temporary fix
      corrValuesDataGrid = mapCorValueToParamInCopula(corValues = trueCorValues, family = familyCap)
      extractThetas = corrValuesDataGrid$predTheta
      corStarValues = rep(NA, length(extractThetas)); names(corStarValues) = trueCorValues
      corLocalValues = rep(NA, length(extractThetas))
      for(theta_i in 1:length(extractThetas)){
        thi = extractThetas[theta_i]
        if(thi != "0"){
          trueCorTMP = trueCorrStarSimulated(family = familyCap, theta = as.numeric(thi), independentCensoring = FALSE, censProbForDependCens = (1-round(as.numeric(ck), 1)), simNum = simNumForRhoStar)
          corStarValues[theta_i] = trueCorTMP["PSR population value"]
          corLocalValues[theta_i] = NA
        }
      }
      for(mi in methods){
        for (ki in knotsVec){
          # awesomeArray[, , , mi, ni, as.character(ki), as.character(ck)]
          meanSampleCor = apply(awesomeArray["est", , , mi, ni, as.character(ki), as.character(ck)], 1, mean)
          smallerArray["bias", , mi, ni, as.character(ki), as.character(ck)] = meanSampleCor - trueCorValues
          smallerArray["biasStar", , mi, ni, as.character(ki), as.character(ck)] = meanSampleCor - corStarValues
          subArray = awesomeArray[c("true", "lower", "upper"), , , mi, ni, as.character(ki), as.character(ck)]
          simNumInAwesome = length(dimnames(awesomeArray)$simNum)
          for (sim_i in 1:simNumInAwesome){
            smallerArray["coverage", , mi, ni, as.character(ki), as.character(ck)] = smallerArray["coverage", , mi, ni, as.character(ki), as.character(ck)] + as.numeric((subArray["true", , sim_i] >= subArray["lower", , sim_i]) & (subArray["true", , sim_i] <= subArray["upper", , sim_i]))
            
            smallerArray["coverageStar", , mi, ni, as.character(ki), as.character(ck)] = smallerArray["coverageStar", , mi, ni, as.character(ki), as.character(ck)] + as.numeric((corStarValues >= subArray["lower", , sim_i]) & (corStarValues <= subArray["upper", , sim_i]))
          }
          smallerArray["coverage", , mi, ni, as.character(ki), as.character(ck)] = smallerArray["coverage", , mi, ni, as.character(ki), as.character(ck)]/simNumInAwesome
          smallerArray["coverageStar", , mi, ni, as.character(ki), as.character(ck)] = smallerArray["coverageStar", , mi, ni, as.character(ki), as.character(ck)]/simNumInAwesome
        }
      }
    }
  }
  save(smallerArray, file = paste0(fileNameStarter, ".rda"), compress = TRUE)
}


############################################################################
### find rhoStar (the actual estimate of PSR cor) using simulations ############################################################################
trueCorrStarSimulated = function(family, theta, censoringProb1 = 0, censoringProb2 = 0, independentCensoring = FALSE, censProbForDependCens = 0, simNum = 100000){
  set.seed(1)
  res12 = as.data.frame(rclayton(nSim = simNum, thetaPar = theta, lambda1 = 1, lambda2 = 1, censoringProb1 = censoringProb1, censoringProb2 = censoringProb2, independentCensoring = independentCensoring, censProbForDependCens = censProbForDependCens, family = family))
  psr1 = computePSR(res12$t1, res12$delta1)$PSR
  psr2 = computePSR(res12$t2, res12$delta2)$PSR
  returnVal = cor(psr1, psr2)
  names(returnVal) = "PSR population value"
  returnVal
}


###################################################################
###################################################################
######### compute PSR from censored data
###################################################################
###################################################################
computePSR = function(timeForKM, deltaForKM, time = timeForKM, delta = deltaForKM){
  ### computes probability scale residuals (PSR)
  ### "timeForKM", "deltaForKM" is the data for computing empirical CDF.
  ### "time" is time to event or censoring for which to compute 
  ### "delta" corresponds to "time". It is 1 if not censored and 0 if censored.
  ### by default we can compute PSR for the same values as were use to construct CDF
  
  getCDFValues = function(dataCDF, times){
    ### returns CDF for x and CDF for x- (see definition of PSR)
    ### where x is parameter "times"
    ### if you wonder why the function is implemented with "whlie" loops
    ### it is because I wanted to go through the values of CDF and times
    ### in one loop. The complexity of this loop is O(n) which is most efficient
    if(nrow(dataCDF) != length(unique(dataCDF$time))) stop("CDF has to have unique times")
    
    dataCDF = dataCDF[order(dataCDF$CDF),]
    times = times[order(times)]
    res = data.frame(x = times, CDFx = rep(NA, length(times)), CDFxMinus=NA)
    timesCount = nrow(res)
    cdfCount = nrow(dataCDF)
    while(timesCount >= 1){
      while(dataCDF$time[cdfCount] > times[timesCount] & cdfCount>0){
        cdfCount = cdfCount - 1
      }
      if(cdfCount == 0){
        res$CDFx[timesCount] = 0
        res$CDFxMinus[timesCount] = 0
      }
      if(cdfCount == 1){
        res$CDFx[timesCount] = dataCDF$CDF[cdfCount]
        res$CDFxMinus[timesCount] = 0
      }else{
        res$CDFx[timesCount] = dataCDF$CDF[cdfCount]
        res$CDFxMinus[timesCount] = dataCDF$CDF[cdfCount-1]
      }
      timesCount = timesCount - 1 
    }
    res
  }
  
  oldOrderOfTimeData = data.frame(order = 1:length(time), time = time, delta = delta)
  oldOrderOfTimeData = oldOrderOfTimeData[order(oldOrderOfTimeData$time),]
  
  dataWithKM = myOwnKM(time = time, delta = delta)
  
  if(FALSE){
    ############################# old way
    #KM <- survfit(Surv(time, delta) ~ 1, type="kaplan-meier", conf.type="log")
    shrinkCDF = unique(dataWithKM[, c("time", "CDF")])
    shrinkCDF = shrinkCDF[order(shrinkCDF$time), ]
    cdfValues = getCDFValues(dataCDF = shrinkCDF, times=time)
    if(any(oldOrderOfTimeData$itme != cdfValues$time)) stop("Something's wrong with time ordering")
    cdfValues = cbind(cdfValues, oldOrderOfTimeData[c("order", "delta")])
    cdfValues$PSR = cdfValues$CDFx - cdfValues$delta*(1 - cdfValues$CDFxMinus)
    ############# we make sure to return everything in the ORIGINAL order (the same order)
    ############# as variables time and delta were supplied
    cdfValues = cdfValues[order(cdfValues$order),]
    PSR = cdfValues$PSR
  }else{
    ############################# NEW way
    PSR = dataWithKM$CDF - dataWithKM$delta*(1 - dataWithKM$CDF_M)
  }
  list(PSR = PSR, KM = dataWithKM)
}


###################################################
### generates clayton distr. with package copula
###################################################
rclayton = function(nSim, thetaPar=0.1432, lambda1=1, lambda2=1, censoringProb1 = 0, censoringProb2 = 0, independentCensoring = TRUE, censProbForDependCens = 0, family = "clayton", reverse = TRUE){
  ### Generates Clayton family with marginal exponential distributions
  ### with rate 1 and two respectively. The resulting Kendall's tau is about 0.1
  ### according Shih and Loius
  
  ### Allows two types of censoring: independent (b/w C1 and C2) and dependent (C1 = C2)
  
  ### this was my attempt to generate Clayton's family using conditional distr.
  ### did not work out because for some reason the marginal of t1 was not exp with rate=1
  # if(thetaPar<=1) stop("Argument thetaPar should be greater than one")
  # t2 = rexp(n = nSim, rate = lambda2)
  # u = runif(n = nSim)
  # t1 = rep(NA, nSim)
  # condSurv = ( (1-pexp(t1, rate=lambda1))^(1-thetaPar) + (1-pexp(t2, rate=lambda2))^(1-thetaPar)  - 1)^(thetaPar/(1-thetaPar))*lambda2*exp(lambda2*t2*thetaPar)
  # logExpression  = 1 + ((1-u)/(exp(lambda2*t2*thetaPar)*lambda2))^((1-thetaPar)/thetaPar)  -  (1-pexp(t2, rate=lambda2))^(1-thetaPar)
  # t1 = - (1/((1-thetaPar)*lambda1)) *log(logExpression)
  
  ### we generate this family using the package of  Belzile, Genest, McNeil, Neslehova
  family = tolower(family)
  if(family == "clayton"){
    if(thetaPar == 0){  ###### trying to decide which one is better
      t1t2 = rarchi(n = nSim, "clayton", d = 2, theta = 0)
    }else{
      cop12 = archmCopula(family = "clayton", param = thetaPar)
      t1t2 <- rCopula(nSim, cop12)
    }
  }else{
    if(family != "frank") stop("Only Clayton's and Frank's families were implemented\n")
    if(thetaPar == 0){
      t1t2 = rarchi(n = nSim, "clayton", d = 2, theta = 0)
    }else{
      cop12 = archmCopula(family = "frank", param = thetaPar)
      t1t2 <- rCopula(nSim, cop12)
    }
  }
  if(reverse){ 
    ### if reverse = TRUE then the increasing U will translate
    ###  into increasing X and therefore is the correct way of simulating
    t1 = -log(1 - t1t2[,1])/lambda1
    t2 = -log(1 - t1t2[,2])/lambda2
  }else{
    t1 = -log(t1t2[,1])/lambda1
    t2 = -log(t1t2[,2])/lambda2
  }
  delta1 = delta2 = rep(1, nSim)
  if(independentCensoring){
    if(censoringProb1 > 0 & censoringProb1 < 1){
      beta1 = lambda1/censoringProb1 - lambda1
      censoringTime1 = rexp(n = nSim, rate = 1/beta1)
      delta1 = as.numeric(t1 <= censoringTime1) 
      t1[delta1==0] = censoringTime1[delta1==0]
    }
    if(censoringProb2 > 0 & censoringProb2 < 1){
      beta2 = lambda2/censoringProb2 - lambda2
      censoringTime2 = rexp(n = nSim, rate = 1/beta2)
      delta2 = as.numeric(t2 <= censoringTime2) 
      t2[delta2==0] = censoringTime2[delta2==0]
    }
    ### below we assume fixed censoring: censoring all after a certain percentile:
    ### the percentile  is supplied as -(proportion of censoring):
    ### if it equals -0.3 then after 70th-% percentile we censore everyone:
    if(censoringProb1 < 0 & censoringProb1 > -1){
      censoringTime1 = rep(qexp(1 + censoringProb1, rate = lambda1), length(delta1))
      delta1 = as.numeric(t1 <= censoringTime1) 
      t1[delta1==0] = censoringTime1[delta1==0]
    }
    if(censoringProb2 < 0 & censoringProb2 > -1){
      censoringTime2 = rep(qexp(1 + censoringProb2, rate = lambda2), length(delta2))
      delta2 = as.numeric(t2 <= censoringTime2) 
      t2[delta2==0] = censoringTime2[delta2==0]
    }
  }else{
    #cat("\n <<<censProbForDependCens  ", censProbForDependCens, " >>>\n")
    if(censProbForDependCens > 0 & censProbForDependCens < 1){
      beta3 = 1/censProbForDependCens - 1
      censoringTime3 = rexp(n = nSim, rate = 1/beta3)
      delta1 = as.numeric(t1 <= censoringTime3) 
      delta2 = as.numeric(t2 <= censoringTime3) 
      t1[delta1==0] = censoringTime3[delta1==0]
      t2[delta2==0] = censoringTime3[delta2==0]  
    }
  }
  
  res = cbind(t1, delta1, t2, delta2)
  colnames(res) = c("t1", "delta1", "t2", "delta2")
  res
}

###################################################################
###################################################################
######### computes KM (the one from survival package was hard to use) ###################################################################
###################################################################
myOwnKM = function(time, delta, returnToOriginalOrder = TRUE){
  ### returnToOriginalOrder = TRUE - the order of time values is the same as in original data
  ###            it also returns the delta (event indicator)
  ### returnToOriginalOrder = FALSE - the time is unique and ordered and the resulting data
  ###            frame does not contain delta (event indicator).
  uniqueAndOrderedTime = unique(time)[order(unique(time))]
  if(TRUE){ 
    ###-----------------------------------------------
    ### computes KM (the one from survival package was hard to use)
    ### the one from survival package returned rounded times values,
    ###    which would not allow to perform precise PSR calculations
    #dataKM = data.frame(time=c(1, 2, 3, 5, 2, 6, 6, 6), delta=c(1, 0, 0, 1, 1, 0, 1, 1))
    # dataKM = data.frame(time=time, delta=delta)
    # dataKM = dataKM[order(dataKM$time),]
    # #KM <- survfit(Surv(dataKM$time, dataKM$delta) ~ 1, type="kaplan-meier", conf.type="log")
    # nEvents = tapply(dataKM$delta, dataKM$time, sum)
    # nDrop = tapply(dataKM$delta, dataKM$time, length)
    # atRisk = c(length(dataKM$time), length(dataKM$time) - cumsum(nDrop))[1:length(nDrop)]
    # probForEachTime = (1-nEvents/atRisk)
    nEvents = tapply(delta, time, sum)
    nDrop = tapply(delta, time, length)
    atRisk = c(length(time), length(time) - cumsum(nDrop))[1:length(nDrop)]
    probForEachTime = (1-nEvents/atRisk)
    dataKM = data.frame(time=uniqueAndOrderedTime, nEvents = nEvents, atRisk = atRisk, KM = cumprod(probForEachTime))
    #fit <- survfit(Surv(time, delta) ~ 1)  ### after some time this line should be removed
    ### together with the next line (temporary check for my code against R's)
    # if(length(dataKM$KM) == length(fit$surv)){  ### sometimes myKM is longer than survfit
    #   if( max(abs(dataKM$KM - fit$surv))>0.0001) {
    #     problemData = list(uniqueAndOrderedTime = uniqueAndOrderedTime, df = data.frame(time = time, delta = delta), fit = fit, dataKM = dataKM)
    #     fileName = gsub("[ :]", "_", paste("./RESULTS/problemData", Sys.time(),".rda"))
    #     save(problemData, file = fileName, compress = TRUE)
    #   }
    # }
  }
  # else{
  #   ### not using this for now b/c it gave me trouble on ACCRE
  #   ### it rounded up some of the time values, which caused problems when "stretching" KM
  #   ### over time...
  #   fit <- survfit(Surv(time, delta) ~ 1)
  #   if(length(uniqueAndOrderedTime) != length(fit$n.event)){
  #     problemData = list(uniqueAndOrderedTime = uniqueAndOrderedTime, fitN = fit$n.event, df = data.frame(time = time, delta = delta), survFit = fit)
  #     save(problemData, file = "./RESULTS/problemData.rda", compress = TRUE)
  #     stop("Saved problematic data\n")
  #   }
  #   dataKM = data.frame(time=uniqueAndOrderedTime, nEvents = fit$n.event, atRisk = fit$n.risk, KM = fit$surv)
  # }
  dataKM$CDF = 1 - dataKM$KM
  # cat(nrow(dataKM), "\n")
  # cat(paste(dataKM$CDF, collapse = ", "), "\n")
  dataKM$CDF_M = c(0, dataKM$CDF[1:(nrow(dataKM)-1)])
  if(returnToOriginalOrder){
    rownames(dataKM) = uniqueAndOrderedTime
    ### let's order the output according to the original order of time
    dataKM = dataKM[as.character(time),]
    dataKM$delta = delta
  }
  rownames(dataKM) = 1:nrow(dataKM)
  dataKM
}


###################################################################
######### unadjusted Rho PSRs for right-censored data
###################################################################
unadjusted.CorPSRs <- function(X, Y, deltaX, deltaY){
  l1 = length(X)
  l2 = length(Y)
  l3 = length(deltaX)
  l4 = length(deltaY)
  if (l1 != l2 | l1 != l3 | l1 != l4) 
    stop("Arguments 'X', 'Y', 'deltaX', and 'deltaY' are numeric vectors of equal length.\n'")
  
  resX = computePSRs(time = X, delta = deltaX)
  resY = computePSRs(time = Y, delta = deltaY)
  
  xz.d2l.dtheta.dtheta = diag(length(X), length(unique(X[deltaX == 1])))
  yz.d2l.dtheta.dtheta = diag(length(Y), length(unique(Y[deltaY == 1])))
  
  bundle1 = make_estimating_equations_stute(myKaplanMeier = resX$KM)
  bundle2 = make_estimating_equations_stute(myKaplanMeier = resY$KM)
  
  psrWithStute = corTS.survreg.modified(xresid = resX$PSRs, yresid = resY$PSRs,
                                        xz.dl.dtheta = t(bundle1$dl.dtheta), yz.dl.dtheta = t(bundle2$dl.dtheta),
                                        xz.d2l.dtheta.dtheta = xz.d2l.dtheta.dtheta, yz.d2l.dtheta.dtheta = yz.d2l.dtheta.dtheta,
                                        dxresid.dthetax = bundle1$dpresid.dtheta, dyresid.dthetay = bundle2$dpresid.dtheta,
                                        fisher=FALSE, inverseA = FALSE)
  
  res = c(est = psrWithStute$TS, strerr = sqrt(psrWithStute$varTS/length(X)), P = psrWithStute$pvalTS, lower.CI = psrWithStute$CI_TS[1], upper.CI = psrWithStute$CI_TS[2])
  res
}

############################################### estimating equations proposed 
############################################### by Stute's [1995]
############################################### and explained in Shepherd's [2007]
make_estimating_equations_stute = function(myKaplanMeier){
  orderedTime = myKaplanMeier$time[order(myKaplanMeier$time)]
  orderedDelta =  myKaplanMeier$delta[order(myKaplanMeier$time)]
  orderedUniqueKM = unique(myKaplanMeier$KM[order(myKaplanMeier$time)][orderedDelta == 1])
  
  uniqueTime = unique(myKaplanMeier$time)
  orderedUniqueTime = uniqueTime[order(uniqueTime)]
  orderedUniqueFailureTime = unique(orderedTime[orderedDelta == 1])
  
  N = length(orderedTime)
  uniqueN = length(orderedUniqueTime)
  uniqueFailureN = length(orderedUniqueFailureTime)
  
  orderMatrix = matrix(c(1:N, (1:N)[order(myKaplanMeier$time)]), ncol=2)
  colnames(orderMatrix) = c("new", "old")
  originalOrder = orderMatrix[order(orderMatrix[,2]), 1]
  
  dpresid.dtheta = matrix(0, nrow = uniqueFailureN, ncol = N)
  colnames(dpresid.dtheta) = orderedTime
  rownames(dpresid.dtheta) = orderedUniqueFailureTime
  condition1 = (orderedTime == orderedUniqueFailureTime[1]) | ((orderedTime <= orderedUniqueFailureTime[1]) & (orderedDelta == 0))
  dpresid.dtheta[1, (1:N)[condition1]] = 1
  
  phi_ji = matrix(0, nrow = uniqueFailureN, ncol = N)
  colnames(phi_ji) = orderedTime
  rownames(phi_ji) = orderedUniqueFailureTime
  indices = (1:N)[orderedTime <= orderedUniqueFailureTime[1]]
  phi_ji[1, indices] = 1
  
  for (j in 2:uniqueFailureN){
    condition2 = (orderedTime == orderedUniqueFailureTime[j] & (orderedDelta == 1))
    condition3 = ((orderedTime > orderedUniqueFailureTime[j-1]) & (orderedTime <= orderedUniqueFailureTime[j]) & (orderedDelta == 0))
    dpresid.dtheta[j, (1:N)[condition2 | condition3]] = 1
    dpresid.dtheta[j-1, (1:N)[condition2]] = 1
    indices = (1:N)[orderedTime <= orderedUniqueFailureTime[j]]
    phi_ji[j, indices] = 1
  }
  condition4 = (orderedTime > orderedUniqueFailureTime[uniqueFailureN] & (orderedDelta == 0))
  dpresid.dtheta[uniqueFailureN, (1:N)[condition4]] = 1
  orderedDeltaMatr = matrix(orderedDelta, nrow = uniqueFailureN, ncol = N, byrow=TRUE)
  
  H = H0 = H1 = rep(NA, N)
  for (i in 1:N){
    H[i] = mean(orderedTime <= orderedTime[i])
    H0[i] = mean((orderedTime <= orderedTime[i])*(1 - orderedDelta))
    H1[i] = mean((orderedTime <= orderedTime[i])*orderedDelta)
  }
  
  H0dv = diff(c(0, H0))
  H1dw = diff(c(0, H1))
  H1dwPerY = H0dvPerY = rep(0, N)
  H0dvPerY[orderedTime * (1 - orderedDelta) == orderedTime] = H0dv[orderedTime * (1 - orderedDelta) == orderedTime]
  H1dwPerY[orderedTime * orderedDelta == orderedTime] = H1dw[orderedTime * orderedDelta == orderedTime]
  
  ################ making sure that 1-H is never zero: Bryan's suggestion
  Hadj = H
  Hadj[Hadj == 1] = .999999
  
  multiplier = 1/(1 - Hadj) 
  multiplier[H==1] = 0   ### this porbably fixes the the problem of division by 0 anyway
  
  gamma0 = exp(cumsum(c(0, (H0dvPerY * multiplier)  )))[1:N]
  Vji = gamma_j2 = gamma_j1 = phi_ji*0
  vValue = wValue  = orderedTime
  
  for(j in 1:nrow(gamma_j1)){    
    for(i in 1:ncol(gamma_j1)){  
      indicatorForGamma1 = as.numeric((orderedTime[i] < wValue) & phi_ji[j, ]) 
      gamma_j1[j, i] =  multiplier[i]  *  sum(indicatorForGamma1  *  gamma0  *  H1dwPerY)
    }
  } 
  
  for(j in 1:nrow(gamma_j1)){    
    for(i in 1:ncol(gamma_j1)){  
      indicatorForGamma2 = as.numeric((vValue < orderedTime[i]))
      gamma_j2[j, i] =   sum( multiplier * indicatorForGamma2 * gamma_j1[j, ] * H0dvPerY)
      Vji[j, i] = phi_ji[j, i] * gamma0[i] * orderedDelta[i] + gamma_j1[j, i] * (1-orderedDelta[i]) - gamma_j2[j, i]
    }
  }   
  res = Vji[, originalOrder] - (1 - matrix(orderedUniqueKM, nrow = uniqueFailureN, ncol = N))
  
  list(dl.dtheta = res, dpresid.dtheta = dpresid.dtheta[, originalOrder])
}

###################################################################
######### compute PSRs for censored data
###################################################################
computePSRs = function(time, delta){
  dataWithKM = myOwnKM(time = time, delta = delta)
  PSRs = dataWithKM$CDF - dataWithKM$delta*(1 - dataWithKM$CDF_M)
  list(PSRs = PSRs, KM = dataWithKM)
}

partial.corPSRs <- function(modX, modY, likelihood = "partial"){
  # res = c(est = psrWithStute$TS, strerr = sqrt(psrWithStute$varTS/length(X)), P = psrWithStute$pvalTS, lower.CI = psrWithStute$CI_TS[1], upper.CI = psrWithStute$CI_TS[2])
  # res
  #
  ### partial
  ##################################### identify the type of model
  modelX = as.character(modX$call)[1]
  modelY = as.character(modY$call)[1]
  if(modelX != modelY){
    stop("Arguments 'modX' and 'modY' should have the same type of model: parametric or Cox regression.")
  }
  
  ###### SURVREG
  if(modelX == "survreg"){
    scoreRes1 = mestimation.survregPSRs(object = modX, inverseA = TRUE)
    scoreRes2 = mestimation.survregPSRs(object = modY, inverseA = TRUE)
    corXY = corTS.survreg.modified(xresid = scoreRes1$presid, yresid = scoreRes2$presid,
                                   xz.dl.dtheta = t(scoreRes1$dl.dtheta),
                                   yz.dl.dtheta = t(scoreRes2$dl.dtheta),
                                   xz.d2l.dtheta.dtheta = scoreRes1$d2l.dtheta.dtheta,
                                   yz.d2l.dtheta.dtheta = scoreRes2$d2l.dtheta.dtheta,
                                   dxresid.dthetax = scoreRes1$dpresid.dtheta,
                                   dyresid.dthetay = scoreRes2$dpresid.dtheta, inverseA = TRUE)
  }
  
  ###### COX
  if(any(modelX %in% c("cph", "coxph"))){
    likelDict = c("partial" = "old", "full" = "full")
    if(is.na(likelDict[likelihood])){
      stop("Argument 'likelihood' can be set to either 'partial' or 'full'")
    }
    corXY = mestimation.coxph(object1 = modX, object2 = modY, method = likelDict[likelihood])
  }
  res = c(est = corXY$TS, strerr = sqrt(corXY$varTS/modX$n), P = corXY$pvalTS, lower.CI = corXY$CI_TS[1], upper.CI = corXY$CI_TS[2])
  res
  
}

############################################### OLD M-estimation for Cox regression
mestimation.coxph_Old <- function(object, ...){
  ### for non-time dependent covariates
  
  stopifnot(require("sandwich"))
  #stopifnot(require("PResiduals"))
  
  ################### for objectects of type survreg.scores
  ################### returns all the necessary stuff for M-estimation of varience
  ################### of PSR correlation
  ### mod - is the fitted survival model, for which we need M-estimation
  
  ### the following function (from package sandwich) is helpful for finding dl.dtheta
  ###     I yet need to figure out where it gets it the stuff probably from 
  ###     residuals but I don't know where
  ###     ask Cole how to see the source code for:    estfun.survreg
  
  # modSum = summary(object)
  ### the body of function estfun.survreg()
  
  ######################################## dl.dtheta (matrix psi)
  if(FALSE){ #### the following is the body of estfun() that part of sandwich library
    stopifnot(require("survival"))
    xmat <- model.matrix(object)
    wts <- weights(object)
    if(is.null(wts)) wts <- 1
    res <- residuals(object, type = "matrix")
    rval <- as.vector(res[,"dg"]) * wts * xmat
    if(NROW(object$var) > length(coef(object))) {
      rval <- cbind(rval, res[,"ds"])
      colnames(rval)[NCOL(rval)] <- "Log(scale)"
    }
    dl.dtheta = t(rval)
  }
  
  dl.dtheta = t(estfun(object))   ### similar to mean(.) function from sandwich library
  
  ######################################## d2l.dtheta.dtheta (matrix A)
  d2l.dtheta.dtheta = solve(object$var)
  
  ######################################## prob. scale residuals (PSR)
  presid = presid.coxph(object)
  
  ######################################## PSR derivatives by theta (dpresid.dtheta)
  Y = object$y[,1]
  nSubj = length(Y)
  delta = object$y[,2]
  orderedY = Y[order(Y)]
  orderedDelta = delta[order(Y)]
  
  X1 = model.matrix(object)
  if(dim(X1)[2] == 1){
    orderedX1 = as.matrix(X1[order(Y),], ncol = 1)
  }else{
    orderedX1 = X1[order(Y),]
  }
  Xbeta =  X1 %*% matrix(object$coefficients, nrow = length(object$coefficients))
  expXbeta = exp(Xbeta)
  orderedExpXbeta = expXbeta[order(Y)]
  
  nParam = dim(object$var)[1]
  dpresid.dtheta = matrix(NA, nrow = nParam, ncol = nSubj)
  
  if(FALSE){
    ######################################## PSR derivatives by beta: old way
    martResid <- residuals(object, type="martingale")
    #multiplier = (1+delta)*exp(martResid - delta + Xbeta)
    multiplier = -(1+delta)*(martResid - delta)*exp(martResid - delta)
    if(nParam>=1){
      for (i in 1:nParam){
        dpresid.dtheta[i, ] = multiplier * X1[, i]
      }
    }
  }else{
    ######################################## PSR derivatives by beta  
    H0 = basehaz(object, centered=FALSE)
    H0_minus = hazardStretcher_take_two(H0, orderedY)
    H0_minus$delta = orderedDelta
    uniqueBlHaz = H0_minus[H0_minus$delta == 1,]
    uniqueBlHaz = unique(uniqueBlHaz[, c("hazard", "time")])
    timeDiff = diff(c(0, uniqueBlHaz[, "time"]))
    lambdas = diff(c(0, uniqueBlHaz[, "hazard"]))/timeDiff
    element1 = exp(-H0_minus[, "hazard"] * orderedExpXbeta)
    element2 = exp(-H0_minus[, "Hminus"] * orderedExpXbeta)
    if(nParam>=1){
      for (j in 1:nParam){
        dpresid.dtheta[j, ] = (element1*H0_minus[, "hazard"] + orderedDelta * element2 * H0_minus[, "Hminus"]) * orderedExpXbeta * orderedX1[, j]
      }
    }
    
    # H0 = basehaz(object, centered=FALSE)
    # dpresid.dtheta = matrix(0, nrow = nParam, ncol = nSubj)
    # element1 = exp(-H0[,1] * orderedExpXbeta)
    # element2 = exp(-c(0, H0[1:(nrow(H0)-1),1]) * orderedExpXbeta)
    # element1[H0[,1] == 0] = 1
    # element2[c(0, H0[1:(nrow(H0)-1),1]) == 0] = 1  
    # ######################################## PSR derivatives by beta
    # if(nParam>=1){
    #   for (j in 1:nParam){
    #     dpresid.dtheta[j, ] = (element1*H0[,1] + orderedDelta * element2*c(0, H0[1:(nrow(H0)-1),1])) * orderedExpXbeta * orderedX1[, j]
    #   }
    # }
    
    orderMatrix = matrix(c(1:nSubj, (1:nSubj)[order(Y)]), ncol=2)
    colnames(orderMatrix) = c("new", "old")
    originalOrder = orderMatrix[order(orderMatrix[,2]), 1]
    dpresid.dtheta = dpresid.dtheta[, originalOrder]
  }
  
  res = list(dl.dtheta = dl.dtheta,
             d2l.dtheta.dtheta = d2l.dtheta.dtheta,
             resid = NULL,
             dresid.dtheta = NULL,
             presid = presid,
             presid.k= NULL,
             dpresid.dtheta = dpresid.dtheta,
             dpresid.dtheta.k = NULL)
  res
}

############################################### M-estimation for Cox using full likelihood
mestimation.coxph.full.likelihood <- function(object, time, continuous = TRUE, ...){
  ### in general, argument time can be extracted from object: time = object$y[, 1]
  ### but there were numerical issue with it 
  ################################# M-estimation for Cox  
  Y = object$y[,1]
  nSubj = length(Y)
  delta = object$y[,2]
  resid <- residuals(object, type="martingale")
  orderedY = Y[order(Y)]
  orderedDelta = delta[order(Y)]
  orderedUniqueFailureTimes = unique(orderedY[orderedDelta == 1])
  orderedFailureTimesWithZero = c(0, orderedUniqueFailureTimes)
  nHazardParam = length(orderedUniqueFailureTimes)
  orderedMartRes = resid[order(Y)]
  
  orderMatrix = matrix(c(1:nSubj, (1:nSubj)[order(Y)]), ncol=2)
  colnames(orderMatrix) = c("new", "old")
  originalOrder = orderMatrix[order(orderMatrix[,2]), 1]
  
  X1 = model.matrix(object)
  if(dim(X1)[2] == 1){
    orderedX1 = as.matrix(X1[order(Y),], ncol = 1)
  }else{
    orderedX1 = X1[order(Y),]
  }
  Xbeta =  X1 %*% matrix(object$coefficients, nrow = length(object$coefficients))
  expXbeta = exp(Xbeta)
  orderedExpXbeta = expXbeta[order(Y)]
  
  nParam = dim(object$var)[1]
  if(is.null(nParam)){
    nParam = 1
  }
  
  H0 = basehaz(object, centered=FALSE)
  # uniqueBlHaz = H0[, ]
  # timeDiff = diff(c(0, uniqueBlHaz[, "time"]))
  # lambdas = diff(c(0, uniqueBlHaz[, "hazard"]))/timeDiff
  ### the following code is for "discrete" time,  but it turns out we need it for everyting
  # #if(!continuous){
  #   H0$delta = 1
  #   H0_minus = hazardStretcher(H0)
  #   H0_minus$delta = orderedDelta
  #   H0 = H0_minus[, c("hazard", "time")]
  # #}
  
  #H0_minus = hazardStretcher_take_two(H0, orderedY)
  H0_minus = hazardStretcher_take_two(H0, Y)
  H0_minus$delta = orderedDelta
  uniqueBlHaz = H0_minus[H0_minus$delta == 1,]
  uniqueBlHaz = unique(uniqueBlHaz[, c("hazard", "time")])
  timeDiff = diff(c(0, uniqueBlHaz[, "time"]))
  lambdas = diff(c(0, uniqueBlHaz[, "hazard"]))/timeDiff
  
  ######################################## dl.dtheta
  ######################################## dl.dtheta
  ######################################## dl.dtheta
  dl.dtheta = matrix(NA, nrow = nParam + nHazardParam, ncol = nSubj)
  ######################################## dl.dtheta   for beta:
  for (j in 1:nParam){
    dl.dtheta[j, ] = (orderedDelta - H0_minus[,"hazard"]*orderedExpXbeta)*orderedX1[,j]
  }  
  ######################################## dl.dtheta   for lambda:
  dl.dlambda = matrix(0, nrow=nHazardParam, ncol = nSubj)
  for (k in 1:nHazardParam){
    hadEventAtT = (orderedDelta == 1 & orderedY == orderedFailureTimesWithZero[k+1])
    hadEventAfterT = (orderedY > orderedFailureTimesWithZero[k+1] | (orderedDelta == 0 & orderedY == orderedFailureTimesWithZero[k+1]))
    dl.dlambda[k, hadEventAtT] = orderedDelta[hadEventAtT]/lambdas[k] - timeDiff[k] * orderedExpXbeta[hadEventAtT]
    dl.dlambda[k, hadEventAfterT] = -timeDiff[k] * orderedExpXbeta[hadEventAfterT]
  }
  ######################################## merge dl.dtheta for betand and for lambda:
  dl.dtheta[(nParam+1):(nParam+nHazardParam), ] = dl.dlambda
  
  ######################################## d2l.dtheta.dtheta (matrix A)
  ######################################## d2l.dtheta.dtheta (matrix A)
  ######################################## d2l.dtheta.dtheta (matrix A)
  d2l.dtheta.dtheta = diag(nParam + nHazardParam)
  ######################################## d2l.dtheta.dtheta  for beta/beta:
  for (i in 1:nParam){
    for (j in 1:i){
      element = H0_minus[,"hazard"] * orderedExpXbeta * orderedX1[,i] * orderedX1[,j]
      d2l.dtheta.dtheta[i, j] = d2l.dtheta.dtheta[j, i] = -sum(element)
    }
  }  
  ######################################## d2l.dtheta.dtheta  for lambda/lambda:
  for (k in (1:nHazardParam)){
    hadEventAtT = (orderedY == orderedFailureTimesWithZero[k+1])
    ############################################
    element = orderedDelta/(lambdas[k]^2)
    d2l.dtheta.dtheta[nParam+k, nParam+k] = -sum(element[hadEventAtT])
  }
  
  ######################################## d2l.dtheta.dtheta  for beta/lambda:
  for (j in 1:nParam){
    for (k in 1:nHazardParam){
      greaterOrEqual = (orderedY >= orderedFailureTimesWithZero[k+1])
      element = timeDiff[k] * orderedExpXbeta * orderedX1[, j]
      d2l.dtheta.dtheta[j, nParam + k] = d2l.dtheta.dtheta[nParam + k, j] = -sum(element[greaterOrEqual])
    }
  }  
  
  ######################################## prob. scale residuals (PSR)
  presid = presid.coxph(object, continuous = continuous)
  orderedPresid = presid[order(Y)]         #### ORDERING BY Y
  
  ######################################## PSR derivatives by beta and lambda_k
  ######################################## PSR derivatives by beta and lambda_k
  dpresid.dtheta = matrix(0, nrow = nParam + nHazardParam, ncol = nSubj)
  if(continuous){
    ######################################## PSR derivatives by beta  
    expression1 = (1 + orderedDelta)*exp(orderedMartRes - orderedDelta)*(orderedDelta-orderedMartRes)
    dpresid.dtheta[1:nParam, ] = matrix(expression1, nrow = nParam, ncol = nSubj, byrow = TRUE) * t(orderedX1)
    
    ######################################## PSR derivatives by lambda_k
    expression2 = (1 + orderedDelta)*exp(orderedMartRes - orderedDelta) * orderedExpXbeta
    for(k in 1:nHazardParam){
      condGreaterOrEqual = (orderedY >= orderedFailureTimesWithZero[k+1])
      dpresid.dtheta[nParam + k, condGreaterOrEqual]  = timeDiff[k] * expression2[condGreaterOrEqual]
    }
  }else{
    element1 = exp(-H0_minus[, "hazard"] * orderedExpXbeta)
    # element2 = exp(-c(0, H0[1:(nrow(H0)-1),1]) * orderedExpXbeta)
    element2 = exp(-H0_minus[, "Hminus"] * orderedExpXbeta)
    # element1[H0[,1] == 0] = 1
    # element2[c(0, H0[1:(nrow(H0)-1),1]) == 0] = 1
    element3 = element1 * orderedExpXbeta
    element4 = element2 * orderedExpXbeta
    
    ######################################## PSR derivatives by beta  
    if(nParam>=1){
      for (j in 1:nParam){
        dpresid.dtheta[j, ] = (element1*H0_minus[, "hazard"] + orderedDelta * element2 * H0_minus[, "Hminus"]) * orderedExpXbeta * orderedX1[, j]
      }
    }
    
    ######################################## PSR derivatives by lambda_k
    for(k in 1:nHazardParam){
      # condGreater = (orderedY > orderedFailureTimesWithZero[k+1])
      # condEqual = (orderedY == orderedFailureTimesWithZero[k+1])
      # dpresid.dtheta[nParam + k, condGreater]  = (timeDiff[k] * (element3 + orderedDelta * element4) )[condGreater]
      # dpresid.dtheta[nParam + k, condEqual]  = (timeDiff[k] * element3)[condEqual]
      condGreaterOrEqual = (orderedY >= orderedFailureTimesWithZero[k+1])
      dpresid.dtheta[nParam + k, condGreaterOrEqual]  = (timeDiff[k] * (element3 + orderedDelta * element4) )[condGreaterOrEqual]
    }
  }
  
  res = list(dl.dtheta = dl.dtheta[, originalOrder],
             d2l.dtheta.dtheta = d2l.dtheta.dtheta,
             resid = NULL,
             dresid.dtheta = NULL,
             presid = orderedPresid[originalOrder],
             presid.k= NULL,
             dpresid.dtheta = dpresid.dtheta[, originalOrder],
             dpresid.dtheta.k = NULL)
  res
}

############################################### estimating equations for  
############################################### parametric regressions
mestimation.survregPSRs <- function(object, inverseA = TRUE, beta = NULL, ...){
  stopifnot(require("sandwich"))
  ################### for objectects of type survreg.scores
  ################### returns the M-estimation stuff
  ### object - is the fitted survival model
  
  if(class(object$y) == "matrix"){
    Y = exp(object$y[,1])
    logY = object$y[,1]
  }
  if(class(object$y) == "Surv"){
    Y = object$y[,1]    
    logY = log(object$y[,1])
  }
  
  ### similar to mean(.) function from sandwich library
  dl.dtheta = t(sandwich::estfun(object))
  
  ######################################## d2l.dtheta.dtheta (matrix A)
  if(inverseA){
    d2l.dtheta.dtheta = object$var
  }else{
    d2l.dtheta.dtheta = solve(object$var)
  }
  
  ######################################## prob. scale residuals (PSR)
  ### When presidTMP.survreg() is merged with presid.survreg() in PResiduals
  ### Replace presid = presidTMP.survreg(object) with presid = presid(object)
  #presid = presid(object)
  presid = presidTMP.survreg(object)
  
  ######################################## prepare stuff for
  ######################################## PSR derivatives by theta (dpresid.dtheta)
  X1 = model.matrix(object)
  if(is.null(beta)){
    Xbeta =  object$linear.predictors
  }else{
    Xbeta =  X1 %*% matrix(beta, ncol=1)
  }
  nSubj = length(Xbeta)
  nParam = dim(object$var)[1]
  delta = object$y[,2]
  dpresid.dtheta = matrix(NA, nrow = nParam, ncol = nSubj)
  
  switch(object$dist,
         ######################################## exponential
         exponential = {
           #Y = exp(object$y[,1])
           dpresid.dtheta[1,] = -(1+delta)*Y*exp(-exp(-Xbeta)*Y - Xbeta)
           if(nParam>=2){
             for (i in 2:nParam){
               dpresid.dtheta[i, ] = dpresid.dtheta[1,] * X1[, i]
             }
           }
         },
         
         ######################################## weibull
         weibull = {
           scale = object$scale
           #Y = exp(object$y[,1])
           dpresid.dtheta[1,] = -(1+delta)*(Y^(1/scale))*exp(-(Y*exp(-Xbeta))^(1/scale) - Xbeta/scale)/scale
           if(nParam-1 >= 2){
             for (i in 2:(nParam-1)){
               dpresid.dtheta[i, ] = dpresid.dtheta[1,] * X1[, i]
             }
           }
           dpresid.dtheta[nParam, ] = -(1+delta)*(log(Y) - Xbeta)*(Y^(1/scale))*exp(-(Y*exp(-Xbeta))^(1/scale) - Xbeta/scale - log(scale))
         },
         
         ######################################## log-logistic
         loglogistic = {
           #Y = exp(object$y[,1])
           gamma = (1/object$scale)
           monsterTerm = (Y^gamma)*exp(-Xbeta*gamma)
           dpresid.dtheta[1,] = -gamma*(1+delta)*monsterTerm/((1+monsterTerm)^2)
           if(nParam-1 >= 2){
             for (i in 2:(nParam-1)){
               dpresid.dtheta[i, ] = dpresid.dtheta[1,] * X1[, i]
             }
           }
           dpresid.dtheta[nParam, ] = -dpresid.dtheta[1,] * (log(Y) - Xbeta)
         },
         
         ######################################## log-normal
         lognormal = {
           # logY = object$y[,1]
           scale = summary(object)$scale
           dpresid.dtheta[1,] = -(1+delta)*dnorm(logY, mean = Xbeta, sd = scale)
           if(nParam-1 >= 2){
             for (i in 2:(nParam-1)){
               dpresid.dtheta[i, ] = dpresid.dtheta[1,] * X1[, i]
             }
           }
           dpresid.dtheta[nParam, ] = (1+delta)* (-(logY - Xbeta)*dnorm(logY, mean = Xbeta, sd = scale)  +  pnorm(logY, mean = Xbeta, sd = scale)  -  plnorm(exp(logY), meanlog=Xbeta, sdlog = scale, log.p=FALSE))
         },
         ### Stop when unknown distribution:
         stop("Unhandled dist", object$dist)
  ) ### end of switch
  
  res = list(dl.dtheta = dl.dtheta,
             d2l.dtheta.dtheta = d2l.dtheta.dtheta,
             resid = NULL,
             dresid.dtheta = NULL,
             presid = presid,
             presid.k= NULL,
             dpresid.dtheta = dpresid.dtheta,
             dpresid.dtheta.k = NULL)
}

presidTMP.survreg <- function(object, ...){
  ### Eventually this function has to be merged with an existing function in PResiduals
  ###   which is called presid.survreg()
  ### The current version of presidTMP.survreg() in Presiduals does not support log-logistic and log-normal distributions
  
  delta <- object$y[,2]
  time <- object$y[,1]
  
  ### In older R version for exp, weilbull, loglogistic, and lognormal
  ###     object$y  was  a matrix and was equal to log(time) instead of time
  if(class(object$y) == "matrix"){
    time = exp(object$y[,1])
  }
  
  switch(object$dist,
         weibull = {
           prob <- pweibull(time, shape=1/summary(object)$scale,
                            scale=exp(object$linear.predictors),
                            lower.tail=TRUE, log.p=FALSE)
           prob + delta*(prob - 1)
         },
         
         exponential = {
           ### should time be exp(time)?  I am pretty sure about this.
           # prob <- pexp(time, rate=1/exp(object$linear.predictors),
           prob <- pexp(time, rate=1/exp(object$linear.predictors),
                        lower.tail=TRUE, log.p=FALSE)
           prob + delta*(prob - 1)
         },
         
         gaussian = {
           prob <- pnorm(time, mean=object$linear.predictors,
                         sd=summary(object)$scale, lower.tail=TRUE,
                         log.p=FALSE)
           prob + delta*(prob - 1)
         },
         
         logistic = {
           prob <- plogis(time, location=object$linear.predictors,
                          scale=summary(object)$scale, lower.tail=TRUE,
                          log.p=FALSE)
           prob + delta*(prob - 1)
         },
         
         ######### Svetlana's update
         loglogistic = {
           gamma = (1/object$scale)
           monsterTerm = (time^gamma)*exp(-object$linear.predictors*gamma)
           ### presid = 1 - (1+delta)*1/(1+monsterTerm)   ### PSR using survival probability
           prob = 1 - 1/(1+monsterTerm)
           prob + delta*(prob - 1)
         },
         
         lognormal = {
           ### should time be exp(time)?
           # prob <- plnorm(exp(time), meanlog=object$linear.predictors,
           prob <- plnorm(time, meanlog=object$linear.predictors,
                          sdlog=summary(object)$scale, lower.tail=TRUE,
                          log.p=FALSE)
           prob + delta*(prob - 1)
         },
         stop("Unhandled dist", object$dist))
}

##############################################################################################
############################################### compute conditional PSRs correlation
##############################################################################################
conditional.corPSRs = function(modX, modY, z, newZ, numKnots = 0){
  ##################################### identify the type of model
  modelX = as.character(modX$call)[1]
  modelY = as.character(modY$call)[1]
  if(modelX != modelY){
    stop("Arguments 'modX' and 'modY' should have the same type of model: coxph().")
  }
  if(modelX != "coxph"){
    stop("Only 'coxph()' with full likelihood esitmating equations is currenlty implemented.")
  }
  ###### COXPH
  scoreResX = mestimation.coxph.full.likelihood(object = modX, continuous = FALSE)
  scoreResY = mestimation.coxph.full.likelihood(object = modY, continuous = FALSE)
  
  resultXY <- prodTS.new.SKE(x.resid = scoreResX$presid, y.resid = scoreResY$presid,
                             z = z, newZ = newZ, numKnots = numKnots,
                             xres2.method="model", yres2.method="model",
                             xz.dl.dtheta = t(scoreResX$dl.dtheta),
                             yz.dl.dtheta = t(scoreResY$dl.dtheta),
                             xz.d2l.dtheta.dtheta = scoreResX$d2l.dtheta.dtheta,
                             yz.d2l.dtheta.dtheta = scoreResY$d2l.dtheta.dtheta,
                             dxresid.dtheta = scoreResX$dpresid.dtheta,
                             dyresid.dtheta = scoreResY$dpresid.dtheta, debugMode = FALSE)
  
  res = resultXY$pointEstAndCIs
  names(res) = c("est", "lower.CI", "upper.CI")
  res
  
}


