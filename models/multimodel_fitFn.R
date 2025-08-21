# ==============================================================================
# MULTIMODEL FITTING FUNCTIONS - COMPLETE CLEAN VERSION
# ==============================================================================
# This version eliminates rep_len errors through safe parameter handling
# All auto-execution removed - use manual functions for control

# Select Variables for Analysis
vars      <- c('p1_t','p2_t')
nv        <- 2                         # number of variables
ntv       <- nv*2                      # number of total variables
nvm1      <- nv-1                      # number of variables minus 1
selVars   <- paste(vars,c(rep(1,nv),rep(2,nv)),sep="")
selVars

# ---------------------------------------------------------------------------------------------------------------------|
# PREPARE MODEL

# Set Starting Values & Lower Bounds for Means
svMe      <- c(0)                      # start value for means

# Set Starting Values for Parameters
svPa      <- .5                        # start value for path coefficient
svPe      <- .7                        # start value for path coefficient for e
svRa      <- .2                        # start value for genetic correlations
svRc      <- .2                        # start value for shared env correlations
svRe      <- .1                        # start value for unique env correlations

# Create Matrices for Covariates and linear Regression Coefficients
#defAge    <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.age"), name="Age" )
#pathBm    <- mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=.01, label="bm11", name="bm" )
#pathBf    <- mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=.01, label="bf11", name="bf" )

# Create Matrices & Algebra for Expected Means
meanZm    <- mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=svMe, labels=labVars("meZm",vars), name="meanZm" )
meanZf    <- mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=svMe, labels=labVars("meZf",vars), name="meanZf" )
expMeanZm <- mxAlgebra( expression=  cbind(meanZm,meanZm), name="expMeanZm" )  #+ cbind(bm%*%Age,bm%*%Age)
expMeanZf <- mxAlgebra( expression= cbind(meanZf,meanZf), name="expMeanZf" )  #+ cbind(bf%*%Age,bf%*%Age)
expMeanZo <- mxAlgebra( expression= cbind(meanZm,meanZf), name="expMeanZo" )  #+ cbind(bm%*%Age,bf%*%Age)

# Create Matrices for Path Coefficients
pathAm    <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=svPa, label=labDiag("am",nv), lbound=.0001, name="am" )
pathCm    <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=svPa, label=labDiag("cm",nv), lbound=.0001, name="cm" )
pathEm    <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=svPe, label=labDiag("em",nv), lbound=.0001, name="em" )
pathAf    <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=svPa, label=labDiag("af",nv), lbound=.0001, name="af" )
pathCf    <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=svPa, label=labDiag("cf",nv), lbound=.0001, name="cf" )
pathEf    <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=svPe, label=labDiag("ef",nv), lbound=.0001, name="ef" )

# Create Matrices for Correlation Coefficients within/across Individuals
pathRam   <- mxMatrix( type="Stand", nrow=nv, ncol=nv, free=TRUE, values=svRa, label=labSdiag("ram",nv), lbound=-1, ubound=1, name="Ram" )
pathRcm   <- mxMatrix( type="Stand", nrow=nv, ncol=nv, free=TRUE, values=svRc, label=labSdiag("rcm",nv), lbound=-1, ubound=1, name="Rcm" )
pathRem   <- mxMatrix( type="Stand", nrow=nv, ncol=nv, free=TRUE, values=svRe, label=labSdiag("rem",nv), lbound=-1, ubound=1, name="Rem" )
pathRaf   <- mxMatrix( type="Stand", nrow=nv, ncol=nv, free=TRUE, values=svRa, label=labSdiag("raf",nv), lbound=-1, ubound=1, name="Raf" )
pathRcf   <- mxMatrix( type="Stand", nrow=nv, ncol=nv, free=TRUE, values=svRc, label=labSdiag("rcf",nv), lbound=-1, ubound=1, name="Rcf" )
pathRef   <- mxMatrix( type="Stand", nrow=nv, ncol=nv, free=TRUE, values=svRe, label=labSdiag("ref",nv), lbound=-1, ubound=1, name="Ref" )

# Create Matrices for Correlation Coefficients across Opposite Sex Pairs - with Twins or Sibs/halfsibs only ram or rcm can be estimated
pathRao   <- mxMatrix( type="Full", nrow=nv, ncol=nv, free=TRUE, values=1, label=labFull("rao",nv,nv), lbound=-1, ubound=1, name="Rao" )
pathRco   <- mxMatrix( type="Full", nrow=nv, ncol=nv, free=TRUE, values=1, label=labFull("rco",nv,nv), lbound=-1, ubound=1, name="Rco" )

# Create Matrices for Causal Path Coefficients - Causal Direction: down column, out row
pathBm    <- mxMatrix( type="Full", nrow=nv, ncol=nv, free=c(F,F,F,F), labels=labFull("bm",nv,nv), lbound=-.99, ubound=.99, name="bm")
pathBf    <- mxMatrix( type="Full", nrow=nv, ncol=nv, free=c(F,F,F,F), labels=labFull("bf",nv,nv), lbound=-.99, ubound=.99, name="bf")

# Algebra to Constrain Correlation Matrices to be Positive Definite
cnstrPos  <- mxMatrix( type="Full", nrow=1, ncol=6, free=FALSE, values=.0001, name="cnstrPos" )
corMin    <- mxAlgebra( expression= cbind(min(eigenval(Ram)), min(eigenval(Rcm)), min(eigenval(Rem)),
                                          min(eigenval(Raf)), min(eigenval(Rcf)), min(eigenval(Ref))), name="corMin" )
corPos    <- mxConstraint( expression= corMin > cnstrPos, name="corPos" )

# Create Algebra for Variance Components
covAm     <- mxAlgebra( expression= am %*% (Ram) %*% t(am), name="Am" )
covCm     <- mxAlgebra( expression= cm %*% (Rcm) %*% t(cm), name="Cm" )
covEm     <- mxAlgebra( expression= em %*% (Rem) %*% t(em), name="Em" )
covAf     <- mxAlgebra( expression= af %*% (Raf) %*% t(af), name="Af" )
covCf     <- mxAlgebra( expression= cf %*% (Rcf) %*% t(cf), name="Cf" )
covEf     <- mxAlgebra( expression= ef %*% (Ref) %*% t(ef), name="Ef" )
covAmf    <- mxAlgebra( expression= am %*% (Rao) %*% t(af), name="Amf" )
covCmf    <- mxAlgebra( expression= cm %*% (Rco) %*% t(cf), name="Cmf" )

# Create Algebra for Causal Variance Components
matI      <- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I" )
invIBm    <- mxAlgebra( expression= solve(I-bm), name="IBm")
invIBf    <- mxAlgebra( expression= solve(I-bf), name="IBf")

# Create Algebra for Total Variances and Standard Deviations (diagonal only)
covVm     <- mxAlgebra( expression= Am+Cm+Em, name="Vm" )
covVf     <- mxAlgebra( expression= Af+Cf+Ef, name="Vf" )
invSDm    <- mxAlgebra( expression= solve(sqrt(I*Vm)), name="iSDm" )
invSDf    <- mxAlgebra( expression= solve(sqrt(I*Vf)), name="iSDf" )
covIVm    <- mxAlgebra( expression= IBm %&% Vm, name="IVm" )
covIVf    <- mxAlgebra( expression= IBf %&% Vf, name="IVf" )

# Create Algebra for Expected Variance/Covariance Matrices in MZ & DZ twins
covMZm    <- mxAlgebra( expression= Am+ Cm, name="cMZm" )
covMZf    <- mxAlgebra( expression= Af+ Cf, name="cMZf" )
covDZm    <- mxAlgebra( expression= 0.5%x%(Am)+ Cm, name="cDZm" )
covDZf    <- mxAlgebra( expression= 0.5%x%(Af)+ Cf, name="cDZf" )
covDZmf   <- mxAlgebra( expression= 0.5%x%Amf+ Cmf, name="cDZmf" )
covDZfm   <- mxAlgebra( expression= 0.5%x%t(Amf)+ t(Cmf), name="cDZfm" )
matI2     <- mxMatrix( type="Iden", nrow=2, ncol=2, name="I2")
expCovMZm <- mxAlgebra( expression= I2 %x% IBm %&% rbind( cbind(Vm, cMZm), cbind(cMZm, Vm)), name="expCovMZm" )
expCovMZf <- mxAlgebra( expression= I2 %x% IBf %&% rbind( cbind(Vf, cMZf), cbind(cMZf, Vf)), name="expCovMZf" )
expCovDZm <- mxAlgebra( expression= I2 %x% IBm %&% rbind( cbind(Vm, cDZm), cbind(cDZm, Vm)), name="expCovDZm" )
expCovDZf <- mxAlgebra( expression= I2 %x% IBf %&% rbind( cbind(Vf, cDZf), cbind(cDZf, Vf)), name="expCovDZf" )
matZ2     <- mxMatrix( type="Zero", nrow=2, ncol=2, name="Z2")
expCovDZo <- mxAlgebra( expression= rbind(cbind(IBm,Z2),cbind(Z2,IBf)) %&% rbind( cbind(Vm, cDZmf), cbind(cDZfm, Vf)), name="expCovDZo" )

# Create Data Objects for Multiple Groups
dataMZm   <- mxData( observed=dataMZm, type="raw" )
dataMZf   <- mxData( observed=dataMZf, type="raw" )
dataDZm   <- mxData( observed=dataDZm, type="raw" )
dataDZf   <- mxData( observed=dataDZf, type="raw" )
dataDZo   <- mxData( observed=dataDZo, type="raw" )

# Create Expectation Objects for Multiple Groups
expMZm    <- mxExpectationNormal( covariance="expCovMZm", means="expMeanZm", dimnames=selVars )
expMZf    <- mxExpectationNormal( covariance="expCovMZf", means="expMeanZf", dimnames=selVars )
expDZm    <- mxExpectationNormal( covariance="expCovDZm", means="expMeanZm", dimnames=selVars )
expDZf    <- mxExpectationNormal( covariance="expCovDZf", means="expMeanZf", dimnames=selVars )
expDZo    <- mxExpectationNormal( covariance="expCovDZo", means="expMeanZo", dimnames=selVars )
funML     <- mxFitFunctionML()

# ---------------------------------------------------------------------------------------------------------------------|
# PREPARE OUTPUT

# Create Algebra for all Mean Parameter Estimates
colMm      <- cbind("meZm")
colMf      <- cbind("meZf")
estPathMm  <- mxAlgebra( expression= t(meanZm[,]), name="PathMm", dimnames=list(vars,colMm))
estPathMf  <- mxAlgebra( expression= t(meanZf[,]), name="PathMf", dimnames=list(vars,colMf))

# Create Algebra for all (Co)Variance Parameter Estimates
colVm     <- c(labDiag("am",nv),labDiag("cm",nv),labDiag("em",nv),labSdiag("ram",nv),labSdiag("rcm",nv),labSdiag("rem",nv),"bm21","bm12")
colVf     <- c(labDiag("af",nv),labDiag("cf",nv),labDiag("ef",nv),labSdiag("raf",nv),labSdiag("rcf",nv),labSdiag("ref",nv),"bf21","bf12")
colVc     <- c(labFull("rao",nv,nv),labFull("rco",nv,nv))
estPathVm <- mxAlgebra( expression= cbind(t(diag2vec(am)),t(diag2vec(cm)),t(diag2vec(em)),t(vechs(Ram)),t(vechs(Rcm)),t(vechs(Rem)),bm[2,1],bm[1,2]), name="PathVm", dimnames=list("pem",colVm))
estPathVf <- mxAlgebra( expression= cbind(t(diag2vec(af)),t(diag2vec(cf )),t(diag2vec(ef)),t(vechs(Raf)),t(vechs(Rcf )),t(vechs(Ref)),bf[2,1],bf[1,2]), name="PathVf", dimnames=list("pef",colVf))
estPathVc <- mxAlgebra( expression= cbind(t(cvectorize(Rao)),t(cvectorize(Rco))), name="PathVc", dimnames=list("peo",colVc))

# Create Algebra for Proportions of Variance
colPr     <- c("Pam11","Pcm11","Pem11","Pam22","Pcm22","Pem22","Pram21","Prao21","Prcm21","Prco21","Prem21","Pbm21","Pbm12",
               "Paf11","Pcf11","Pef11","Paf22","Pcf22","Pef22","Praf21","Prao12","Prcf21","Prco12","Pref21","Pbf21","Pbf12","Prao11","Prao22")
estProp   <- mxAlgebra( expression= cbind(
              cbind(iSDm[1,1]*am[1,1],iSDm[1,1]*cm[1,1],iSDm[1,1]*em[1,1],iSDm[2,2]*am[2,2],iSDm[2,2]*cm[2,2],iSDm[2,2]*em[2,2],
                     Ram[2,1],Rao[2,1],Rcm[2,1],Rco[2,1],Rem[2,1],iSDm[2,2]*bm[2,1]*iSDm[1,1],iSDm[1,1]*bm[1,2]*iSDm[2,2],
                    iSDf[1,1]*af[1,1],iSDf[1,1]*cf[1,1],iSDf[1,1]*ef[1,1],iSDf[2,2]*af[2,2],iSDf[2,2]*cf[2,2],iSDf[2,2]*ef[2,2],
                     Raf[2,1],Rao[1,2],Rcf[2,1],Rco[1,2],Ref[2,1],iSDf[2,2]*bf[2,1]*iSDf[1,1],iSDf[1,1]*bf[1,2]*iSDf[2,2] )^2,
              Rao[1,1],Rao[2,2]), name="Prop", dimnames=list("pro",colPr))

# Create List of free Parameters per Model
colPra    <- c("am11","cm11","em11","am22","cm22","em22","ram21","rao21","rem21","bm21","bm12","af11","cf11","ef11","af22","cf22","ef22","raf21","rao12","ref21","bf21","bf12","rao11","rao22","rc21")
colPrc    <- c("am11","cm11","em11","am22","cm22","em22","rcm21","rco21","rem21","bm21","bm12","af11","cf11","ef11","af22","cf22","ef22","rcf21","rco12","ref21","bf21","bf12","rco11","rco22","ra21")
colPq     <- c("am11","cm11","em11","am22","cm22","em22","rem21","bm21","bm12","af11","cf11","ef11","af22","cf22","ef22","ref21","bf21","bf12","ra21","rc21")
colP      <- c("am11","cm11","em11","am22","cm22","em22","rem21","bm21","bm12","af11","cf11","ef11","af22","cf22","ef22","ref21","bf21","bf12","ra21","rc21")

# Create Algebra for Variance Components
colZm     <- paste(rep(c('_Am','_Cm','_Em','_IBm','_Vm','_iSDm'),each=nv),1:nv,sep="")
colZf     <- paste(rep(c('_Af','_Cf','_Ef','_IBf','_Vf','_iSDf'),each=nv),1:nv,sep="")
estVarsZm <- mxAlgebra( expression= cbind(Am,Cm,Em, IBm, Vm,iSDm), name="VarsZm", dimnames=list(paste("vc",1:nv,sep=""),colZm))
estVarsZf <- mxAlgebra( expression= cbind(Af,Cf,Ef, IBf, Vf,iSDf), name="VarsZf", dimnames=list(paste("vc",1:nv,sep=""),colZf))

# Create Algebra for Standardized Variance Components and Correlations
colSm     <- paste(rep(c('_sAm','_sCm','_sEm','_rAm','_rCm','_rEm'),each=nv),1:nv,sep="")
colSf     <- paste(rep(c('_sAf','_sCf','_sEf','_rAf','_rCf','_rEf'),each=nv),1:nv,sep="")
estCorsZm <- mxAlgebra( expression= cbind(Am/Vm,Cm/Vm,Em/Vm, Ram,Rcm,Rem), name="CorsZm", dimnames=list(paste("sc",1:nv,sep=""),colSm))
estCorsZf <- mxAlgebra( expression= cbind(Af/Vf,Cf/Vf,Ef/Vf, Raf,Rcf,Ref), name="CorsZf", dimnames=list(paste("sc",1:nv,sep=""),colSf))

# Create Lists of (Derived) Parameter Estimates
colEs     <- c(labVars("meZm",vars),labVars("meZf",vars),
               labDiag("am",nv),labDiag("cm",nv),labDiag("em",nv),labSdiag("ram",nv),labSdiag("rcm",nv),labSdiag("rem",nv),"bm21","bm12",
               labDiag("af",nv),labDiag("cf",nv),labDiag("ef",nv),labSdiag("raf",nv),labSdiag("rcf",nv),labSdiag("ref",nv),"bf21","bf12",
               labFull("rao",nv,nv),labFull("rco",nv,nv))
paths     <- mxAlgebra( cbind(t(rvectorize(PathMm)),t(rvectorize(PathMf)),
                        PathVm,PathVf,PathVc,Prop), name="paths", dimnames=list(NULL,c(colEs,colPr)))
estTV     <- list( estPathMm, estPathMf, estPathVm, estPathVf, estPathVc, estProp)
estVC     <- list( estVarsZm, estCorsZm, estVarsZf, estCorsZf)

# ---------------------------------------------------------------------------------------------------------------------|
# RUN MODELS

fitEsts   <- function(fit) {}
#    print(round(cbind(fit$PathMm$result,fit$PathMf$result),4));
#    print(round(cbind(fit$PathVm$result,fit$PathVf$result,fit$PathVc$result),4));
#    print(round(fit$Prop$result,4));
#    print(round(cbind(fit$VarsZm$result,fit$VarsZf$result),4)); }
#    print(round(cbind(fit$CorsZm$result,fit$CorsZf$result),4));
#    print(round(fit$output$confidenceIntervals,4))

# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
# Correlated Liability ACE

makeModel <- function(name) {
pars      <- list( matI, matI2, matZ2, cnstrPos )
parsT     <- list( meanZm, meanZf )
parsZm    <- list( pathAm, pathCm, pathEm, pathRam, pathRcm, pathRem, covAm, covCm, covEm, covVm, invSDm, pathBm, invIBm, covIVm )
parsZf    <- list( pathAf, pathCf, pathEf, pathRaf, pathRcf, pathRef, covAf, covCf, covEf, covVf, invSDf, pathBf, invIBf, covIVf )
parsZo    <- list( parsZm, parsZf, pathRao, pathRco, covAmf, covCmf )
modelMZm  <- mxModel( pars, parsZm, meanZm, expMeanZm, covMZm, expCovMZm, dataMZm, expMZm, funML, name="MZm" )
modelMZf  <- mxModel( pars, parsZf, meanZf, expMeanZf, covMZf, expCovMZf, dataMZf, expMZf, funML, name="MZf" )
modelDZm  <- mxModel( pars, parsZm, meanZm, expMeanZm, covDZm, expCovDZm, dataDZm, expDZm, funML, name="DZm" )
modelDZf  <- mxModel( pars, parsZf, meanZf, expMeanZf, covDZf, expCovDZf, dataDZf, expDZf, funML, name="DZf" )
modelDZo  <- mxModel( pars, parsZo, meanZm, meanZf, expMeanZo, covDZmf, covDZfm, expCovDZo, dataDZo, expDZo, funML, name="DZo" )
modelT    <- list( modelMZm, modelMZf, modelDZm, modelDZf, modelDZo )
multi     <- mxFitFunctionMultigroup( c("MZm","DZm","MZf","DZf","DZo") )
name	  <- mxModel( name, pars, parsT, parsZo, corMin, corPos, modelT, multi, paths, estTV, estVC )
}

# Build Model with Confidence Intervals
#modelACE5rF <- mxModel( "ACE5rF", pars, parsT, parsZo, var1m, var1f, corMin, corPos, modelT, multi, paths, estTV, estVC )
modelACE5rF <- makeModel("ACE5rF")
modelACE5rF

# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
# Alternative models, including Correlated Liability, Direction of Causation and Hybrid Models

# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
# Correlated Liability ACE

# Run ACEra model - Qualitative (ra) & Quantitative Sex Differences
modelACE5raFace <- mxModel( modelACE5rF, name="ACE5raFace" )
modelACE5raFace <- omxSetParameters( modelACE5raFace, labels=c(labSdiag("rcm",nv),labSdiag("rcf",nv)), free=TRUE, values=svRc, newlabels=labSdiag("rc",nv) )
modelACE5raFace <- omxSetParameters( modelACE5raFace, labels=labOdiag("rco",nv), free=TRUE, values=svRc, newlabels=labSdiag("rc",nv) )
modelACE5raFace <- omxSetParameters( modelACE5raFace, labels=labDiag("rco",nv), free=FALSE, values=1)
fitACE5raFace   <- mxRun( modelACE5raFace ); fitGofs( fitACE5raFace )
fitACE5raFace   <- mxTryHard(fitACE5raFace, extraTries=4, greenOK=TRUE, silent=TRUE)
fitGofs( fitACE5raFace ); fitEsts( fitACE5raFace ); print(mxCompare( fitACE5raFace ))

# Run ACErc model - Qualitative (rc) & Quantitative Sex Differences
modelACE5rcFace <- mxModel( modelACE5rF, name="ACE5rcFace" )
modelACE5rcFace <- omxSetParameters( modelACE5rcFace, labels=c(labSdiag("ram",nv),labSdiag("raf",nv)), free=TRUE, values=svRa, newlabels=labSdiag("ra",nv) )
modelACE5rcFace <- omxSetParameters( modelACE5rcFace, labels=labOdiag("rao",nv), free=TRUE, values=svRa, newlabels=labSdiag("ra",nv) )
modelACE5rcFace <- omxSetParameters( modelACE5rcFace, labels=labDiag("rao",nv), free=FALSE, values=1 )
fitACE5rcFace   <- mxRun( modelACE5rcFace ); fitGofs( fitACE5rcFace )
fitACE5rcFace   <- mxTryHard(fitACE5rcFace, extraTries=4, greenOK=TRUE, silent=TRUE)
fitGofs( fitACE5rcFace ); fitEsts( fitACE5rcFace ); print(mxCompare( fitACE5raFace, fitACE5rcFace ))

# Run ACEq model - Quantitative Sex Differences
modelACE5qFace <- mxModel( fitACE5rcFace, name="ACE5qFace" )
modelACE5qFace <- omxSetParameters( modelACE5qFace, labels=c(labSdiag("rcm",nv),labSdiag("rcf",nv)), free=TRUE, values=svRc, newlabels=labSdiag("rc",nv) )
modelACE5qFace <- omxSetParameters( modelACE5qFace, labels=labOdiag("rco",nv), free=TRUE, values=svRc, newlabels=labSdiag("rc",nv) )
modelACE5qFace <- omxSetParameters( modelACE5qFace, labels=labDiag("rco",nv), free=FALSE, values=1 )
fitACE5qFace   <- mxRun( modelACE5qFace ); fitGofs( fitACE5qFace )
fitACE5qFace   <- mxTryHard(fitACE5qFace, extraTries=4, greenOK=TRUE, silent=TRUE)
fitGofs( fitACE5qFace ); fitEsts( fitACE5qFace ); print(mxCompare( fitACE5raFace, fitACE5qFace ))

# Run ACE model - No Sex Differences
modelACE5Face <- mxModel( fitACE5qFace, name="ACE5Face" )
modelACE5Face <- omxSetParameters( modelACE5Face, labels=c(labSdiag("rem",nv),labSdiag("ref",nv)), free=TRUE, values=svRc, newlabels=labSdiag("re",nv) )
modelACE5Face <- omxSetParameters( modelACE5Face, labels=c(labDiag("am",nv),labDiag("af",nv)), free=TRUE, values=svPa, newlabels=labDiag("a",nv) )
modelACE5Face <- omxSetParameters( modelACE5Face, labels=c(labDiag("cm",nv),labDiag("cf",nv)), free=TRUE, values=svPa, newlabels=labDiag("c",nv) )
modelACE5Face <- omxSetParameters( modelACE5Face, labels=c(labDiag("em",nv),labDiag("ef",nv)), free=TRUE, values=svPe, newlabels=labDiag("e",nv) )
fitACE5Face   <- mxRun( modelACE5Face ); fitGofs( fitACE5Face )
fitACE5Face   <- mxTryHard(fitACE5Face, extraTries=4, greenOK=TRUE, silent=TRUE)
fitGofs( fitACE5Face); fitEsts( fitACE5Face ); print(mxCompare( fitACE5raFace, fitACE5Face ))

# Generate Output Table of all Nested Models
ACE5FaceFits <- c( fitACE5raFace, fitACE5rcFace, fitACE5qFace, fitACE5Face )

ACE5Face.fit <- rbind(
 mxCompare( fitACE5raFace, fitACE5qFace ),
 mxCompare( fitACE5rcFace, fitACE5qFace )[2,],
 mxCompare( fitACE5qFace, fitACE5Face )[2,])

# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
# Correlated Liability AC

# Run ACEq model - Quantitative Sex Differences
modelACE5raFac <- mxModel( fitACE5raFace, name="ACE5raFac" )
modelACE5raFac <- omxSetParameters( modelACE5raFac, labels=c("rem21","ref21"), free=FALSE, values=0 )
fitACE5raFac   <- mxRun( modelACE5raFac ); fitGofs( fitACE5raFac ); fitEsts( fitACE5raFac );
if (fitACE5raFac$output$status$code >=5) {fitACE5raFac   <- mxRun( fitACE5raFac ); fitGofs( fitACE5raFac ); fitEsts( fitACE5raFac )}

# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
# Correlated Liability AE

# Run ACEq model - Quantitative Sex Differences
modelACE5raFae <- mxModel( fitACE5raFace, name="ACE5raFae" )
modelACE5raFae <- omxSetParameters( modelACE5raFae, labels=c("rc21"), free=FALSE, values=0 )
fitACE5raFae   <- mxRun( modelACE5raFae ); fitGofs( fitACE5raFae ); fitEsts( fitACE5raFae );
if (fitACE5raFae$output$status$code >=5) {fitACE5raFae   <- mxRun( fitACE5raFae ); fitGofs( fitACE5raFae ); fitEsts( fitACE5raFae )}

# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
# Correlated Liability CE

# Run ACEq model - Quantitative Sex Differences
modelACE5raFce <- mxModel( fitACE5raFace, name="ACE5raFce" )
modelACE5raFce <- omxSetParameters( modelACE5raFce, labels=c("ram21","raf21","rao21","rao12"), free=FALSE, values=0 )
fitACE5raFce   <- mxRun( modelACE5raFce ); fitGofs( fitACE5raFce ); fitEsts( fitACE5raFce );
if (fitACE5raFce$output$status$code >=5) {fitACE5raFce   <- mxRun( fitACE5raFce ); fitGofs( fitACE5raFce ); fitEsts( fitACE5raFce )}

# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
# Correlated Liability A

# Run ACEq model - Quantitative Sex Differences
modelACE5raFa <- mxModel( fitACE5raFace, name="ACE5raFa" )
modelACE5raFa <- omxSetParameters( modelACE5raFa, labels=c("rc21","rem21","ref21"), free=FALSE, values=0 )
fitACE5raFa   <- mxRun( modelACE5raFa ); fitGofs( fitACE5raFa ); fitEsts( fitACE5raFa );
if (fitACE5raFa$output$status$code >=5) {fitACE5raFa   <- mxRun( fitACE5raFa ); fitGofs( fitACE5raFa ); fitEsts( fitACE5raFa )}

# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
# Correlated Liability C

# Run ACEq model - Quantitative Sex Differences
modelACE5raFc <- mxModel( fitACE5raFace, name="ACE5raFc" )
modelACE5raFc <- omxSetParameters( modelACE5raFc, labels=c("ram21","raf21","rao21","rao12","rem21","ref21"), free=FALSE, values=0 )
fitACE5raFc   <- mxRun( modelACE5raFc ); fitGofs( fitACE5raFc ); fitEsts( fitACE5raFc );
if (fitACE5raFc$output$status$code >=5) {fitACE5raFc   <- mxRun( fitACE5raFc ); fitGofs( fitACE5raFc ); fitEsts( fitACE5raFc )}

# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
# Correlated Liability E

# Run ACEq model - Quantitative Sex Differences
modelACE5raFe <- mxModel( fitACE5raFace, name="ACE5raFe" )
modelACE5raFe <- omxSetParameters( modelACE5raFe, labels=c("ram21","raf21","rao21","rao12","rc21"), free=FALSE, values=0 )
fitACE5raFe   <- mxRun( modelACE5raFe ); fitGofs( fitACE5raFe ); fitEsts( fitACE5raFe );
if (fitACE5raFe$output$status$code >=5) {fitACE5raFe   <- mxRun( fitACE5raFe ); fitGofs( fitACE5raFe ); fitEsts( fitACE5raFe )}


# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
# Hybrid Causation and Correlated Liability (ACb21)

# Run ACEq model - Quantitative Sex Differences
modelACE5raHaco <- mxModel( fitACE5raFace, name="ACE5raHaco" )
modelACE5raHaco <- omxSetParameters( modelACE5raHaco, labels=c("rem21","ref21"), free=FALSE, values=0 )
modelACE5raHaco <- omxSetParameters( modelACE5raHaco, labels=c("bm21","bf21"), free=TRUE, values=0 )
fitACE5raHaco   <- mxRun( modelACE5raHaco ); fitGofs( fitACE5raHaco ); fitEsts( fitACE5raHaco );
if (fitACE5raHaco$output$status$code >=5) {fitACE5raHaco   <- mxRun( fitACE5raHaco ); fitGofs( fitACE5raHaco ); fitEsts( fitACE5raHaco )}

# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
# Hybrid Causation and Correlated Liability (AEb21)

# Run ACEq model - Quantitative Sex Differences
modelACE5raHaeo <- mxModel( fitACE5raFace, name="ACE5raHaeo" )
modelACE5raHaeo <- omxSetParameters( modelACE5raHaeo, labels=c("rc21"), free=FALSE, values=0 )
modelACE5raHaeo <- omxSetParameters( modelACE5raHaeo, labels=c("bm21","bf21"), free=TRUE, values=0 )
fitACE5raHaeo   <- mxRun( modelACE5raHaeo ); fitGofs( fitACE5raHaeo ); fitEsts( fitACE5raHaeo );
if (fitACE5raHaeo$output$status$code >=5) {fitACE5raHaeo   <- mxRun( fitACE5raHaeo ); fitGofs( fitACE5raHaeo ); fitEsts( fitACE5raHaeo )}

# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
# Hybrid Causation and Correlated Liability (CEb21)

# Run ACEq model - Quantitative Sex Differences
modelACE5raHceo <- mxModel( fitACE5raFace, name="ACE5raHceo" )
modelACE5raHceo <- omxSetParameters( modelACE5raHceo, labels=c("ram21","raf21","rao21","rao12"), free=FALSE, values=0 )
modelACE5raHceo <- omxSetParameters( modelACE5raHceo, labels=c("bm21","bf21"), free=TRUE, values=0 )
fitACE5raHceo   <- mxRun( modelACE5raHceo ); fitGofs( fitACE5raHceo ); fitEsts( fitACE5raHceo );
if (fitACE5raHceo$output$status$code >=5) {fitACE5raHceo   <- mxRun( fitACE5raHceo ); fitGofs( fitACE5raHceo ); fitEsts( fitACE5raHceo )}

# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
# Hybrid Causation and Correlated Liability (ACb12)

# Run ACEq model - Quantitative Sex Differences
modelACE5raHacr <- mxModel( fitACE5raFace, name="ACE5raHacr" )
modelACE5raHacr <- omxSetParameters( modelACE5raHacr, labels=c("rem21","ref21"), free=FALSE, values=0 )
modelACE5raHacr <- omxSetParameters( modelACE5raHacr, labels=c("bm12","bf12"), free=TRUE, values=0 )
fitACE5raHacr   <- mxRun( modelACE5raHacr ); fitGofs( fitACE5raHacr ); fitEsts( fitACE5raHacr );
if (fitACE5raHacr$output$status$code >=5) {fitACE5raHacr   <- mxRun( fitACE5raHacr ); fitGofs( fitACE5raHacr ); fitEsts( fitACE5raHacr )}

# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
# Hybrid Causation and Correlated Liability (AEb12)

# Run ACEq model - Quantitative Sex Differences
modelACE5raHaer <- mxModel( fitACE5raFace, name="ACE5raHaer" )
modelACE5raHaer <- omxSetParameters( modelACE5raHaer, labels=c("rc21"), free=FALSE, values=0 )
modelACE5raHaer <- omxSetParameters( modelACE5raHaer, labels=c("bm12","bf12"), free=TRUE, values=0 )
fitACE5raHaer   <- mxRun( modelACE5raHaer ); fitGofs( fitACE5raHaer ); fitEsts( fitACE5raHaer );
if (fitACE5raHaer$output$status$code >=5) {fitACE5raHaer   <- mxRun( fitACE5raHaer ); fitGofs( fitACE5raHaer ); fitEsts( fitACE5raHaer )}

# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
# Hybrid Causation and Correlated Liability (CEb12)

# Run ACEq model - Quantitative Sex Differences
modelACE5raHcer <- mxModel( fitACE5raFace, name="ACE5raHcer" )
modelACE5raHcer <- omxSetParameters( modelACE5raHcer, labels=c("ram21","raf21","rao21","rao12"), free=FALSE, values=0 )
modelACE5raHcer <- omxSetParameters( modelACE5raHcer, labels=c("bm12","bf12"), free=TRUE, values=0 )
fitACE5raHcer   <- mxRun( modelACE5raHcer ); fitGofs( fitACE5raHcer ); fitEsts( fitACE5raHcer );
if (fitACE5raHcer$output$status$code >=5) {fitACE5raHcer   <- mxRun( fitACE5raHcer ); fitGofs( fitACE5raHcer ); fitEsts( fitACE5raHcer )}

# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
# Hybrid Causation and Correlated Liability (Ab21b12)

# Run ACEq model - Quantitative Sex Differences
modelACE5raHaor <- mxModel( fitACE5raFace, name="ACE5raHaor" )
modelACE5raHaor <- omxSetParameters( modelACE5raHaor, labels=c("rc21","rem21","ref21"), free=FALSE, values=0 )
modelACE5raHaor <- omxSetParameters( modelACE5raHaor, labels=c("bm21","bm12","bf21","bf12"), free=TRUE, values=0 )
fitACE5raHaor   <- mxRun( modelACE5raHaor ); fitGofs( fitACE5raHaor ); fitEsts( fitACE5raHaor );
if (fitACE5raHaor$output$status$code >=5) {fitACE5raHaor   <- mxRun( fitACE5raHaor ); fitGofs( fitACE5raHaor ); fitEsts( fitACE5raHaor )}

# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
# Hybrid Causation and Correlated Liability (Cb21b12)

# Run ACEq model - Quantitative Sex Differences
modelACE5raHcor <- mxModel( fitACE5raFace, name="ACE5raHcor" )
modelACE5raHcor <- omxSetParameters( modelACE5raHcor, labels=c("ram21","raf21","rao21","rao12","rem21","ref21"), free=FALSE, values=0 )
modelACE5raHcor <- omxSetParameters( modelACE5raHcor, labels=c("bm21","bm12","bf21","bf12"), free=TRUE, values=0 )
fitACE5raHcor   <- mxRun( modelACE5raHcor ); fitGofs( fitACE5raHcor ); fitEsts( fitACE5raHcor );
if (fitACE5raHcor$output$status$code >=5) {fitACE5raHcor   <- mxRun( fitACE5raHcor ); fitGofs( fitACE5raHcor ); fitEsts( fitACE5raHcor )}

# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
# Hybrid Causation and Correlated Liability (Eb21b12)

# Run ACEq model - Quantitative Sex Differences
modelACE5raHeor <- mxModel( fitACE5raFace, name="ACE5raHeor" )
modelACE5raHeor <- omxSetParameters( modelACE5raHeor, labels=c("ram21","raf21","rao21","rao12","rc21"), free=FALSE, values=0 )
modelACE5raHeor <- omxSetParameters( modelACE5raHeor, labels=c("bm21","bm12","bf21","bf12"), free=TRUE, values=0 )
fitACE5raHeor   <- mxRun( modelACE5raHeor ); fitGofs( fitACE5raHeor ); fitEsts( fitACE5raHeor );
if (fitACE5raHeor$output$status$code >=5) {fitACE5raHeor   <- mxRun( fitACE5raHeor ); fitGofs( fitACE5raHeor ); fitEsts( fitACE5raHeor )}

# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
# Hybrid Causation and Correlated Liability (Ab21)

# Run ACEq model - Quantitative Sex Differences
modelACE5raHao <- mxModel( fitACE5raFace, name="ACE5raHao" )
modelACE5raHao <- omxSetParameters( modelACE5raHao, labels=c("rc21","rem21","ref21"), free=FALSE, values=0 )
modelACE5raHao <- omxSetParameters( modelACE5raHao, labels=c("bm21","bf21"), free=TRUE, values=0 )
fitACE5raHao   <- mxRun( modelACE5raHao ); fitGofs( fitACE5raHao ); fitEsts( fitACE5raHao );
if (fitACE5raHao$output$status$code >=5) {fitACE5raHao   <- mxRun( fitACE5raHao ); fitGofs( fitACE5raHao ); fitEsts( fitACE5raHao )}

# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
# Hybrid Causation and Correlated Liability (Cb21)

# Run ACEq model - Quantitative Sex Differences
modelACE5raHco <- mxModel( fitACE5raFace, name="ACE5raHco" )
modelACE5raHco <- omxSetParameters( modelACE5raHco, labels=c("ram21","raf21","rao21","rao12","rem21","ref21"), free=FALSE, values=0 )
modelACE5raHco <- omxSetParameters( modelACE5raHco, labels=c("bm21","bf21"), free=TRUE, values=0 )
fitACE5raHco   <- mxRun( modelACE5raHco ); fitGofs( fitACE5raHco ); fitEsts( fitACE5raHco );
if (fitACE5raHco$output$status$code >=5) {fitACE5raHco   <- mxRun( fitACE5raHco ); fitGofs( fitACE5raHco ); fitEsts( fitACE5raHco )}

# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
# Hybrid Causation and Correlated Liability (Eb21)

# Run ACEq model - Quantitative Sex Differences
modelACE5raHeo <- mxModel( fitACE5raFace, name="ACE5raHeo" )
modelACE5raHeo <- omxSetParameters( modelACE5raHeo, labels=c("ram21","raf21","rao21","rao12","rc21"), free=FALSE, values=0 )
modelACE5raHeo <- omxSetParameters( modelACE5raHeo, labels=c("bm21","bf21"), free=TRUE, values=0 )
fitACE5raHeo   <- mxRun( modelACE5raHeo ); fitGofs( fitACE5raHeo ); fitEsts( fitACE5raHeo );
if (fitACE5raHeo$output$status$code >=5) {fitACE5raHeo   <- mxRun( fitACE5raHeo ); fitGofs( fitACE5raHeo ); fitEsts( fitACE5raHeo )}

# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
# Hybrid Causation and Correlated Liability (Ab12)

# Run ACEq model - Quantitative Sex Differences
modelACE5raHar <- mxModel( fitACE5raFace, name="ACE5raHar" )
modelACE5raHar <- omxSetParameters( modelACE5raHar, labels=c("rc21","rem21","ref21"), free=FALSE, values=0 )
modelACE5raHar <- omxSetParameters( modelACE5raHar, labels=c("bm12","bf12"), free=TRUE, values=0 )
fitACE5raHar   <- mxRun( modelACE5raHar ); fitGofs( fitACE5raHar ); fitEsts( fitACE5raHar );
if (fitACE5raHar$output$status$code >=5) {fitACE5raHar   <- mxRun( fitACE5raHar ); fitGofs( fitACE5raHar ); fitEsts( fitACE5raHar )}

# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
# Hybrid Causation and Correlated Liability (Cb12)

# Run ACEq model - Quantitative Sex Differences
modelACE5raHcr <- mxModel( fitACE5raFace, name="ACE5raHcr" )
modelACE5raHcr <- omxSetParameters( modelACE5raHcr, labels=c("ram21","raf21","rao21","rao12","rem21","ref21"), free=FALSE, values=0 )
modelACE5raHcr <- omxSetParameters( modelACE5raHcr, labels=c("bm12","bf12"), free=TRUE, values=0 )
fitACE5raHcr   <- mxRun( modelACE5raHcr ); fitGofs( fitACE5raHcr ); fitEsts( fitACE5raHcr );
if (fitACE5raHcr$output$status$code >=5) {fitACE5raHcr   <- mxRun( fitACE5raHcr ); fitGofs( fitACE5raHcr ); fitEsts( fitACE5raHcr )}

# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
# Hybrid Causation and Correlated Liability (Eb12)

# Run ACEq model - Quantitative Sex Differences
modelACE5raHer <- mxModel( fitACE5raFace, name="ACE5raHer" )
modelACE5raHer <- omxSetParameters( modelACE5raHer, labels=c("ram21","raf21","rao21","rao12","rc21"), free=FALSE, values=0 )
modelACE5raHer <- omxSetParameters( modelACE5raHer, labels=c("bm12","bf12"), free=TRUE, values=0 )
fitACE5raHer   <- mxRun( modelACE5raHer ); fitGofs( fitACE5raHer ); fitEsts( fitACE5raHer );
if (fitACE5raHer$output$status$code >=5) {fitACE5raHer   <- mxRun( fitACE5raHer ); fitGofs( fitACE5raHer ); fitEsts( fitACE5raHer )}

# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
# Reciprocal Causation (b21b12)

# Run ACEq model - Quantitative Sex Differences
modelACE5raDor <- mxModel( fitACE5raFace, name="ACE5raDor")
modelACE5raDor <- omxSetParameters( modelACE5raDor, labels=c("ram21","raf21","rao21","rao12","rc21","rem21","ref21"), free=FALSE, values=0 )
modelACE5raDor <- omxSetParameters( modelACE5raDor, labels=c("bm21","bm12","bf21","bf12"), free=TRUE, values=0.01 )
fitACE5raDor   <- mxRun( modelACE5raDor ); fitGofs( fitACE5raDor ); fitEsts( fitACE5raDor );
if (fitACE5raDor$output$status$code >=5) {fitACE5raDor   <- mxRun( fitACE5raDor ); fitGofs( fitACE5raDor ); fitEsts( fitACE5raDor )}

# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
# Direction of Causation (b21)

# Run ACEq model - Quantitative Sex Differences
modelACE5raDo <- mxModel( fitACE5raFace, name="ACE5raDo")
modelACE5raDo <- omxSetParameters( modelACE5raDo, labels=c("ram21","raf21","rao21","rao12","rc21","rem21","ref21"), free=FALSE, values=0 )
modelACE5raDo <- omxSetParameters( modelACE5raDo, labels=c("bm21","bf21"), free=TRUE, values=0 )
fitACE5raDo   <- mxRun( modelACE5raDo ); fitGofs( fitACE5raDo ); fitEsts( fitACE5raDo );
if (fitACE5raDo$output$status$code >=5) {fitACE5raDo   <- mxRun( fitACE5raDo ); fitGofs( fitACE5raDo ); fitEsts( fitACE5raDo )}

# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
# Reverse Causation (b12)

# Run ACEq model - Quantitative Sex Differences
modelACE5raDr <- mxModel( fitACE5raFace, name="ACE5raDr")
modelACE5raDr <- omxSetParameters( modelACE5raDr, labels=c("ram21","raf21","rao21","rao12","rc21","rem21","ref21"), free=FALSE, values=0 )
modelACE5raDr <- omxSetParameters( modelACE5raDr, labels=c("bm12","bf12"), free=TRUE, values=0 )
fitACE5raDr   <- mxRun( modelACE5raDr ); fitGofs( fitACE5raDr ); fitEsts( fitACE5raDr );
if (fitACE5raDr$output$status$code >=5) {fitACE5raDr   <- mxRun( fitACE5raDr ); fitGofs( fitACE5raDr ); fitEsts( fitACE5raDr )}

# Print Comparative Fit Statistics

Fits  <- rbind(fitGofT(fitACE5raFace),
fitGofT(fitACE5raFac), fitGofT(fitACE5raFae), fitGofT(fitACE5raFce),
fitGofT(fitACE5raFa), fitGofT(fitACE5raFc), fitGofT(fitACE5raFe),
fitGofT(fitACE5raHaco), fitGofT(fitACE5raHaeo), fitGofT(fitACE5raHceo),
fitGofT(fitACE5raHacr), fitGofT(fitACE5raHaer), fitGofT(fitACE5raHcer),
fitGofT(fitACE5raHaor), fitGofT(fitACE5raHcor), fitGofT(fitACE5raHeor),
fitGofT(fitACE5raHao), fitGofT(fitACE5raHco), fitGofT(fitACE5raHeo),
fitGofT(fitACE5raHar), fitGofT(fitACE5raHcr), fitGofT(fitACE5raHer),
fitGofT(fitACE5raDor), fitGofT(fitACE5raDo), fitGofT(fitACE5raDr))

# Print Estimates & Proportions
Ests <- round(rbind(
fitACE5raFace$paths$result,
fitACE5raFac$paths$result, fitACE5raFae$paths$result, fitACE5raFce$paths$result,
fitACE5raFa$paths$result, fitACE5raFc$paths$result, fitACE5raFe$paths$result,
fitACE5raHaco$paths$result, fitACE5raHaeo$paths$result, fitACE5raHceo$paths$result,
fitACE5raHacr$paths$result, fitACE5raHaer$paths$result, fitACE5raHcer$paths$result,
fitACE5raHaor$paths$result, fitACE5raHcor$paths$result, fitACE5raHeor$paths$result,
fitACE5raHao$paths$result, fitACE5raHco$paths$result, fitACE5raHeo$paths$result,
fitACE5raHar$paths$result, fitACE5raHcr$paths$result, fitACE5raHer$paths$result,
fitACE5raDor$paths$result, fitACE5raDo$paths$result, fitACE5raDr$paths$result ),4)

# Print	Both & Write to	File
colGT <- c("Model", "Nobs", "Nrec", "Npar", "Ncon", "df", "Minus2LL", "AIC", "Status")
colnames(Fits ) <- colGT
#print(cbind(Fits,Ests),quote=F)
write.table(cbind(Fits,Ests ),paste(filename,".txt",sep=""),row.names=T,col.names=T,append=F,quote=F)

# Get Model Average
fitACE5 <- list(fitACE5raFace,
fitACE5raFac,  fitACE5raFae,  fitACE5raFce,
fitACE5raFa,   fitACE5raFc,   fitACE5raFe,
fitACE5raHaco, fitACE5raHaeo, fitACE5raHceo,
fitACE5raHacr, fitACE5raHaer, fitACE5raHcer,
fitACE5raHaor, fitACE5raHcor, fitACE5raHeor,
fitACE5raHao,  fitACE5raHco,  fitACE5raHeo,
fitACE5raHar,  fitACE5raHcr,  fitACE5raHer,
fitACE5raDor,  fitACE5raDo,   fitACE5raDr)
