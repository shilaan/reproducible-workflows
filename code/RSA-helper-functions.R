## This code is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
## 
## For a copy of the GNU General Public License, 
## see <http://www.gnu.org/licenses/>.

## (c) 2018 Sarah Humberg, Felix Schönbrodt. 

## This source code accompanies the following paper:
## Humberg, S., Dufner, M., Schönbrodt, F. D., Geukes, K., Hutteman, R., Küfner, A. C. P., van Zalk, M. H. W., Denissen, J. J. A., Nestler, S., & Back, M. D. (2018). Is accurate, positive, or inflated self-perception most advantageous for psychological adjustment? A competitive test of key hypotheses. Retrieved from osf.io/9w3bh


##############################
### USEFUL FUNCTIONS
##############################

# function to delete leading zero in decimal number
# easysub <- function(x){
#           sub("0.", ".", x)
# }

# function to create nice looking correlation matrix of data frame
# corcons <- function(data){
#     
#           #Compute correlation matrix
#           require(Hmisc)
#           data <- as.matrix(data)
#           correlation_matrix <- rcorr(data, type="pearson")
#           R <- correlation_matrix$r # Matrix of correlation coefficients
#           
#           # trunctuate the correlation matrix to two decimals and format nicely
#           R <- format(round(cbind(rep(-1.11, ncol(data)), R), 2))[,-1]
# 
#           # define row names and column names
#           rownames(R) <- colnames(data)
#           colnames(R) <- colnames(data) # paste0(colnames(x), "")
#           
#           # remove upper triangle of correlation matrix
#           R <- as.matrix(R)
#           R[upper.tri(R, diag = T)] <- ""
#           R <- R[-1,]
#           R <- R[,-ncol(R)]
#                    
#           R <- as.data.frame(R)
#           R <- apply(R, 1:2, easysub)
#           
#           # print matrix
#           R
# }


# define function irm ("find redundant models") that identifies models with essentially the same log-likelihood as a simpler model
# the input AIC/AICc table must be sorted by AIC (smalles AIC in first row)
frm <- function(aictable){
          
          # define list with all sets of models which are nested in each other, respectively
          nestinglist <- list(c("full", "SSQDposCneg", "SQDpos", "null"),
                              c("full", "IApos", "discrpos", "onlyxpos", "null"),
                              c("full", "IApos", "xandypos", "onlyxpos", "null"),
                              c("full", "IApos", "xandypos", "onlyypos", "null"),
                              c("full", "IApos", "discrneg", "onlyypos", "null"),
                              c("full", "IApos", "discrneg", "onlyxneg", "null"),
                              c("full", "IApos", "xandyneg", "onlyxneg", "null"),
                              c("full", "onlyx2pos", "onlyxpos", "null"),
                              c("full", "onlyx2pos", "onlyxneg", "null"),
                              c("full", "onlyy2pos", "onlyypos", "null"))
          
          # create empty vector to be filled with names of models we might want to remove
          toremove <- c()
          
          # go through all rows of the table, starting with the second row
          for (i in 2:nrow(aictable)){
                    # compare this row i to each of the previous rows k
                    for (k in 1:(i-1)){
                              # Are the log-likelihoods of the two models essentially the same?
                              similar_LL <- abs(aictable[i,"LL"] - aictable[k,"LL"]) < 1
                              
                              # Is model k nested in model i? --> test whether there is a vector in the nesting-list which contains both of the two models i and k
                              modelnames <- c(as.character(aictable[i,"Modnames"]), as.character(aictable[k,"Modnames"]))
                              matches <- c()
                              
                              for (j in 1:length(nestinglist)){
                                        matches = c(matches,sum(modelnames %in% nestinglist[[j]]))
                              }
                              
                              nested <- any(matches==2)          
                              
                              # If both conditions are true (similar LL and nested), add the name of model i to the list of models we might want to remove
                              if(similar_LL & nested){
                                        toremove <- c(toremove, as.character(aictable[i,"Modnames"]))
                              }
                    }
          }
          
          toremove <- as.vector(unique(toremove))
          
          # Print warning if models were detected that we might want to remove
          # if(length(toremove) > 0){
          #           warning(paste0("There were nested models with log-likelihood differences < 1. You might want to remove the following models from the set: ", "c(", paste0(toremove, collapse=", "), ")"))
          # }
          
          res <- list(toremove = toremove)
          
}


##############################
### PREPARE MAIN ANALYSES
##############################

# define function eam ("estimate all models") that estimates all models in our initial model set and prepares respective output table

# options in this function that do not explain themselves:
# variables Vector of variables the analysis is based on, must be ordered as follows: c("dependent variable","first predictor","second predictor")
# modelset Vector of model names that should be estimated and compared. Per default, only the full model is estimated.
# nr Optional string, for example number of analysis. Will be pasted into the file name of the output tables.
# controlvariables Optional.

eam <- function(variables, 
                data, 
                controlvariables="", 
                startvaluesfull="", 
                modelset = c("full"), 
                nr = "", 
                se="robust", 
                estimator="MLR", 
                missing="fiml", 
                outliers="liberal"){
          
          # rename dataframe, for consistency with RSA() function
          df <- data
          
          # set all result objects to NULL as default
          sem_xandypos <- NULL
          sem_onlyxpos <- NULL
          sem_xandyneg <- NULL
          sem_onlyxneg <- NULL
          sem_discrpos <- NULL
          sem_discrneg <- NULL
          sem_SQDpos <- NULL
          sem_SSQDposCneg <- NULL
          sem_onlyypos <- NULL
          sem_onlyx2pos <- NULL
          sem_onlyy2pos <- NULL
          sem_IApos <- NULL
          sem_null <- NULL
          sem_full <- NULL
          
          # extract variable names
          DV <- variables[1]
          IV1 <- variables[2]
          IV2 <- variables[3]
          IV12 <- paste0(IV1, "2")
          IV_IA <- paste0(IV1, "_", IV2)
          IV22 <- paste0(IV2, "2")
          
          # add higher order variables
          df[, IV12] <- df[, IV1]^2
          df[, IV22] <- df[, IV2]^2
          df[, IV_IA] <- df[, IV1]*df[, IV2]
          
          # run polynomial regresion as a OLS linear model
          f <- paste0(DV, " ~ ", paste("1", IV1, IV2, IV12, IV_IA, IV22, sep=" + "), controlvariables) 
          lm.full <- lm(f, data=df, na.action=na.exclude)
          
          # test full model for significance
          LM <- summary(lm.full)
          r2_full <- LM$adj.r.squared
          F <- LM$fstatistic
          p_full <- 1-pf(F[1], F[2], F[3])
          
          
          ## print warning if full model does not explain a significant amount of variance and re-define model set:
          if (p_full > 0.10){
                    warning(paste0("Do not interpret results of model comparisons, the full model is not significant: Adjusted R2 = ", round(r2_full,4), " with p = ", round(p_full,4)))
                    
                    ## build output object
                    res <- list(lm_full = lm.full,
                                r2_full = r2_full,
                                p_full = p_full
                    )         
                    
          }
          
          
          ## proceed if full model does explain significant amount of variance:
          if (p_full < 0.10){
                     
                    # outlier treatment: identify and save positions of outliers
                    df$out <- FALSE

                    # liberal outlier treatment
                    if (outliers == "liberal"){
                              inf <- influence.measures(lm.full)
                              df$out <- apply(inf$is.inf[, c("dffit", "cook.d", "hat")], 1, sum) == 3
                    }
                    
                    # more conservative outlier treatment (Cohen, Cohen, West, & Aiken, 2003)
                    if (outliers == "conservative"){
                              inf <- data.frame(influence.measures(lm.full)$infmat)
                              
                              # number of predictors and observations
                              k <- length(coef(lm.full))
                              n <- nobs(lm.full)
                              
                              # flag outliers
                              dffit.flag <- abs(inf$dffit) > 2*sqrt((k+1)/n)        
                              dfbetas <- inf[, grepl("dfb.", colnames(inf), fixed=TRUE)]
                              dfbetas.flag <- abs(dfbetas) > 2/sqrt(n)
                              
                              flags <- cbind(dffit.flag, dfbetas.flag)
                              
                              df$out <- apply(flags, 1, any)
                    }
                    

                    n.out <- sum(na.omit(df$out) == TRUE)
                    if (n.out>0) {
                              warning(paste("Removed", n.out, "multivariate outlier(s). Outliers are in row(s):", paste(which(df$out == TRUE) , collapse=", ")))
                    }
                    
                    # keep cases with missing value on df$out (missings will be replaced with missing option)
                    df$out[is.na(df$out)] <- FALSE
                    
                    
                    ## estimate all models that are contained in the vector modelset
                    
                    # define full polynomial model
                    poly <- paste(paste0(DV, " ~ 1 + b1*", IV1, " + b2*", IV2, " + b3*", IV12, " + b4*", IV_IA, " + b5*", IV22, controlvariables), 
                                  sep="\n")
								  
		# define RSA parameters to be added to each model
                    a1_to_a4 <- paste("a1 := b1+b2",
                                      "a2 := b3+b4+b5",
                                      "a3 := b1-b2",
                                      "a4 := b3-b4+b5",
                                      sep="\n")
                    
                    # positive main (PSV) effect model
                    if ("xandypos" %in% modelset){
                              xandypos <- paste(poly,
                                                    "b1>0",
                                                    "b2>0", # implementation of inequality constraint b2>=0 (see Rosseel, 2012)
                                                    "b3==0",
                                                    "b4==0",
                                                    "b5==0",
                                                    a1_to_a4,
                                                    sep="\n")
                              sem_xandypos <- sem(xandypos, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing)
                    }
                    
                    # only positive effect of S model
                    if ("onlyxpos" %in% modelset){
                              onlyxpos <- paste(poly,
                                                "b1>0",
                                                "b2==0", 
                                                "b3==0",
                                                "b4==0",
                                                "b5==0",
                                                a1_to_a4,
                                                sep="\n")
                              sem_onlyxpos <- sem(onlyxpos, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing)
                    }          
                    
                    # negative main (PSV) effect model
                    if ("xandyneg" %in% modelset){
                              xandyneg <- paste(poly,
                                                    "b1<0",
                                                    "b2<0", # implementation of inequality constraint b2>=0 (see Rosseel, 2012)
                                                    "b3==0",
                                                    "b4==0",
                                                    "b5==0",
                                                    a1_to_a4,
                                                    sep="\n")
                              sem_xandyneg <- sem(xandyneg, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing)
                    }
                    
                    # only negative effect of S model
                    if ("onlyxneg" %in% modelset){
                              onlyxneg <- paste(poly,
                                                "b1<0",
                                                "b2==0", 
                                                "b3==0",
                                                "b4==0",
                                                "b5==0",
                                                a1_to_a4,
                                                sep="\n")
                              sem_onlyxneg <- sem(onlyxneg, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing)
                    }          
                    
                    # positive discrepancy (SE) effect model
                    if ("discrpos" %in% modelset){
                              discrpos <- paste(poly,
                                                "b1>0",
                                                "b2<0", 
                                                "b3==0",
                                                "b4==0",
                                                "b5==0",
                                                a1_to_a4,
                                                # compute strength of the SE effect following Humberg et al. (in press)
                                                "abs := abs(a3) - abs(a1)",
                                                sep="\n")
                              sem_discrpos <- sem(discrpos, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing)
                    }
                    
                    # negative discrepancy (SE) effect model
                    if ("discrneg" %in% modelset){
                              discrneg <- paste(poly,
                                                "b1<0",
                                                "b2>0", 
                                                "b3==0",
                                                "b4==0",
                                                "b5==0",
                                                a1_to_a4,
                                                # compute strength of the SE effect following Humberg et al. (2018b)
                                                "abs := abs(a3) - abs(a1)",
                                                sep="\n")
                              sem_discrneg <- sem(discrneg, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing)
                    }
                    
                    # positive agreement effect model
                    if ("SQDpos" %in% modelset){
                              SQDpos <- paste(poly,
                                              "b1==0",
                                              "b2==0",
                                              "b3<0",
                                              "b3==b5",
                                              "b3+b4+b5==0",
                                              sep="\n")
                              sem_SQDpos <- sem(SQDpos, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing)
                    }  
                    
                    # positive optimal margin effect model with positive optimal amount of S-R
                    if ("SSQDposCneg" %in% modelset){
                              SSQDposCneg <- paste(poly,
                                                   "b1==-b2",
                                                   "b3<0",
                                                   "b3==b5",
                                                   "b3+b4+b5==0",
                                                   #
                                                   # the surface is maximal for S + C = R, with C as defined below (notation as in Schönbrodt, 2016)
                                                   # to constrain the "optimal amount of S-R" to be positive (surface highest for positve discrepancy S-R), 
                                                   # we need to constrain C to be negative. 
                                                   # Because b3<0, the following inequality is equivalent to C := (b1-b2)/(4*b3) < 0.
                                                   "b1-b2 > 0", 
                                                   #
                                                   # compute C, such that surface is highest for S + C = R:
                                                   "C := (b1-b2) / (4*b3)",
                                                   sep="\n")
                              sem_SSQDposCneg <- sem(SSQDposCneg, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing)
                    }
                    
                    # positive R effect model
                    if ("onlyypos" %in% modelset){
                              onlyypos <- paste(poly,
                                                "b1==0",
                                                "b2>0", 
                                                "b3==0",
                                                "b4==0",
                                                "b5==0",
                                                a1_to_a4,
                                                sep="\n")
                              sem_onlyypos <- sem(onlyypos, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing)
                    }
                    
                    # curvilinear effect of S model
                    if ("onlyx2pos" %in% modelset){
                              onlyx2pos <- paste(poly,
                                                 "b2==0",
                                                 "b3<0",
                                                 "b4==0",
                                                 "b5==0",
                                                 a1_to_a4,
                                                 # compute S value of vertex of the parabola:
                                                 "V := -b1/(2*b3)",
                                                 sep="\n")
                              sem_onlyx2pos <- sem(onlyx2pos, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing)
                    }
                    
                    # curvilinear effect of R model
                    if ("onlyy2pos" %in% modelset){
                              onlyy2pos <- paste(poly,
                                                 "b1==0",
                                                 "b3==0",
                                                 "b4==0",
                                                 "b5<0",
                                                 a1_to_a4,
                                                 # compute R value of vertex of the parabola:
                                                 "V := -b2/(2*b5)",
                                                 sep="\n")
                              sem_onlyy2pos <- sem(onlyy2pos, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing)
                    }
                    
                    # positive interaction term model
                    if ("IApos" %in% modelset){
                              IApos <- paste(poly,
                                             "b3==0",
                                             "b5==0",
                                             "b4>0",
                                             a1_to_a4,
                                             sep="\n")
                              sem_IApos <- sem(IApos, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing)
                    }
                    
                    # null model
                    if ("null" %in% modelset){
                              null <- paste(poly,
                                            "b1==0",
                                            "b2==0",
                                            "b3==0",
                                            "b4==0",
                                            "b5==0",
                                            a1_to_a4,
                                            sep="\n")
                              sem_null <- sem(null, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing)
                    }
                    
                    # full polynomial model of second degree
                    if ("full" %in% modelset){
                              # full <- paste(poly,
                              #               a1_to_a4,
                              #               sep="\n")
                              full <- paste(paste0(DV, " ~ 1 + b1*", IV1, " + b2*", IV2, " + b3*", IV12, " + b4*", IV_IA, " + b5*", IV22, controlvariables, startvaluesfull), 
                                            a1_to_a4,
                                            sep="\n")
                              sem_full <- sem(full, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing)
                    }
                    
                    
                    ## Build result objects of interest
                    
                    # list of all estimated models
                    maximal_modellist <- list(xandypos=sem_xandypos, 
                                              onlyxpos = sem_onlyxpos,
                                              xandyneg=sem_xandyneg, 
                                              onlyxneg = sem_onlyxneg,
                                              discrpos=sem_discrpos,
                                              discrneg=sem_discrneg,
                                              SQDpos=sem_SQDpos,
                                              SSQDposCneg=sem_SSQDposCneg,
                                              onlyypos=sem_onlyypos,
                                              onlyx2pos = sem_onlyx2pos,
                                              onlyy2pos = sem_onlyy2pos,
                                              IApos = sem_IApos,
                                              null = sem_null,
                                              full = sem_full
                    )
                    
                    modellist <- maximal_modellist[!sapply(maximal_modellist,is.null)]
                    
					
					
                    ## Prepare output table including coefficient estimates of all models in modelset, as well as further parameters of interest
                    coeftab <- data.frame(names(modellist))
                    names(coeftab)[1] <- "Modnames"
                    coeftab[,paste0("b",1:5)] <- NA
                    for (model in as.vector(coeftab$Modnames)){
                              coeftab[coeftab$Modnames==model,paste0("b",1:5)] <- round(t(dplyr::select(dplyr::filter(parameterEstimates(modellist[[model]]),
							label=="b1"|label=="b2"|label=="b3"|label=="b4"|label=="b5"
                              ),est)),3)
                    }
                    
                    
                    # extend table by further parameters of interest for specific models
                    
                    # R^2 of the full model
                    coeftab[,"R2"] <- ""
                    coeftab[coeftab$Modnames=="full","R2"] <- paste0(round(r2_full, 4), " (p=", round(p_full,4),")")
                    
                    # strength of the SE effect for the discrepancy effect models
                    discrepancy_models <- c("discrpos","discrneg")
                    if (any(discrepancy_models %in% names(modellist))){ 
                              
                              # identify respective models in fitted models
                              which <- which(discrepancy_models %in% names(modellist))
                              models <- discrepancy_models[which]
                              
                              # write their shift values into coeftab
                              coeftab[,"abs"] <- ""
                              
                              for (model in models){
                                        coeftab[coeftab$Modnames==model,"abs"] <- round(t(dplyr::select(dplyr::filter(parameterEstimates(modellist[[model]]),label=="abs"),est)),2)
                              }
                    }

                    # shift parameter C for the optimal margin model
                    shifted_models <- c("SSQDposCneg")
                    if (any(shifted_models %in% names(modellist))){ 
                              
                              # identify respective models in fitted models
                              which <- which(shifted_models %in% names(modellist))
                              models <- shifted_models[which]
                              
                              # write their shift values into coeftab
                              coeftab[,"Shift.C"] <- ""
                              
                              for (model in models){
                                        coeftab[coeftab$Modnames==model,"Shift.C"] <- round(t(dplyr::select(dplyr::filter(parameterEstimates(modellist[[model]]),label=="C"),est)),2)
                              }
                    }
                    
                    # position of the vertex of the parabola for the curvilinear effect of S/R models
                    curvi_models <- c("onlyx2pos","onlyy2pos")
                    if (any(curvi_models %in% names(modellist))){ 
                              
                              # identify respective models in fitted models
                              which <- which(curvi_models %in% names(modellist))
                              models <- curvi_models[which]
                              
                              # write their vertex values into coeftab
                              coeftab[,"Vertex"] <- ""
                              
                              for (model in models){
                                        coeftab[coeftab$Modnames==model,"Vertex"] <- round(t(dplyr::select(dplyr::filter(parameterEstimates(modellist[[model]]),label=="V"),est)),2)
                              }
                    }
                    
                    
                    # find out how many percent of the values of S/R lie "left" of the vertex of onlyx2pos/onlyy2pos
                    curvi_x_models <- c("onlyx2pos")
                    curvi_y_models <- c("onlyy2pos")
                    
                    if (any(curvi_models %in% names(modellist))){ 
                              coeftab[,"left.of.V"] <- ""

                              if (any(curvi_x_models %in% names(modellist))){ 
                                        
                                        # identify respective models in fitted models
                                        which <- which(curvi_x_models %in% names(modellist))
                                        models <- curvi_x_models[which]
                                        
                                        # write the percentage of x values left of V into coeftab
                                        for (model in models){
                                                  V <- dplyr::select(dplyr::filter(parameterEstimates(modellist[[model]]),label=="V"),est)
                                                  left_of_V <- round(ecdf(as.matrix(df[df$out==FALSE, IV1]))(V)*100,0)
                                                  coeftab[coeftab$Modnames==model,"left.of.V"] <- round(left_of_V,2)
                                        }
                              }
                              
                              if (any(curvi_y_models %in% names(modellist))){ 
                                        
                                        # identify respective models in fitted models
                                        which <- which(curvi_y_models %in% names(modellist))
                                        models <- curvi_y_models[which]
                                        
                                        # write the percentage of y values left of V into coeftab
                                        for (model in models){
                                                  V <- dplyr::select(dplyr::filter(parameterEstimates(modellist[[model]]),label=="V"),est)
                                                  left_of_V <- round(ecdf(as.matrix(df[df$out==FALSE, IV2]))(V)*100,0)
                                                  coeftab[coeftab$Modnames==model,"left.of.V"] <- round(left_of_V,2)
                                        }
                              }
                    
                    }
                    
					
                    # find out how many percent of the (S,R) points lie "left" of the highest line 
                    # of the models that are shifted perpendicular to the S=R line 
                    # <-> percent of (S,R) left of the line where S+C = R 
		# <-> how many percent of (S,R) points satisfy S+C < R <-> S-R < -C
                    shifted_nonrotated_models <- c("SSQDposCneg")

                    if (any(shifted_nonrotated_models %in% names(modellist))){ 
                              
                              coeftab[,"left.of.ridge"] <- ""
                              
                              # identify respective models in fitted models
                              which <- which(shifted_nonrotated_models %in% names(modellist))
                              models <- shifted_nonrotated_models[which]
                              
                              # write the percentage of S values left of V into coeftab
                              for (model in models){
                                        C <- as.numeric(dplyr::select(dplyr::filter(parameterEstimates(modellist[[model]]),label=="C"),est))
                                        D <- as.matrix(df[df$out==FALSE, IV2] - df[df$out==FALSE, IV1])
                                        Cut_at_Ridge <- cut(D, breaks=c(-Inf, -C, Inf), labels=c("S < R Direction of Ridge", "S > R Direction of Ridge"))
                                        pointpositions_ridge <- round(prop.table(table(Cut_at_Ridge))*100, 2)
                                        coeftab[coeftab$Modnames==model,"left.of.ridge"] <- round(pointpositions_ridge[1],2)
                              }
                    }
                    
                    
                    ## AICc table
                    
                    # build AICctable with AICcmodavg package, delete variables we do not need
                    aictab <- as.data.frame(AICcmodavg::aictab(modellist, modnames=names(modellist)))
                    aictab$ModelLik <- NULL
                    
                    # round all variables in AICc table to reasonable number of digits
                    aictab$AICc <- round(aictab$AICc,2)
                    aictab$Delta_AICc <- round(aictab$Delta_AICc,2)
                    aictab$AICcWt <- round(aictab$AICcWt,3)
                    aictab$LL <- round(aictab$LL,2)
                    aictab$Cum.Wt <- round(aictab$Cum.Wt,3)
                    
                    # rename Akaike weight variable
                    names(aictab)[names(aictab)=="AICcWt"] <- "w"
                    
                    # # print warning if Log-Likelihood is the same for two models in the set
                    # if (any(duplicated(aictab$LL))) {
                    #           dup <- which(duplicated(aictab$LL))
                    #           dupvalues <- unique(aictab$LL[dup])
                    #           duplist <- list()
                    #           for (k in 1:length(dupvalues)){
                    #                     duplist[[k]] <- paste(aictab[aictab$LL==dupvalues[k], "Modnames"], collapse=", ")
                    #           }
                    #           warning(paste0("These models have the same Log-Likelihood: ", duplist, ".\n"))
                    # }
                    
                    # identify models with essentially the same log-likelihood as a simpler model and print warning if there are any
                    frm <- frm(aictab)
                    
                    # add variable with less technical model names to AICc table
                    aictab$nicenames <- NA
                    if ("xandypos" %in% modelset){aictab[aictab$Modnames=="xandypos",]$nicenames <- "Beneficial PSV and Ability"}
                    if ("onlyxpos" %in% modelset){aictab[aictab$Modnames=="onlyxpos",]$nicenames <- "Beneficial PSV Only"}
                    if ("xandyneg" %in% modelset){aictab[aictab$Modnames=="xandyneg",]$nicenames <- "Detrimental PSV and Ability"}
                    if ("onlyxneg" %in% modelset){aictab[aictab$Modnames=="onlyxneg",]$nicenames <- "Detrimental PSV Only"}
                    if ("discrpos" %in% modelset){aictab[aictab$Modnames=="discrpos",]$nicenames <- "Beneficial SE"}
                    if ("discrneg" %in% modelset){aictab[aictab$Modnames=="discrneg",]$nicenames <- "Detrimental SE"}
                    if ("SQDpos" %in% modelset){aictab[aictab$Modnames=="SQDpos",]$nicenames <- "Self-Knowledge"}
                    if ("SSQDposCneg" %in% modelset){aictab[aictab$Modnames=="SSQDposCneg",]$nicenames <- "Optimal Margin"}
                    if ("onlyypos" %in% modelset){aictab[aictab$Modnames=="onlyypos",]$nicenames <- "Beneficial Ability Only"}
                    if ("onlyx2pos" %in% modelset){aictab[aictab$Modnames=="onlyx2pos",]$nicenames <- "Curvilinear PSV"}
                    if ("onlyy2pos" %in% modelset){aictab[aictab$Modnames=="onlyy2pos",]$nicenames <- "Curvilinear Ability"}
                    if ("IApos" %in% modelset){aictab[aictab$Modnames=="IApos",]$nicenames <- "Interaction"}
                    if ("null" %in% modelset){aictab[aictab$Modnames=="null",]$nicenames <- "Null model"}
                    if ("full" %in% modelset){aictab[aictab$Modnames=="full",]$nicenames <- "Full model"}
                    
                    # arrange variable order of aictab 
                    aictab <- aictab[,c(1,8,2:7)]
                    
                    # identify models belonging to the confidence set (first confnr rows of aictab) and indicate them by a 1 in a new confset variable
                    ifelse (any(aictab$Cum.Wt < 0.95),
                            confnr <- max(which(aictab$Cum.Wt < 0.95)) + 1,
                            confnr <- 0)
                    aictab$confset <- NA
                    aictab[1:confnr, "confset"] <- 1
                    aictab$confset[is.na(aictab$confset)] <- 0
                    
		# combine AICctable with table of coefficients
                    aiccoeftab <- merge(aictab, coeftab, by="Modnames")
                    aiccoeftab <- aiccoeftab[match(aictab$Modnames, aiccoeftab$Modnames),]
                    			
                    # If needed for nicer display of AIC table: Rearrange rows by model order in the manuscript and select columns of interest
                    # aiccoeftab$position <- 1:nrow(aiccoeftab)
                    # aiccoeftab_paper <- aiccoeftab[match(modelset, aiccoeftab$Modnames),]
                    # aiccoeftab_paper <- select(aiccoeftab_paper, nicenames, b1, b2, b3, b4, b5, LL, K, AICc, Delta_AICc, w, confset)
                    
                    
                    ## BIC table (for Reviewer 3)
                    
                    # build BIC table with AICcmodavg package
                    bictab <- as.data.frame(AICcmodavg::bictab(modellist, modnames=names(modellist)))
                    bictab$ModelLik <- NULL
                    
                    # round all variables to 3 digits
                    bictab[,3:7] <- round(bictab[,3:7],3)
                    
                    # add variable with less technical model names to AICc table
                    bictab$nicenames <- NA
                    if ("xandypos" %in% modelset){bictab[bictab$Modnames=="xandypos",]$nicenames <- "Beneficial PSV and Ability"}
                    if ("onlyxpos" %in% modelset){bictab[bictab$Modnames=="onlyxpos",]$nicenames <- "Beneficial PSV Only"}
                    if ("xandyneg" %in% modelset){bictab[bictab$Modnames=="xandyneg",]$nicenames <- "Detrimental PSV and Ability"}
                    if ("onlyxneg" %in% modelset){bictab[bictab$Modnames=="onlyxneg",]$nicenames <- "Detrimental PSV Only"}
                    if ("discrpos" %in% modelset){bictab[bictab$Modnames=="discrpos",]$nicenames <- "Beneficial SE"}
                    if ("discrneg" %in% modelset){bictab[bictab$Modnames=="discrneg",]$nicenames <- "Detrimental SE"}
                    if ("SQDpos" %in% modelset){bictab[bictab$Modnames=="SQDpos",]$nicenames <- "Self-Knowledge"}
                    if ("SSQDposCneg" %in% modelset){bictab[bictab$Modnames=="SSQDposCneg",]$nicenames <- "Optimal Margin"}
                    if ("onlyypos" %in% modelset){bictab[bictab$Modnames=="onlyypos",]$nicenames <- "Beneficial Ability Only"}
                    if ("onlyx2pos" %in% modelset){bictab[bictab$Modnames=="onlyx2pos",]$nicenames <- "Curvilinear PSV"}
                    if ("onlyy2pos" %in% modelset){bictab[bictab$Modnames=="onlyy2pos",]$nicenames <- "Curvilinear Ability"}
                    if ("IApos" %in% modelset){bictab[bictab$Modnames=="IApos",]$nicenames <- "Interaction"}
                    if ("null" %in% modelset){bictab[bictab$Modnames=="null",]$nicenames <- "Null model"}
                    if ("full" %in% modelset){bictab[bictab$Modnames=="full",]$nicenames <- "Full model"}
                    
                    # arrange variable order of bictab 
                    bictab <- bictab[,c(1,8,2:7)]
		
                    			
		## Define automatic output
                    
                    # print variables that the analysis is based on
                    #print(variables)
					
                    ## write output table into respective file
                    #dir.create("Result_tables", showWarnings = FALSE)
                    if (nr==""){filename <- paste(DV,IV1,IV2, sep="_")}
                    if (nr!=""){filename <- paste(DV,IV1,IV2,nr, sep="_")}
                    
		# Save result table
                    #write.table(aictab, file=paste0("Result_tables/aictab_", filename,".dat"), sep="\t", row.names=FALSE)
                    # write.csv(
                    #   aictab, 
                    #   file=paste0("pilot-pool/output/aictab_", nr,".csv"),
                    #   row.names = FALSE
                    #   )
                    #write.csv(
                    #  aiccoeftab,
                    #  file=paste0("pilot-pool/output/RSA_", nr,".csv"),
                    #  row.names = FALSE
                    #)
                    #write.table(aiccoeftab[,-which(names(aiccoeftab)=="position")], file=paste0("Result_tables/aiccoeftab_", filename,".dat"), sep="\t", row.names=FALSE)
                    # write.table(aiccoeftab_paper, file=paste0("Result_tables/aiccoeftab_paper_", filename,".dat"), sep="\t", row.names=FALSE)
                    
                    # Print table including AICc variables and coefficients
                    #print(aiccoeftab)
                    
                    ## build output object
                    res <- list(models = modellist,
                                coeftab = coeftab,
			  aictab = aictab,
			  aiccoeftab = aiccoeftab,
			  bictab = bictab,
			  influence.measures = inf,
			  toremove = frm$toremove,
			  # aiccoeftab_paper = aiccoeftab_paper,
                                poly = poly,
                                r2_full = r2_full,
                                p_full = p_full,
                                data_used = df[df$out==FALSE, ],
                                variables = variables,
                                filename = filename,
                                outliers = which(df$out == TRUE)
                    )         
                    
          } # end of if-loop for "full model significant" condition
          
          attr(res, "class") <- "eam"
          
          return(res)
}


######################################
## Note to myself:
## some ways to inspect the result 
## object of the eam function
######################################
#
# # plot specific model (the plot is automatically saved in the folder Result_plots when executing this function)
# plotmodel(eamobject, model="full")
# 
# # eamobject$models is a list of all fitted models, which can be inspected using lavaan-functions
# # For example, show parameters for specific model with:
# par <- parameterEstimates(eamobject$models[["full"]])
# 
# # show table of estimated coefficients for all estimated models
# eamobject$coeftab
# 
# # show polynomial of full model
# eamobject$poly
# 
# # show R^2 of full model and its p-value
# eamobject$r2_full
# eamobject$p_full
# 
# # show removed outliers
# eamobject$outliers
# 
# # inspect data that was used for the analyses (excluding outliers)
# # e.g.:
# nrow(eamobject$data_used)
# 
# # show file name of the table(s) saved during the analyses
# eamobject$filename


######################################
### PREPARE PLOTTING OF SURFACES 
######################################

## create directory for plots
#dir.create("Result_plots", showWarnings = FALSE)

# define selection of axes, projections, and parameters that should be plotted for a specific model
# (depending on the choice of model, different lines are helpful for easier interpretation of the RSA plot)
helpful_axes_project <- function(model){
          axes <- c()
          project <- c()
          param <- FALSE
          
          # select reasonable axes and projection vector for each model
          if(model == "xandypos"){axes <- project <- c()}
          if(model == "onlyxpos"){axes <- project <- c()}
          if(model == "xandyneg"){axes <- project <- c()}
          if(model == "onlyxneg"){axes <- project <- c()}
          if(model == "discrpos"){axes <- project <- c("LOIC")}
          if(model == "discrneg"){axes <- project <- c("LOIC")}
          if(model == "SQDpos"){axes <- project <- c("LOC")}
          if(model == "SSQDposCneg"){axes <- project <- c("LOC","PA1")}
          if(model == "onlyypos"){axes <- project <- c()}
          if(model == "onlyx2pos"){axes <- project <- c()}
          if(model == "onlyy2pos"){axes <- project <- c()}
          if(model == "IApos"){axes <- project <- c("LOC")}
          if(model == "full"){axes <- project <- c("LOC","LOIC","PA1","PA2")
          param <- TRUE}
          
          # build results object
          res <- list(axes = axes, 
                      project = project,
                      param = param)
          return(res)
}

# define function that extracts the regression coefficients of a pre-defined model 
# from the output of the estimate-all-models function, puts them into the RSA plotting function
# and saves and shows the resulting plot simultaneously

# plotmodel <- function(eam, 
#                       type="3d", model="full",
#                       xlim=NULL, ylim=NULL, zlim=NULL, 
#                       xlab=var[2], ylab=var[3], zlab=var[1], b0 = 0,
#                       axes=helpful_axes_project(model)$axes, 
#                       project=helpful_axes_project(model)$project, 
#                       #param=helpful_axes_project(model)$param,
#                       param = TRUE,
#                       main = model,
#                       nr = "",
#                       coefs = FALSE,
#                       legend=TRUE){
#           
#           # extract regression coefficients and variable names of the respective model
#           par <- dplyr::filter(eam$coeftab, Modnames==model)
#           var <- eam$variables
#           
#           # build plot of specified model
#           plot <- plotRSA(x=par$b1, y=par$b2, x2=par$b3, xy=par$b4, y2=par$b5,
#                           b0 = b0,
#                           type=type, main = main, 
#                           xlim=xlim, ylim=ylim, zlim=zlim, 
#                           xlab=xlab, ylab=ylab, zlab=zlab,
#                           cex.axesLabel=1, 
#                           axes=axes, 
#                           project=project, 
#                           param=param,
#                           legend=legend,
#                           coefs = coefs
#                           )
#           
#           # show plot
#           print(plot)
#           
#           # prepare saving of plots 
#           dir.create("Result_plots", showWarnings = FALSE)
#           filename <- paste0("Result_plots/plot_",eam$filename,"_",model,"_",nr,".png")
# 
#           # save plot
#           png(filename=filename,units="in",width=5,height=5,res=300,pointsize=18)
#           print(plot)
#           dev.off()
# }

################################
## Literature cited in this file
################################

# Bollen, K. A., & Jackman, R. W. (1985). Regression diagnostics: An expository treatment of outliers and influential cases. Sociological Methods & Research, 13(4), 510–542. doi:10.1177/0049124185013004004
# Humberg, S., Dufner, M., Schönbrodt, F. D., Geukes, K., Hutteman, R., van Zalk, M. H. W., Denissen, J. J. A., Nestler, S., & Back, M. D. (2018b). Enhanced versus simply positive: A new condition-based regression analysis to disentangle effects of self-enhancement from effects of positivity of self-view. Journal of Personality and Social Psychology, 114(2), 303-322. doi:10.1037/pspp0000134. Preprint available at osf.io/smmh7
# Schönbrodt, F. D. (2016). RSA: An R package for response surface analysis. Retrieved from https://cran.r-project.org/package=RSA
# Rosseel, Y. (2012). lavaan: An R package for structural equation modeling. Journal of Statistical Software, 48(2), 1–36.

