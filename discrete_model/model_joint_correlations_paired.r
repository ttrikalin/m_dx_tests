# this is the modelfile for separate meta-analysis 
# for humerus and femur data 
# the two meta-analyses are run in the same file to allow easy 
# construction of results 
#

model {
    
#### Both tests, complete crosstabs for TPR and FPR: group 1
    for(k1 in 1:K.G1) {
        # observational part of the model
        # x and p order is 11 10 01 00 
        # first humerus, second femur 

        #                  1          2           3           4        5         6
        # eta.xi order is eta1 (hum) eta2 (fem)  eta3 (jtpr) xi1 (hum) xi2 (fem) xi3 (jfpr) 
        #                  
        x.down.cells.G1[k1,1:4] ~ dmulti(p.down.cells.G1[k1,1:4], N.down.G1[k1])
        x.nodown.cells.G1[k1,1:4] ~ dmulti(p.nodown.cells.G1[k1, 1:4], N.nodown.G1[k1])


       
        # p11.DOWN <-ilogit(eta3)
        p.down.cells.G1[k1,1] <- ilogit(eta.xi.G1[k1,3])
        # Sensitivity for test 1 (humerus)
        p.down.1dot.G1[k1] <- ilogit(eta.xi.G1[k1,1])
        # Sensitivity for test 2 (femur)
        p.down.dot1.G1[k1] <- ilogit(eta.xi.G1[k1,2])

        # remaining counts (down): p10, p01, p00
        p.down.cells.G1[k1,2] <-  p.down.1dot.G1[k1] - p.down.cells.G1[k1,1]
        p.down.cells.G1[k1,3] <-  p.down.dot1.G1[k1] - p.down.cells.G1[k1,1]
        p.down.cells.G1[k1,4] <- 1.0 - p.down.cells.G1[k1,1] - p.down.cells.G1[k1,2] -
                                 p.down.cells.G1[k1,3]


        
        # p11.NODOWN
        p.nodown.cells.G1[k1,1] <- ilogit(eta.xi.G1[k1,6])
        # FPR for test 1 (humerus)
        p.nodown.1dot.G1[k1] <- ilogit(eta.xi.G1[k1,4])
        # FPR for test 2 (femur)
        p.nodown.dot1.G1[k1] <- ilogit(eta.xi.G1[k1,5])

        # remaining counts (nodown): p10, p01, p00
        p.nodown.cells.G1[k1,2] <-  p.nodown.1dot.G1[k1] - p.nodown.cells.G1[k1,1]
        p.nodown.cells.G1[k1,3] <-  p.nodown.dot1.G1[k1] - p.nodown.cells.G1[k1,1]
        p.nodown.cells.G1[k1,4] <- 1.0 - p.nodown.cells.G1[k1,1] - p.nodown.cells.G1[k1,2] -
                                 p.nodown.cells.G1[k1,3]
        
        # structural part of the model 
        eta.xi.G1[k1, 1:6] ~dmnorm(Eta.Xi[1:6], Omega[1:6,1:6])
        
    }
    



    # prior for mean -- independent vague priors suffice 
    #                -- otherwise, mvnormal with vague Wishart 
    Eta.Xi[1] ~ dnorm(mu0[1], prec.mu0[1])
    Eta.Xi[2] ~ dnorm(mu0[2], prec.mu0[2])
    Eta.Xi[3] ~ dnorm(mu0[3], prec.mu0[3])
    Eta.Xi[4] ~ dnorm(mu0[4], prec.mu0[4])
    Eta.Xi[5] ~ dnorm(mu0[5], prec.mu0[5])
    Eta.Xi[6] ~ dnorm(mu0[6], prec.mu0[6])

    # S is a diagonal matrix of sd's
    S[1,1] ~ dunif(s.lower[1], s.upper[1])  # or use inverse gamma
    S[2,2] ~ dunif(s.lower[2], s.upper[2])
    S[3,3] ~ dunif(s.lower[3], s.upper[3])
    S[4,4] ~ dunif(s.lower[4], s.upper[4])
    S[5,5] ~ dunif(s.lower[5], s.upper[5])
    S[6,6] ~ dunif(s.lower[6], s.upper[6])

    for(r1 in 1:5) {
        for(c1 in (r1+1):6) {
            S[r1,c1] <- 0
            S[c1,r1] <- 0
        }
    }
    
    # LL' = U'U is the cholesky decomposition of the correlation matrix (L=U')
    L[1,1] <- 1
    L[1,2] <- 0
    L[1,3] <- 0
    L[1,4] <- 0
    L[1,5] <- 0
    L[1,6] <- 0


    L[2,1] <- cos(phi.21)
    L[2,2] <- sin(phi.21)
    L[2,3] <- 0
    L[2,4] <- 0
    L[2,5] <- 0
    L[2,6] <- 0

    L[3,1] <- cos(phi.31)
    L[3,2] <- cos(phi.32)*sin(phi.31)
    L[3,3] <- sin(phi.32)*sin(phi.31)
    L[3,4] <- 0
    L[3,5] <- 0
    L[3,6] <- 0
    
    L[4,1] <- cos(phi.41)
    L[4,2] <- cos(phi.42)*sin(phi.41)
    L[4,3] <- cos(phi.43)*sin(phi.42)*sin(phi.41)
    L[4,4] <- sin(phi.43)*sin(phi.42)*sin(phi.41)
    L[4,5] <- 0
    L[4,6] <- 0

    L[5,1] <- cos(phi.51)
    L[5,2] <- cos(phi.52)*sin(phi.51)
    L[5,3] <- cos(phi.53)*sin(phi.52)*sin(phi.51)
    L[5,4] <- cos(phi.54)*sin(phi.53)*sin(phi.52)*sin(phi.51)
    L[5,5] <- sin(phi.54)*sin(phi.53)*sin(phi.52)*sin(phi.51)
    L[5,6] <- 0

    L[6,1] <- cos(phi.61)
    L[6,2] <- cos(phi.62)*sin(phi.61)
    L[6,3] <- cos(phi.63)*sin(phi.62)*sin(phi.61)
    L[6,4] <- cos(phi.64)*sin(phi.63)*sin(phi.62)*sin(phi.61)
    L[6,5] <- cos(phi.65)*sin(phi.64)*sin(phi.63)*sin(phi.62)*sin(phi.61)
    L[6,6] <- sin(phi.65)*sin(phi.64)*sin(phi.63)*sin(phi.62)*sin(phi.61)


    phi.21 ~ dunif(chol.r.lower[1], chol.r.upper[1]) 
    phi.31 ~ dunif(chol.r.lower[2], chol.r.upper[2]) 
    phi.32 ~ dunif(chol.r.lower[3], chol.r.upper[3]) 
    phi.41 ~ dunif(chol.r.lower[4], chol.r.upper[4]) 
    phi.42 ~ dunif(chol.r.lower[5], chol.r.upper[5]) 
    phi.43 ~ dunif(chol.r.lower[6], chol.r.upper[6])
    phi.51 ~ dunif(chol.r.lower[7], chol.r.upper[7]) 
    phi.52 ~ dunif(chol.r.lower[8], chol.r.upper[8]) 
    phi.53 ~ dunif(chol.r.lower[9], chol.r.upper[9]) 
    phi.54 ~ dunif(chol.r.lower[10], chol.r.upper[10]) 
    phi.61 ~ dunif(chol.r.lower[11], chol.r.upper[11]) 
    phi.62 ~ dunif(chol.r.lower[12], chol.r.upper[12]) 
    phi.63 ~ dunif(chol.r.lower[13], chol.r.upper[13]) 
    phi.64 ~ dunif(chol.r.lower[14], chol.r.upper[14]) 
    phi.65 ~ dunif(chol.r.lower[15], chol.r.upper[15]) 


    R <- L %*% t(L)
    Tau <- S %*% R %*% S
    Omega <- inverse(Tau) 
        

    
    # Data to report 
    tpr.hum <- ilogit(Eta.Xi[1])
    tpr.fem <- ilogit(Eta.Xi[2])
    jtpr <- ilogit(Eta.Xi[3])
    fpr.hum <- ilogit(Eta.Xi[4])    
    fpr.fem <- ilogit(Eta.Xi[5])
    jfpr <- ilogit(Eta.Xi[6])

    diff.tpr <- tpr.hum - tpr.fem
    diff.fpr <- fpr.hum - fpr.fem
    
    or.tpr <- exp(Eta.Xi[1]-Eta.Xi[2])
    or.fpr <- exp(Eta.Xi[4]-Eta.Xi[5])

    diff.eta <- Eta.Xi[1] - Eta.Xi[2] 
    diff.xi <- Eta.Xi[4] - Eta.Xi[5] 
}
