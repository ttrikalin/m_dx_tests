// Research code for paper  
// "Methods for the Joint Meta-Analysis of Multiple Tests" 
// TA Trikalinos, DC Hoaglin, KM Small, N Terrin, CH Schmid
// 
// Normal approximation models (see appendix to paper)
//

//
// Fit 2 variants of random effects with REML. They differ 
// in the way the T matrix is structured (random effect variance)
//
// case A: unstructured 21 parameters for 2 tests)
// case B: common sigma_eta and sigma_xi, sigma_eta*, sigma_xi* 
//         and 9 correlations (13 parameters).
//		   Case B (structured matrix) does not converge in the example
//
// For case A, we show code with Ian White's mvmeta which uses a 
// data-augmentation scheme to handle values that are missing completely
// at random, as well as custom research code developed to avoid the data
// augmentation scheme. 

clear
insheet using paired_data.txt // has ready-calculated empirical logit-TPR,
						      // logit-FPR and covariances. 
// b1: \eta_1, logit(TPR), humerus
// b2: \eta_2, logit(TPR), femur
// b3: \eta_*, logit(JTPR)
// b4: \xi_1, logit(FPR), humerus
// b5: \xi_2, logit(FPR), femur
// b6: \xi_*, logit(JFPR)
// sij: the respective covariance matrix elements, calculated using the 
// formulas in the appendix  

local p 6   // the length of the parameter vector (needed below)

// Hand-correct singular covariance matrices 
// Adding epsilon to the diagonal is the same as adding epsilon 
// to all eigenvalues 
// See the appendix of TA Trikalinos, DC Hoaglin and CH Schmid, Stat Med 2014 
// "An empirical comparison of univariate and multivariate meta-analysis"
// for as=n algorithm for choosing the correction factor for 
// reproducible reseach code

local epsilon = 0.5 
replace s11 = s11 + `epsilon' if inlist(_n, 4, 7)
replace s22 = s22 + `epsilon' if inlist(_n, 4, 7)
replace s33 = s33 + `epsilon' if inlist(_n, 4, 7)


///////////////////////////////////////////////////////////
// Case A - unstructured random effects (paired design)  //
///////////////////////////////////////////////////////////

// mvmeta code - install Ian White's mvmeta 
mvmeta b s


// research code -- does not use the data augmentation scheme 
qui count 
local n = r(N)

global restricted = 1 // do REML rather than ML 
// These globals and matrices are needed to pass information to the likelihood 
// optimization programs. 

global n `n'
global ymat  "ymat"
global Smat  "Smat"
global p `p'

forval i=1/`n' {
    mat ymat`i' = J(1, `p', 0)  
    forval j =1/`p' {
        qui summ b`j' in `i' , meanonly      
        mat ymat`i'[1,`j'] = r(mean) 
    }
    if (matmissing(ymat`i')) {
        noi di in yellow "Study `i' has missing values:"
        noi mat li ymat`i'
    }

    mat Smat`i' = J(`p', `p' , 0)
    forval k =1/`p' {
        forval m =`k'/`p' {
            qui summ  s`k'`m' in `i' , meanonly 
            mat Smat`i'[`k', `m'] = r(mean)
            mat Smat`i'[`m', `k'] = r(mean)
        }
    }   
}



// To obtain starting values, fit a fixed effects model with GLS. 
if (0==0) {  
    mat Wsum = J(`p', `p', 0)
        mat Wysum = J(1, `p', 0)

    // This works out to be the same as the method by 
    // Gleser and Olkin referenced in the paper.
    // this strategy will not work OK for the REML calculations!
    // If we havemissing data we can calculate them by padding with 0's

        forval k=1/`n' {
                mat W`k' = syminv(Smat`k')
                mat Wysum = Wysum + ymat`k' * W`k'
                mat Wsum = Wsum + W`k'
        }

        mat covbetaFixedGLS = syminv(Wsum)
        mat betaFixedGLS = Wysum * covbetaFixedGLS
    mat betaFixed = betaFixedGLS 

}

// Assume that study_K has missing values. Its likelihood contribution is different
// than that of the other studies.  
// Here we replace the missing values 
// with *structural* zeros.  
// The REML programs will not utilize the zeros; they will work only with the 
// non zero entries 

count 
forval K=1/`r(N)' {
    forval i=1/`p' {
        if (el(ymat`K',1,`i')>=.)  mat ymat`K'[1,`i']=0
        forval j=1/`p' {
            if (el(Smat`K',`i',`j')>=.)  mat Smat`K'[`i',`j']=0
        }
    }
}


// All matrices must be non-singular. This is guaranteed because they are diagonal!
// These are the three examples of random effects that are discussed in the paper

if (1==1) {
    // case A: unstructured 21 parameters for 2 tests)

    program drop _all 
    global be_verbose = 1  // force ll program to identify itself - avoid blunders 

    ml model d0 ll_caseA_miss (mu: b1 b2 b3 b4 b5 b6 , nocons) ///
        (S11:) (S12:) (S13:) (S14:) (S15:) (S16:) ///
               (S22:) (S23:) (S24:) (S25:) (S26:) ///
                      (S33:) (S34:) (S35:) (S36:) ///
                             (S44:) (S45:) (S46:) ///
							        (S55:) (S56:) /// 
							               (S66:), obs(`n')  collinear 

    ml search S11: -10 10 S12: -10 10 S13: -10 10 S14: -10 10 S15: -10 10 S16: -10 10 ///
                          S22: -10 10 S23: -10 10 S24: -10 10 S25: -10 10 S26: -10 10 ///
                                      S33: -10 10 S34: -10 10 S35: -10 10 S36: -10 10 ///
                                                  S44: -10 10 S45: -10 10 S46: -10 10 ///          
                                                              S55: -10 10 S56: -10 10 ///    
															              S66: -10 10

    mat b0 = [betaFixed, 1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,1,0,1] 
    ml init b0 , copy skip
    //ml check

    ml max , difficult iterate(30) ltolerance(1e-7)
    est save multi_estimates, replace
    est store A
    mat betaA = BETA
    mat TauA = T 
    mat covbetaA = COVBETA
}



// The same analyses, treating the studies as parallel designs 

if (1==1) {
    // case Anonpaired: unstructured 21 parameters for 2 tests)
	global Smat "Smat_indep"
	forval i=1/`n' {
		mat Smat_indep`i' = diag(vecdiag(Smat`i')) 
	}
    program drop _all 

    global be_verbose = 1  // force ll program to identify itself - avoid blunders 


    ml model d0 ll_caseA_miss (mu: b1 b2 b3 b4 b5 b6 , nocons) ///
        (S11:) (S12:) (S13:) (S14:) (S15:) (S16:) ///
               (S22:) (S23:) (S24:) (S25:) (S26:) ///
                      (S33:) (S34:) (S35:) (S36:) ///
                             (S44:) (S45:) (S46:) ///
							        (S55:) (S56:) ///
							               (S66:), obs(`n')  collinear 

    ml search S11: -10 10 S12: -10 10 S13: -10 10 S14: -10 10 S15: -10 10 S16: -10 10 ///
                          S22: -10 10 S23: -10 10 S24: -10 10 S25: -10 10 S26: -10 10 ///
                                      S33: -10 10 S34: -10 10 S35: -10 10 S36: -10 10 ///
                                                  S44: -10 10 S45: -10 10 S46: -10 10 ///          
                                                              S55: -10 10 S56: -10 10 ///    
															              S66: -10 10

    mat b0 = [betaFixed, 1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,1,0,1] 
    ml init b0 , copy skip
    //ml check

    ml max , difficult iterate(30) ltolerance(1e-7)
    est save multi_estimates, append
    est store Anonpaired
    mat betaAnonpaired = BETA
    mat TauAnonpaired = T 
    mat covbetaAnonpaired = COVBETA

	// restore the covariance matrix
	global Smat "Smat"
}



///////////////////////////////////////////////////////////
// Case B - structured random effects (paired design)    //
///////////////////////////////////////////////////////////

/* 
// Case B does not converge in the current example 
// kept in the appendix code for completeness 
// needs re-parameterization of the covariance matrix to ensure 
// robust convergence -- TO DO. 

if (10==1) {
    // case B: common sigma for all eta and xi, and 4 correlations  (5 parameters)
    // parameters passed in the following order:
    // log of sigma 
    // tanh^-1 of rho_eta rho_eta_xi rho_star_eta_xi rho_xi
    program drop _all 
    global be_verbose = 1  // force ll program to identify itself - avoid blunders 

    ml model d0 ll_caseB_miss (mu: b1 b2 b3 b4 b5 b6 , nocons) ///
        (sigma_eta:) (sigma_xi:) (sigma_eta_star:) (sigma_xi_star:) ///
        (rho_eta:) (rho_xi:) (rho_star:) (rho_eta_xi:) ///
        (varrho_eta_xi:) (rho_eta_star:) (rho_xi_star:)   ///
		(varrho_eta_star:)  (varrho_xi_star:)   ///
        , obs(`n')  collinear 



    ml search sigma_eta: -3 0 sigma_xi: -3 0 ///
	         sigma_eta_star: -3 0 sigma_xi_star: -3 0 ///
             rho_eta: -1 2 rho_xi: -1 2 rho_star: -1 2 rho_eta_xi: -1 2 ///
             varrho_eta_xi: -1 2 rho_eta_star: -1 2  rho_xi_star: -1 2   ///
	         varrho_eta_star: -1 2  varrho_xi_star: -1 2

    
    mat b0 = [betaFixed, -0.259, -1.096, -1.149, -1.335, ///
	0.401, ///
	0.149, ///
	0.244, ///
	0.504, /// 
	-0.259, ///
	1.175, ///
	1.027, ///
	0.312, ///
	0.041] 
    ml init b0 , copy skip
    //ml check

    ml max , difficult iterate(30) ltolerance(1e-5)
    est save multi_estimates, append
    est store B
    mat betaB = BETA
    mat TauB = T 
    mat covbetaB = COVBETA
}

*/

