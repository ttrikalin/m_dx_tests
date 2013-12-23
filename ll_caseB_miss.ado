*! version 0.5 15aug2012
*! thomas_trikalinos@brown.edu 
// case B in the report -- hard coded for 2 tests only !!!

prog def ll_caseB_miss

if ($be_verbose) noi di as text _n "THIS HANDLES ONLY 2 TESTS AND THE STRUCTURE T_B in the report"
if ($be_verbose) noi di as text    "                                   [can parse missing values]" _n 

global be_verbose = 0

args todo b  lnf

local y $ymat
local S $Smat
local n $n
//local p $p
local p=6     // hard coded!!!

local epsilon 1e-6  // for back transforming correlations

tempname BETA 
tempname C C_eta C_xi C_offdiag C_upper C_lower SD T W dev minus2ll Wsum ll
tempname X P Wsum_miss

// get input arguments  
mat `BETA' = J(1, `p', 0)
forval i =1/`p' {
	mat `BETA'[1, `i']=el(`b',1,`i')
}
local k `p'

// variance passed as log-transformed value

foreach v in sigma_eta sigma_xi sigma_eta_star sigma_xi_star {	
	local ++k
	scalar `v' = exp(el(`b', 1, `k'))
	if (`v' >= .) {
		scalar `v' = `epsilon'
	}
}

// correlations passed as tanh transformed values
foreach v in rho_eta rho_xi rho_star rho_eta_xi ///
             varrho_eta_xi rho_eta_star rho_xi_star   ///
             varrho_eta_star varrho_xi_star {

	local ++k 
	scalar `v'1 = el(`b', 1, `k')	 
	local width = 0.50
	local start = 0.47
	scalar `v' = `width'*invlogit(`v'1) + `start'
	
	if (`v' >=. ) {
	        if (sign(`v'1)==-1) {
	                scalar `v' = `start' + `epsilon'
	        }
	        if (sign(`v'1)==1) {
	                scalar `v' = `start' + `width' - `epsilon'
	        }
	}
}



//  Now for the correlation matrix assuming that you get the eta's first 
//  and then the xi's in the same order


mat `T' = J(`p', `p', 0)

// diagonals 
mat `T'[1,1] = sigma_eta^2
mat `T'[2,2] = sigma_eta^2
mat `T'[3,3] = sigma_eta_star^2
mat `T'[4,4] = sigma_xi^2
mat `T'[5,5] = sigma_xi^2
mat `T'[6,6] = sigma_xi_star^2

// Block 1: T_eta
mat `T'[1,2] = rho_eta*sigma_eta^2
mat `T'[1,3] = rho_eta_star*sigma_eta*sigma_eta_star
mat `T'[2,3] = rho_eta_star*sigma_eta*sigma_eta_star

// Block 2: T_eta_xi
mat `T'[1,4] = rho_eta_xi*sigma_eta*sigma_xi
mat `T'[1,5] = varrho_eta_xi*sigma_eta*sigma_xi
mat `T'[1,6] = varrho_eta_star*sigma_eta*sigma_xi_star


mat `T'[2,4] = varrho_eta_xi*sigma_eta*sigma_xi
mat `T'[2,5] = rho_eta_xi*sigma_eta*sigma_xi
mat `T'[2,6] = varrho_eta_star*sigma_eta*sigma_xi_star

mat `T'[3,4] = varrho_xi_star*sigma_eta_star*sigma_xi
mat `T'[3,5] = varrho_xi_star*sigma_eta_star*sigma_xi
mat `T'[3,6] = rho_star*sigma_eta_star*sigma_xi_star


// Block 4: T_xi
mat `T'[4,5] = rho_xi*sigma_xi^2
mat `T'[4,6] = rho_xi_star*sigma_xi*sigma_xi_star
mat `T'[5,6] = rho_xi_star*sigma_xi*sigma_xi_star

// Mirror below the diagonal


forval r=2/6 {
	forval c=1/`=`r'-1' {
		mat `T'[`r',`c'] = `T'[`c',`r']
	}
}

mat `Wsum' = J(`p',`p',0)
scalar `ll'= 0

forvalues i = 1/`n' {

	local hasmissing = (diag0cnt(`S'`i')>0)
	if (`hasmissing' == 0) {
		cap mat `W' = invsym(`S'`i'+`T')
		if  _rc {
			noi di as error "Problem in study `i'"
			//exit -1
			
			mat tocorrect = `S'`i'+`T'
			correctmemat, matname(tocorrect)
			cap mat `W' = invsym(r(M))
			if _rc {
				noi di as error "Problem in study `i' after correction"
				exit -1
			}
  		}
		else {
			mat `dev' = `y'`i'-`BETA'
			mat `minus2ll' = `p'*log(2*_pi) - log(det(`W')) + `dev' * `W' * `dev''
			mat `Wsum' = `Wsum' + `W'
		}
	}
	if (`hasmissing'==1) {
		// get the permutation matrix to isolate the non-missing rows/columns
		gimmepermmat, matname(`S'`i')

		local nonmissing = `=colsof(`S'`i')' - diag0cnt(`S'`i')

		mat `P' = r(P)
		mat `X' = `P'*(`S'`i'+`T')*`P''
		mat `X' = `X'[1..`nonmissing', 1..`nonmissing']
		mat `W' = syminv(`X')
		mat `dev' = (`y'`i'-`BETA')*`P''
		mat `dev' = `dev'[1,1..`nonmissing']
		mat `minus2ll' = `=colsof(`S'`i')'*log(2*_pi) - log(det(`W')) + `dev' * `W' * `dev''

		// pad W with 0's for the missing rows/columns ma
		mat  `Wsum_miss' = J(`p',`p', 0)
		forval j=1/`nonmissing' {
			forval k=1/`nonmissing' {
				mat `Wsum_miss'[`j', `k'] = `W'[`j', `k']
			}
		}

		// rearrange to match the row/column order of the pxp Wsum
		// and add to complete Wsum 
		mat `Wsum_miss' = `P''*`Wsum_miss'*`P'
		mat `Wsum' = `Wsum' + `Wsum_miss'
	}
	scalar `ll' = `ll' - el(`minus2ll',1,1)/2
}
if ($restricted == 1 ) {
	// of there is a study with a missing outcome, the constant part is not correct! - but this does not matter!!!
	scalar `ll' = `ll' - log(det(`Wsum'))/2 + `p'*log(2*_pi)/2
}


scalar `lnf' = `ll'

mat BETA = `BETA'
mat T = `T'
mat COVBETA = invsym(`Wsum')
end
