capture program drop XTEWreg
program define XTEWreg, eclass
	version 12
	//syntax varlist(min=2 numeric) [if] [in] , MAXDeg(integer) [MISmeasured(integer 1)  METHod(string) PANmethod(string) BXint(numlist) HAScons NOPRN VCE(string) OPTim(integer 2) CENTmom(string)]
	syntax varlist(min=2 numeric) [if] [in] , MAXDeg(integer) [MISmeasured(integer 1)  METHod(string) PANmethod(string) BXint(numlist) CENTmom(string) HAScons NOPRN]
	marksample touse

	local optim = 2

	quietly count if `touse'
	if `r(N)' == 0 error 2000
	
	if (`mismeasured'<1) {
		di as err "Option <mismeasured> is required to be at least 1. Aborting."
		exit 197
	}
	if (`maxdeg'<3) {
		di as err "Option <maxdeg> is required to be at least 3. Aborting."
		exit 197
	}
	if (`optim'!=1 & `optim'!=2) {
		di as err "Option <optim> is required to be either 1(DFP) or 2(GN). Aborting."
		exit 197
	}
	
	local vce = strupper("`vce'")
	if ("`vce'"=="" | "`vce'"=="TWOSTEP") 	local vcemet = 0
	else if ("`vce'"=="OLS") 				local vcemet = 1
	else {
		di as err "Unknown value provided for VCE. Aborting."
		exit 197
	}

	local centmom = strupper("`centmom'")
	if ("`centmom'"=="" | "`centmom'"=="RESET") local centmet = 0
	else if ("`centmom'"=="SET") 				local centmet = 1
	else if ("`centmom'"=="USE") 				local centmet = 2
	else {
		di as err "Unknown value provided for CentMom. Aborting."
		exit 197
	}
	

	// separate varlist
	gettoken depvar varlist: varlist
	local wc = wordcount("`varlist'")
	if (`wc'<`mismeasured') {
		di as err "Less regressors provided than specified in option <mismeasurd>. Aborting."
		exit 197
	}
	local misindep = ""
	forvalues i=1/`mismeasured' {
		gettoken mis varlist: varlist
		local misindep = "`misindep' `mis'"
	}
	local indep = "`varlist'"

	// verify bxint
	if ("`bxint'"!="" & ///
		wordcount("`bxint'")/`mismeasured'!=round(wordcount("`bxint'")/`mismeasured')) {
		di as err "If BXint is specified, it must have c*<mismeasured> members for some integer c. Aborting."
		exit 197
	}
	if "`bxint'"=="" local bxint = "0"

	// parse method
	local method = strupper("`method'")
	if ("`method'"=="" | "`method'"=="CML") {
		local met = 0
	}
	else if ("`method'"=="MOM" | "`method'"=="GMM") {
		local met = 1
	}
	else {
		di as err "Unknown value provided for METHOD. Aborting."
		exit 197
	}
	
	// parse panmethod
	local panmethod = strupper("`panmethod'")
	if ("`panmethod'"=="" | "`panmethod'"=="CLS") {
		local clustmet = 0
	}
	else if ("`panmethod'"=="CMD") {
		local clustmet = -1
	}
	else if ("`panmethod'"=="IDN") {
		local clustmet = 1
	}
	else if ("`panmethod'"=="NON") {
		local clustmet = 2
	}
	else {
		di as err "Unknown value provided for PANMETHOD. Aborting."
		exit 197
	}

	qui: xtset
	local idname `r(panelvar)'
	local tmname `r(timevar)'

	if "`idname'" == "" {
		di "Warning: xtset was not run before this [XT] command. Use xtset <panel variable>. Continuing without clustering."
	}
	if "`tmname'" == "" & `clustmet'==-1 {
		di as error "Error: xtset was not used to set a time variable for panmethod(CMD) to use. Aborting.\n"
		exit(197)
	}
	
	tempname b V cst
	mata: doEW("`depvar'", "`misindep'", "`indep'", "`idname'", "`tmname'", `met', `clustmet', `maxdeg', "`bxint'", 0, "`hascons'", `optim',  `centmet', "`touse'")
	local cst = cond("`hascons'" == "", "_cons", "")
	local vnames `misindep' `cst' `indep'
	mat `b' = r(beta)
	mat `V' = r(VCmat)
	matname `V' `vnames'
	matname `b' `vnames', c(.)
	local N = r(N)
	ereturn post `b' `V', depname(`depvar') obs(`N') esample(`touse')

	ereturn scalar rho = r(rho)
	ereturn matrix tau = tau
	ereturn matrix serr  = serr
	ereturn scalar SErho = r(SErho)
	ereturn matrix SEtau = SEtau
	ereturn matrix vcrhotau  = vcrhotau
	ereturn matrix w  = w
	ereturn local  bxint = "`bxint'"
	ereturn local  vcetype = "`vce'"
	ereturn local  panmethod = "`panmethod'"
	ereturn scalar Jstat = r(Jstat)
	ereturn scalar Jval = r(Jval)
	ereturn scalar dfree = r(dfree)
	ereturn scalar obj = r(obj)
	ereturn matrix IDval = IDval
	ereturn matrix IDstat = IDstat
	

	if ("`noprn'"=="") {
		display _newline "`method'`maxdeg'(`mismeasured') EIV results" _col(65) "N = " %9.0g e(N)
		display _col(65) "Rho^2 = " %5.3f r(rho)
		display _col(65) "       (" %5.3f r(SErho) ")"
		ereturn display
		forvalues i=1/`mismeasured' {
			display "Tau" `i' "^2: " %5.3f el(e(tau),`i',1) " (" %5.3f el(e(SEtau),`i',1) ")"
		}
		if ( e(Jstat) > 0 ) {
			display "Sargan-Hansen J statistic: " %7.3f e(Jstat) ///
					"   (p=" %5.3f e(Jval) ", d=" e(dfree) ")"
		}
	}

end



version 12
mata:

mata clear

// return value struct
struct stats {
	real 	colvector 	beta
	real 	matrix 		VCmat
	real 	colvector 	serr
	real 	scalar 		N
	real 	scalar 		Jstat
	real 	scalar 		Jval
	real 	scalar 		dfree
	real 	scalar 		rho
	real 	colvector	tau
	real 	scalar 		SErho
	real 	colvector	SEtau
	real 	matrix 		vcrhotau
	real 	matrix 		inflncXZ
	real 	matrix 		inflncRT
	real 	matrix 		w
	real 	scalar 		obj
	real 	colvector 	IDval
	real 	colvector 	IDstat
}
// end struct stats


///////////////////////////////////////////////////////////////////////////////////////
// The MATA entrypoint - mostly communicates with stata and calls doCLS/doCMD
///////////////////////////////////////////////////////////////////////////////////////
void doEW(string scalar depvar, string scalar misindep,	string scalar indepvars, 	///
		string scalar idname, string scalar tmname,	numeric scalar met, 			///
		numeric scalar clustmet, numeric scalar maxdeg, string scalar bXint, 		///
		numeric scalar vcemet, string scalar hascons, numeric scalar optmet, 		///
		real scalar centmet, string scalar touse)
{
	st_view(y, ., depvar, touse)
	n  = rows(y)
	st_view(x=., ., tokens(misindep), touse)
	st_view(z=., ., tokens(indepvars), touse)
	
	if (idname=="") id = J(n,1,1)
	else st_view(id, ., idname, touse)
	
	if (tmname=="") tm = J(n,1,1)
	else st_view(tm, ., tmname, touse)

	bXint1 = strtoreal(tokens(bXint))
	
	//verify no missing in y,x,z,id,tm. we don't deal with missing values well.
	if (sum(y:==.)+sum(x:==.)+sum(z:==.)+sum(id:==.)+sum(tm:==.)>0) 
	{
		errprintf("Command does not support missing values. Aborting.\n")
		exit(error(197))
	}
	
	if (hascons=="")
	{
		con = J(n,1,1)
		z = (con,z)
	}
	
	// make sure everything is sorted on id, will be crucial later
	bigmat 	= (id,tm,y,x,z)
	bigmat 	= sort(bigmat,(1,2,3))
	id 		= bigmat[,1]
	tm 		= bigmat[,2]
	y 		= bigmat[,3]
	x 		= bigmat[,4..(cols(x)+4-1)]
	z 		= bigmat[,(cols(x)+4)..cols(bigmat)]
	
	struct stats 	scalar 		retval
	class EWreg 	scalar 		ew
	
	ew = EWreg()

	if (clustmet==-1) { //do classical minimum distance
		retval = ew.doCMD(id, tm, y, x, z, met, maxdeg, bXint1, vcemet, optmet)
	}
	else { // do clustered weighting matrix
		retval = ew.doCLS(id, y, x, z, met, clustmet, maxdeg, bXint1, vcemet, optmet, centmet)
	}

	// Return values
	st_matrix("r(beta)",(retval.beta)')
	st_matrix("r(VCmat)",retval.VCmat)
	st_numscalar("r(N)",n)
	st_matrix("serr",retval.serr)
	st_numscalar("r(Jstat)",retval.Jstat)
	st_numscalar("r(Jval)",retval.Jval)
	st_numscalar("r(dfree)",retval.dfree)
	st_numscalar("r(rho)",retval.rho)
	st_numscalar("r(SErho)",retval.SErho)
	st_matrix("tau",retval.tau)
	st_matrix("SEtau",retval.SEtau)
	st_matrix("vcrhotau",retval.vcrhotau)
	st_matrix("w",retval.w)
	st_numscalar("r(obj)",retval.obj)
	st_matrix("IDval",retval.IDval)
	st_matrix("IDstat",retval.IDstat)
}
// end function doEW

end
