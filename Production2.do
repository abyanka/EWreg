/* General initial setup */
qui {
clear all
set matsize 10000

global plc = 0
if ${plc} == 0 {
	global workfold = "C:\Users\finguy\SkyDrive\Documents\Toni\\"
}
if ${plc} == 1 {
	global workfold = "D:\SkyDrive\Documents\Toni\\"
}
global scriptfold = "Scripts\"
global resfold = "Results\"
global datfold = "Data\"

cd ${workfold}
adopath + "${workfold}${scriptfold}"

/*
mata: 
mata clear

do "${scriptfold}Symbolic.mata"
do "${scriptfold}Partitions.mata"
do "${scriptfold}EWproblem.mata"
do "${scriptfold}EWopt.mata"
do "${scriptfold}EWdata.mata"
do "${scriptfold}EWreg.mata"

mata: 
mata mlib create lEW , dir("${workfold}${scriptfold}") replace
mata mlib add lEW Symbolic()
mata mlib add lEW sereqn()
mata mlib add lEW Partitions()
mata mlib add lEW EWproblem()
mata mlib add lEW serprob()
mata mlib add lEW EWopt()
mata mlib add lEW EWdata()
mata mlib add lEW EWreg()
mata mlib add lEW i_crit()
mata mlib add lEW stats()
mata mlib add lEW mt()

mata mlib index
end
*/
}
// compile block

use "${datfold}EJW1.dta", clear

//qui: xtset gvkey
//XTEWreg2 lever mtb tangib logsales oi , mis(2) maxd(3) meth(MOM) pan(CLS) optim(2)

qui: tsset gvkey fyear

//XTEWreg1 ik q cfk , meth(CML4) pan(CLS) 
//XTEWreg1 lever mtb tangib logsales oi , meth(CML3) pan(CLS) two
//disp e(rho) " " e(SErho) " " e(tau1) " " e(SEtau1) " " e(tau2) " " e(SEtau2) " " ///
//	e(Jstat) " " e(Jval)

XTEWreg0 ik q cfk , meth(CML4) pan(CLS) 

//XTEWreg ik q cfk , meth(GMM3)



/************************************/
/* Prepare non-bootstrap versions	*/
/************************************/
local t1 	= "ik#q#cfk#CML#CLS#5#6"
local t2 	= "ik#q#cfk#MOM#CLS#5#6"
local t3 	= "lever#mtb tangib#logsales oi#CML#CLS#8#6"
local t4 	= "lever#mtb tangib#logsales oi#MOM#CLS#8#6"
local t5 	= "lever#mtb tangib logsales#oi#CML#CLS#9#6"
local t6 	= "lever#mtb tangib logsales#oi#MOM#CLS#9#6"

forvalues tst = 5/6 {
	//if `tst'!=1 {
	//	continue
	//}
	tokenize "`t`tst''" , p("#")
	local dep 	= "`1'"
	local mis 	= "`3'"
	local zvar 	= "`5'"
	local meth 	= "`7'"
	local pan  	= "`9'"
	local rcnt 	= "`11'"
	local ccnt 	= "`13'"

	matrix results = J(`rcnt',`ccnt',0)
	matrix sd_results = results
	
	local nx = wordcount("`mis'")
	local nz = wordcount("`zvar'")
	local nv = `nx' + `nz'
	

	forvalues i=1/`ccnt' {
		local i1 = `i'+2
		XTEWreg2 `dep' `mis' `zvar', meth("`meth'") pan("`pan'") maxd(`i1') mis(`nx') optim(2)
		mat tmp1 = e(b)
		mat tmp2 = e(serr)'
		
		local iv = 1

		forvalues j=1/`nx' {
			qui: mat results[`iv',`i'] = el("tmp1",1,`j')
			qui: mat sd_results[`iv',`i'] = el("tmp2",1,`j')
			local iv = `iv' + 1
		}

		forvalues j=1/`nz' {
			local j1 = `nx' + `j' + 1
			qui: mat results[`iv',`i'] = el("tmp1",1,`j1')
			qui: mat sd_results[`iv',`i'] = el("tmp2",1,`j1')
			local iv = `iv' + 1
		}

		qui: mat results[`iv',`i'] = `e(rho)'
		qui: mat sd_results[`iv',`i'] = `e(SErho)'
		local iv = `iv' + 1

		forvalues j=1/`nx' {
			local tau = el(e(tau),`j',1)
			local SEtau = el(e(SEtau),`j',1)
			qui: mat results[`iv',`i'] = `tau'
			qui: mat sd_results[`iv',`i'] = `SEtau'
			local iv = `iv' + 1
		}

		qui: mat results[`iv',`i'] = `e(Jstat)'
		qui: mat sd_results[`iv',`i'] = `e(Jval)'
	}
	
	mat list results
	mat list sd_results
	outtable2 using ${resfold}P1_`tst', /*sig("B")*/ mat(results) nobox replace f(%6.3f) format2(%5.3f) mat2(sd_results)
}


bootstrap _b e(rho) el(e(tau),1,1) el(e(tau),2,1), rep(50) seed(1337) : XTEWreg2 lever mtb tangib logsales oi, meth(CML) pan(CLS) maxd(4) mis(2) optim(2)


/********************************/
/* Prepare bootstrap versions	*/
/********************************/
local t1 	= "ik#q#cfk#CML#CLS#4#1"
local t2 	= "ik#q#cfk#MOM#CLS#4#4"
local t3 	= "lever#mtb tangib#logsales oi#CML#CLS#7#4"
local t4 	= "lever#mtb tangib#logsales oi#MOM#CLS#7#4"

forvalues tst = 1/1 {
	//if `tst'!=1 {
	//	continue
	//}
	tokenize "`t`tst''" , p("#")
	local dep 	= "`1'"
	local mis 	= "`3'"
	local zvar 	= "`5'"
	local meth 	= "`7'"
	local pan  	= "`9'"
	local rcnt 	= "`11'"
	local ccnt 	= "`13'"

	matrix results = J(`rcnt',`ccnt',0)
	matrix sd_results = results
	
	local nx = wordcount("`mis'")
	local nz = wordcount("`zvar'")
	local nv = `nx' + `nz'
	
	local bstr = "bootstrap _b e(rho)"
	forvalues j=1/`nx' {
		local bstr = "`bstr' el(e(tau),`j',1)"
	}
	local bstr = "`bstr', rep(50) seed(1337) : "
	

	forvalues i=1/`ccnt' {
		local i1 = `i'+2
		`bstr' XTEWreg2 `dep' `mis' `zvar', meth("`meth'") pan("`pan'") maxd(`i1') mis(`nx') optim(2)
		mat tmp1 = e(b_bs)
		mat tmp2 = e(se)
		
		local iv = 1

		forvalues j=1/`nx' {
			qui: mat results[`iv',`i'] = el("tmp1",1,`j')
			qui: mat sd_results[`iv',`i'] = el("tmp2",1,`j')
			local iv = `iv' + 1
		}

		forvalues j=1/`nz' {
			local j1 = `nx' + `j' + 1
			qui: mat results[`iv',`i'] = el("tmp1",1,`j1')
			qui: mat sd_results[`iv',`i'] = el("tmp2",1,`j1')
			local iv = `iv' + 1
		}

		qui: mat results[`iv',`i'] = `e(rho)'
		qui: mat sd_results[`iv',`i'] = `e(SErho)'
		local iv = `iv' + 1

		forvalues j=1/`nx' {
			local j1 = `nv' + `j' + 2
			qui: mat results[`iv',`i'] = el("tmp1",1,`j1')
			qui: mat sd_results[`iv',`i'] = el("tmp1",1,`j1')
			local iv = `iv' + 1
		}
	}
	
	mat list results
	mat list sd_results
	outtable2 using ${resfold}P2_`tst', /*sig("B")*/ mat(results) nobox replace f(%6.3f) format2(%5.3f) mat2(sd_results)
}




/*


// Example on how to do bootstrap with centered moment conditions
XTEWreg2 lever mtb tangib logsales oi , mis(2) maxd(4) meth(CML) pan(CLS) optim(2) cent(set)
bootstrap _b e(rho) el(e(tau),1,1) el(e(tau),2,1), rep(50) seed(1337) :  ///
	XTEWreg2 lever mtb tangib logsales oi , mis(2) maxd(4) meth(CML) pan(CLS) optim(2) cent(use)
bootstrap _b e(rho) el(e(tau),1,1) el(e(tau),2,1), rep(50) seed(1337) :  ///
	XTEWreg2 lever mtb tangib logsales oi , mis(2) maxd(4) meth(CML) pan(CLS) optim(2)

	
timer clear
timer on 1
XTEWreg2 lever mtb tangib logsales oi, meth(CML) pan(CLS) maxd(8) mis(3) optim(2)
timer off 1
//mata: EWSAVEDprb = NULL
timer on 2
XTEWreg2 lever mtb tangib logsales oi, meth(CML) pan(CLS) maxd(8) mis(3) optim(2)
timer off 2
timer list
*/

/*

foreach tp in MOM CML {
	forvalues maxdeg = 3/6 {
		XTEWreg2 ik q cfk, meth(`tp') pan(CLS) maxd(`maxdeg') mis(1) optim(2)
	}
}

foreach tp in MOM CML {
	forvalues maxdeg = 3/6 {
		XTEWreg2 lever mtb tangib logsales oi, meth(`tp') pan(CLS) maxd(`maxdeg') mis(2) optim(2)
	}
}

foreach tp in /*MOM*/ CML {
	forvalues maxdeg = 3/6 {
		XTEWreg2 lever mtb tangib logsales oi, meth(`tp') pan(CLS) maxd(`maxdeg') mis(3) optim(2)
	}
}

foreach tp in /*MOM*/ CML {
	forvalues maxdeg = 3/6 {
		XTEWreg2 lever mtb tangib logsales oi, meth(`tp') pan(CLS) maxd(`maxdeg') mis(4) optim(2)
	}
}



XTEWreg2 lever mtb tangib logsales oi, meth(CML) pan(CMD) maxd(5) mis(2) optim(2) //vce(ols)
//bootstrap _b e(rho) el(e(tau),1,1) el(e(tau),2,1), rep(100) seed(1337) : XTEWreg2 lever mtb tangib logsales oi, meth(CML) pan(CLS) maxd(3) mis(2) optim(2) //vce(ols)

qui {
/*
//XTEWreg1 lever mtb tangib logsales oi, meth(CML4) pan(CLS) two 
//disp `e(rho)' " " `e(SErho)' " " `e(tau1)' " " `e(SEtau1)' " " `e(tau2)' " " `e(SEtau2)' 
XTEWreg2 lever mtb tangib logsales oi, meth(CML) pan(CLS) maxd(4) mis(2) optim(2) //vce(ols)
local tau1 = el(e(tau),1,1)
local SEtau1 = el(e(SEtau),1,1)
local tau2 = el(e(tau),2,1)
local SEtau2 = el(e(SEtau),2,1)
disp `e(rho)' " " `e(SErho)' " " `tau1' " " `SEtau1' " " `tau2' " " `SEtau2' 
*/

//xtset, clear

/*
bootstrap _b, rep(20) seed(1337) /*cluster(gvkey)*/ : regress lever mtb
XTEWreg2 lever mtb tangib logsales oi, meth(CML) pan(CLS) maxd(4) mis(2) optim(2) //vce(ols)
bootstrap _b, rep(100) seed(1337) : XTEWreg2 lever mtb tangib logsales oi, meth(CML) pan(CLS) maxd(4) mis(2) optim(2) //vce(ols)
*/

//XTEWreg1 lever mtb tangib logsales oi, meth(MOM5) pan(CLS) two 
//XTEWreg3 lever mtb tangib logsales oi, meth(MOM) pan(CLS) maxd(4) mis(2) optim(2) //vce(ols)
}

//mata mosave Symbolic() , dir("${workfold}${scriptfold}") replace
//mata mosave Partitions() , dir("${workfold}${scriptfold}") replace


*/
