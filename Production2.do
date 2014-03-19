// General initial setup
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
}



// Compile an mlib
qui {
// /*
mata: mata clear

do "${scriptfold}Symbolic.mata"
do "${scriptfold}Partitions.mata"
do "${scriptfold}EWproblem.mata"
do "${scriptfold}EWopt.mata"
do "${scriptfold}EWdata.mata"
do "${scriptfold}EWreg.mata"

mata: 
mata mlib create lEW , dir("${workfold}${scriptfold}") replace
mata mlib add lEW Symbolic()
mata mlib add lEW Partitions()
mata mlib add lEW EWproblem()
mata mlib add lEW EWopt()
mata mlib add lEW EWdata()
mata mlib add lEW EWreg()
mata mlib add lEW i_crit()
mata mlib add lEW stats()
mata mlib add lEW mt()

mata mlib index
end
// */
}



// FE and time limits
qui {
/*
insheet using "${datfold}EJWdata.csv" , clear
ren v1 gvkey
ren v2 fyear
ren v3 ik
ren v4 q
ren v5 cfk
ren v6 lever
ren v7 mtb
ren v8 tangib
ren v9 logsales
ren v10 oi
drop v11

keep if fyear>=1970 & fyear<=1970+42-1

bysort gvkey: egen tot = count(gvkey)
drop if tot==1
drop tot
foreach vr in ik q cfk lever mtb tangib logsales oi {
	by gvkey: egen m`vr' = mean(`vr')
	replace `vr' = `vr'-m`vr' //if tot>1
	drop m`vr'
}
foreach vr in ik q cfk lever mtb tangib logsales oi {
	regress `vr' i.fyear
	predict tmp , resid
	replace `vr' = tmp
	drop tmp
}

// support Gauss code as well
gen intercep = 1

saveold "${datfold}EJW1.dta", replace
*/
}



// Start work
{
use "${datfold}EJW1.dta", clear

qui: xtset gvkey

//XTEWreg ik q cfk , maxd(3) mis(1)
//XTEWreg ik q cfk , maxd(4) mis(1)
//XTEWreg ik q cfk , maxd(5) mis(1)
//XTEWreg lever mtb tangib logsales oi , maxd(3) mis(2) meth(MOM)
//XTEWreg lever mtb tangib logsales oi , maxd(4) mis(2) meth(MOM)
XTEWreg lever mtb tangib logsales oi , maxd(5) mis(2) meth(MOM)

/*
XTEWreg ik q cfk , maxd(5) cent(set)
bootstrap _b e(rho) el(e(tau),1,1), rep(100) seed(1337) cluster(gvkey) : ///
	XTEWreg ik q cfk , maxd(5) cent(use)
*/
}
