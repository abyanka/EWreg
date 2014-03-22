/* General initial setup */
clear all

global workfold = "C:\Users\finguy\SkyDrive\Documents\Toni\\"

cd ${workfold}
adopath + "${workfold}"

use "EJW1.dta", clear

qui: xtset gvkey

XTEWreg lever mtb tangib logsales oi , maxd(5) mis(2)

//XTEWreg ik q cfk , maxd(5) cent(set)
//bootstrap _b e(rho) el(e(tau),1,1), rep(100) seed(1337) cluster(gvkey) : ///
//	XTEWreg ik q cfk , maxd(5) cent(use)
