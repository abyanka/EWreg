*! version 1.0.0 class EWopt

/* 
Author: Robert Parham, University of Rochester

This work is licensed under the Creative Commons Attribution 4.0 International License. 
To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.
*/

version 11

mata:

class EWopt {

	public:
		real			scalar		n
		real			scalar		nx
		real			scalar		maxdeg
		real			scalar		nz
		real			scalar		nt
		real			scalar		nt2
		real			scalar		neq
		real			scalar		neq2
		real			scalar		N1neq
		real			scalar		ncml
		real			scalar		nk

		real			scalar 		met
		real			scalar 		clustmet
		real			scalar 		vcmet
		real			scalar 		optmet
		real			scalar 		centmet

		real			scalar		nbXint

		string			scalar		fname
		
		real 			scalar		CML_maxiter
		real 			scalar		GN_maxiter
		real 			scalar		GN_maxsqz
		real 			scalar		tol

		void 						einit()
		void 						setprb()
		real 			scalar 		doMOM()
}



// initialize all options
void EWopt::einit(real scalar n, real scalar nx, real scalar maxdeg, real scalar nz, 		///
				real scalar met, real scalar clustmet, real scalar vcmet, 					///
				real scalar optmet, real scalar centmet)
{
	this.n 			= n
	this.nx 		= nx
	this.maxdeg 	= maxdeg
	this.nz			= nz
	this.met		= met
	this.clustmet	= clustmet
	this.vcmet		= vcmet
	this.optmet		= optmet
	this.fname 		= "EWcache_" + strofreal(nx) + "_" + strofreal(maxdeg) + ".data"
	this.CML_maxiter = 500
	this.GN_maxiter = 2999
	this.GN_maxsqz	= 240
	this.tol		= 1e-9
	this.centmet	= centmet
	
}
// end einit


// sets model dependant options
void EWopt::setprb(class EWproblem scalar prb)
{
	neq 	= rows(prb.Drhs)
	nt 		= cols(prb.Drhs)
	neq2 	= rows(prb.mom2)
	nt2 	= prb.mom2[1,1].getN()
	ncml	= rows(prb.cml)
	N1neq	= rows(prb.N1lhs)
	nk		= cols(prb.yidx)
}
// end setprb


// is the method to use MOM?
real scalar EWopt::doMOM()
{
	assert (n!=.)
	
	return (met==1)
}
// end doMOM

end
