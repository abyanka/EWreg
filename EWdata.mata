/* 
Author: Robert Parham, University of Rochester

This work is licensed under the Creative Commons Attribution 4.0 International License. 
To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.
*/

version 12

mata:

// Create and hold all data pertaining to the EW problem
class EWdata {
	private:
		real 			colvector	find_beta()

	public:
		real			matrix		yx		//[n,nx+1]
		real			matrix		z		//[n,nz]
		real			matrix		mom		//[n,neq]
		real			matrix		emom	//[neq,1]
		real			matrix		zmom	//[neq,nz]
		real			matrix		N1mom	//[n,N1neq]
		real			matrix		N1emom	//[N1neq,1]
		real			matrix		N1zmom	//[N1neq,nz]
		real 			colvector	muy		//[nz,1], projection of y on z, y_i=z_i*muy
		real 			matrix		mux		//[nz,nx], projection of x on z, x_i=z_i*mux
		real			matrix		inEzz	//[nz,nz], cholinv((z'*z)/n)
		real			matrix		zyx		//[n,nz*(nx+1)], yx(_d) dot-multiplied by z's columns 

		real			rowvector	ni		//[1,??], number of clustered obs

		real			colvector	fCent	//[neq/nk,1], moments for centered bootstrap

		real 			colvector	ky		//[nk,1]
		real 			matrix		kx		//[nk,nx]
		real 			matrix		Dky		//[nk,neq]
		real 			matrix		Dkx		//[nk,nx*neq]
		
		real 			rowvector	bXint	//[1,c*nx]
		real 			colvector	bX0_Gry	//[nx,1]
		real 			colvector	bX0_OLS	//[nx,1]
		
		real 			colvector	bX		//[nx,1]
		real 			colvector	bZ		//[nz,1]
		real 			scalar		rho		//[1,1]
		real 			colvector	tau		//[nx,1]
		
		void 						einit()
		real			colvector	get_beta()
}
// end EWdata


// find starting beta based on Eq26 of EW2002
real colvector EWdata::find_beta(class EWopt scalar opt, class EWproblem scalar prb)
{
	nr = comb(opt.nx+1,2)
	V = J(nr,1,0)
	W1 = J(opt.nx,opt.nx,0)
	W2 = J(nr-opt.nx,opt.nx,0)

	// fill V
	pos = sum(1..opt.nx)+opt.nx+1 + 1*(opt.maxdeg>=5) + 1
	wpos = pos
	V[1..opt.nx,1] = emom[(pos)..(pos+opt.nx-1),1]
	vpos = opt.nx + 1
	pos = pos + opt.nx

	for (i=opt.nx-1;i>=1;i--) {
		pos++
		V[(vpos)..(vpos+i-1),1] = emom[(pos)..(pos+i-1),1]
		vpos = vpos+i
		pos = pos+i
	}
	
	// fill W1
	pos = wpos + opt.nx
	for (i=1;i<=opt.nx;i++) {
		for (j=i;j<=opt.nx;j++) {
			W1[i,j] = emom[pos,1]			// E(y*x_i*x_j)
			W1[j,i] = W1[i,j] 				// use symmetry of W1 to save a few steps
			pos++
		}
	}
	
	// fill W2
	if (opt.nx>1) {
		W2 = prb.getW2(opt.nx, nr-opt.nx, opt.neq)
		for (i=1;i<=nr-opt.nx;i++) {
			W2[i,.] = emom[W2[i,.],1]'
		}
		W = W1\W2
	}
	else {
		W = W1
	}

	// verify identification, find beta and return it
	pW = pinv(W,rnk)		//this should be pinv, not cholinv!
	if (rnk<opt.nx) {
		errprintf("Third moment matrix not full-rank. System not identified. Aborting.\n")
		exit(error(459))
	}
	beta = pW*V

	return (beta)
}
// end find_beta


// initialize data
void EWdata::einit(class EWopt scalar opt, class EWproblem scalar prb, transmorphic colvector id,	///
					real colvector y, real matrix x, real matrix z1, real rowvector bXint)
{
	// save z
	z = z1
	
	// partial out z from x and y
	invzz = cholinv(quadcross(z,z))
	mux = invzz*quadcross(z,x)
	muy = invzz*quadcross(z,y)
	inEzz = opt.n*invzz
	yx = (y-z*muy,x-z*mux)
	
	// generate ezmom
	mom = J(opt.n, opt.neq, 0)
	emom = J(opt.neq, 1, 0)
	zmom = J(opt.neq, opt.nz, 0)
	
	for (i=1;i<=opt.neq;i++) {
		mom[.,i] = prb.lhs[i,1].evalmat(yx)
	}
	
	emom = mean(mom)'
	
	for (i=1;i<=opt.nz;i++) {
		zmom[.,i] = mean(mom:*(J(1,opt.neq,1)#z[.,i]))'
	}

	// generate N1ezmom
	N1mom = J(opt.n, opt.N1neq, 0)
	N1emom = J(opt.N1neq, 1, 0)
	N1zmom = J(opt.N1neq, opt.nz, 0)
	
	for (i=1;i<=opt.N1neq;i++) {
		N1mom[.,i] = prb.N1lhs[i,1].evalmat(yx)
	}
	
	N1emom = mean(N1mom)'
	
	for (i=1;i<=opt.nz;i++) {
		N1zmom[.,i] = mean(N1mom:*(J(1,opt.N1neq,1)#z[.,i]))'
	}

	// generate zyx
	zyx = J(opt.n,opt.nz*(opt.nx+1),0)
	for (i=1;i<=(opt.nx+1);i++) {
		cur_beg = (i-1)*opt.nz+1
		cur_end = i*opt.nz
		zyx[.,cur_beg..cur_end] = (J(1,opt.nz,1)#yx[.,i]):*z
	}
	
	// build the ni for clustering (as per section 2.7 of EJW)
	uid = uniqrows(id)
	ni = J(1,rows(uid),0)
	j=1
	for (i=1;i<=opt.n;i++)
	{
		if (id[i,1]==uid[j,1]) continue
		ni[1,j]=i-1
		j=j+1
	}
	ni[1,j]=opt.n
	ni = (0,ni)
	assert (j==rows(uid))
	
	// count the implied bXint iterations, and update opt
	this.bXint = bXint
	if (cols(bXint)==1 & bXint[1,1] == 0) {
		this.bXint = J(1,opt.nx,0) , J(1,opt.nx,0)
	}
	opt.nbXint = (cols(this.bXint)/opt.nx)
	assert(opt.nbXint==round(opt.nbXint))
	
	// verify and put fCent
	external real colvector EWSAVEDfCent
	if (opt.centmet == 0 | opt.centmet == 1 | EWSAVEDfCent==NULL) {
		EWSAVEDfCent = J(0,1,0)
		fCent = EWSAVEDfCent
	}
	else { //opt.centmet == 2, i.e. Use
		valid_fCent = opt.doMOM()*(rows(EWSAVEDfCent)==opt.neq) + ///
					  (!opt.doMOM())*(rows(EWSAVEDfCent)==opt.nk)
		
		if (valid_fCent) 	fCent = EWSAVEDfCent
		else				fCent = J(0,1,0)
	} // if fCent is set, then it is ready to be used.
	
	// do cumulants
	cml = J(opt.neq, 1, 0)
	for (i=1;i<=opt.neq;i++) {
		cml[i,1] = prb.cml[i,1].eval(emom)
	}
	// split cml into ky and kxj
	ky = cml[prb.yidx,1]
	kx = J(opt.nk,0,0)
	for (i=1;i<=opt.nx;i++) {
		kx = kx, cml[prb.xidx[i,.],1]
	}
	
	// do Dcumulants
	Dcml = J(opt.neq, opt.neq, 0)
	for (i=1;i<=opt.neq;i++) {
		for (j=1;j<=opt.neq;j++) {
			Dcml[i,j] = prb.Dcml[i,j].eval(emom)
		}
	}
	// split Dcml into Dky and Dkxj
	Dky = Dcml[prb.yidx,.]
	Dkx = J(opt.nk,0,0)
	for (i=1;i<=opt.nx;i++) {
		Dkx = Dkx, Dcml[prb.xidx[i,.],.]
	}

	// create bX0_Gry - the Eq26 least-squares guess for beta
	bX0_Gry = find_beta(opt, prb)
	
	// create bX0_OLS - the simple ols estimate
	bX0_OLS = qrsolve((x,z), y)
	bX0_OLS = bX0_OLS[1..opt.nx,1]
	
	// initialize bX, bZ and tau, for good measure
	bX = 	J(opt.nx,1,0)
	bZ = 	J(opt.nz,1,0)
	tau = 	J(opt.nx,1,0)

}
// end einit


// return i-th beta to attempt (possibly generating a candidate)
real colvector EWdata::get_beta(numeric scalar i, class EWopt scalar opt)
{
	retval = bXint[1,((i-1)*opt.nx+1)..(i*opt.nx)]
	iszero = sum(retval:!=0)==0
	isdbl  = 0
	if (i>1) {
		retval1 = bXint[1,((i-2)*opt.nx+1)..((i-1)*opt.nx)]
		isdbl = sum(retval1:!=0)==0
	}
	if (iszero) {
		retval = bX0_Gry
	}
	if (isdbl) {
		retval = bX0_OLS
	}
	
	return (retval)
}
// end get_beta

end
