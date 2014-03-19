/* 
Author: Robert Parham, University of Rochester

This work is licensed under the Creative Commons Attribution 4.0 International License. 
To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.
*/

version 12

mata:

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


// a struct to allow vectors of matrices
struct mt
{
	transmorphic matrix inmat
}
// end struct mt


// an EW regression class
class EWreg {
	private:
		pointer(class EWproblem) scalar 		pprb

		class EWdata 	scalar 		dta
		class EWopt 	scalar 		opt
		
		real 			matrix 		ff				//[n,neq]
		real 			matrix 		omega
		real 			matrix 		w
		real 			colvector 	t
		real			scalar		obj

		real 			colvector 	fsave
		
		real			matrix		optw()
		real			matrix		getomega()
		real 			colvector 	dogmm()
		struct stats	scalar		getret()
		real 			colvector	EWgmm1()
		real 			colvector	EWgmm2()
		real 			colvector	sqeez()
		real 			matrix		dogeeRT()
		
		real 			scalar		mutual_cnt()
		real			matrix		CMDshare()
		real 			colvector	findpos()
		real 			rowvector 	cmd()
		
	public:
		struct stats	scalar		doCLS()
		struct stats	scalar		doCMD()

		real 			colvector	deff()
		real 			matrix		gradMOM()
		pointer(class EWproblem) scalar 		get_pprb()
		
}
// end EWreg


// a public interface for pprb, for i_crit. A hack, but we'll survive.
pointer(class EWproblem) scalar EWreg::get_pprb()
{
	return (pprb)
}
// end get_pprb


// define the f = gi(mu_hat)-gbar(mu_hat)+Gbar(mu_hat)*KSImui vector for the optimal weight matrix
real matrix EWreg::optw()
{
	nemom = (J(opt.n,1,1)#dta.emom')
	retval = dta.mom-nemom
	
	// do standard error adjustment
	if (opt.vcmet==0) {
	
		ei=J(opt.n,opt.neq,0)

		// define G(mu) - maybe easier with Symbolic, but i kinda gave up here.
		G = J(opt.neq,opt.nz*(opt.nx+1),0)
		for (eq=1;eq<=opt.neq;eq++) {
			for (j=1;j<=opt.nx+1;j++) {
				d = pprb->rmat[eq,j]

				if (d>0) {
					mr = pprb->rmat[eq,.]
					mr[1,j] = mr[1,j]-1
					m = sum((colsum(pprb->rmat':==mr'):==opt.nx+1):*(1..opt.neq))
					if (m==0) {
						m = sum((colsum(pprb->N1rmat':==mr'):==opt.nx+1):*(1..opt.N1neq))
						if (m==0) continue

						G[eq,((j-1)*opt.nz+1)..(j*opt.nz)] = -1:*d:*dta.N1zmom[m,.]
						continue
					}
					G[eq,((j-1)*opt.nz+1)..(j*opt.nz)] = -1:*d:*dta.zmom[m,.]
				}
			}
		}
		
		// multiply G, Q^-1 (inEee), Rmu (zyx)
		ei = (G*(I(opt.nx+1)#dta.inEzz)*dta.zyx')'
		
		// add correction into retval
		retval = retval + ei
	}
	
	// do numerical adjustment - reweigthing f, so that optw will be well-conditioned
	if (opt.doMOM()) retval = retval:/nemom
	
	return (retval)	
}
// end optw


// ff->omega, take care of clustering
real matrix EWreg::getomega(real matrix ff1)
{
	// 0: clustered
	// 1: identity
	// 2: regular ((ff'*ff)/n)

	if (opt.clustmet==1) {
		return (I(cols(ff1)))
	}
	if (opt.clustmet==2) {
		return ((ff1'*ff1)/rows(ff1))
	}
	if (opt.clustmet==0) {
		retval = J(cols(ff1),cols(ff1),0)
	}
	else {
		errprintf("getomega recieved unexpected clustmet value. Please submit bug report.\n")
		exit(error(198))
	}
		
	for (i=2;i<=cols(dta.ni);i++)
	{
		hi = colsum(ff1[dta.ni[1,i-1]+1..dta.ni[1,i],.])
        retval = retval + hi'*hi
    }
    retval = retval:/rows(ff1)
	
	return (retval)
}
// end getomega


// compute gradient matrix of moment problem
real matrix EWreg::gradMOM(class Symbolic matrix systm, real colvector t, real scalar wgt)
{
	dvdc = J(rows(systm), cols(systm), 0)
	for (eq=1;eq<=rows(systm);eq++) {
		for (it=1;it<=cols(systm);it++) {
			dvdc[eq,it] = systm[eq,it].eval(t)
		}
	}
	
	if (wgt) dvdc = dvdc:/(J(1,opt.nt,1)#dta.emom)
	
	return (-1:*dvdc)
}
// end gradMOM


// define the f vector (distance between data moments and constructed moments)
real colvector EWreg::deff(real colvector t)
{
	// f'wf is the GMM objective function to minimize to find t
	// for an example defining f, see page 783 of EW2002
	
	f=J(opt.neq, 1, 0)
	
	for (eq=1;eq<=opt.neq;eq++) {
		f[eq,1] = (dta.emom[eq,1] - pprb->rhs[eq,1].eval(t))
	}
	
	if (opt.doMOM()) f = f:/dta.emom
	
	if (rows(dta.fCent)>0) f = f - dta.fCent

	return (f)
}
// end deff


// criterion function for EWgmm1. Provide the obj func and the gradient thereof
void i_crit(real scalar todo, real rowvector b, real scalar crit, real rowvector g, real matrix H)
{
	class EWreg scalar EWobj
	
	p = findexternal("EWexternal")
	EWobj  	= *((*p)[1])
	w 		= *((*p)[2])

	f=EWobj.deff(b')
	crit = f'*w*f
	
	if (todo==1) {
		gr = EWobj.gradMOM(EWobj.get_pprb()->Drhs, b', 1)
		gwg = gr'*w*gr
		gwf = gr'*w*f
		g = (cholinv(gwg)*gwf)'
	}
}
// end function i_crit


// compute squeezes
real colvector EWreg::sqeez(real colvector s_t, real colvector s_dt, real scalar obj1, real matrix w)
{
	// compare the values of the objective function at the
	// points c+s*dc and c+0.5*s*dc, with dc = the proposed change in the vector
	// of parameters, and with step length s initially at 1. s is halved until
	// minus the objective function stops declining.
	
	s_t1=s_t-s_dt
	s_lm=1/2
	s_itr=1

	s_f1=deff(s_t1)
	lc1 = s_f1'*w*s_f1

	while (s_itr<=opt.GN_maxsqz)
	{
		s_t2=s_t-s_lm*s_dt
		s_f2=deff(s_t2)
		lc2 = s_f2'*w*s_f2

		if (lc1 <= lc2 && lc1 <= obj1) {
			return (s_t1)
		}
		else {
			s_t1=s_t2
			s_lm=s_lm/2
			lc1=lc2
			s_itr=s_itr+1
		}
	}
	
	return(s_t2)
}
// end sqeez


// do gauss-newton
real colvector EWreg::EWgmm2(real matrix w, real colvector t)
{
	dt=1	// Initialize the step length.
	
	for (iter=1; iter<=opt.GN_maxiter & norm(dt,.)>=opt.tol ; iter++) 
	{
		// find current objective
		f = deff(t)
		g = gradMOM(pprb->Drhs, t, 1)
		obj1 = f'*w*f
		
		// use the GAUSS-NEWTON method to compute optimal full step dt
		gwg= g'*w*g
		gwf= g'*w*f
		dt= cholinv(gwg)*gwf
		
		if (opt.GN_maxsqz > 0) 	t_new=sqeez(t,dt,obj1,w)
		else 					t_new=t-dt

		dt=t_new-t	// Update variables for the next iteration, also better for numerical stability
		t=t_new
	}
	
	if (iter>opt.GN_maxiter)
	{
		printf("Gauss-Newton stopped after maximum iterations. \n")
		return (0\t)
	}

	return (1\t)
}
// end EWgmm2


// do mata internal optimization (BFGS)
real colvector EWreg::EWgmm1(real matrix w, real colvector t)
{
	pointer(pointer() vector) scalar p 
	p = crexternal("EWexternal")
	(*p) = (&this, &w)
	
	S=optimize_init()
	optimize_init_evaluator(S, &i_crit())
	optimize_init_which(S,"min")
	optimize_init_tracelevel(S,"none")
	optimize_init_evaluatortype(S,"d1")
	optimize_init_technique(S, "bfgs")
	optimize_init_conv_warning(S, "off")
	optimize_init_params(S,t')
	if (err=_optimize(S)) {
		printf("%s\n", optimize_result_errortext(S))
		printf("Using value from Gauss-Newton. Suspect problem is ill-condtioned. \n")
		t = J(1,0,0)
	}
	else {
		t=optimize_result_params(S)
	}

	rmexternal("EWexternal")	
	return(t')
}
// end EWgmm1


// do optimization based GMM computation
real colvector EWreg::dogmm(real matrix w, real colvector t)
{
	if (opt.optmet == 1) {
		//printf("Using BFGS optimizer. \n")
		t1 = EWgmm1(w,t)
		if (rows(t1)!=0) return (t1)
		// if here it means BFGS failed
		printf("BFGS failed. Attempting Gauss-Newton optimizer. \n")
	}

	// if here it means BFGS failed or optmet==2
	//printf("Using Gauss-Newton optimizer. \n")
	t1 = EWgmm2(w,t)
	
	if (opt.optmet == 2 & t1[1,1]==0) {
		//GN was first, we ran it, and it stopped after maxiter
		printf("Gauss-Newton failed. Attempting BFGS optimizer. \n")
		t2 = EWgmm1(w,t1[2..rows(t1),1])
		if (rows(t2)==0) {
			printf("BFGS failed to find a better solution than Gauss-Newton.")
			return (t1[2..rows(t1),1])
		}
		return (t2)
	}
	
	return (t1[2..rows(t1),1])
}
// end dogmm


// calculate the gradient matrix going with bigphiRT
real matrix EWreg::dogeeRT(real matrix sigz, real scalar denomy, real scalar numery, 			///
					real colvector denomx, real colvector numerx, real matrix sigeta, 			///
					numeric rowvector eta2p)
{

	//			muyx				sigz				beta		etaj				eps	  u
	nrows = (opt.nx+1)*opt.nz + opt.nz*(opt.nz+1)/2 + opt.nx + opt.nx*(opt.nx+1)/2 + opt.nx + 1

	// first column rho2, next nx columns tau2j
	gee = J(nrows,(1+opt.nx),0)
	
	/*** start with rho2 column (column 1) ***/
	
	// derivatives wrt muy
	counter=1
	for (i=1;i<=opt.nz;i++) {
		gee[counter,1]= (2*dta.muy'*sigz[,i])/denomy - numery*(2*dta.muy'*sigz[,i])/(denomy^2)
		counter++
	}

	// derivatives wrt sigz
	counter=(opt.nx+1)*opt.nz + 1
	for (i=1;i<=opt.nz;i++) {
		for (j=i;j<=opt.nz;j++) {
			x = ((i!=j)*1+1)*dta.muy[i,1]*dta.muy[j,1]
			gee[counter,1] = x/denomy - numery*x/(denomy^2)
			counter++
		}
	}
	
	// derivatives wrt the members of t we care about
	// first, betaj
	counter=(opt.nx+1)*opt.nz + opt.nz*(opt.nz+1)/2 + 1
	for (i=1;i<=opt.nx;i++) {
		gee[counter,1] = (2*dta.bX'*sigeta[.,i])/denomy - numery*(2*dta.bX'*sigeta[,i])/(denomy^2)
		counter++
	}

	// Now E(etaj1*etaj2) for all combs
	for (i=1;i<=opt.nx;i++) {
		for (j=i;j<=opt.nx;j++) {
			x = ((i!=j)*1+1)*dta.bX[i,1]*dta.bX[j,1]
			gee[counter,1] = x/denomy - numery*x/(denomy^2)
			counter++
		}
	}
	
	// E(u^2) is lonely
	//			muyx				sigz				 beta		etaj				eps
	counter = (opt.nx+1)*opt.nz + opt.nz*(opt.nz+1)/2 + opt.nx + opt.nx*(opt.nx+1)/2 + opt.nx + 1
	gee[counter,1] = -numery/(denomy^2)
	
	
	/*** now for the tau columns ***/

	for (jt=1;jt<=opt.nx;jt++) {
		// derivatives wrt mux
		counter = opt.nz*jt + 1
		for (i=1;i<=opt.nz;i++) {
			gee[counter,1+jt]= (2*dta.mux[.,jt]'*sigz[,i])/denomx[jt,1] - numerx[jt,1]*(2*dta.mux[.,jt]'*sigz[,i])/(denomx[jt,1]^2)
			counter++
		}

		// derivatives wrt sigz
		counter=(opt.nx+1)*opt.nz + 1
		for (i=1;i<=opt.nz;i++) {
			for (j=i;j<=opt.nz;j++) {
				x = ((i!=j)*1+1)*dta.mux[i,jt]*dta.mux[j,jt]
				gee[counter,1+jt] = x/denomx[jt,1] - numerx[jt,1]*x/(denomx[jt,1]^2)
				counter++
			}
		}
		
		// derivatives wrt the members of t we care about
		// E(eta^2)j
		counter = (opt.nx+1)*opt.nz + opt.nz*(opt.nz+1)/2 + opt.nx + eta2p[1,jt]
		gee[counter,1+jt] = 1/denomx[jt,1] - numerx[jt,1]/(denomx[jt,1]^2)
		// E(eps^2)j
		counter = (opt.nx+1)*opt.nz + opt.nz*(opt.nz+1)/2 + opt.nx + opt.nx*(opt.nx+1)/2 + jt
		gee[counter,1+jt] = - numerx[jt,1]/(denomx[jt,1]^2)
	}
	
	return (gee)
}
// end dogeeRT


// do post-estimation and create struct stats return value
struct stats scalar EWreg::getret()
{
	struct stats scalar retval
	retval = stats()

	
	/*** The basic stuff ***/
	
	retval.N			= opt.n
	retval.w			= w
	retval.IDval  		= J(1,1,-1)
	retval.IDstat 		= J(1,1,-1)

	
	/*** Save centered moments for next time ***/

	external real colvector EWSAVEDfCent
	if (opt.centmet==1)		EWSAVEDfCent = fsave
	
	
	/*** Sargan-Hansen J-stats ***/
	
	retval.obj		= obj
	retval.Jstat 	= obj*opt.n
	if (opt.doMOM())	retval.dfree = opt.neq-opt.nt
	else 				retval.dfree = opt.nk-opt.nx
	retval.Jval  	= 1-chi2(retval.dfree,retval.Jstat)

	
	/*** inflncX + inflncT ***/
	
	if (opt.doMOM()) {
		g 		= gradMOM(pprb->Drhs, t, 1)			//[neq,nt]
		nvcX 	= cholinv(g'*w*g)
		vcX 	= nvcX/opt.n
		inflncT = (nvcX*g'*w*ff')'  				//[n,nt], NOTE - Minus removed!
		inflncX = inflncT[.,1..opt.nx]				//[n,nx]
		inflncT = inflncT[.,(opt.nx+1)..opt.nt]		//[n,nt-nx]
	}
	else {
		D 		= dta.Dky - dta.Dkx*(dta.bX#I(opt.neq))			//[nk,neq]
		nvcX 	= cholinv(dta.kx'*w*dta.kx)
		vcX 	= nvcX/opt.n
		inflncX = (nvcX*dta.kx'*w*D*ff')'						//[n,nx]
		t		= pprb->resolve_eq(pprb->mom2, dta.bX, opt.neq2, opt.nx, opt.nt2, dta.emom)	//[nt2,1]
		g		= gradMOM(pprb->Dmom2, t, 0)					//[neq2,nt2]
		tmp		= I(opt.nt2)
		g 		= tmp[1..opt.nx,.]\g							//[neq2+nx,nt2] == [nt2,nt2]
		ff2		= (inflncX, ff[.,1..opt.neq2])					//[n,nt2]
		w2		= cholinv(getomega(ff2))						//[nt2,nt2]
		inflncT = (cholinv(g'*w2*g)*g'*w2*ff2')' 				//[n,nt2], NOTE - Minus removed!
		inflncT = inflncT[.,(opt.nx+1)..opt.nt2]				//[n,nt2-nx] == [n,neq2]
	}
	

	// redefine a local version of nt and ntvec to fit both MOM and CML
	if (opt.doMOM()) {
		nt = opt.nt
		ntvec = pprb->ntvec
	}
	else {
		nt = opt.nt2
		ntvec = pprb->ntvec2
	}
	
	
	/*** bZ ***/
	
	dta.bZ = dta.muy-dta.mux*dta.bX

	
	/*** inflncZ ***/
	
	geeZ = J((opt.nz*(opt.nx+1)+opt.nx), opt.nz, 0)
	geeZ[1..opt.nz,.] = -I(opt.nz)
	for (j=1;j<=opt.nx;j++) {
		geeZ[(opt.nz*j+1)..(opt.nz*(j+1)),.] = dta.bX[j,1]*I(opt.nz)
		geeZ[(opt.nz*(opt.nx+1)+j),.] = dta.mux[.,j]'
	}

	phimuyx = ((I(opt.nx+1)#dta.inEzz)*dta.zyx')'
	bigphiZ = phimuyx, inflncX
	avarZ = getomega(bigphiZ):/opt.n
	vcZ = geeZ'*avarZ*geeZ
	seZ = sqrt(diagonal(vcZ))
	inflncZ = bigphiZ*geeZ


	/*** rho and tau ***/

	ez = mean(dta.z)
	sigz = crossdev(dta.z, ez, dta.z, ez) :/ opt.n

	// build sigeta symmetric matrix [nx,nx]
	sigeta = J(opt.nx,opt.nx,0)
	eta2 = J(1,opt.nx,0)
	curr = ntvec[1,1] + 1
	for (i=1;i<=opt.nx;i++) {
		// also grab E(etaji^2] rowvector [1,nx] of addresses in t - only etaj*etaj
		eta2[1,i] = curr
		for (j=i;j<=opt.nx;j++) {
			sigeta[i,j] = t[curr]
			sigeta[j,i] = sigeta[i,j]
			curr++
		}
	}

	// grab E(ui^2] address in t
	u2 = sum(ntvec[1..3,1]) + 1
	
	// grab E(epsji^2] rowvector [1,nx] of addresses in t
	base_eps = sum(ntvec[1..2,1])
	eps2 = ((base_eps+1)..(base_eps+opt.nx))

	// grab E(etaji^2] rowvector [1,nx*(nx+1)/2] of addresses in t - all etaj1*etaj2 combs
	eta2_all=((ntvec[1,1] + 1)..(ntvec[1,1] + opt.nx*(opt.nx+1)/2))

	// do rho
	numery=dta.muy'*sigz*dta.muy + dta.bX'*sigeta*dta.bX
	denomy=(numery+t[u2,1])
	dta.rho=numery/denomy

	// do tau
	numerx=diagonal(dta.mux'*sigz*dta.mux) :+ t[eta2,1]
	denomx=(numerx:+t[eps2,1])
	dta.tau=numerx:/denomx

	
	/*** SE rho and tau ***/

	vecsigz=J(opt.nz*(opt.nz+1)/2,1,0)
	phiz=J(opt.n,opt.nz*(opt.nz+1)/2,0)

	counter=1
	for (i=1;i<=opt.nz;i++) {
		for (j=i;j<=opt.nz;j++) {
			vecsigz[counter,1] = sigz[i,j]
			phiz[.,counter]=(dta.z[,i]-J(opt.n,1,1)#ez[1,i]):*(dta.z[,j]-J(opt.n,1,1)#ez[1,j])
			counter++
		}
	}
	phiz=phiz:-J(opt.n,1,1)#vecsigz'
	
	// make the influence functions for rhotau
	bigphiRT = (phimuyx, phiz, inflncX, inflncT[.,eta2_all:-opt.nx], inflncT[.,eps2:-opt.nx], inflncT[.,u2:-opt.nx])
	avarRT = getomega(bigphiRT):/opt.n
	
	
	// build geeRT
	eta2p = eta2 :- ntvec[1,1]
	geeRT = dogeeRT(sigz, denomy, numery, denomx, numerx, sigeta, eta2p)

	
	/*** finish retval ***/
	
	retval.inflncXZ 	= inflncX, inflncZ
	retval.inflncRT 	= bigphiRT*geeRT
	retval.vcrhotau		= geeRT'*avarRT*geeRT
	retval.SErho		= sqrt(retval.vcrhotau[1,1])
	retval.SEtau		= sqrt(diagonal(retval.vcrhotau[2..(opt.nx+1),2..(opt.nx+1)]))
	retval.rho 			= dta.rho
	retval.tau 			= dta.tau
	retval.beta 		= dta.bX\dta.bZ
	retval.VCmat 		= getomega(retval.inflncXZ):/opt.n
	retval.serr 		= sqrt(diagonal(retval.VCmat))
	
	return (retval)
}
// end getret


// do J mismeasured regressors using moments or cumulants :)
struct stats scalar EWreg::doCLS(transmorphic colvector id, real colvector y, real matrix x, 	///
							real matrix z, numeric scalar met, numeric scalar clustmet, 		///
							numeric scalar maxdeg, real rowvector bXint, numeric scalar vcemet, ///
							numeric scalar optmet, real scalar centmet)
{
	opt = EWopt()
	opt.einit(rows(id), cols(x), maxdeg, cols(z), met, clustmet, vcemet, optmet, centmet)

	external class EWproblem scalar EWSAVEDprb
	if (EWSAVEDprb!=NULL) {
		if (cols(EWSAVEDprb.rocmat) != maxdeg | rows(EWSAVEDprb.xidx) != cols(x)) {
			printf("Problem structure different from last executed. Rebuilding problem.\n")
			EWSAVEDprb = EWproblem()
			EWSAVEDprb.einit(opt.maxdeg, opt.nx)
		}
	}
	else {
		EWSAVEDprb = EWproblem()
		EWSAVEDprb.einit(opt.maxdeg, opt.nx)
	}
	pprb = &EWSAVEDprb
	opt.setprb(*pprb)
	
	dta = EWdata()
	dta.einit(opt, *pprb, id, y ,x, z, bXint)
	
	// do identification tests
	// RPG - skip id tests for now

	
	// calculate optimal inverse weighting matrix, omega
	ff = optw()
	omega = getomega(ff)
	
	if (opt.doMOM()) {
		// do moments
		objsave = J(opt.nbXint,1,0)
		t = J(opt.nt,1,0)
		
		for (rep=1;rep<=opt.nbXint+1;rep++) {
			if (rep<=opt.nbXint) {
				bXinit = dta.get_beta(rep,opt)
			}
			else {
				if (opt.nbXint>1) {
					// restore best inital value
					minindex(objsave, 1, ind, where)
					bXinit = dta.get_beta(ind[1,1],opt)
				}
				else {
					// no point in restoring - only one value tested
					break
				}
			}
			
			t = pprb->resolve_eq(pprb->rhs, bXinit, opt.neq, opt.nx, opt.nt, dta.emom)
			w = cholinv(omega)
			t = dogmm(w, t)
			
			// save objective value for usage in final loop and reporting
			f = deff(t)
			obj = f'*w*f
			if (rep<=opt.nbXint) objsave[rep,1]=obj
		}
		dta.bX = t[1..opt.nx,1]
		fsave = f
	}
	else {
		// do cumulants
		
		bX = J(opt.nx,1,0)
		diffW = opt.tol + 1
		
		if (rows(dta.fCent)>0) dta.ky = dta.ky - dta.fCent

		for (iter = 0; iter<=opt.CML_maxiter & diffW > opt.tol ; iter++)
		{
			// D as defined in prop 1 of EJW2013
			D = dta.Dky - dta.Dkx*(bX#I(opt.neq))
			w = cholinv(D*omega*D')					//[nk,nk]
			iKWK = cholinv(dta.kx'*w*dta.kx)		//[nx,nx]
			bX1 = iKWK*dta.kx'*w*dta.ky				//[nx,1]
			diffW = norm(bX - bX1)
			bX = bX1
		}
		
		if (diffW > opt.tol) printf("Reached maxiter in CML. Continuing nevertheless.\n")
		dta.bX = bX
		fsave = (dta.ky-dta.kx*dta.bX)
		obj = fsave'*w*fsave
	}

	return (getret())
}
// end doCLS


// count mutual appearances in pan1 and pan2
real scalar EWreg::mutual_cnt(transmorphic colvector pan1, transmorphic colvector pan2)
{
	mutual = 0
	p1 = 1
	p2 = 1
	
	while(p1<=rows(pan1) && p2<=rows(pan2))
	{
		if(pan1[p1,1]<pan2[p2,1])
			p1=p1+1
		else if (pan1[p1,1]>pan2[p2,1])
			p2=p2+1
		else
		{
			p1=p1+1
			p2=p2+1
			mutual = mutual+1
		}
	}
	
	return (mutual)
}
// end mutual_cnt


// sharefrac hold nc/ni as defined in eq A4 of EW2011 RFS
real matrix EWreg::CMDshare(transmorphic colvector id, transmorphic colvector tm, 		///
					transmorphic colvector periods, transmorphic colvector pans)
{
	bigmat 	= (id,tm)
	bigmat 	= sort(bigmat,(2,1))
	id 		= bigmat[,1]
	tm 		= bigmat[,2]
	
	n_per   = rows(periods)
	n_panel = rows(pans)

	struct mt rowvector used
	used  = J(1, n_per, mt())
	
	ibeg = 1
	j = 1
	for (i=1;i<=rows(tm);i++)
	{
		if (tm[i,1]==periods[j,1]) continue
		used[1,j].inmat = id[ibeg..(i-1),1]
		ibeg=i
		j=j+1
	}
	used[1,j].inmat = id[ibeg..rows(id),1]
	
	if (j<n_per)
	{
		errprintf("Internal error in CMDshare. Please submit bug report.\n")
		exit(error(198))
	}
	
	sharefrac = J(n_per,n_per,0)
	for (t1=1;t1<=n_per;t1++)
	{
		n_pan = rows(used[1,t1].inmat)
		for (t2=1;t2<=n_per;t2++)
		{
			// is the entity in the panel in both time periods?
			sharefrac[t1,t2] = mutual_cnt(used[1,t1].inmat,used[1,t2].inmat)/n_pan
		}
	}
	
	return(sharefrac)
}
// end sharefrac


// find beginning position of each of ids in allids (assume both are sorted) - same as find ni
real colvector EWreg::findpos(transmorphic colvector ids, transmorphic colvector allids)
{
	pos = J(rows(ids),1,0)
	
	p_allids = 1
	
	for (p_ids=1;p_ids<=rows(ids);p_ids++)
	{
		while(allids[p_allids,1]!=ids[p_ids,1])
		{
			p_allids = p_allids+1
			if (p_allids > rows(allids))
			{
				errprintf("Internal error in findpos - please submit bug report\n")
				exit(error(198))
			}
		}
		pos[p_ids,1]=p_allids
	}
	
	return (pos)
}
//end findpos


// subroutine to compute classical minimum distance estimator
real rowvector EWreg::cmd(real rowvector csave, real matrix isave, real scalar n_period,		///
					real matrix sharefrac)
{
	// csave [1,S] - coefficients to be merged
	// isave [N,S] - influence functions of coefficients
	// retval = [theta, stderr, cmdtest, cmdval]

	ninfl = (isave:!=0)
	divideby = ninfl'*ninfl
	w        = pinv(isave'*isave:*sharefrac:*sharefrac':/(divideby:^2))
	g        = J(n_period,1,1)
	theta    = pinv(g'*w*g)*g'*w*csave'		// weighted average csave
	f        = csave' :- theta
	stderr   = sqrt(pinv(g'*w*g))
	cmdtest  = f'*w*f						// Sargan test value for the mini-gmm in here
	cmdval   = 1-chi2(n_period-1,cmdtest)	// Sargan p-val

	retval = (theta, stderr, cmdtest, cmdval)
	return (retval)
}
// end cmd


// wrap doCLS with a classical minimum distance estimator
struct stats scalar EWreg::doCMD(transmorphic colvector id, transmorphic colvector tm,			///
							real colvector y, real matrix x, real matrix z, real scalar met, 	///
							real scalar maxdeg, real rowvector bXint, real scalar vcemet, 		///
							real scalar optmet)
{
	n = rows(y)
	nx = cols(x)
	nz = cols(z)

	periods = uniqrows(tm)
	n_per   = rows(periods)
	pans    = uniqrows(id)
	n_panel = rows(pans)

	struct stats scalar ret
	struct stats scalar tmp
	ret = stats()
	tmp = stats()

	printf("\nDoing CMD post-estimation, please wait.\n")

	xz_csave = J(nx+nz,n_per,0)
	xz_isave = J(n_panel,(nx+nz)*n_per,0)
	rt_csave = J(nx+1,n_per,0)
	rt_isave = J(n_panel,(nx+1)*n_per,0)

	for (it=1;it<=n_per;it++)
	{
		//printf("Now at period %f of %f : %f\n",t,n_per,periods[t,1])
		us = (tm:==periods[it,1])
		tmp = doCLS(select(id,us), select(y,us), select(x,us), select(z,us), met, 2, maxdeg, bXint, vcemet, optmet, 0)
		
		xz_csave[.,it] = tmp.beta
		rt_csave[.,it] = tmp.rho\tmp.tau

		pos = findpos(select(id,us), pans)
		xz_isave[pos,((it-1)*(nx+nz)+1)..(it*(nx+nz))] = tmp.inflncXZ
		rt_isave[pos,((it-1)*(nx+1)+1)..(it*(nx+1))]   = tmp.inflncRT
	}
	
	sharefrac = CMDshare(id, tm, periods, pans)
	
	ret.beta 			= J(nx+nz,1,0)
	ret.VCmat 			= J(nx+nz,nx+nz,0)
	ret.tau 			= J(nx,1,0)
	ret.SEtau 			= J(nx,1,0)

	// do XZ
	for (i=1;i<=nx+nz;i++) {
		tmp1 = cmd(xz_csave[i,.], xz_isave[.,((i-1)*n_per+1)..(i*n_per)], n_per, sharefrac)
		ret.beta[i,1]  = tmp1[1,1]
		ret.VCmat[i,i] = tmp1[1,2]
	}
	
	// do R
	tmp1 = cmd(rt_csave[1,.], rt_isave[.,1..n_per], n_per, sharefrac)
	ret.rho   = tmp1[1,1]
	ret.SErho = tmp1[1,2]
	
	// do T
	for (i=2;i<=nx+1;i++) {
		tmp1 = cmd(rt_csave[i,.], rt_isave[.,((i-1)*n_per+1)..(i*n_per)], n_per, sharefrac)
		ret.tau[i-1,1]   = tmp1[1,1]
		ret.SEtau[i-1,1] = tmp1[1,2]
	}

	ret.serr 			= sqrt(diagonal(ret.VCmat))
	ret.N 				= n
	ret.Jstat 			= -1
	ret.Jval 			= -1
	ret.dfree 			= -1
	ret.vcrhotau 		= J(1,1,-1)
	ret.w 				= J(1,1,-1)
	ret.obj 			= -1
	ret.IDval 			= J(1,1,-1)
	ret.IDstat 			= J(1,1,-1)

	return (ret)
}
// end doCMD

end
