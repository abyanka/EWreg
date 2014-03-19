/* 
Author: Robert Parham, University of Rochester

This work is licensed under the Creative Commons Attribution 4.0 International License. 
To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.
*/

version 12

mata:

// a problem definition class
class EWproblem {
	private:
		void						get_lhs()
		void						get_N1lhs()
		void 						get_rhs()
		void 						get_mom2()
		void 						get_cml()
		void						get_D()
		void 						get_idx()
		
		real 			rowvector 	combrep()
		real 			rowvector 	get_neqvec()
		real 			colvector 	get_tvec()
		real 			rowvector 	next_ind()
		real 			matrix 		get_rocmat()
		real 			rowvector 	occur2count()
		real 			matrix 		get_rmat()
		real 			scalar 		colmul()
		real 			scalar 		facto()
		real 			scalar 		get_avk()
		real 			rowvector	get_locvk()		
		real 			colvector	getexploc()
		real 			colvector	resolve()

	public:
		class Symbolic 	colvector	lhs				//[neq,1]
		class Symbolic 	matrix		Dlhs			//[neq,nx1]
		class Symbolic 	colvector	N1lhs			//[N1neq,1]
		class Symbolic 	colvector	rhs				//[neq,1]
		class Symbolic 	matrix		Drhs			//[neq,nt]
		class Symbolic 	colvector	mom2			//[neq2,1]
		class Symbolic 	matrix		Dmom2			//[neq2,nt2]
		class Symbolic 	colvector	cml				//[neq,1]
		class Symbolic 	matrix		Dcml			//[neq,neq]
		
		real			matrix		rmat			//[neq,nx1]
		real			matrix		rocmat			//[neq,maxdeg]
		real			matrix		N1rmat			//[N1neq,nx1]
		real			matrix		N1rocmat		//[N1neq,maxdeg]
		real			rowvector	neqvec			//[1,maxdeg-1]
		real			colvector	ntvec			//[4,1]
		real			colvector	ntvec2			//[4,1]
		real			rowvector	tmap			//[1,nt] - transforms nt to nt2 vec

		real			rowvector	yidx			//[1,nk]
		real			matrix		xidx			//[nx,nk]

		void 						einit()
		real 			matrix		getW2()
		real			colvector	resolve_eq()
}
// end EWproblem


// choose k from n with repetitions
real rowvector EWproblem::combrep(real rowvector n, real rowvector k)
{
	n1 = n :+ k :-1
	return (comb(n1,k))
}
// end combrep


// Return vector of equation numbers per each degree - returns [1,maxdeg-1] rowvector
real rowvector EWproblem::get_neqvec(real scalar maxdeg, real scalar nvars, real scalar correct)
{
	if (maxdeg==0) return (J(1,1,1))
	if (maxdeg==1) return (J(1,1,nvars))
	n = J(1,maxdeg-1,nvars)
	k = (2..maxdeg)
	degs = combrep(n,k)
	degs[1,cols(degs)] = degs[1,cols(degs)] - nvars*(maxdeg>=3)*(correct>0)
	tmp = (cols(degs)-1+1*(maxdeg==2)) //avoid subscript errors when maxdeg==2
	degs[1,tmp] = degs[1,tmp] - nvars*(maxdeg>=4)*(correct>0)
	return (degs)
}
// end get_neqvec


// count the stacked elements (including optional correction) for vector t
real colvector EWproblem::get_tvec(real scalar maxdeg, real scalar nx, real scalar correct)
{
	retval = J(4,1,0)
	//the betas
	retval[1,1] = nx
	//the E(etas)
	retval[2,1] = rowsum(get_neqvec(maxdeg,nx,0))
	//the E(epsilon^i), 2..maxdeg per xi (full)
	retval[3,1] = (maxdeg-1-1*(correct>0)-1*(correct>0 & maxdeg>=4))*nx
	//the E(u^i)
	retval[4,1] = (maxdeg-1-1*(correct>0)-1*(correct>0 & maxdeg>=4))
	
	return (retval)
}
// end function get_tcnt


// advance the indices denoting which moment combination to do next
real rowvector EWproblem::next_ind(real rowvector indx, real scalar nx1, real scalar correct)
{
	deg = cols(indx)
	retindx = indx
	for (b_d=deg;b_d>=1;b_d--) {
		if (retindx[1,b_d]+1>nx1) continue
		retindx[1,b_d] = retindx[1,b_d]+1
		// here indx is correct up to position b_d, now take care of fwd positions
		for (f_d=b_d+1;f_d<=deg;f_d++) {
			retindx[1,f_d] = retindx[1,f_d-1]
		}
		if (retindx[1,deg] == retindx[1,1] & correct>0) {
			retindx[1,deg] = retindx[1,deg] + 1
		}
		return (retindx)
	}
	// We're done
	return (J(1,0,0))
}
// end function next_ind


// create r (occurances) vectors for all eqs
real matrix EWproblem::get_rocmat(real scalar maxdeg, real scalar nx1, real scalar neqfull, 	///
								real rowvector neqvecfull)
{
	retval = J(neqfull,maxdeg,0)
	//printf("maxdeg: %g, nx1: %g, neqfull: %g\n",maxdeg, nx1, neqfull)
	
	pos=1
	for (deg=2;deg<=maxdeg;deg++) {
		indx = J(1,deg,1) // start with y^deg
		nmom = neqvecfull[1,deg-1]
		for (j=1;j<=nmom;j++) {
			retval[pos,1..deg] = indx
			pos++
			indx = next_ind(indx,nx1,0)
		}
	}
	assert(pos==neqfull+1)
	
	return(retval)
}
// end get_rocmat


// turn occur rep to count rep
real rowvector EWproblem::occur2count(real rowvector oc, real scalar nvar)
{
	retval = J(1,nvar,0)
	
	for (j=1;j<=nvar;j++) {
		retval[1,j] = sum(oc[1,.]:==j)
	}
	
	return (retval)
}
// end occur2count


// create r (counts) vectors for all eqs
real matrix EWproblem::get_rmat(real scalar nx1, real scalar neqfull)
{
	retval = J(neqfull,nx1,0)
	
	for (i=1;i<=neqfull;i++) {
		retval[i,.] = occur2count(rocmat[i,.], nx1)
	}
	
	return(retval)
}
// end get_rmat


// return the multiplication of all members of a vector
real scalar EWproblem::colmul(real colvector vec)
{
	// better to do this in a bi-section method for numerical stability, at some point
	total = 1
	for (i=1;i<=rows(vec);i++) {
		total=total*vec[i,1]
	}
	
	return (total)
}
// end function colmul


// numer!/mul(denom[i]!)
real scalar EWproblem::facto(real scalar numer, real rowvector denom)
{
	tmp = factorial(numer)
	return (tmp[1,1]/colmul(factorial(denom')))
}
//end function facto


// generates the constant part a[v,k] of a term. zeroes it if term has vanishing moments multiplied
real scalar EWproblem::get_avk(real rowvector r, real rowvector v, real rowvector k,		///
					real scalar maxdeg)
{
	// take care of vanishing u
	if (v[1,1]==1 | (v[1,1]==(maxdeg-1) & maxdeg!=3)) return (0)
	
	// take care of vanishing eta
	vk = v[1,2..(cols(v))]:+k
	if (sum(vk:>0)==1 & (sum(vk)==1 /*| (sum(vk)==(maxdeg-1) & maxdeg!=3)*/ )) return (0)
	
	a = 1
	a = a*facto(r[1,1],v)
	
	for (j=1;j<=(cols(r)-1);j++) {
		// take care of vanishing epsilon
		if (r[1,(j+1)]-k[1,j]==1 | (r[1,(j+1)]-k[1,j]==maxdeg-1 & maxdeg!=3)) return (0)
		a = a*facto(r[1,(j+1)],(k[1,j],(r[1,(j+1)]-k[1,j])))
	}
	
	return (a)
}
//end get_avk


// locvk is a list of powers of each variable in the term
real rowvector EWproblem::get_locvk(real rowvector r, real rowvector v, real rowvector k,	///
					real scalar maxdeg, real scalar nx, real scalar nt)
{
	
	retval = J(1,nt,0)
	
	base_beta	= 1
	base_eta 	= ntvec[1,1] + 1
	base_eps 	= sum(ntvec[1..2,1]) + 1
	base_u		= sum(ntvec[1..3,1]) + 1

	// u^v0
	if (v[1,1]!=0) {
		retval[1,(base_u+v[1,1]-2)] = 1
	}

	for (j=1;j<=nx;j++) {
		rj = r[1,(j+1)]
		kj = k[1,j]
		vj = v[1,(j+1)]

		// epsj^(rj-kj)
		if (rj-kj!=0) {
			retval[1,(base_eps+nx*(rj-kj-2)+(j-1))] = 1
		}

		// betaj^vj
		if (vj!=0) {
			retval[1,j] = vj
		}
	}	
	
	// v and k define a specific eta combination
	vk = v[1,2..(cols(v))]:+k
	indx = J(1,maxdeg,0)
	for (i=0;i<=nx;i++) {
		indx = next_ind(indx, nx, 0)
	}
	for (i=1;i<=ntvec[2,1];i++) {
		if (sum(occur2count(indx,nx):==vk)==nx) {
			retval[1,(base_eta + (i-1))] = 1
			break
		}
		indx = next_ind(indx, nx, 0)
	}

	return (retval)
}
// end function get_locvk


// take a subexp and return loc vector
real colvector EWproblem::getexploc(real colvector exp, real scalar nx1, real scalar neq)
{

	retval = J(neq,1,0)
	ibeg = 1
	while (1) {
		for (iend=ibeg ; exp[iend+1,1]!=-1 ; iend++) {}
		rtmp = occur2count(exp[ibeg..iend,1]',nx1)
		tmp = (colsum(rmat':==rtmp'):==nx1):*(1..neq)
		assert (sum(tmp:!=0)==1)
		pos = sum(tmp)
		retval[pos,1] = retval[pos,1]+1
		if (iend+2>rows(exp)) break
		ibeg = iend+2
		if (exp[ibeg,1]==-1) break
	}
	return (retval)
}	
// end getexploc


// build all lhs equations (EQ7 EW2002)
void EWproblem::get_lhs(real scalar nx, real scalar neq)
{
	lhs = J(neq,1,Symbolic())

	for (i=1;i<=neq;i++) {
		lhs[i,1].einit(nx+1)		// lhs equation variables are y and the xj
		lhs[i,1].addterm(1,rmat[i,.]')
	}
}
// end get_lhs


// build all N1lhs equations (y^(maxdeg-1), xj^(maxdeg-1), for SE correction)
void EWproblem::get_N1lhs(real scalar nx, real scalar N1neq)
{
	N1lhs = J(N1neq,1,Symbolic())

	for (i=1;i<=N1neq;i++) {
		N1lhs[i,1].einit(nx+1)		// lhs equation variables are y and the xj
		N1lhs[i,1].addterm(1,N1rmat[i,.]')
	}
}
// end get_N1lhs


// build all rhs equations (EQ7 EW2002)
void EWproblem::get_rhs(real scalar maxdeg, real scalar nx, real scalar neq, real scalar nt)
{
	rhs = J(neq,1,Symbolic())
	
	for (i=1;i<=neq;i++) {
		rhs[i,1].einit(nt)		// rhs equation variables are t
		r = rmat[i,.]
		
		vvec = get_neqvec(r[1,1],nx+1,0)
		nv = vvec[1,cols(vvec)]
		max_nv = rowsum(vvec)+(nx+1)+1 //+(nx+1)+1 is for deg==0,1
		m = sum(r)
		max_nk = rowsum(get_neqvec(m,nx,0))+nx+1 //+nx+1 is for deg==0,1
		

		curr_k = J(1,maxdeg,0)
		for (ik=1;ik<=max_nk;ik++) {
			count_k = occur2count(curr_k,nx)

			// if curr_k not admissable, continue
			if (sum(count_k)>m | sum(count_k:>r[1,2..(nx+1)])>0) {
				curr_k = next_ind(curr_k, nx, 0)
				continue
			}

			curr_v = J(1,maxdeg,0)
			for (iv=1;iv<=max_nv;iv++) {
				count_v = occur2count(curr_v,nx+1)
				
				// if curr_v not admissable, continue
				if (sum(count_v)!=r[1,1]) {
					curr_v = next_ind(curr_v, nx+1, 0)
					continue
				}
				
				//get_a zeroes vanishing terms
				a = get_avk(r, count_v, count_k, maxdeg)
				
				if (a>0) {
					locvk = get_locvk(r, count_v, count_k, maxdeg, nx, nt)
					rhs[i,1].addterm(a,locvk')
				}
				
				curr_v = next_ind(curr_v, nx+1, 0)
			}
			
			curr_k = next_ind(curr_k, nx, 0)
		}
	}
}
// end get_rhs


// build equations for mom2 problem
void EWproblem::get_mom2(real scalar nx, real scalar nt, real colvector ntvec)
{
	neta2 = combrep(nx,2)
	neq2 = rowsum(get_neqvec(2,(nx+1),0))
	ntvec2 = (nx,neta2,nx,1)'
	nt2 = colsum(ntvec2)
	
	mom2 = J(neq2,1,Symbolic())
	
	pos2 = ntvec[1,1] + 1
	pos3 = ntvec[1,1] + ntvec[2,1] + 1
	pos4 = ntvec[1,1] + ntvec[2,1] + ntvec[3,1] + 1
	
	tmap = J(1,nt,0)
	tmap[1,1..nx] = J(1,nx,1)										// all beta
	tmap[1,pos2..(pos2+neta2-1)] = J(1,neta2,1)						// E(etai*etaj)
	tmap[1,pos3..(pos3+nx-1)] = J(1,nx,1)							// E(ej^2)
	tmap[1,pos4] = 1												// E(u^2)

	for (i=1;i<=neq2;i++) {
		mom2[i,1].einit(nt2)
		a = rhs[i,1].getA(0)
		locmat = rhs[i,1].getloc(0)
		
		locmat_nt2 = select(locmat, tmap')
		assert (sum(select(locmat, (!tmap)'))==0)
		
		for (j=1;j<=rhs[i,1].getT();j++)
		{
			mom2[i,1].addterm(a[j],locmat_nt2[.,j])
		}
		//mom2[i,1].print(1)
	}
}
// end get_mom2


// build equations for cumulants
void EWproblem::get_cml(real scalar maxdeg, real scalar nx, real scalar neq, real scalar nt)
{
	class Partitions scalar part
	cml = J(neq,1,Symbolic())
	
	deg = rowsum(rmat)

	for (eq=1;eq<=neq;eq++) {
		cml[eq,1].einit(neq)
		
		part = Partitions()
		part.einit(deg[eq,1])
		
		for (i=1;i<=part.getP();i++) {
			pt = part.getnext()							// [deg,1] part template
			if (sum(pt:==1)!=0) continue				// E(xj^1) = E(y^1) = 0
			q = sum(pt:>0)
			a = ((-1)^(q-1))*factorial(q-1)
			expmat = part.expand(pt)					// [2*deg,?] expanded exp
			for (j=1;j<=cols(expmat);j++) {
				subexp = part.sub(expmat[.,j],rocmat[eq,1..deg[eq,1]]')
				locv = getexploc(subexp, nx+1, neq)
				cml[eq,1].addterm(a,locv)
			}
		}
	}
}
// end get_cml


// build indices for choosing cumulants
void EWproblem::get_idx(real scalar maxdeg, real scalar nx, real scalar nk,			///
				real scalar neq, real rowvector neqvecfull)
{
	yidx = J(1,nk,0)
	xidx = J(nx,nk,0)
	
	pos = 1
	for (d=2;d<=maxdeg-1;d++) {
		idx = J(1,d,1)										// start with y^d
		idx = next_ind(idx,nx+1,1)							// and advance one
		for (i=1; i<=(neqvecfull[1,d-1]-(nx+1)); i++) {
			ry = occur2count((1,idx), nx+1)
			yidx[1,pos] = sum((colsum(rmat':==ry'):==nx+1):*(1..neq))
			for (j=1; j<=nx; j++) {
				rx = occur2count((j+1,idx), nx+1)
				xidx[j,pos] = sum((colsum(rmat':==rx'):==nx+1):*(1..neq))
			}
			pos++
			idx = next_ind(idx, nx+1, 1)
		}
	}
}
// end get_idx


// build equations for gradient matrix of system
void EWproblem::get_D(class Symbolic colvector systm, class Symbolic matrix Dsys)
{
	neq   = rows(systm)
	nvars = systm[1,1].getN()
	
	class Symbolic colvector tmp
	
	for (i=1;i<=neq;i++) {
		tmp = systm[i,1].getgrad()	//tmp[nvars,1]
		
		for (j=1;j<=nvars;j++) {
			Dsys[i,j].einit(nvars)
			Dsys[i,j].addeq(tmp[j,1])
		}
	}
}
// end get_D


// initialize problem
void EWproblem::einit(real scalar maxdeg, real scalar nx)
{
	correct = 1
	neqvec = get_neqvec(maxdeg, (nx+1), correct)
	neqvecfull = get_neqvec(maxdeg, (nx+1), 0)
	ntvec = get_tvec(maxdeg, nx, correct)
	neq = rowsum(neqvec)
	neqfull = rowsum(neqvecfull)
	nt = colsum(ntvec)
	nk = rowsum(neqvecfull[1,1..cols(neqvecfull)-1]:-(nx+1))
	
	rocmat = get_rocmat(maxdeg, (nx+1), neqfull, neqvecfull)
	rmat = get_rmat((nx+1), neqfull)
	
	N1sel = ((rowsum(rmat:!=0):==1):*rowsum(rmat):==(maxdeg-1))
	N1rocmat = select(rocmat, N1sel)
	N1rmat = select(rmat, N1sel)
	N1neq = rows(N1rmat)
	
	if (correct) {
		corr_sel = (rowsum(rmat:!=0):==1):*rowsum(rmat)
		corr_sel = (corr_sel:==(maxdeg-1)):*(maxdeg>=4) + (corr_sel:==maxdeg):*(maxdeg>=3)
		rocmat = select(rocmat, !corr_sel)
		rmat = select(rmat, !corr_sel)
	}
	assert (rows(rmat)==neq)
	
	get_lhs(nx, neq)
	get_N1lhs(nx, N1neq)
	get_rhs(maxdeg, nx, neq, nt)
	get_mom2(nx, nt, ntvec)
	get_cml(maxdeg, nx, neq, nt)
	
	// RPG - Dlhs currently unused, so don't waste time
	//Dlhs  = J(rows(lhs),lhs[1,1].getN(),Symbolic())
	Drhs  = J(rows(rhs),rhs[1,1].getN(),Symbolic())
	Dmom2 = J(rows(mom2),mom2[1,1].getN(),Symbolic())
	Dcml  = J(rows(cml),cml[1,1].getN(),Symbolic())

	//get_D(lhs, Dlhs)
	get_D(rhs, Drhs)
	get_D(mom2, Dmom2)
	get_D(cml, Dcml)
	
	get_idx(maxdeg, nx, nk, neq, neqvecfull)
	
	return
}
// end einit


// return the lower part of the W matrix of Eq26
real matrix EWproblem::getW2(real scalar nx, real scalar nr, real scalar neq)
{
	retval = J(nr,nx,1)
	
	idx = (2,3)
	for (i=1;i<=nr;i++) {
		for (j=1;j<=nx;j++) {
			r = idx , j+1
			rc = occur2count(r,nx+1)
			retval[i,j] = sum((colsum(rmat':==rc'):==nx+1):*(1..neq))
		}
		idx = next_ind(idx,nx+1,1)
	}
	
	return(retval)	
}
// end getW2


// resolve a general set of equations by iterated substitution
real colvector EWproblem::resolve(real colvector retval, real colvector LHS, 
						class Symbolic colvector RHS, real scalar neq, real scalar nt)
{
	// repeat until all of t is resolved or too many iterations
	for (iter=1;(sum(retval:==.)>0 & iter<=10);iter++) {
		for (eq=1;eq<=neq;eq++) {
			// make sure everything in retval is substituted into the equation
			cst = RHS[eq,1].partialeval(retval)
			// deduct constant from LHS
			LHS[eq,1] = LHS[eq,1] - cst
			// for 1-term equations, divide by a
			if (RHS[eq,1].getT()==1) {
				LHS[eq,1] = LHS[eq,1]/RHS[eq,1].getA(1)
				locv = RHS[eq,1].getloc(1)
				RHS[eq,1].replaceterm(1, 1, locv)
				
				// for 1-member 1-term equations, resolve if needed, discard anyway
				if (sum(locv:!=0)==1) {
					loc = sum((locv:!=0):*(1..nt)')
					if (retval[loc,1]==.) retval[loc,1] = LHS[eq,1]:^(1/sum(locv))
					RHS[eq,1].delterm(1)
				}
			}
		}
	}
	
	if (sum(retval:==.)>0) {
		printf("Iterated substitution failed. Please submit bug report.")
		exit(197)
	}
	
	return (retval)
}
// end resolve


// resolve a set of equations
real colvector EWproblem::resolve_eq(class Symbolic colvector systm, real colvector start,		///
						real scalar neq, real scalar nx, real scalar nt, real colvector emom)
{
	retval = J(nt,1,.)
	retval[1..nx,1] = start
	
	class Symbolic colvector RHS
	RHS = J(neq,1,Symbolic())
	LHS = J(neq,1,0)
	for (i=1;i<=neq;i++) {
		RHS[i,1].einit(nt)
		RHS[i,1].addeq(systm[i,1])
		LHS[i,1] = emom[i,1]
	}

	return (resolve(retval, LHS, RHS, neq, nt))	
}
// end resolve_eq

end
