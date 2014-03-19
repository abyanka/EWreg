/* 
Author: Robert Parham, University of Rochester

This work is licensed under the Creative Commons Attribution 4.0 International License. 
To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.
*/

version 12

mata:

// represent a symbolic equation
class Symbolic {
	private:
		real matrix loc					// [N,T] - each column is a location vector for a term
										// loc holds the power of each variable in the term [N,1]
		numeric rowvector A				// [1,T]
		real scalar N					// number of variables in the multinomial
		real scalar T					// number of terms in the equation
		
		numeric scalar 				colmul()
		numeric colvector 			matmul()
		transmorphic matrix			del_mem()
		
	public:
		void 						einit()
		real scalar 				getN()
		real scalar 				getT()
		numeric rowvector 			getA()
		real matrix					getloc()
		void 						print()
		void 						addterm()
		void 						delterm()
		void 						replaceterm()
		numeric scalar 				eval()
		numeric colvector 			evalmat()
		numeric scalar 				partialeval()
		class Symbolic colvector	getgrad()
		class Symbolic matrix		gethess()
		void 						addeq()
		void 						muleq()
}
// end class


// return the multiplication of all members of a column vector
numeric scalar Symbolic::colmul(numeric colvector vec)
{
	// better to do this in a bi-section method for numerical stability, at some point
	total = 1
	for (i=1;i<=rows(vec);i++) {
		total=total*vec[i,1]
	}
	
	return (total)
}
// end function colmul


// return the multiplication of all rows of a matrix
numeric colvector Symbolic::matmul(numeric matrix mat)
{
	// better to do this in a bi-section method for numerical stability, at some point
	total = J(rows(mat),1,1)
	
	for (i=1;i<=cols(mat);i++) {
		total = total:*mat[.,i]
	}
	
	return (total)
}
// end function colmul


// remove a member at any position in a rowvector
transmorphic matrix Symbolic::del_mem(transmorphic matrix mt, real scalar to_del) {
	assert(to_del>=1 & to_del<=cols(mt))
	
	sel = J(1,cols(mt),1)
	sel[1,to_del] = 0
	return (select(mt,sel))
}
// end function mem_del


// initialize the equation, sets T=0 and sets N
void Symbolic::einit(real scalar n)
{
	assert (N==round(N) & N>0)
	T = 0
	N = n
	loc = J(N,0,0)
	A = J(1,0,0)
}
// end einit


// return N value
real scalar Symbolic::getN()
{
	return (N)
}
// end getN


// return T value
real scalar Symbolic::getT()
{
	return (T)
}
// end getT


// return the a of specific term or all of them
numeric rowvector Symbolic::getA(real scalar tr)
{
	assert (T>0 & tr==round(tr) & tr>=0 & tr<=T)
	if (tr==0) return (A)
	return (A[1,tr])
}
// end getA


// return the loc of specific term or all of them
real matrix Symbolic::getloc(real scalar tr)
{
	assert (T>0 & tr==round(tr) & tr>=0 & tr<=T)
	if (tr==0) return (loc)
	return (loc[.,tr])
}
// end getloc


// print equation
void Symbolic::print(real scalar verbosity)
{
	if (verbosity==0) return
	
	if (T==0) printf("[empty]")
	
	for (tr=1;tr<=T;tr++) {
		printf ("%g*",A[1,tr])
		for (j=1;j<=N;j++) {
			po = loc[j,tr]
			if (po != 0) {
				printf ("[%g]^%g",j,po)
				if (j<N) printf("*")
			}
		}
		if (tr<T) printf(" + ")
	}
	printf ("\n")
}
// end print


// add a term to equation - might be first term!
void Symbolic::addterm(numeric scalar a, real colvector locv)
{
	assert (a!=0 & rows(locv)==N)

	A = A, a
	loc = loc, locv
	T++
}
// end addterm


// replace a term in equation
void Symbolic::replaceterm(real scalar tr, numeric scalar a, real colvector locv)
{
	assert (a!=0 & rows(locv)==N & T>0 & tr==round(tr) & tr>0 & tr<=T)

	A[1,tr] = a
	loc[.,tr] = locv
}
// end replaceterm


// delete a term from equation - might leave it empty!
void Symbolic::delterm(real scalar tr)
{
	assert (T>0 & tr==round(tr) & tr>0 & tr<=T)
	
	A = del_mem(A,tr)
	loc = del_mem(loc,tr)
	T--
}
// end delterm


// evaluate equation at V=vars
numeric scalar Symbolic::eval(numeric colvector vars)
{
	assert (rows(vars)==N)
	
	numeric scalar retval 
	retval = 0
	
	if (T==0 | T==.) return (retval)

	for (tr=1;tr<=T;tr++) {
		retval = retval + A[1,tr]*colmul(vars:^loc[.,tr])
	}
	
	return (retval)
}
// end eval


// evaluate equation at V=vars
numeric colvector Symbolic::evalmat(numeric matrix vars)
{
	assert (cols(vars)==N)

	numeric colvector retval
	retval = J(rows(vars),1,0)
	
	if (T==0 | T==.) return (retval)
	
	for (tr=1;tr<=T;tr++) {
		retval = retval :+ A[1,tr]*matmul(vars:^loc[.,tr]')
	}
	
	return (retval)
}
// end eval


// evaluate equation at a partial vector. return the constant term resulting or 0
numeric scalar Symbolic::partialeval(numeric colvector vars)
{
	assert (rows(vars)==N)
	//vars
	
	numeric scalar cst
	cst = 0
	if (T==0 | T==.) return (cst)

	todel = J(1,0,0)
	for (tr=1;tr<=T;tr++) {
		// mark variables to eval
		evlmask = (vars:!=.):*(loc[.,tr]:!=0)
		
		if (sum(evlmask)!=0) {
			// mul into a
			A[1,tr] = A[1,tr]*colmul(select(vars,evlmask):^select(loc[.,tr],evlmask))
			// zero these variables in loc
			loc[select(1..N,evlmask'),tr] = select(J(1,N,0),evlmask')'
		}
	
		// build cst and mark empty terms
		if (sum(loc[.,tr]:!=0)==0) {
			cst = cst + A[1,tr]
			todel = todel , tr
		}
	}
	
	// remove empty terms
	for (i=0;i<cols(todel);i++) {
		delterm(todel[1,i+1]-i)
	}
	
	return (cst)
}
// end eval


// add two equations
void Symbolic::addeq(class Symbolic scalar toadd)
{
	assert (N==toadd.getN())
	
	if (toadd.getT()==0) return
	
	for (tr_add=1;tr_add<=toadd.getT();tr_add++){
		found = 0
		if (T>0) {
			for (tr=1;tr<=T;tr++){
				if (loc[.,tr]==toadd.getloc(tr_add)){
					A[1,tr] = A[1,tr] + toadd.getA(tr_add)
					found = 1
					break
				}
			}
		}
		if (found==0) {
			addterm(toadd.getA(tr_add), toadd.getloc(tr_add))
		}
	}
}
// end addeq


// multiply two equations
void Symbolic::muleq(class Symbolic scalar tomul)
{
	assert (T>0 & tomul.getT()>0 & N==tomul.getN())
	
	class Symbolic scalar retval
	retval.einit(N)
	class Symbolic scalar empty
	empty.einit(N)
	
	
	for (tr_mul=1;tr_mul<=tomul.getT();tr_mul++){
		for (tr=1;tr<=T;tr++){
			t_a = A[1,tr] * tomul.getA(tr_mul)
			t_loc = loc[.,tr] + tomul.getloc(tr_mul)
			retval.addterm(t_a, t_loc)
		}
	}
	
	// this bit has the side effect of getting rid of repeat members of retval
	empty.addeq(retval)

	T = empty.getT()
	A = empty.getA(0)
	loc = empty.getloc(0)
	
}
// end muleq


// get gradient equations for EQ
class Symbolic colvector Symbolic::getgrad()
{
	assert (T>0)
		
	class Symbolic colvector retval
	retval = J(N,1,Symbolic())
	for (iv=1;iv<=N;iv++) {
		retval[iv,1].einit(N)
	}
	
	for (tr=1;tr<=T;tr++) {
		for (iv=1;iv<=N;iv++) {
			if (loc[iv,tr] == 0) continue
			t_a = A[1,tr] * loc[iv,tr]
			t_loc = loc[.,tr]
			t_loc[iv,1] = t_loc[iv,1]-1
			retval[iv,1].addterm(t_a,t_loc)
		}
	}
	
	return (retval)
}
// end getgrad


// get hessian equations for EQ
class Symbolic matrix Symbolic::gethess()
{
	assert (T>0)

	class Symbolic matrix retval
	retval = J(N,N,Symbolic())
	for (jv=1;jv<=N;jv++){
		for (iv=1;iv<=N;iv++) {
			retval[iv,jv].einit(N)
		}
	}

	class Symbolic colvector grad
	grad = getgrad()

	class Symbolic colvector tmp
	
	for (jv=1;jv<=N;jv++){
		tmp = grad[jv,1].getgrad()
		for (iv=1;iv<=N;iv++){
			retval[iv,jv].addeq(tmp[iv,1]) 
		}
	}
	
	return (retval)
}
// end gethess

end
