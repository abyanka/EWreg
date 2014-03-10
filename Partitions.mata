*! version 1.0.0 class Partitions

/* 
Author: Robert Parham, University of Rochester

This work is licensed under the Creative Commons Attribution 4.0 International License. 
To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.

Several functions are adapted from Hankin RKS (2005). “Additive integer partitions in R.”, 
Journal of Statistical Software, Code Snippets, 16(1)
*/


version 11

mata:

// define Partitions class
class Partitions {
	private:
		real scalar D
		real scalar P
		real matrix PMAT
		static real scalar CURR
		
		real rowvector 				getpartvec()
		real colvector 				nextpart()
		real colvector				nextexp()
		
	public:
		void 						einit()
		real scalar 				getP()
		real colvector 				getnext()
		real matrix 				expand()
		real colvector				sub()
}
// end class


// return a number of parts per cumulant degree vector (e.g. f(4)=15)
real rowvector Partitions::getpartvec()
{
	/* These functions are adapted from Hankin RKS (2005). “Additive integer partitions in R.”, 
	Journal of Statistical Software, Code Snippets, 16(1) */
	
	p = J(1,D,0) 
	p[1,1] = 1
	p[1,2] = 1
	for (i=3 ; i<=D ; i++)
	{
		// first do r = m(3m+1)/2
		s = 1  /* "s" for "sign" */
		f = 5  /* f is first difference  */
		r = 2  /* initial value (viz m(3m+1)/2 for m=1) */

		p[1,i] = 0
		while(i-r >= 1)
		{
			p[1,i] = p[1,i] + s*p[1,i-r]
			r = r+f
			f = f+3  /* 3 is the second difference */
			/* change sign of s */
			s = s*(-1)
		}
		/* now do r = m(3m-1)/2 */
		s = 1
		f = 4 /*first difference now 4 */
		r = 1 /* initial value (viz m(3m-1)/2 for m=1) */
		while(i-r >= 1)
		{
			p[1,i] = p[1,i] + s*p[1,i-r];
			r = r+f
			f = f+3  /* 3 is the second difference */
			s = s*(-1)
		}
	}
	
	return (p)
}
// end get_partvec


// build the next parts vector given the current one
real colvector Partitions::nextpart(real scalar ind)
{
	x = PMAT[.,ind]'
	
	a=1
	while(x[1,a] > 0) {
		a++
	} 
	a-- 					/* a: pos of last nonzero */
	
	b=a
	while(x[1,b] == 1){
		b--
	}

	if(x[1,a]>1){ 			/* if last nonzero number >1 */
		x[1,a] = x[1,a]-1 	/* subtract one from it */
		x[1,a+1] = 1 		/* and put a 1 next to it */
		return (x')
	}
	
	n = a-b   				/* n: number of 1s*/
	x[1,b] = x[1,b]-1 		/* decrement final nonzero digit (perforce >1) */
	yy = x[1,b]				/* and prepare to replicate it as many times as possible */
	n++;					/* n is now number of 1s plus 1 (thus "n" is to
							 * to be distributed as "evenly" as possible
							 * after x[a]) */ 
	j = b
	while(n >= yy){			/* distribute n among elements x[b] onwards */
		j++
		x[1,j] = yy			/* by repeatedly adding yy elements until no longer */
		n = n-yy       		/* possible */ 
	}
	if(n > 0){
		j++
		x[1,j] = n			/* add the remainder to x */
	}
	while(j < a){			/* set remaining elements to 0 */
		j++
		x[1,j] = 0
	}
	
	return (x')
}
// end nextpart


// take partition-specific exp, "advance by one", return. return J(1,1,-1) if last exp
real colvector Partitions::nextexp(real colvector exp, real colvector begvec, real colvector pt)
{
	/* if only one group, return */
	if (rows(begvec)==1) return (J(1,1,-1))
	
	
	/* try advancing tail */
	res = nextexp(exp[begvec[2,1]..rows(exp),1],begvec[2..rows(begvec),1]:-(begvec[2,1]-1), pt[2..rows(pt),1])
	if (res[1,1]!=-1) return (exp[1..(begvec[1,1]+pt[1,1]),1] \ res)
	
	
	/* advance first group, initialize tail */
	uid = sort(uniqrows(exp),1)
	uid = uid[2..rows(uid)] 			// get rid of -1
	
	// last cell of first group
	j = begvec[1,1]+pt[1,1]-1
	
	// find cell to increase in first block
	for(pos=rows(uid);j>0;j--,pos--) {
		if (exp[j]!=uid[pos]) break		
	}
	// no cell left to increase in first block
	if (j==0) return (J(1,1,-1))

	// increase cell
	for(;exp[j]!=uid[pos];pos--){}		
	pos++
	exp[j]=uid[pos]
	
	// remove first block members upto j from uid
	good = J(rows(uid),1,1)
	for(i=j;i>0;i--){
		for(;exp[i]!=uid[pos];pos--){}
		good[pos] = 0
	}
	uid = select(uid,good)
	
	// first block must be index-increasing
	if (exp[j+1]!=-1) {
		for (;uid[pos]<=exp[j];pos++) {}
		j++
		good = J(rows(uid),1,1)
		for (;exp[j]!=-1;j++) {
			exp[j]=uid[pos]
			good[pos]=0
			pos++
		}
		uid = select(uid,good)
	}

	// spread what's left orderly across non(-1) members of exp
	for (i=1;i<=rows(uid);i++) {
		if (exp[j+i]==-1) j++
		exp[j+i] = uid[i]
	}
	
	return (exp)
}
// end next_exp


// initialize a partition
void Partitions::einit(real scalar deg)
{
	assert (deg>=2)
	
	D = deg
	pvec = getpartvec()
	P = pvec[1,D]
	PMAT = J(D,P,0)
	CURR = 0
	PMAT[1,1] = D
	for(i=2 ; i<=P ; i++){
		PMAT[.,i] = nextpart(i-1)
	}
}
// end init


// return the number of partitions
real scalar Partitions::getP()
{
	assert (D!=.)
	return (P)
}
// end getP


// return next partition
real colvector Partitions::getnext()
{
	assert (D!=. & CURR <= P)
	CURR++
	return (PMAT[.,CURR])
}
// end getnext


// return a matrix whose columns are raw exp vectors (implement [k] in tensor notation)
real matrix Partitions::expand(real colvector pt)
{
	D2 = 2*D
	ngrp = sum(pt:!=0)
	retval = J(D2,1,-1)
	logeq = 0

	// build block beginning vector for pt
	begvec = J(ngrp,1,1)
	for (i=2;i<=ngrp;i++) {
		begvec[i,1] = begvec[i-1,1]+pt[i-1,1]+1
	}

	// fill the first exp
	for (i=1;i<=ngrp;i++) {
		curgrp = (begvec[i,1]-i+1)..(sum(pt[1..i,1]))
		curidx = (begvec[i,1])..(begvec[i,1]+pt[i,1]-1)
		retval[curidx,1] = curgrp'
	}

	// fill the rest of them
	while (1) {
		tmp = nextexp(retval[.,cols(retval)], begvec, pt)
		if (tmp[1,1]==-1) break
		retval = retval, tmp
	}
	
	// verify all generated exps
	good = J(1,cols(retval),1)
	for (ie=1;ie<=cols(retval);ie++) {
		
		// verify exp obeys the index order condition for same sized blocks
		keepexp = 1
		for (i=1;i<ngrp & keepexp==1;i++) {
			for(j=i+1;j<=ngrp;j++) {
				if (pt[i]!=pt[j]) continue
				if (retval[begvec[i,1],ie] > retval[begvec[j,1],ie]) {
					keepexp = 0
					break
				}
			}
		}
		if (!keepexp) good[1,ie] = 0
	}
	
	retval = select(retval,good)

	return (retval)
}
// end expand


// get a mapping vector and an exp, sub map into exp
real colvector Partitions::sub(real colvector exp, real colvector map) {
	assert(max(exp)==rows(map))
	retval = J(rows(exp),1,-1)
	for (i=1;i<=rows(exp);i++) {
		if (exp[i] == -1) continue
		retval[i] = map[exp[i]]
	}
	
	return (retval)
}
// end sub

end
