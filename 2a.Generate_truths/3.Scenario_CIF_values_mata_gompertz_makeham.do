


 
clear all  
include "$N/1b.Scenarios_folder/Scenario10_new.do"

 
mata:

CIF1_mat=J(5,11,.)
CIF2_mat=J(5,11,.)

CIF1_mat[1,11]=50
CIF1_mat[2,11]=60
CIF1_mat[3,11]=70
CIF1_mat[4,11]=80
CIF1_mat[5,11]=90

CIF2_mat[1,11]=50
CIF2_mat[2,11]=60
CIF2_mat[3,11]=70
CIF2_mat[4,11]=80
CIF2_mat[5,11]=90

alpha = strtoreal(st_global("alpha"))
lambda = strtoreal(st_global("lambda"))
gamma = strtoreal(st_global("gamma"))


lambda1 = strtoreal(st_global("lambda1"))
gamma1 = strtoreal(st_global("gamma1"))
lambda2 = strtoreal(st_global("lambda2"))
gamma2 = strtoreal(st_global("gamma2"))
pmix = strtoreal(st_global("pmix"))

bsex_c2 = strtoreal(st_global("bsex_c2"))
bsex_c2_tvc_1 = strtoreal(st_global("bsex_c2_tvc_1"))
bsex_c2_tvc_2 = strtoreal(st_global("bsex_c2_tvc_2"))
bsex_c1 = strtoreal(st_global("bsex_c1"))
b1 = strtoreal(st_global("b1"))
b2 = strtoreal(st_global("b2"))

// gq() - function to calculae gauss-legendre nodes and weights

//Define function not giving direct results unless you ask for them by typing them at the end of the function
// weights matrix, nodes matrix, scalar n
void gq(real matrix weights, real matrix nodes, real scalar n)
{

//Define objects like C: il, j , muzero will be scalars while a,A,vec will be matrices
real scalar i1, j, muzero, b
real matrix a, A, vec 

//il will range from 1 to n-1 by a step of 1
//muzero must be a starting value of sorts
//a is a zero vector (1 row, n-1 columns)  
//b is a simple function of il. It will be a vector of n-1 columns. weight function
//A will be a diagonal matrix with the diagonal elements be those of the vector a (a.k.a zeros) at the start.
//However with the loop, A's j diagonal element will be the b[j] result of funtion b. 
 
        i1 = range(1,n-1,1)'
        muzero = 2
        a = J(1,n,0)
        b = i1:/sqrt(4 :* i1:^2 :- 1)   
        A= diag(a)
        for(j=1;j<=n-1;j++){
                A[j,j+1] = b[j]
                A[j+1,j] = b[j]
        }       
		
//pragma unset X suppresses the warning message note: variable X may be used before set

                                pragma unset vec
								
//eigensystem(A, X, L, rcond, nobalance) calculates eigenvectors and eigenvalues of a general, real or complex, square matrix A. Eigenvectors are returned in X and eigenvalues in L. 
//in symeigensystem the difference is that A is Hermitian-symmetric
                                symeigensystem(A,vec,nodes)
								
//weights is a vector. It is the result of the product of the vector vec^2 with muzero and then it gets ordered by nodes matrix???
//Nodes is a matrix of the given nodes ordered by them?


        weights = (vec[1,]:^2:*muzero)'
        weights = weights[order(nodes',1)]
        nodes = nodes'[order(nodes',1)']
 
}
 
// function to derive hazard function for other cause mortality
real matrix hazard2(t,age,alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2) 
{
    return(
	((alpha :+lambda:*exp(gamma:*(age:+t))):*exp(bsex_c2 :+ bsex_c2_tvc_1:*(age:+t) :+ bsex_c2_tvc_2:*(age:+t):^2)  ) )
	
}


hazard2((1..10),80,alpha,lambda,gamma,bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2)

// function to derive survival function for other cause mortality
real matrix survival2(t, age, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2, n) 
{
	gq(wt=.,nodes=.,n)
    enter = 0
    nodes_t = 0.5:*((age:+t):-enter):*nodes' :+ 0.5:*((age:+t):+enter)
    
    ht_nodes = hazard2(nodes_t, 0, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2 )

    return(exp(-0.5:*((age:+t):-enter):*quadrowsum(wt':*ht_nodes)))
}

survival2(80,0,alpha,lambda,gamma,bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,70)

survival2(10,80,alpha,lambda,gamma,bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,70)

// function to derive survival function for cancer mortality
real matrix survival1(t, age, lambda1, gamma1, lambda2, gamma2, pmix, bsex_c1, b1, b2) 
{
	return((pmix :*exp(-lambda1:*t:^gamma1):+(1:-pmix):*exp(-lambda2:*t:^gamma2 )):^ 
	        exp((b1:*(age:-65):+(b2:*(age:-65):^2):+bsex_c1)))

}

survival1(10, 80, lambda1, gamma1, lambda2, gamma2, pmix, bsex_c1, b1, b2) 

// function to derive hazard function for cancer mortality
real matrix hazard1(t, age, lambda1, gamma1, lambda2, gamma2, pmix, bsex_c1, b1, b2) 
{
	return((((lambda1:*gamma1:*t:^(gamma1:-1):*pmix:*exp(-lambda1:*t:^gamma1)):+(lambda2:*gamma2:*t:^(gamma2:-1)):* 
	       (1:-pmix):*exp(-lambda2:*t:^gamma2)):/(pmix :*exp(-lambda1 :*t:^gamma1):+ 
		   (1:-pmix):*exp(-lambda2:*t:^gamma2))):*exp((b1:*(age:-65):+ b2:*((age:-65):^2):+ bsex_c1) ))

}

hazard1(10, 80, lambda1, gamma1, lambda2, gamma2, pmix, bsex_c1, b1, b2) 


// function to derive CIF for cancer mortality (CIF1)
real matrix CIF1(t, age, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, n) {
	gq(wt=.,nodes=.,n)
	enter = 0
	nodes_t = 0.5:*(t:-enter):*nodes' :+ 0.5:*(t:+enter)
 
	nodes_t
// create empty matrix
	CIF1t_nodes = J(1,n,.)
    CIF1t_nodes            
// loop over nodes
// evaluation of each node requires numerical integration
	for(i=1;i<=n;i++) {
		tn = nodes_t[1,i]
		CIF1t_nodes[1,i] =  survival1(tn, age, lambda1, gamma1, lambda2, gamma2, pmix, bsex_c1, b1, b2) :*
							survival2(tn, age, alpha,lambda,gamma,bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,n) :*
							hazard1(tn, age, lambda1, gamma1, lambda2, gamma2, pmix, bsex_c1, b1, b2) :/
							survival2(age, 0, alpha,lambda,gamma,bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,n)
	}
    return(0.5:*(t:-enter):*quadrowsum(wt':*CIF1t_nodes))
} 
 
 CIF1(10, 80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)

  // function to derive CIF for other cause mortality (CIF2)
real matrix CIF2(t, age, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, n) {
	gq(wt=.,nodes=.,n)
	enter = 0
	nodes_t = 0.5:*(t:-enter):*nodes' :+ 0.5:*(t:+enter)
 
	nodes_t
// create empty matrix
	CIF2t_nodes = J(1,n,.)
    CIF2t_nodes            
// loop over nodes
// evaluation of each node requires numerical integration
	for(i=1;i<=n;i++) {
		tn = nodes_t[1,i]
		CIF2t_nodes[1,i] =  survival1(tn, age, lambda1, gamma1, lambda2, gamma2, pmix, bsex_c1, b1, b2) :*
							survival2(tn, age, alpha,lambda,gamma,bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,n) :*
							hazard2(  tn, age, alpha,lambda,gamma,bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2) :/
							survival2(age,  0, alpha,lambda,gamma,bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,n)
	}
    return(0.5:*(t:-enter):*quadrowsum(wt':*CIF2t_nodes))
} 
 
 CIF2(10,80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)

   //Create the matrices that contain the true CIF1 for 50,60,70,80,90 years of age at diagnosis
 CIF1_mat[1,1]=  CIF1(1,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,2]=  CIF1(2,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,3]=  CIF1(3,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,4]=  CIF1(4,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,5]=  CIF1(5,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,6]=  CIF1(6,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,7]=  CIF1(7,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,8]=  CIF1(8,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,9]=  CIF1(9,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,10]= CIF1(10, 50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 CIF1_mat[2,1]=  CIF1(1,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,2]=  CIF1(2,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,3]=  CIF1(3,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,4]=  CIF1(4,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,5]=  CIF1(5,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,6]=  CIF1(6,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,7]=  CIF1(7,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,8]=  CIF1(8,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,9]=  CIF1(9,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,10]= CIF1(10, 60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 CIF1_mat[3,1]=  CIF1(1,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,2]=  CIF1(2,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,3]=  CIF1(3,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,4]=  CIF1(4,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,5]=  CIF1(5,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,6]=  CIF1(6,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,7]=  CIF1(7,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,8]=  CIF1(8,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,9]=  CIF1(9,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,10]= CIF1(10, 70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 CIF1_mat[4,1]=  CIF1(1,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,2]=  CIF1(2,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,3]=  CIF1(3,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,4]=  CIF1(4,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,5]=  CIF1(5,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,6]=  CIF1(6,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,7]=  CIF1(7,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,8]=  CIF1(8,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,9]=  CIF1(9,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,10]= CIF1(10, 80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 CIF1_mat[5,1]=  CIF1(1,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,2]=  CIF1(2,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,3]=  CIF1(3,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,4]=  CIF1(4,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,5]=  CIF1(5,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,6]=  CIF1(6,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,7]=  CIF1(7,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,8]=  CIF1(8,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,9]=  CIF1(9,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,10]= CIF1(10, 90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 st_matrix("CIF1_mat", CIF1_mat)
 
  //Create the matrices that contain the true CIF2 for 50,60,70,80,90 years of age at diagnosis
 CIF2_mat[1,1]=  CIF2(1,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,2]=  CIF2(2,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,3]=  CIF2(3,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,4]=  CIF2(4,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,5]=  CIF2(5,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,6]=  CIF2(6,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,7]=  CIF2(7,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,8]=  CIF2(8,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,9]=  CIF2(9,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,10]= CIF2(10, 50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 CIF2_mat[2,1]=  CIF2(1,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,2]=  CIF2(2,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,3]=  CIF2(3,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,4]=  CIF2(4,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,5]=  CIF2(5,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,6]=  CIF2(6,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,7]=  CIF2(7,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,8]=  CIF2(8,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,9]=  CIF2(9,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,10]= CIF2(10, 60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 CIF2_mat[3,1]=  CIF2(1,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,2]=  CIF2(2,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,3]=  CIF2(3,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,4]=  CIF2(4,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,5]=  CIF2(5,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,6]=  CIF2(6,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,7]=  CIF2(7,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,8]=  CIF2(8,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,9]=  CIF2(9,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,10]= CIF2(10, 70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 CIF2_mat[4,1]=  CIF2(1,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,2]=  CIF2(2,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,3]=  CIF2(3,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,4]=  CIF2(4,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,5]=  CIF2(5,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,6]=  CIF2(6,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,7]=  CIF2(7,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,8]=  CIF2(8,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,9]=  CIF2(9,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,10]= CIF2(10, 80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 CIF2_mat[5,1]=  CIF2(1,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,2]=  CIF2(2,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,3]=  CIF2(3,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,4]=  CIF2(4,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,5]=  CIF2(5,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,6]=  CIF2(6,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,7]=  CIF2(7,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,8]=  CIF2(8,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,9]=  CIF2(9,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,10]= CIF2(10, 90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 st_matrix("CIF2_mat", CIF2_mat)

 
end

mat li CIF1_mat
mat li CIF2_mat

cd "$N\2b.Truths_folder" 

svmat double CIF1_mat
svmat double CIF2_mat

 //matsave CIF1_mat, replace p("CIF1_mat") saving replace
 //matsave CIF2_mat, replace p("CIF2_mat") saving replace
gen scenario=10

rename CIF1_mat1  CIF1_true_f1_at1
rename CIF1_mat2  CIF1_true_f1_at2
rename CIF1_mat3  CIF1_true_f1_at3
rename CIF1_mat4  CIF1_true_f1_at4
rename CIF1_mat5  CIF1_true_f1_at5
rename CIF1_mat6  CIF1_true_f1_at6
rename CIF1_mat7  CIF1_true_f1_at7
rename CIF1_mat8  CIF1_true_f1_at8
rename CIF1_mat9  CIF1_true_f1_at9
rename CIF1_mat10  CIF1_true_f1_at10


rename CIF2_mat1  CIF2_true_f1_at1
rename CIF2_mat2  CIF2_true_f1_at2
rename CIF2_mat3  CIF2_true_f1_at3
rename CIF2_mat4  CIF2_true_f1_at4
rename CIF2_mat5  CIF2_true_f1_at5
rename CIF2_mat6  CIF2_true_f1_at6
rename CIF2_mat7  CIF2_true_f1_at7
rename CIF2_mat8  CIF2_true_f1_at8
rename CIF2_mat9  CIF2_true_f1_at9
rename CIF2_mat10  CIF2_true_f1_at10


save "$N\2b.Truths_folder/True_CIFs_scenario10", replace

rename CIF1_mat11 givenage
drop CIF2_mat11

gen id=_n

//Reshape true CIF values from wide to long
reshape long CIF1_true_f1_at CIF2_true_f1_at , i(id) j(year)


//Save true CIF values as datasets
save "$N\2b.Truths_folder/True_CIFs_scenario10", replace

/********************************************************************/
 
clear all  
include "$N/1b.Scenarios_folder/Scenario12_new.do"

 
mata:

CIF1_mat=J(5,11,.)
CIF2_mat=J(5,11,.)

CIF1_mat[1,11]=50
CIF1_mat[2,11]=60
CIF1_mat[3,11]=70
CIF1_mat[4,11]=80
CIF1_mat[5,11]=90

CIF2_mat[1,11]=50
CIF2_mat[2,11]=60
CIF2_mat[3,11]=70
CIF2_mat[4,11]=80
CIF2_mat[5,11]=90

alpha = strtoreal(st_global("alpha"))
lambda = strtoreal(st_global("lambda"))
gamma = strtoreal(st_global("gamma"))


lambda1 = strtoreal(st_global("lambda1"))
gamma1 = strtoreal(st_global("gamma1"))
lambda2 = strtoreal(st_global("lambda2"))
gamma2 = strtoreal(st_global("gamma2"))
pmix = strtoreal(st_global("pmix"))

bsex_c2 = strtoreal(st_global("bsex_c2"))
bsex_c2_tvc_1 = strtoreal(st_global("bsex_c2_tvc_1"))
bsex_c2_tvc_2 = strtoreal(st_global("bsex_c2_tvc_2"))
bsex_c1 = strtoreal(st_global("bsex_c1"))
b1 = strtoreal(st_global("b1"))
b2 = strtoreal(st_global("b2"))

// gq() - function to calculae gauss-legendre nodes and weights

//Define function not giving direct results unless you ask for them by typing them at the end of the function
// weights matrix, nodes matrix, scalar n
void gq(real matrix weights, real matrix nodes, real scalar n)
{

//Define objects like C: il, j , muzero will be scalars while a,A,vec will be matrices
real scalar i1, j, muzero, b
real matrix a, A, vec 

//il will range from 1 to n-1 by a step of 1
//muzero must be a starting value of sorts
//a is a zero vector (1 row, n-1 columns)  
//b is a simple function of il. It will be a vector of n-1 columns. weight function
//A will be a diagonal matrix with the diagonal elements be those of the vector a (a.k.a zeros) at the start.
//However with the loop, A's j diagonal element will be the b[j] result of funtion b. 
 
        i1 = range(1,n-1,1)'
        muzero = 2
        a = J(1,n,0)
        b = i1:/sqrt(4 :* i1:^2 :- 1)   
        A= diag(a)
        for(j=1;j<=n-1;j++){
                A[j,j+1] = b[j]
                A[j+1,j] = b[j]
        }       
		
//pragma unset X suppresses the warning message note: variable X may be used before set

                                pragma unset vec
								
//eigensystem(A, X, L, rcond, nobalance) calculates eigenvectors and eigenvalues of a general, real or complex, square matrix A. Eigenvectors are returned in X and eigenvalues in L. 
//in symeigensystem the difference is that A is Hermitian-symmetric
                                symeigensystem(A,vec,nodes)
								
//weights is a vector. It is the result of the product of the vector vec^2 with muzero and then it gets ordered by nodes matrix???
//Nodes is a matrix of the given nodes ordered by them?


        weights = (vec[1,]:^2:*muzero)'
        weights = weights[order(nodes',1)]
        nodes = nodes'[order(nodes',1)']
 
}
 
// function to derive hazard function for other cause mortality
real matrix hazard2(t,age,alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2) 
{
    return(
	((alpha :+lambda:*exp(gamma:*(age:+t))):*exp(bsex_c2 :+ bsex_c2_tvc_1:*(age:+t) :+ bsex_c2_tvc_2:*(age:+t):^2)  ) )
	
}


hazard2((1..10),80,alpha,lambda,gamma,bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2)

// function to derive survival function for other cause mortality
real matrix survival2(t, age, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2, n) 
{
	gq(wt=.,nodes=.,n)
    enter = 0
    nodes_t = 0.5:*((age:+t):-enter):*nodes' :+ 0.5:*((age:+t):+enter)
    
    ht_nodes = hazard2(nodes_t, 0, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2 )

    return(exp(-0.5:*((age:+t):-enter):*quadrowsum(wt':*ht_nodes)))
}

survival2(80,0,alpha,lambda,gamma,bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,70)

survival2(10,80,alpha,lambda,gamma,bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,70)

// function to derive survival function for cancer mortality
real matrix survival1(t, age, lambda1, gamma1, lambda2, gamma2, pmix, bsex_c1, b1, b2) 
{
	return((pmix :*exp(-lambda1:*t:^gamma1):+(1:-pmix):*exp(-lambda2:*t:^gamma2 )):^ 
	        exp((b1:*(age:-65):+(b2:*(age:-65):^2):+bsex_c1)))

}

survival1(10, 80, lambda1, gamma1, lambda2, gamma2, pmix, bsex_c1, b1, b2) 

// function to derive hazard function for cancer mortality
real matrix hazard1(t, age, lambda1, gamma1, lambda2, gamma2, pmix, bsex_c1, b1, b2) 
{
	return((((lambda1:*gamma1:*t:^(gamma1:-1):*pmix:*exp(-lambda1:*t:^gamma1)):+(lambda2:*gamma2:*t:^(gamma2:-1)):* 
	       (1:-pmix):*exp(-lambda2:*t:^gamma2)):/(pmix :*exp(-lambda1 :*t:^gamma1):+ 
		   (1:-pmix):*exp(-lambda2:*t:^gamma2))):*exp((b1:*(age:-65):+ b2:*((age:-65):^2):+ bsex_c1) ))

}

hazard1(10, 80, lambda1, gamma1, lambda2, gamma2, pmix, bsex_c1, b1, b2) 


// function to derive CIF for cancer mortality (CIF1)
real matrix CIF1(t, age, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, n) {
	gq(wt=.,nodes=.,n)
	enter = 0
	nodes_t = 0.5:*(t:-enter):*nodes' :+ 0.5:*(t:+enter)
 
	nodes_t
// create empty matrix
	CIF1t_nodes = J(1,n,.)
    CIF1t_nodes            
// loop over nodes
// evaluation of each node requires numerical integration
	for(i=1;i<=n;i++) {
		tn = nodes_t[1,i]
		CIF1t_nodes[1,i] =  survival1(tn, age, lambda1, gamma1, lambda2, gamma2, pmix, bsex_c1, b1, b2) :*
							survival2(tn, age, alpha,lambda,gamma,bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,n) :*
							hazard1(tn, age, lambda1, gamma1, lambda2, gamma2, pmix, bsex_c1, b1, b2) :/
							survival2(age, 0, alpha,lambda,gamma,bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,n)
	}
    return(0.5:*(t:-enter):*quadrowsum(wt':*CIF1t_nodes))
} 
 
 CIF1(10, 80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)

  // function to derive CIF for other cause mortality (CIF2)
real matrix CIF2(t, age, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, n) {
	gq(wt=.,nodes=.,n)
	enter = 0
	nodes_t = 0.5:*(t:-enter):*nodes' :+ 0.5:*(t:+enter)
 
	nodes_t
// create empty matrix
	CIF2t_nodes = J(1,n,.)
    CIF2t_nodes            
// loop over nodes
// evaluation of each node requires numerical integration
	for(i=1;i<=n;i++) {
		tn = nodes_t[1,i]
		CIF2t_nodes[1,i] =  survival1(tn, age, lambda1, gamma1, lambda2, gamma2, pmix, bsex_c1, b1, b2) :*
							survival2(tn, age, alpha,lambda,gamma,bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,n) :*
							hazard2(  tn, age, alpha,lambda,gamma,bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2) :/
							survival2(age,  0, alpha,lambda,gamma,bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,n)
	}
    return(0.5:*(t:-enter):*quadrowsum(wt':*CIF2t_nodes))
} 
 
 CIF2(10,80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)

   //Create the matrices that contain the true CIF1 for 50,60,70,80,90 years of age at diagnosis
 CIF1_mat[1,1]=  CIF1(1,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,2]=  CIF1(2,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,3]=  CIF1(3,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,4]=  CIF1(4,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,5]=  CIF1(5,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,6]=  CIF1(6,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,7]=  CIF1(7,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,8]=  CIF1(8,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,9]=  CIF1(9,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,10]= CIF1(10, 50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 CIF1_mat[2,1]=  CIF1(1,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,2]=  CIF1(2,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,3]=  CIF1(3,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,4]=  CIF1(4,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,5]=  CIF1(5,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,6]=  CIF1(6,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,7]=  CIF1(7,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,8]=  CIF1(8,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,9]=  CIF1(9,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,10]= CIF1(10, 60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 CIF1_mat[3,1]=  CIF1(1,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,2]=  CIF1(2,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,3]=  CIF1(3,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,4]=  CIF1(4,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,5]=  CIF1(5,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,6]=  CIF1(6,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,7]=  CIF1(7,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,8]=  CIF1(8,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,9]=  CIF1(9,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,10]= CIF1(10, 70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 CIF1_mat[4,1]=  CIF1(1,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,2]=  CIF1(2,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,3]=  CIF1(3,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,4]=  CIF1(4,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,5]=  CIF1(5,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,6]=  CIF1(6,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,7]=  CIF1(7,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,8]=  CIF1(8,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,9]=  CIF1(9,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,10]= CIF1(10, 80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 CIF1_mat[5,1]=  CIF1(1,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,2]=  CIF1(2,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,3]=  CIF1(3,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,4]=  CIF1(4,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,5]=  CIF1(5,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,6]=  CIF1(6,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,7]=  CIF1(7,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,8]=  CIF1(8,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,9]=  CIF1(9,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,10]= CIF1(10, 90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 st_matrix("CIF1_mat", CIF1_mat)
 
  //Create the matrices that contain the true CIF2 for 50,60,70,80,90 years of age at diagnosis
 CIF2_mat[1,1]=  CIF2(1,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,2]=  CIF2(2,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,3]=  CIF2(3,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,4]=  CIF2(4,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,5]=  CIF2(5,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,6]=  CIF2(6,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,7]=  CIF2(7,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,8]=  CIF2(8,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,9]=  CIF2(9,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,10]= CIF2(10, 50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 CIF2_mat[2,1]=  CIF2(1,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,2]=  CIF2(2,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,3]=  CIF2(3,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,4]=  CIF2(4,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,5]=  CIF2(5,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,6]=  CIF2(6,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,7]=  CIF2(7,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,8]=  CIF2(8,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,9]=  CIF2(9,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,10]= CIF2(10, 60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 CIF2_mat[3,1]=  CIF2(1,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,2]=  CIF2(2,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,3]=  CIF2(3,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,4]=  CIF2(4,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,5]=  CIF2(5,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,6]=  CIF2(6,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,7]=  CIF2(7,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,8]=  CIF2(8,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,9]=  CIF2(9,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,10]= CIF2(10, 70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 CIF2_mat[4,1]=  CIF2(1,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,2]=  CIF2(2,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,3]=  CIF2(3,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,4]=  CIF2(4,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,5]=  CIF2(5,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,6]=  CIF2(6,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,7]=  CIF2(7,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,8]=  CIF2(8,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,9]=  CIF2(9,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,10]= CIF2(10, 80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 CIF2_mat[5,1]=  CIF2(1,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,2]=  CIF2(2,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,3]=  CIF2(3,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,4]=  CIF2(4,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,5]=  CIF2(5,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,6]=  CIF2(6,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,7]=  CIF2(7,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,8]=  CIF2(8,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,9]=  CIF2(9,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,10]= CIF2(10, 90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 st_matrix("CIF2_mat", CIF2_mat)

 
end

mat li CIF1_mat
mat li CIF2_mat

cd "$N\2b.Truths_folder" 

svmat double CIF1_mat
svmat double CIF2_mat

 //matsave CIF1_mat, replace p("CIF1_mat") saving replace
 //matsave CIF2_mat, replace p("CIF2_mat") saving replace
gen scenario=12

rename CIF1_mat1  CIF1_true_f1_at1
rename CIF1_mat2  CIF1_true_f1_at2
rename CIF1_mat3  CIF1_true_f1_at3
rename CIF1_mat4  CIF1_true_f1_at4
rename CIF1_mat5  CIF1_true_f1_at5
rename CIF1_mat6  CIF1_true_f1_at6
rename CIF1_mat7  CIF1_true_f1_at7
rename CIF1_mat8  CIF1_true_f1_at8
rename CIF1_mat9  CIF1_true_f1_at9
rename CIF1_mat10  CIF1_true_f1_at10


rename CIF2_mat1  CIF2_true_f1_at1
rename CIF2_mat2  CIF2_true_f1_at2
rename CIF2_mat3  CIF2_true_f1_at3
rename CIF2_mat4  CIF2_true_f1_at4
rename CIF2_mat5  CIF2_true_f1_at5
rename CIF2_mat6  CIF2_true_f1_at6
rename CIF2_mat7  CIF2_true_f1_at7
rename CIF2_mat8  CIF2_true_f1_at8
rename CIF2_mat9  CIF2_true_f1_at9
rename CIF2_mat10  CIF2_true_f1_at10


save "$N\2b.Truths_folder/True_CIFs_scenario12", replace

rename CIF1_mat11 givenage
drop CIF2_mat11

gen id=_n

//Reshape true CIF values from wide to long
reshape long CIF1_true_f1_at CIF2_true_f1_at , i(id) j(year)


//Save true CIF values as datasets
save "$N\2b.Truths_folder/True_CIFs_scenario12", replace

use "$N\2b.Truths_folder/True_CIFs_scenario12", replace


/***************************************************************/



 
clear all  
include "$N/1b.Scenarios_folder/Scenario22_new.do"

 
mata:

CIF1_mat=J(5,11,.)
CIF2_mat=J(5,11,.)

CIF1_mat[1,11]=50
CIF1_mat[2,11]=60
CIF1_mat[3,11]=70
CIF1_mat[4,11]=80
CIF1_mat[5,11]=90

CIF2_mat[1,11]=50
CIF2_mat[2,11]=60
CIF2_mat[3,11]=70
CIF2_mat[4,11]=80
CIF2_mat[5,11]=90

alpha = strtoreal(st_global("alpha"))
lambda = strtoreal(st_global("lambda"))
gamma = strtoreal(st_global("gamma"))


lambda1 = strtoreal(st_global("lambda1"))
gamma1 = strtoreal(st_global("gamma1"))
lambda2 = strtoreal(st_global("lambda2"))
gamma2 = strtoreal(st_global("gamma2"))
pmix = strtoreal(st_global("pmix"))

bsex_c2 = strtoreal(st_global("bsex_c2"))
bsex_c2_tvc_1 = strtoreal(st_global("bsex_c2_tvc_1"))
bsex_c2_tvc_2 = strtoreal(st_global("bsex_c2_tvc_2"))
bsex_c1 = strtoreal(st_global("bsex_c1"))
b1 = strtoreal(st_global("b1"))
b2 = strtoreal(st_global("b2"))

// gq() - function to calculae gauss-legendre nodes and weights

//Define function not giving direct results unless you ask for them by typing them at the end of the function
// weights matrix, nodes matrix, scalar n
void gq(real matrix weights, real matrix nodes, real scalar n)
{

//Define objects like C: il, j , muzero will be scalars while a,A,vec will be matrices
real scalar i1, j, muzero, b
real matrix a, A, vec 

//il will range from 1 to n-1 by a step of 1
//muzero must be a starting value of sorts
//a is a zero vector (1 row, n-1 columns)  
//b is a simple function of il. It will be a vector of n-1 columns. weight function
//A will be a diagonal matrix with the diagonal elements be those of the vector a (a.k.a zeros) at the start.
//However with the loop, A's j diagonal element will be the b[j] result of funtion b. 
 
        i1 = range(1,n-1,1)'
        muzero = 2
        a = J(1,n,0)
        b = i1:/sqrt(4 :* i1:^2 :- 1)   
        A= diag(a)
        for(j=1;j<=n-1;j++){
                A[j,j+1] = b[j]
                A[j+1,j] = b[j]
        }       
		
//pragma unset X suppresses the warning message note: variable X may be used before set

                                pragma unset vec
								
//eigensystem(A, X, L, rcond, nobalance) calculates eigenvectors and eigenvalues of a general, real or complex, square matrix A. Eigenvectors are returned in X and eigenvalues in L. 
//in symeigensystem the difference is that A is Hermitian-symmetric
                                symeigensystem(A,vec,nodes)
								
//weights is a vector. It is the result of the product of the vector vec^2 with muzero and then it gets ordered by nodes matrix???
//Nodes is a matrix of the given nodes ordered by them?


        weights = (vec[1,]:^2:*muzero)'
        weights = weights[order(nodes',1)]
        nodes = nodes'[order(nodes',1)']
 
}
 
// function to derive hazard function for other cause mortality
real matrix hazard2(t,age,alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2) 
{
    return(
	((alpha :+lambda:*exp(gamma:*(age:+t))):*exp(bsex_c2 :+ bsex_c2_tvc_1:*(age:+t) :+ bsex_c2_tvc_2:*(age:+t):^2)  ) )
	
}


hazard2((1..10),80,alpha,lambda,gamma,bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2)

// function to derive survival function for other cause mortality
real matrix survival2(t, age, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2, n) 
{
	gq(wt=.,nodes=.,n)
    enter = 0
    nodes_t = 0.5:*((age:+t):-enter):*nodes' :+ 0.5:*((age:+t):+enter)
    
    ht_nodes = hazard2(nodes_t, 0, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2 )

    return(exp(-0.5:*((age:+t):-enter):*quadrowsum(wt':*ht_nodes)))
}

survival2(80,0,alpha,lambda,gamma,bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,70)

survival2(10,80,alpha,lambda,gamma,bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,70)

// function to derive survival function for cancer mortality
real matrix survival1(t, age, lambda1, gamma1, lambda2, gamma2, pmix, bsex_c1, b1, b2) 
{
	return((pmix :*exp(-lambda1:*t:^gamma1):+(1:-pmix):*exp(-lambda2:*t:^gamma2 )):^ 
	        exp((b1:*(age:-65):+(b2:*(age:-65):^2):+bsex_c1)))

}

survival1(10, 80, lambda1, gamma1, lambda2, gamma2, pmix, bsex_c1, b1, b2) 

// function to derive hazard function for cancer mortality
real matrix hazard1(t, age, lambda1, gamma1, lambda2, gamma2, pmix, bsex_c1, b1, b2) 
{
	return((((lambda1:*gamma1:*t:^(gamma1:-1):*pmix:*exp(-lambda1:*t:^gamma1)):+(lambda2:*gamma2:*t:^(gamma2:-1)):* 
	       (1:-pmix):*exp(-lambda2:*t:^gamma2)):/(pmix :*exp(-lambda1 :*t:^gamma1):+ 
		   (1:-pmix):*exp(-lambda2:*t:^gamma2))):*exp((b1:*(age:-65):+ b2:*((age:-65):^2):+ bsex_c1) ))

}

hazard1(10, 80, lambda1, gamma1, lambda2, gamma2, pmix, bsex_c1, b1, b2) 


// function to derive CIF for cancer mortality (CIF1)
real matrix CIF1(t, age, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, n) {
	gq(wt=.,nodes=.,n)
	enter = 0
	nodes_t = 0.5:*(t:-enter):*nodes' :+ 0.5:*(t:+enter)
 
	nodes_t
// create empty matrix
	CIF1t_nodes = J(1,n,.)
    CIF1t_nodes            
// loop over nodes
// evaluation of each node requires numerical integration
	for(i=1;i<=n;i++) {
		tn = nodes_t[1,i]
		CIF1t_nodes[1,i] =  survival1(tn, age, lambda1, gamma1, lambda2, gamma2, pmix, bsex_c1, b1, b2) :*
							survival2(tn, age, alpha,lambda,gamma,bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,n) :*
							hazard1(tn, age, lambda1, gamma1, lambda2, gamma2, pmix, bsex_c1, b1, b2) :/
							survival2(age, 0, alpha,lambda,gamma,bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,n)
	}
    return(0.5:*(t:-enter):*quadrowsum(wt':*CIF1t_nodes))
} 
 
 CIF1(10, 80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)

  // function to derive CIF for other cause mortality (CIF2)
real matrix CIF2(t, age, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, n) {
	gq(wt=.,nodes=.,n)
	enter = 0
	nodes_t = 0.5:*(t:-enter):*nodes' :+ 0.5:*(t:+enter)
 
	nodes_t
// create empty matrix
	CIF2t_nodes = J(1,n,.)
    CIF2t_nodes            
// loop over nodes
// evaluation of each node requires numerical integration
	for(i=1;i<=n;i++) {
		tn = nodes_t[1,i]
		CIF2t_nodes[1,i] =  survival1(tn, age, lambda1, gamma1, lambda2, gamma2, pmix, bsex_c1, b1, b2) :*
							survival2(tn, age, alpha,lambda,gamma,bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,n) :*
							hazard2(  tn, age, alpha,lambda,gamma,bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2) :/
							survival2(age,  0, alpha,lambda,gamma,bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,n)
	}
    return(0.5:*(t:-enter):*quadrowsum(wt':*CIF2t_nodes))
} 
 
 CIF2(10,80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)

   //Create the matrices that contain the true CIF1 for 50,60,70,80,90 years of age at diagnosis
 CIF1_mat[1,1]=  CIF1(1,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,2]=  CIF1(2,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,3]=  CIF1(3,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,4]=  CIF1(4,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,5]=  CIF1(5,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,6]=  CIF1(6,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,7]=  CIF1(7,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,8]=  CIF1(8,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,9]=  CIF1(9,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,10]= CIF1(10, 50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 CIF1_mat[2,1]=  CIF1(1,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,2]=  CIF1(2,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,3]=  CIF1(3,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,4]=  CIF1(4,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,5]=  CIF1(5,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,6]=  CIF1(6,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,7]=  CIF1(7,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,8]=  CIF1(8,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,9]=  CIF1(9,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,10]= CIF1(10, 60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 CIF1_mat[3,1]=  CIF1(1,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,2]=  CIF1(2,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,3]=  CIF1(3,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,4]=  CIF1(4,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,5]=  CIF1(5,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,6]=  CIF1(6,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,7]=  CIF1(7,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,8]=  CIF1(8,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,9]=  CIF1(9,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,10]= CIF1(10, 70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 CIF1_mat[4,1]=  CIF1(1,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,2]=  CIF1(2,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,3]=  CIF1(3,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,4]=  CIF1(4,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,5]=  CIF1(5,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,6]=  CIF1(6,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,7]=  CIF1(7,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,8]=  CIF1(8,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,9]=  CIF1(9,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,10]= CIF1(10, 80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 CIF1_mat[5,1]=  CIF1(1,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,2]=  CIF1(2,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,3]=  CIF1(3,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,4]=  CIF1(4,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,5]=  CIF1(5,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,6]=  CIF1(6,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,7]=  CIF1(7,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,8]=  CIF1(8,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,9]=  CIF1(9,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,10]= CIF1(10, 90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 st_matrix("CIF1_mat", CIF1_mat)
 
  //Create the matrices that contain the true CIF2 for 50,60,70,80,90 years of age at diagnosis
 CIF2_mat[1,1]=  CIF2(1,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,2]=  CIF2(2,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,3]=  CIF2(3,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,4]=  CIF2(4,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,5]=  CIF2(5,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,6]=  CIF2(6,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,7]=  CIF2(7,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,8]=  CIF2(8,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,9]=  CIF2(9,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,10]= CIF2(10, 50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 CIF2_mat[2,1]=  CIF2(1,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,2]=  CIF2(2,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,3]=  CIF2(3,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,4]=  CIF2(4,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,5]=  CIF2(5,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,6]=  CIF2(6,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,7]=  CIF2(7,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,8]=  CIF2(8,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,9]=  CIF2(9,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,10]= CIF2(10, 60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 CIF2_mat[3,1]=  CIF2(1,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,2]=  CIF2(2,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,3]=  CIF2(3,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,4]=  CIF2(4,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,5]=  CIF2(5,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,6]=  CIF2(6,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,7]=  CIF2(7,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,8]=  CIF2(8,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,9]=  CIF2(9,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,10]= CIF2(10, 70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 CIF2_mat[4,1]=  CIF2(1,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,2]=  CIF2(2,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,3]=  CIF2(3,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,4]=  CIF2(4,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,5]=  CIF2(5,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,6]=  CIF2(6,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,7]=  CIF2(7,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,8]=  CIF2(8,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,9]=  CIF2(9,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,10]= CIF2(10, 80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 CIF2_mat[5,1]=  CIF2(1,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,2]=  CIF2(2,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,3]=  CIF2(3,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,4]=  CIF2(4,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,5]=  CIF2(5,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,6]=  CIF2(6,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,7]=  CIF2(7,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,8]=  CIF2(8,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,9]=  CIF2(9,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,10]= CIF2(10, 90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 st_matrix("CIF2_mat", CIF2_mat)

 
end

mat li CIF1_mat
mat li CIF2_mat

cd "$N\2b.Truths_folder" 

svmat double CIF1_mat
svmat double CIF2_mat

 //matsave CIF1_mat, replace p("CIF1_mat") saving replace
 //matsave CIF2_mat, replace p("CIF2_mat") saving replace
gen scenario=22

rename CIF1_mat1  CIF1_true_f1_at1
rename CIF1_mat2  CIF1_true_f1_at2
rename CIF1_mat3  CIF1_true_f1_at3
rename CIF1_mat4  CIF1_true_f1_at4
rename CIF1_mat5  CIF1_true_f1_at5
rename CIF1_mat6  CIF1_true_f1_at6
rename CIF1_mat7  CIF1_true_f1_at7
rename CIF1_mat8  CIF1_true_f1_at8
rename CIF1_mat9  CIF1_true_f1_at9
rename CIF1_mat10  CIF1_true_f1_at10


rename CIF2_mat1  CIF2_true_f1_at1
rename CIF2_mat2  CIF2_true_f1_at2
rename CIF2_mat3  CIF2_true_f1_at3
rename CIF2_mat4  CIF2_true_f1_at4
rename CIF2_mat5  CIF2_true_f1_at5
rename CIF2_mat6  CIF2_true_f1_at6
rename CIF2_mat7  CIF2_true_f1_at7
rename CIF2_mat8  CIF2_true_f1_at8
rename CIF2_mat9  CIF2_true_f1_at9
rename CIF2_mat10  CIF2_true_f1_at10


save "$N\2b.Truths_folder/True_CIFs_scenario22", replace

rename CIF1_mat11 givenage
drop CIF2_mat11

gen id=_n

//Reshape true CIF values from wide to long
reshape long CIF1_true_f1_at CIF2_true_f1_at , i(id) j(year)

//Save true CIF values as datasets
save "$N\2b.Truths_folder/True_CIFs_scenario22", replace


/*********************************************************************/



 
clear all  
include "$N/1b.Scenarios_folder/Scenario24_new.do"

 
mata:

CIF1_mat=J(5,11,.)
CIF2_mat=J(5,11,.)

CIF1_mat[1,11]=50
CIF1_mat[2,11]=60
CIF1_mat[3,11]=70
CIF1_mat[4,11]=80
CIF1_mat[5,11]=90

CIF2_mat[1,11]=50
CIF2_mat[2,11]=60
CIF2_mat[3,11]=70
CIF2_mat[4,11]=80
CIF2_mat[5,11]=90

alpha = strtoreal(st_global("alpha"))
lambda = strtoreal(st_global("lambda"))
gamma = strtoreal(st_global("gamma"))


lambda1 = strtoreal(st_global("lambda1"))
gamma1 = strtoreal(st_global("gamma1"))
lambda2 = strtoreal(st_global("lambda2"))
gamma2 = strtoreal(st_global("gamma2"))
pmix = strtoreal(st_global("pmix"))

bsex_c2 = strtoreal(st_global("bsex_c2"))
bsex_c2_tvc_1 = strtoreal(st_global("bsex_c2_tvc_1"))
bsex_c2_tvc_2 = strtoreal(st_global("bsex_c2_tvc_2"))
bsex_c1 = strtoreal(st_global("bsex_c1"))
b1 = strtoreal(st_global("b1"))
b2 = strtoreal(st_global("b2"))

// gq() - function to calculae gauss-legendre nodes and weights

//Define function not giving direct results unless you ask for them by typing them at the end of the function
// weights matrix, nodes matrix, scalar n
void gq(real matrix weights, real matrix nodes, real scalar n)
{

//Define objects like C: il, j , muzero will be scalars while a,A,vec will be matrices
real scalar i1, j, muzero, b
real matrix a, A, vec 

//il will range from 1 to n-1 by a step of 1
//muzero must be a starting value of sorts
//a is a zero vector (1 row, n-1 columns)  
//b is a simple function of il. It will be a vector of n-1 columns. weight function
//A will be a diagonal matrix with the diagonal elements be those of the vector a (a.k.a zeros) at the start.
//However with the loop, A's j diagonal element will be the b[j] result of funtion b. 
 
        i1 = range(1,n-1,1)'
        muzero = 2
        a = J(1,n,0)
        b = i1:/sqrt(4 :* i1:^2 :- 1)   
        A= diag(a)
        for(j=1;j<=n-1;j++){
                A[j,j+1] = b[j]
                A[j+1,j] = b[j]
        }       
		
//pragma unset X suppresses the warning message note: variable X may be used before set

                                pragma unset vec
								
//eigensystem(A, X, L, rcond, nobalance) calculates eigenvectors and eigenvalues of a general, real or complex, square matrix A. Eigenvectors are returned in X and eigenvalues in L. 
//in symeigensystem the difference is that A is Hermitian-symmetric
                                symeigensystem(A,vec,nodes)
								
//weights is a vector. It is the result of the product of the vector vec^2 with muzero and then it gets ordered by nodes matrix???
//Nodes is a matrix of the given nodes ordered by them?


        weights = (vec[1,]:^2:*muzero)'
        weights = weights[order(nodes',1)]
        nodes = nodes'[order(nodes',1)']
 
}
 
// function to derive hazard function for other cause mortality
real matrix hazard2(t,age,alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2) 
{
    return(
	((alpha :+lambda:*exp(gamma:*(age:+t))):*exp(bsex_c2 :+ bsex_c2_tvc_1:*(age:+t) :+ bsex_c2_tvc_2:*(age:+t):^2)  ) )
	
}


hazard2((1..10),80,alpha,lambda,gamma,bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2)

// function to derive survival function for other cause mortality
real matrix survival2(t, age, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2, n) 
{
	gq(wt=.,nodes=.,n)
    enter = 0
    nodes_t = 0.5:*((age:+t):-enter):*nodes' :+ 0.5:*((age:+t):+enter)
    
    ht_nodes = hazard2(nodes_t, 0, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2 )

    return(exp(-0.5:*((age:+t):-enter):*quadrowsum(wt':*ht_nodes)))
}

survival2(80,0,alpha,lambda,gamma,bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,70)

survival2(10,80,alpha,lambda,gamma,bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,70)

// function to derive survival function for cancer mortality
real matrix survival1(t, age, lambda1, gamma1, lambda2, gamma2, pmix, bsex_c1, b1, b2) 
{
	return((pmix :*exp(-lambda1:*t:^gamma1):+(1:-pmix):*exp(-lambda2:*t:^gamma2 )):^ 
	        exp((b1:*(age:-65):+(b2:*(age:-65):^2):+bsex_c1)))

}

survival1(10, 80, lambda1, gamma1, lambda2, gamma2, pmix, bsex_c1, b1, b2) 

// function to derive hazard function for cancer mortality
real matrix hazard1(t, age, lambda1, gamma1, lambda2, gamma2, pmix, bsex_c1, b1, b2) 
{
	return((((lambda1:*gamma1:*t:^(gamma1:-1):*pmix:*exp(-lambda1:*t:^gamma1)):+(lambda2:*gamma2:*t:^(gamma2:-1)):* 
	       (1:-pmix):*exp(-lambda2:*t:^gamma2)):/(pmix :*exp(-lambda1 :*t:^gamma1):+ 
		   (1:-pmix):*exp(-lambda2:*t:^gamma2))):*exp((b1:*(age:-65):+ b2:*((age:-65):^2):+ bsex_c1) ))

}

hazard1(10, 80, lambda1, gamma1, lambda2, gamma2, pmix, bsex_c1, b1, b2) 


// function to derive CIF for cancer mortality (CIF1)
real matrix CIF1(t, age, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, n) {
	gq(wt=.,nodes=.,n)
	enter = 0
	nodes_t = 0.5:*(t:-enter):*nodes' :+ 0.5:*(t:+enter)
 
	nodes_t
// create empty matrix
	CIF1t_nodes = J(1,n,.)
    CIF1t_nodes            
// loop over nodes
// evaluation of each node requires numerical integration
	for(i=1;i<=n;i++) {
		tn = nodes_t[1,i]
		CIF1t_nodes[1,i] =  survival1(tn, age, lambda1, gamma1, lambda2, gamma2, pmix, bsex_c1, b1, b2) :*
							survival2(tn, age, alpha,lambda,gamma,bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,n) :*
							hazard1(tn, age, lambda1, gamma1, lambda2, gamma2, pmix, bsex_c1, b1, b2) :/
							survival2(age, 0, alpha,lambda,gamma,bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,n)
	}
    return(0.5:*(t:-enter):*quadrowsum(wt':*CIF1t_nodes))
} 
 
 CIF1(10, 80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)

  // function to derive CIF for other cause mortality (CIF2)
real matrix CIF2(t, age, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, n) {
	gq(wt=.,nodes=.,n)
	enter = 0
	nodes_t = 0.5:*(t:-enter):*nodes' :+ 0.5:*(t:+enter)
 
	nodes_t
// create empty matrix
	CIF2t_nodes = J(1,n,.)
    CIF2t_nodes            
// loop over nodes
// evaluation of each node requires numerical integration
	for(i=1;i<=n;i++) {
		tn = nodes_t[1,i]
		CIF2t_nodes[1,i] =  survival1(tn, age, lambda1, gamma1, lambda2, gamma2, pmix, bsex_c1, b1, b2) :*
							survival2(tn, age, alpha,lambda,gamma,bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,n) :*
							hazard2(  tn, age, alpha,lambda,gamma,bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2) :/
							survival2(age,  0, alpha,lambda,gamma,bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,n)
	}
    return(0.5:*(t:-enter):*quadrowsum(wt':*CIF2t_nodes))
} 
 
 CIF2(10,80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)

   //Create the matrices that contain the true CIF1 for 50,60,70,80,90 years of age at diagnosis
 CIF1_mat[1,1]=  CIF1(1,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,2]=  CIF1(2,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,3]=  CIF1(3,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,4]=  CIF1(4,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,5]=  CIF1(5,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,6]=  CIF1(6,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,7]=  CIF1(7,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,8]=  CIF1(8,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,9]=  CIF1(9,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[1,10]= CIF1(10, 50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 CIF1_mat[2,1]=  CIF1(1,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,2]=  CIF1(2,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,3]=  CIF1(3,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,4]=  CIF1(4,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,5]=  CIF1(5,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,6]=  CIF1(6,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,7]=  CIF1(7,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,8]=  CIF1(8,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,9]=  CIF1(9,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[2,10]= CIF1(10, 60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 CIF1_mat[3,1]=  CIF1(1,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,2]=  CIF1(2,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,3]=  CIF1(3,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,4]=  CIF1(4,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,5]=  CIF1(5,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,6]=  CIF1(6,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,7]=  CIF1(7,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,8]=  CIF1(8,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,9]=  CIF1(9,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[3,10]= CIF1(10, 70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 CIF1_mat[4,1]=  CIF1(1,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,2]=  CIF1(2,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,3]=  CIF1(3,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,4]=  CIF1(4,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,5]=  CIF1(5,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,6]=  CIF1(6,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,7]=  CIF1(7,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,8]=  CIF1(8,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,9]=  CIF1(9,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[4,10]= CIF1(10, 80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 CIF1_mat[5,1]=  CIF1(1,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,2]=  CIF1(2,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,3]=  CIF1(3,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,4]=  CIF1(4,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,5]=  CIF1(5,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,6]=  CIF1(6,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,7]=  CIF1(7,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,8]=  CIF1(8,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,9]=  CIF1(9,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF1_mat[5,10]= CIF1(10, 90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 st_matrix("CIF1_mat", CIF1_mat)
 
  //Create the matrices that contain the true CIF2 for 50,60,70,80,90 years of age at diagnosis
 CIF2_mat[1,1]=  CIF2(1,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,2]=  CIF2(2,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,3]=  CIF2(3,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,4]=  CIF2(4,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,5]=  CIF2(5,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,6]=  CIF2(6,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,7]=  CIF2(7,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,8]=  CIF2(8,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,9]=  CIF2(9,  50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[1,10]= CIF2(10, 50, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 CIF2_mat[2,1]=  CIF2(1,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,2]=  CIF2(2,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,3]=  CIF2(3,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,4]=  CIF2(4,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,5]=  CIF2(5,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,6]=  CIF2(6,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,7]=  CIF2(7,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,8]=  CIF2(8,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,9]=  CIF2(9,  60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[2,10]= CIF2(10, 60, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 CIF2_mat[3,1]=  CIF2(1,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,2]=  CIF2(2,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,3]=  CIF2(3,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,4]=  CIF2(4,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,5]=  CIF2(5,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,6]=  CIF2(6,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,7]=  CIF2(7,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,8]=  CIF2(8,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,9]=  CIF2(9,  70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[3,10]= CIF2(10, 70, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 CIF2_mat[4,1]=  CIF2(1,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,2]=  CIF2(2,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,3]=  CIF2(3,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,4]=  CIF2(4,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,5]=  CIF2(5,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,6]=  CIF2(6,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,7]=  CIF2(7,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,8]=  CIF2(8,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,9]=  CIF2(9,  80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[4,10]= CIF2(10, 80, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 CIF2_mat[5,1]=  CIF2(1,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,2]=  CIF2(2,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,3]=  CIF2(3,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,4]=  CIF2(4,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,5]=  CIF2(5,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,6]=  CIF2(6,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,7]=  CIF2(7,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,8]=  CIF2(8,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,9]=  CIF2(9,  90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 CIF2_mat[5,10]= CIF2(10, 90, alpha,lambda,gamma, bsex_c2, bsex_c2_tvc_1, bsex_c2_tvc_2,lambda1,gamma1,lambda2,gamma2,pmix, bsex_c1,b1, b2, 100)
 
 st_matrix("CIF2_mat", CIF2_mat)

 
end

mat li CIF1_mat
mat li CIF2_mat

cd "$N\2b.Truths_folder" 

svmat double CIF1_mat
svmat double CIF2_mat

 //matsave CIF1_mat, replace p("CIF1_mat") saving replace
 //matsave CIF2_mat, replace p("CIF2_mat") saving replace
gen scenario=24

rename CIF1_mat1  CIF1_true_f1_at1
rename CIF1_mat2  CIF1_true_f1_at2
rename CIF1_mat3  CIF1_true_f1_at3
rename CIF1_mat4  CIF1_true_f1_at4
rename CIF1_mat5  CIF1_true_f1_at5
rename CIF1_mat6  CIF1_true_f1_at6
rename CIF1_mat7  CIF1_true_f1_at7
rename CIF1_mat8  CIF1_true_f1_at8
rename CIF1_mat9  CIF1_true_f1_at9
rename CIF1_mat10  CIF1_true_f1_at10


rename CIF2_mat1  CIF2_true_f1_at1
rename CIF2_mat2  CIF2_true_f1_at2
rename CIF2_mat3  CIF2_true_f1_at3
rename CIF2_mat4  CIF2_true_f1_at4
rename CIF2_mat5  CIF2_true_f1_at5
rename CIF2_mat6  CIF2_true_f1_at6
rename CIF2_mat7  CIF2_true_f1_at7
rename CIF2_mat8  CIF2_true_f1_at8
rename CIF2_mat9  CIF2_true_f1_at9
rename CIF2_mat10  CIF2_true_f1_at10


save "$N\2b.Truths_folder/True_CIFs_scenario24", replace

rename CIF1_mat11 givenage
drop CIF2_mat11

gen id=_n

//Reshape true CIF values from wide to long
reshape long CIF1_true_f1_at CIF2_true_f1_at , i(id) j(year)


//Save true CIF values as datasets
save "$N\2b.Truths_folder/True_CIFs_scenario24", replace


