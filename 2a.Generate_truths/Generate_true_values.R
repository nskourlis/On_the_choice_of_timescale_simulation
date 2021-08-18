
#Set times
times=c(1,2,3,4,5,6,7,8,9,10)

#True parameters for scenario 1
 obs           =  c(1000         )
 lambda1       =  c(.71312397193535    )
 gamma1        =  c(.7665924306563239)
 lambda2       =  c(.0075545896481545)
 gamma2        =  c(.7919159630956031)
 pmix          =  c(.2816078681168775)
 lambda        =  c(1.00204431863e-24)
 gamma         =  c(12.2742061562415)
 lambda        =  c(1.00204431863e-24)
 gamma         =  c(12.2742061562415)
 b1            =  c(-.0030723147292845)
 b2            =  c(.0001359907528383)
 c1            =  c(-.00143)
 c2            =  c(.3429)
 c3            =  c(-18)
 alpha         =  c(.002)
 bsex_c1       =  c(-.0512932943875506)
 bsex_c2       =  c(-.3566749439387324)
 bsex_c2_tvc_1 =  c(0)
 bsex_c2_tvc_2 =  c(0)
 minage        =  c(20)
 mintime_cancer=  c(0.0027)
 maxtime_cancer=  c(10)
 independence  =  c(0)
 prop_sex_other=  c("PH_sex_other")
 base_function =  c("Weibull")
 crit_age      =  c(70)
 tt1=matrix(NA,nrow =1,ncol = 10)
 tt1[1,1]=1;  tt1[1,2]=2;   tt1[1,3]=3;   tt1[1,4]=4;  tt1[1,5]=5;  tt1[1,6]=6;  tt1[1,7]=7;  tt1[1,8]=8;  tt1[1,9]=9;  tt1[1,10]=10
 ageloop=c(50,60,70,80,90)
 female=1
 
 h1_f1=array(NA,dim=c(10,5))
 h1_f0=array(NA,dim=c(10,5))
 s1_f1=array(NA,dim=c(10,5))
 s1_f0=array(NA,dim=c(10,5))
 
 h2_f1=array(NA,dim=c(10,5))
 h2_f0=array(NA,dim=c(10,5))
 s2_f1=array(NA,dim=c(10,5))
 s2_f0=array(NA,dim=c(10,5))
 s2_f1_age=array(NA,dim=c(10,5))
 s2_f0_age=array(NA,dim=c(10,5))
 s2c_f0=array(NA,dim=c(10,5))
 s2c_f1=array(NA,dim=c(10,5))
 CIF_f0=array(NA,dim=c(10,5))
 CIF1_f1=matrix(NA,nrow = ncol(tt1),ncol = 5)
 CIF2_f1=matrix(NA,nrow = ncol(tt1),ncol = 5)
 
 
 for (i in 1:5 ) {
   
   for (j in 1:10 ) {

 h1_f1[j,i]=(((lambda1[1] *gamma1[1]* tt1[1,j]^( gamma1[1] -1 )*pmix[1] * exp(-lambda1[1]  *tt1[1,j]^(gamma1[1] ) ) )+ 
                (lambda2[1]* gamma2[1]* (tt1[1,j]^(gamma2[1]-1) )* 
        (1-pmix[1] )*exp(-lambda2[1] *tt1[1,j]^( gamma2[1]) ))) /(pmix[1] *exp(- lambda1[1] *tt1[1,j]^( gamma1[1] ) ) + 
        (1-pmix[1] )*exp(-lambda2[1] *tt1[,j]^( gamma2[1] ) ))) *exp((b1[1]*ageloop[i]+ (b2[1]*ageloop[i]^2) +  bsex_c1[1]*female ) )  

 
  h1_f0[j,i]=(((lambda1[1] *gamma1[1]* tt1[1,j]^( gamma1[1] -1 )*pmix[1] * exp(-lambda1[1]  *tt1[1,j]^(gamma1[1] ) ) )+ 
                (lambda2[1]* gamma2[1]* (tt1[j]^(gamma2[1]-1) )* 
                (1-pmix[1] )*exp(-lambda2[1] *tt1[1,j]^( gamma2[1]) ))) /(pmix[1] *exp(- lambda1[1] *tt1[1,j]^( gamma1[1] ) ) + 
                (1-pmix[1] )*exp(-lambda2[1] *tt1[1,j]^( gamma2[1] ) ))) *exp((b1[1]*ageloop[i]+ (b2[1]*ageloop[i]^2)) ) 
 
  s1_f1[j,i]=(pmix[1] *exp(-lambda1[1] *tt1[1,j]^gamma1[1] )+(1-pmix[1] )*exp(-lambda2[1] *tt1[j]^gamma2[1] )) ^exp((b1[1]*ageloop[i]+ (b2[1]*ageloop[i]^2) +   bsex_c1[1]*female) ) 
  
  s1_f0[j,i]=(pmix[1] *exp(-lambda1[1] *tt1[1,j]^(gamma1[1]) )+(1-pmix[1] )*exp(-lambda2[1] *tt1[j]^gamma2[1] )) ^exp((b1[1]*ageloop[i]+ (b2[1]*ageloop[i]^2)) ) 
   
  if (base_function[1]=="Weibull") {
    
    h2_f1[j,i]=(alpha[1]+ lambda[1] *gamma[1] *(ageloop[1]+tt1[1,j])^(gamma[1] -1)) *exp(bsex_c2[1]*female + bsex_c2_tvc_1[1]*female*(ageloop[1]+tt1[1,j]) + bsex_c2_tvc_2[1]*female*(ageloop[i]+tt1[1,j])^2  )
    h2_f0[j,i]=(alpha[1]+ lambda[1] *gamma[1] *(ageloop[1]+tt1[1,j])^(gamma[1] -1)) 


    s2_f0[j,i]=(alpha[1]+ lambda[1] *(ageloop[1]+tt1[j])^(gamma[1]) )
    s2_f0_age[j,i]=(alpha[1]+ lambda[1] *(ageloop[1])^(gamma[1]) )
    s2c_f0[j,i]=s2_f0[j,i]/s2_f0_age[j,i]
    
    InnerFunc_s2c_f0= function(s){(alpha[1]+ lambda[1] *(ageloop[1]+s)^(gamma[1]) )/(alpha[1]+ lambda[1] *(ageloop[i])^(gamma[1]) )   }
    
    CIF_f0[j,i]=integrate(InnerFunc_s2c_f0 , 0,tt1[1,j] )$value
    
   }
   
   }
   
 }
 
 
 for (i in 1:5 ) {

 
 InnerFunc_h2= function(x){(alpha[1]+ lambda[1] *gamma[1] *(ageloop[i]+x)^(gamma[1] -1)) *
                        exp(bsex_c2[1]*female + bsex_c2_tvc_1[1]*female*(ageloop[i]+x) + bsex_c2_tvc_2[1]*female*(ageloop[i]+x)^2  ) }
 
 
 InnerFunc_h2_0= function(a0){(alpha[1]+ lambda[1] *gamma[1] *(a0)^(gamma[1] -1)) *
     exp(bsex_c2[1]*female + bsex_c2_tvc_1[1]*female*(a0) + bsex_c2_tvc_2[1]*female*(a0)^2  ) }

  s2_f1[i]=exp(-integrate( InnerFunc_h2_0 , 0,ageloop[i] )$value)
 
 
 InnerIntegral_CIF1=function(y){
   (exp(-(sapply(y, function(z) { integrate(InnerFunc_h2, 0, z)$value } )))/s2_f1[i]) *
     (((lambda1[1] *gamma1[1]*y ^( gamma1[1] -1 )*pmix[1] * exp(-lambda1[1]  *y^(gamma1[1] ) ) )+ 
     (lambda2[1]* gamma2[1]* (y^(gamma2[1]-1) )* 
     (1-pmix[1] )*exp(-lambda2[1] *y^( gamma2[1]) ))) /(pmix[1] *exp(- lambda1[1] *y^( gamma1[1] ) ) + 
     (1-pmix[1] )*exp(-lambda2[1] *y^( gamma2[1] ) ))) *exp((b1[1]*ageloop[i]+ (b2[1]*ageloop[i]^2) +  bsex_c1[1]*female ) )  *
     (pmix[1] *exp(-lambda1[1] *tt1[1,j]^gamma1[1] )+(1-pmix[1] )*exp(-lambda2[1] *y^gamma2[1] )) ^exp((b1[1]*ageloop[i]+ (b2[1]*ageloop[i]^2) +   bsex_c1[1]*female) ) 
 }
 
 
 InnerIntegral_CIF2=function(y){
   (exp(-(sapply(y, function(z) { integrate(InnerFunc_h2, 0, z)$value } )))/s2_f1[i]) *
     (alpha[1]+ lambda[1] *gamma[1] *(ageloop[1]+y)^(gamma[1] -1)) *exp(bsex_c2[1]*female + bsex_c2_tvc_1[1]*female*(ageloop[i]+y) + bsex_c2_tvc_2[1]*female*(ageloop[i]+y)^2  ) *
     (pmix[1] *exp(-lambda1[1] *y^gamma1[1] )+(1-pmix[1] )*exp(-lambda2[1] *y^gamma2[1] )) ^exp((b1[1]*ageloop[i]+ (b2[1]*ageloop[i]^2) +   bsex_c1[1]*female) ) 
 }
 
 #CIF_f1=matrix(NA,nrow = ncol(tt1),ncol = 1)
 
 for(k in 1:ncol(tt1)){
   
   CIF1_f1[k,i]=integrate(InnerIntegral_CIF1 , 0,tt1[1,k] )$value
   CIF2_f1[k,i]=integrate(InnerIntegral_CIF2 , 0,tt1[1,k] )$value
   
 }
 }
 
 CIF1_f1
 CIF2_f1
 