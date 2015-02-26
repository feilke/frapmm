/*
##
##
## Copyright (c) 2015, Martina Feilke
## All rights reserved.
## 
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
## 
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer. 
##     * Redistributions in binary form must reproduce the above
##       copyright notice, this list of conditions and the following
##       disclaimer in the documentation and/or other materials provided
##       with the distribution.
##     * The names of the authors may not be used to endorse or promote
##       products derived from this software without specific prior
##       written permission.
## 
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
## 
##
*/

#include <math.h> 
#include <R.h>
#include <Rmath.h> 
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "zufall.h"

#define min(a1,a2) ((a1) > (a2) ? (a2):(a1))
#define max(a1,a2) ((a1) > (a2) ? (a1):(a2))

// function for tuning
double tune4k1(double what,int acc,int n)
{
	  double newvar = 0.0;
	  if (acc<.3*n){newvar=what*.3;}
      if (acc<.2*n){newvar=what*.1;}
      if (acc<.1*n){newvar=what*.1;}
      if (acc<.02*n){newvar=what*.1;}
      if (acc>.6*n){newvar=what*1.1;}
      if (acc>.7*n){newvar=what*2;}
      if (acc>.8*n){newvar=what*5;}
      if (acc>.9*n){newvar=what*5;}
      if (acc>.99*n){newvar=what*15;}
      if(newvar!=0.0){
      return(newvar);
      }
      else{
      return(what);
      }
}

// main function
void k1_mhgibbs(double* timeall, double* concall, int* nobs,
		int* ttotal, int* J, double* sv, int* size, int* burnin, int* thin, int* tune, int* nobs_cum, double* sigma2b1, double* sigma2beta1,
		double* a0start, double* a1start, double* b1start,
		double* tau2_alpha0, double* tau2_alpha1, double* tau2_beta1, double* sigma2,
		double* a0, double* a1, double* b1, double* alpha0, double* alpha1, double* beta1, int* acc_b1, int* acc_beta1)
{
	
GetRNGstate();

// MH and Gibbs (hybrid algorithm)
	int i;
	double temp;
	int tu_b1;
	int tu_beta1;
	
	double tau2_alpha0_new;
	double tau2_alpha1_new;
	double tau2_beta1_new;
	double sigma2_new;
	
	double a0_new=a0start[0];
	double a1_new=a1start[0];
	double b1_new=b1start[0];
	
	double alpha0j_new[J[0]];
	int i5;
	for(i5=0;i5<J[0];i5++){
		alpha0j_new[i5] = 0.0;
	}
	double alpha1j_new[J[0]];
	int i3;
	for(i3=0;i3<J[0];i3++){
		alpha1j_new[i3] = 0.0;
	}
	double beta1j_new[J[0]];
	int i4;
	for(i4=0;i4<J[0];i4++){
		beta1j_new[i4] = 0.0;
	}
			
	int sample=0;
	int sample_vec=0;
	
	for(i=1; i<size[0]+1; i++) 
	{		
		/*************** GIBBS UPDATES *********************************************************/
		// update of tau2_alpha0
		
		double param5b = sv[2]+(J[0]/2.0); 
				
		double sum5=0;
		int x;
		for(x=0;x<J[0];x++){
			sum5 += pow((alpha0j_new[x]),2);
		}
						
		double param6b = 0.5*sum5+sv[3]; 
				
		double temp5 = RNDGAM(param5b,param6b);	
		tau2_alpha0_new = 1.0/temp5;
				
		/***************************************************************************************/
		// update of tau2_alpha1
		
		double param1b = sv[4]+(J[0]/2.0); 
				
		double sum2=0;
		int k;
		for(k=0;k<J[0];k++){
			sum2 += pow((alpha1j_new[k]),2);
		}
						
		double param2b = 0.5*sum2+sv[5]; 
				
		double temp2 = RNDGAM(param1b,param2b);	
		tau2_alpha1_new = 1.0/temp2;
		
		/***************************************************************************************/
		// update of tau2_beta1
		
		double param1c = sv[6]+(J[0]/2.0);
				
		double sum3=0;
		int l;
		for(l=0;l<J[0];l++){
			sum3 += pow((beta1j_new[l]),2);
		}

		double param2c = 0.5*sum3+sv[7]; 
				
		double temp3 = RNDGAM(param1c,param2c);	
		tau2_beta1_new = 1.0/temp3;
		
		/***************************************************************************************/
		// update of sigma2
		
		double param1d = sv[0] + (ttotal[0]/2.0); 
		
		double sum4=0;
		double sum4help=0;
		int y;
		int z;
		for(y=0;y<J[0];y++){ 
			for(z=nobs_cum[y];z<nobs_cum[y+1];z++){ 
				sum4help = pow((concall[z] - ((1-a0_new-alpha0j_new[y])-(a1_new+alpha1j_new[y])*exp(-(exp(beta1j_new[y]))*(exp(b1_new))*(timeall[z])))),2);
				sum4 = sum4 + sum4help;
			}
		}
		
		double param2d = sv[1] + 0.5*sum4; 
		
		double temp4 = RNDGAM(param1d,param2d);
		sigma2_new = 1.0/temp4;

		/***************************************************************************************/
		// update of a1
		
		double va1=0;
		double va1help=0;
		int p;
		int q;
		for(p=0;p<J[0];p++){ 
			for(q=nobs_cum[p];q<nobs_cum[p+1];q++){
				va1help = (exp(-2*(exp(beta1j_new[p]))*(exp(b1_new))*(timeall[q])))/(sigma2_new);
				va1 = va1 + va1help;
			}
		}
		
		double ma1=0;
		double ma1help=0;
		int a;
		int b;
		for(a=0;a<J[0];a++){ 
			for(b=nobs_cum[a];b<nobs_cum[a+1];b++){
				ma1help = ((-concall[b] + (1-a0_new-alpha0j_new[a]) - (alpha1j_new[a])*exp(-(exp(beta1j_new[a]))*(exp(b1_new))*(timeall[b])))*(exp(-(exp(beta1j_new[a]))*(exp(b1_new))*(timeall[b]))))/(sigma2_new);
				ma1 = ma1 + ma1help;
			}
		}
		
		a1_new = normal(ma1/va1,1/va1);	
				
		/***************************************************************************************/
		// update of a0
				
		double va0=ttotal[0]/sigma2_new;
				
		double ma0=0;
		double ma0help=0;
		int m;
		int r;
		for(m=0;m<J[0];m++){ 
			for(r=nobs_cum[m];r<nobs_cum[m+1];r++){
				ma0help = (1 - concall[r] - alpha0j_new[m] - (a1_new+alpha1j_new[m])*exp(-(exp(beta1j_new[m]))*(exp(b1_new))*(timeall[r])))/(sigma2_new);
				ma0 = ma0 + ma0help;
			}
		}
			
		a0_new = normal(ma0/va0,1/va0);	
			
		/***************************************************************************************/
		// update of alpha1j
		
		int f;
		for(f=0;f<J[0];f++){ 
			double valpha1=0;
			double valpha1help=0;
			double malpha1=0;
			double malpha1help=0;
			
			int g;
			for(g=nobs_cum[f];g<nobs_cum[f+1];g++){ 
				valpha1help = (exp(-2*(exp(beta1j_new[f]))*(exp(b1_new))*(timeall[g])))/(sigma2_new);
				valpha1 = valpha1 + valpha1help;				
			
				malpha1help = ((-concall[g] + (1-a0_new-alpha0j_new[f]) - (a1_new)*exp(-(exp(beta1j_new[f]))*(exp(b1_new))*(timeall[g])))*(exp(-(exp(beta1j_new[f]))*(exp(b1_new))*(timeall[g]))))/(sigma2_new);
				malpha1 = malpha1 + malpha1help;
			}
		alpha1j_new[f] = normal(malpha1/(valpha1+(1/(tau2_alpha1_new))),1/(valpha1+(1/(tau2_alpha1_new))));
		}
		
		/***************************************************************************************/
		// update of alpha0j
				
		int v;
		double valpha0 = 0; 
		
		for(v=0;v<J[0];v++){ 
			
			valpha0 = nobs[v]/sigma2_new;
			double malpha0=0;
			double malpha0help=0;
					
			int t;
			for(t=nobs_cum[v];t<nobs_cum[v+1];t++){ 		
				malpha0help = (1 - concall[t] - a0_new - (a1_new+alpha1j_new[v])*exp(-(exp(beta1j_new[v]))*(exp(b1_new))*(timeall[t])))/(sigma2_new);
				malpha0 = malpha0 + malpha0help;
			}
		alpha0j_new[v] = normal(malpha0/(valpha0+(1/(tau2_alpha0_new))),1/(valpha0+(1/(tau2_alpha0_new))));			
		}			
						
		/*************** MH-UPDATES ************************************************************/
		// update of b1
    
		double b1_vorschlag = normal(b1_new,sigma2b1[0]);
		
		double nump2=0;
		int h1;
		for(h1=0;h1<J[0];h1++){ 
			double nump1=0;
			double nump1help=0;
			int x1;
			for(x1=nobs_cum[h1];x1<nobs_cum[h1+1];x1++){ 
				nump1help = pow((concall[x1]-((1-a0_new-alpha0j_new[h1]) - (a1_new+alpha1j_new[h1])*exp(-(exp(beta1j_new[h1]))*(exp(b1_vorschlag))*(timeall[x1])))),2);
				nump1 = nump1 + nump1help;
			}
		nump2 = nump2 + nump1;
		}
		double num = -(1/(2*(sigma2_new)))*nump2;
		
		double denomp2=0;
		int h2;
		for(h2=0;h2<J[0];h2++){ 
			double denomp1=0;
			double denomp1help=0;
			int x2;
			for(x2=nobs_cum[h2];x2<nobs_cum[h2+1];x2++){ 
				denomp1help = pow((concall[x2]-((1-a0_new-alpha0j_new[h2]) - (a1_new+alpha1j_new[h2])*exp(-(exp(beta1j_new[h2]))*(exp(b1_new))*(timeall[x2])))),2);
				denomp1 = denomp1 + denomp1help;
			}
		denomp2 = denomp2 + denomp1;
		}
		double denom = -(1/(2*(sigma2_new)))*denomp2;
		
		double paccept_b1=0;
		paccept_b1 = min(1,exp(num-denom));
		
		double randunif1 = nulleins();
		if(randunif1 < paccept_b1){
			b1_new = b1_vorschlag;
			acc_b1[0]++;
		}
		else{
			b1_new = b1_new;
		}

  	/***************************************************************************************/
  	// update of beta1j

    int h3;
		for(h3=0;h3<J[0];h3++){ 
			
		double beta1_vorschlag = normal(beta1j_new[h3],sigma2beta1[0]);

    double numpart1 = -log((sqrt((tau2_beta1_new)*2*M_PI))*(exp(beta1_vorschlag))) -(pow((beta1_vorschlag),2))/(2*(tau2_beta1_new));
		
		double numpart2=0;
		double numpart2help=0;
		int x3;
		for(x3=nobs_cum[h3];x3<nobs_cum[h3+1];x3++){ 
			numpart2help = pow((concall[x3]-((1-a0_new-alpha0j_new[h3]) - (a1_new+alpha1j_new[h3])*exp(-(exp(beta1_vorschlag))*(exp(b1_new))*(timeall[x3])))),2);
			numpart2 = numpart2 + numpart2help;
		}
		double num2 = numpart1 - (1/(2*(sigma2_new)))*numpart2;
		
		double denompart1 = -log((sqrt((tau2_beta1_new)*2*M_PI))*(exp(beta1j_new[h3]))) -(pow((beta1j_new[h3]),2))/(2*(tau2_beta1_new));
		
		double denompart2=0;
		double denompart2help=0;
		int x4;
		for(x4=nobs_cum[h3];x4<nobs_cum[h3+1];x4++){ 
			denompart2help = pow((concall[x4]-((1-a0_new-alpha0j_new[h3]) - (a1_new+alpha1j_new[h3])*exp(-(exp(beta1j_new[h3]))*(exp(b1_new))*(timeall[x4])))),2);
			denompart2 = denompart2 + denompart2help;
		}
		double denom2 = denompart1 - (1/(2*(sigma2_new)))*denompart2;
		
		double paccept_beta1=0;
		paccept_beta1 = min(1,exp(num2-denom2)); 
    
		double randunif2 = nulleins(); 
		if(randunif2 < paccept_beta1){
			beta1j_new[h3] = beta1_vorschlag;
			acc_beta1[0]++;
		}
		else{
			beta1j_new[h3] = beta1j_new[h3];
		}
		}
		
    /***************************************************************************************/
		// tuning
		if(i==tune[0])
		{
		// b1 
		tu_b1=0;
		temp=tune4k1(sigma2b1[0],acc_b1[0],tune[0]);

		if (sigma2b1[0]!=temp)
		{
			sigma2b1[0]=temp;
			tu_b1=1;
		}
		
		// beta1
		tu_beta1=0;
		temp=tune4k1(sigma2beta1[0],acc_beta1[0],(tune[0]*J[0]));
		if (sigma2beta1[0]!=temp)
		{
			sigma2beta1[0]=temp;
			tu_beta1=1;
		}

		if (tu_b1!=0 || tu_beta1!=0) 
		{
		   i=1; 
		   sample=0;
		   sample_vec=0;
		   acc_b1[0]=0;
		   acc_beta1[0]=0;
		}
		}
		
    // save results
		if(i>burnin[0] && fmod(i,thin[0])==0){
			tau2_alpha0[sample] = tau2_alpha0_new;
			tau2_alpha1[sample] = tau2_alpha1_new;
			tau2_beta1[sample] = tau2_beta1_new;
			sigma2[sample] = sigma2_new;
			a1[sample] = a1_new;
			a0[sample] = a0_new;
			b1[sample] = b1_new;

			int shelp;
			for(shelp=0;shelp<J[0];shelp++){ 
				alpha1[sample_vec] = alpha1j_new[shelp];
				alpha0[sample_vec] = alpha0j_new[shelp];
				beta1[sample_vec] = beta1j_new[shelp];
				sample_vec++;
			}
			sample++;	
		}			
	}
PutRNGstate();
}