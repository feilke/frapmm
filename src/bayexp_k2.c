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
double tune4k2(double what,int acc,int n)
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
void k2_mhgibbs(double* timeall, double* concall, int* nobs,
		int* ttotal, int* J, double* sv, int* size, int* burnin, int* thin, int* tune, int* nobs_cum, double* sigma2b1, double* sigma2beta1,
		double* sigma2b2, double* sigma2beta2, double* a0start, double* a1start, double* b1start, double* a2start, double* b2start,
		double* tau2_alpha0, double* tau2_alpha1, double* tau2_beta1, double* tau2_alpha2, double* tau2_beta2, double* sigma2,
		double* a0, double* a1, double* b1, double* a2, double* b2, double* alpha0, double* alpha1, double* beta1, double* alpha2, double* beta2, 
		int* acc_b1, int* acc_beta1, int* acc_b2, int* acc_beta2)
{
	
GetRNGstate();

// MH and Gibbs (hybrid algorithm)
	int i;
	double temp;
	int tu_b1;
	int tu_beta1;
	int tu_b2;
	int tu_beta2;

	double tau2_alpha0_new;
	double tau2_alpha1_new;
	double tau2_beta1_new;
	double tau2_alpha2_new;
	double tau2_beta2_new;
	double sigma2_new;
	
	double a0_new=a0start[0];
	double a1_new=a1start[0];
	double b1_new=b1start[0];
	double a2_new=a2start[0];
	double b2_new=b2start[0];
	
	double alpha0j_new[J[0]];
	int i2;
	for(i2=0;i2<J[0];i2++){
		alpha0j_new[i2] = 0.0;
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
	double alpha2j_new[J[0]];
	int i5;
	for(i5=0;i5<J[0];i5++){
		alpha2j_new[i5] = 0.0;
	}
	double beta2j_new[J[0]];
	int i6;
	for(i6=0;i6<J[0];i6++){
		beta2j_new[i6] = 0.0;
	}
		
	int sample=0;
	int sample_vec=0;
	
	for(i=1; i<size[0]+1; i++) 
	{
		/*************** GIBBS UPDATES *********************************************************/
		// update of tau2_alpha0
		
		double param1g = sv[2]+(J[0]/2.0); 
				
		double sum1=0;
		int j;
		for(j=0;j<J[0];j++){
			sum1 += pow((alpha0j_new[j]),2);
		}
						
		double param2g = 0.5*sum1+sv[3]; 
				
		double temp1 = RNDGAM(param1g,param2g);	
		tau2_alpha0_new = 1.0/temp1;

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
		// update of tau2_alpha2
				
		double param1f = sv[8]+(J[0]/2.0);
						
		double sum5=0;
		int m;
		for(m=0;m<J[0];m++){
			sum5 += pow((alpha2j_new[m]),2);
		}
								
		double param2f = 0.5*sum5+sv[9]; 
				
		double temp5 = RNDGAM(param1f,param2f);	
		tau2_alpha2_new = 1.0/temp5;
			
		/***************************************************************************************/
		// update of tau2_beta2
			
		double param1e = sv[10]+(J[0]/2.0); 
						
		double sum6=0;
		int r;
		for(r=0;r<J[0];r++){
			sum6 += pow((beta2j_new[r]),2);
		}

		double param2e = 0.5*sum6+sv[11]; 
						
		double temp6 = RNDGAM(param1e,param2e);	
		tau2_beta2_new = 1.0/temp6;
	
		/***************************************************************************************/
		// update of sigma2
		
		double param1d = sv[0] + (ttotal[0]/2.0); 
		
		double sum4=0;
		double sum4help=0;
		int y;
		int z;
		for(y=0;y<J[0];y++){ 
			for(z=nobs_cum[y];z<nobs_cum[y+1];z++){ 
				sum4help = pow((concall[z] - ((1-a0_new-alpha0j_new[y])-(a1_new+alpha1j_new[y])*exp(-(exp(beta1j_new[y]))*(exp(b1_new))*(timeall[z]))-(a2_new+alpha2j_new[y])*exp(-(exp(beta2j_new[y]))*(exp(b2_new))*(timeall[z])))),2);
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
				ma1help = ((-concall[b] + (1-a0_new-alpha0j_new[a]) - (alpha1j_new[a])*exp(-(exp(beta1j_new[a]))*(exp(b1_new))*(timeall[b])) -  (a2_new+alpha2j_new[a])*exp(-(exp(beta2j_new[a]))*(exp(b2_new))*(timeall[b])))*(exp(-(exp(beta1j_new[a]))*(exp(b1_new))*(timeall[b]))))/(sigma2_new);
				ma1 = ma1 + ma1help;
			}
		}
		
		a1_new = normal(ma1/va1,1/va1);	
    
		/***************************************************************************************/
		// update of a2
		
		double va2=0;
		double va2help=0;
		int w;
		int x;
		for(w=0;w<J[0];w++){ 
			for(x=nobs_cum[w];x<nobs_cum[w+1];x++){
				va2help = (exp(-2*(exp(beta2j_new[w]))*(exp(b2_new))*(timeall[x])))/(sigma2_new);
				va2 = va2 + va2help;
			}
		}
		
		double ma2=0;
		double ma2help=0;
		int u;
		int v;
		for(u=0;u<J[0];u++){ 
			for(v=nobs_cum[u];v<nobs_cum[u+1];v++){
				ma2help = ((-concall[v] + (1-a0_new-alpha0j_new[u]) - (alpha2j_new[u])*exp(-(exp(beta2j_new[u]))*(exp(b2_new))*(timeall[v])) -  (a1_new+alpha1j_new[u])*exp(-(exp(beta1j_new[u]))*(exp(b1_new))*(timeall[v])))*(exp(-(exp(beta2j_new[u]))*(exp(b2_new))*(timeall[v]))))/(sigma2_new);
				ma2 = ma2 + ma2help;
			}
		}
		
		a2_new = normal(ma2/va2,1/va2);	

		/***************************************************************************************/
		// update of a0
					
		double va0=ttotal[0]/sigma2_new;
					
		double ma0=0;
		double ma0help=0;
		int o;
		int o1;
		for(o=0;o<J[0];o++){ 
		for(o1=nobs_cum[o];o1<nobs_cum[o+1];o1++){
				ma0help = (1 - concall[o1] - alpha0j_new[o] - (a1_new+alpha1j_new[o])*exp(-(exp(beta1j_new[o]))*(exp(b1_new))*(timeall[o1])) - (a2_new+alpha2j_new[o])*exp(-(exp(beta2j_new[o]))*(exp(b2_new))*(timeall[o1])))/(sigma2_new);
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
			
				malpha1help = ((-concall[g] + (1-a0_new-alpha0j_new[f]) - (a1_new)*exp(-(exp(beta1j_new[f]))*(exp(b1_new))*(timeall[g])) - (a2_new+alpha2j_new[f])*exp(-(exp(beta2j_new[f]))*(exp(b2_new))*(timeall[g])))*(exp(-(exp(beta1j_new[f]))*(exp(b1_new))*(timeall[g]))))/(sigma2_new);
				malpha1 = malpha1 + malpha1help;
			}
		alpha1j_new[f] = normal(malpha1/(valpha1+(1/(tau2_alpha1_new))),1/(valpha1+(1/(tau2_alpha1_new))));
		}

    /***************************************************************************************/
		// update of alpha2j
		
		int s;
		for(s=0;s<J[0];s++){
			double valpha2=0;
			double valpha2help=0;
			double malpha2=0;
			double malpha2help=0;
			
			int t;
			for(t=nobs_cum[s];t<nobs_cum[s+1];t++){ 
				valpha2help = (exp(-2*(exp(beta2j_new[s]))*(exp(b2_new))*(timeall[t])))/(sigma2_new);
				valpha2 = valpha2 + valpha2help;				
			
				malpha2help = ((-concall[t] + (1-a0_new-alpha0j_new[s]) - (a2_new)*exp(-(exp(beta2j_new[s]))*(exp(b2_new))*(timeall[t])) - (a1_new+alpha1j_new[s])*exp(-(exp(beta1j_new[s]))*(exp(b1_new))*(timeall[t])))*(exp(-(exp(beta2j_new[s]))*(exp(b2_new))*(timeall[t]))))/(sigma2_new);
				malpha2 = malpha2 + malpha2help;
			}
		alpha2j_new[s] = normal(malpha2/(valpha2+(1/(tau2_alpha2_new))),1/(valpha2+(1/(tau2_alpha2_new))));
		}
		
		/***************************************************************************************/
		// update of alpha0j
			
		int c;
		double valpha0 = 0; 
		
		for(c=0;c<J[0];c++){ 
			
			valpha0 = nobs[c]/sigma2_new;
			double malpha0=0;
			double malpha0help=0;
			int d;
			for(d=nobs_cum[c];d<nobs_cum[c+1];d++){
				malpha0help = (1 - concall[d] - a0_new - (a1_new+alpha1j_new[c])*exp(-(exp(beta1j_new[c]))*(exp(b1_new))*(timeall[d])) - (a2_new+alpha2j_new[c])*exp(-(exp(beta2j_new[c]))*(exp(b2_new))*(timeall[d])))/(sigma2_new);
				malpha0 = malpha0 + malpha0help;
			}
		alpha0j_new[c] = normal(malpha0/(valpha0+(1/(tau2_alpha0_new))),1/(valpha0+(1/(tau2_alpha0_new))));	
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
				nump1help = pow((concall[x1]-((1-a0_new-alpha0j_new[h1]) - (a1_new+alpha1j_new[h1])*exp(-(exp(beta1j_new[h1]))*(exp(b1_vorschlag))*(timeall[x1])) - (a2_new+alpha2j_new[h1])*exp(-(exp(beta2j_new[h1]))*(exp(b2_new))*(timeall[x1])) )),2);
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
				denomp1help = pow((concall[x2]-((1-a0_new-alpha0j_new[h2]) - (a1_new+alpha1j_new[h2])*exp(-(exp(beta1j_new[h2]))*(exp(b1_new))*(timeall[x2])) - (a2_new+alpha2j_new[h2])*exp(-(exp(beta2j_new[h2]))*(exp(b2_new))*(timeall[x2])) )),2);
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
			numpart2help = pow((concall[x3]-((1-a0_new-alpha0j_new[h3]) - (a1_new+alpha1j_new[h3])*exp(-(exp(beta1_vorschlag))*(exp(b1_new))*(timeall[x3])) - (a2_new+alpha2j_new[h3])*exp(-(exp(beta2j_new[h3]))*(exp(b2_new))*(timeall[x3])) )),2);
			numpart2 = numpart2 + numpart2help;
		}
		double num2 = numpart1 - (1/(2*(sigma2_new)))*numpart2;
		
		double denompart1 = -log((sqrt((tau2_beta1_new)*2*M_PI))*(exp(beta1j_new[h3]))) -(pow((beta1j_new[h3]),2))/(2*(tau2_beta1_new));
		
		double denompart2=0;
		double denompart2help=0;
		int x4;
		for(x4=nobs_cum[h3];x4<nobs_cum[h3+1];x4++){ 
			denompart2help = pow((concall[x4]-((1-a0_new-alpha0j_new[h3]) - (a1_new+alpha1j_new[h3])*exp(-(exp(beta1j_new[h3]))*(exp(b1_new))*(timeall[x4])) - (a2_new+alpha2j_new[h3])*exp(-(exp(beta2j_new[h3]))*(exp(b2_new))*(timeall[x4])) )),2);
			denompart2 = denompart2 + denompart2help;
		}
		double denom2 = denompart1 - (1/(2*(sigma2_new)))*denompart2;
		
		double paccept_beta1=0;
		paccept_beta1 = min(1,exp(num2-denom2)); // Akzeptanzwkeit für ein beta1j berechnen
				
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
		// update of b2
    
		double b2_vorschlag = normal(b2_new,sigma2b2[0]);
		
		double nump4=0;
		int h4;
		for(h4=0;h4<J[0];h4++){
			double nump3=0;
			double nump3help=0;
			int x5;
			for(x5=nobs_cum[h4];x5<nobs_cum[h4+1];x5++){ 
				nump3help = pow((concall[x5]-((1-a0_new-alpha0j_new[h4]) - (a1_new+alpha1j_new[h4])*exp(-(exp(beta1j_new[h4]))*(exp(b1_new))*(timeall[x5])) - (a2_new+alpha2j_new[h4])*exp(-(exp(beta2j_new[h4]))*(exp(b2_vorschlag))*(timeall[x5])) )),2);
				nump3 = nump3 + nump3help;
			}
		nump4 = nump4 + nump3;
		}
		double nume = -(1/(2*(sigma2_new)))*nump4;
		
		double denomp4=0;
		int h5;
		for(h5=0;h5<J[0];h5++){ 
			double denomp3=0;
			double denomp3help=0;
			int x6;
			for(x6=nobs_cum[h5];x6<nobs_cum[h5+1];x6++){ 
				denomp3help = pow((concall[x6]-((1-a0_new-alpha0j_new[h5]) - (a1_new+alpha1j_new[h5])*exp(-(exp(beta1j_new[h5]))*(exp(b1_new))*(timeall[x6])) - (a2_new+alpha2j_new[h5])*exp(-(exp(beta2j_new[h5]))*(exp(b2_new))*(timeall[x6])) )),2);
				denomp3 = denomp3 + denomp3help;
			}
		denomp4 = denomp4 + denomp3;
		}
		double denome = -(1/(2*(sigma2_new)))*denomp4;
		
		double paccept_b2=0;
		paccept_b2 = min(1,exp(nume-denome));
			
		double randunif3 = nulleins();
		if(randunif3 < paccept_b2){
			b2_new = b2_vorschlag;
			acc_b2[0]++;
		}
		else{
			b2_new = b2_new;
		}

		
		/***************************************************************************************/
		// update of beta2j

    int h6;
		for(h6=0;h6<J[0];h6++){ 
					
		double beta2_vorschlag = normal(beta2j_new[h6],sigma2beta2[0]);
				
		double numpart3 = -log((sqrt((tau2_beta2_new)*2*M_PI))*(exp(beta2j_new[h6]))) -(pow((beta2_vorschlag),2))/(2*(tau2_beta2_new));
			
		double numpart4=0;
		double numpart4help=0;
		int x7;
		for(x7=nobs_cum[h6];x7<nobs_cum[h6+1];x7++){ 
			numpart4help = pow((concall[x7]-((1-a0_new-alpha0j_new[h6]) - (a1_new+alpha1j_new[h6])*exp(-(exp(beta1j_new[h6]))*(exp(b1_new))*(timeall[x7])) - (a2_new+alpha2j_new[h6])*exp(-(exp(beta2_vorschlag))*(exp(b2_new))*(timeall[x7])) )),2);
			numpart4 = numpart4 + numpart4help;
		}
		double num3 = numpart3 - (1/(2*(sigma2_new)))*numpart4;
		
		double denompart3 = -log((sqrt((tau2_beta2_new)*2*M_PI))*(exp(beta2j_new[h6]))) -(pow((beta2j_new[h3]),2))/(2*(tau2_beta2_new));
			
		double denompart4=0;
		double denompart4help=0;
		int x8;
		for(x8=nobs_cum[h6];x8<nobs_cum[h6+1];x8++){ 
			denompart4help = pow((concall[x8]-((1-a0_new-alpha0j_new[h6]) - (a1_new+alpha1j_new[h6])*exp(-(exp(beta1j_new[h6]))*(exp(b1_new))*(timeall[x8])) - (a2_new+alpha2j_new[h6])*exp(-(exp(beta2j_new[h6]))*(exp(b2_new))*(timeall[x8])) )),2);
			denompart4 = denompart4 + denompart4help;
		}
		double denom3 = denompart3 - (1/(2*(sigma2_new)))*denompart4;
		
		double paccept_beta2=0;
		paccept_beta2 = min(1,exp(num3-denom3)); 
    
		double randunif4 = nulleins(); 
		if(randunif4 < paccept_beta2){
			beta2j_new[h6] = beta2_vorschlag;
			acc_beta2[0]++;
		}
		else{
			beta2j_new[h6] = beta2j_new[h6];
		}
		}
		
		/***************************************************************************************/
		// tuning
		if(i==tune[0])
		{
		// b1 
		tu_b1=0;
		temp=tune4k2(sigma2b1[0],acc_b1[0],tune[0]);
		if (sigma2b1[0]!=temp)
		{
			sigma2b1[0]=temp;
			tu_b1=1;
		}
		
		// beta1
		tu_beta1=0;
		temp=tune4k2(sigma2beta1[0],acc_beta1[0],(tune[0]*J[0]));
		if (sigma2beta1[0]!=temp)
		{
			sigma2beta1[0]=temp;
			tu_beta1=1;
		}
		
		// b2 
		tu_b2=0;
		temp=tune4k2(sigma2b2[0],acc_b2[0],tune[0]);
		if (sigma2b2[0]!=temp)
		{
			sigma2b2[0]=temp;
			tu_b2=1;
		}
				
		// beta2
		tu_beta2=0;
		temp=tune4k2(sigma2beta2[0],acc_beta2[0],(tune[0]*J[0]));
		if (sigma2beta2[0]!=temp)
		{
			sigma2beta2[0]=temp;
			tu_beta2=1;
		}

				   
		if (tu_b1!=0 || tu_beta1!=0 || tu_b2!=0 || tu_beta2!=0) 
		{
		   i=1; 
		   sample=0;
		   sample_vec=0;
		   acc_b1[0]=0;
		   acc_beta1[0]=0;
		   acc_b2[0]=0;
		   acc_beta2[0]=0;
		}
		}
		
		
		// save results
		if(i>burnin[0] && fmod(i,thin[0])==0){
			tau2_alpha0[sample] = tau2_alpha0_new;
			tau2_alpha1[sample] = tau2_alpha1_new;
			tau2_beta1[sample] = tau2_beta1_new;
			tau2_alpha2[sample] = tau2_alpha2_new;
			tau2_beta2[sample] = tau2_beta2_new;
			sigma2[sample] = sigma2_new;
			a0[sample] = a0_new;
			a1[sample] = a1_new;
			b1[sample] = b1_new;
			a2[sample] = a2_new;
			b2[sample] = b2_new;
			int shelp;
			for(shelp=0;shelp<J[0];shelp++){ 
				alpha0[sample_vec] = alpha0j_new[shelp];
				alpha1[sample_vec] = alpha1j_new[shelp];
				beta1[sample_vec] = beta1j_new[shelp];
				alpha2[sample_vec] = alpha2j_new[shelp];
				beta2[sample_vec] = beta2j_new[shelp];
				sample_vec++;
			}
			sample++;	
		}	
	}
PutRNGstate();
}