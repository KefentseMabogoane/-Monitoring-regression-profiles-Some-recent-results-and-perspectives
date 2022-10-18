proc iml;
**********************************************************************;
*Number of simulations;
sim = 20000;
rlvec = j(sim,11,.);

*Size of the Phase II sample;
n = 4;

*lambda is the smothing constant;
lambda = 0.2; 
l1 = 1-lambda;

*Process parameters: Case K i.e. distribution parameters are known;
 mu0 = 0;
 sigma0 = 1;

stdev = sigma0;
* df, v and wp to be used under the t, gamma and Weibull distribution-not included here;
df=5;
v=1;
wp=2;

*IC regression parameters;
a0=3;
a1=2;
*Fixed explanatory variable;
xi={2,4,6,8};
xbar=mean(xi);
xp=xi-xbar*j(n,1,1);
b0=3+2*xbar;
b1=2;

*Control limits constants;
Lb0=3.0156; Lb1=3.0109; Le=1.3729;
*Shift, beta0, in the slope;
shift=0;
do beta0=0.000 to 2 by 0.2 ;
shift=shift+1;

delta = 0.00;

do k=1 to sim;

zib0_1 = b0;zib1_1 = b1;zie_1 = log(Sigma0**2);

vecb0={};vecb1={};vece={};
count=0;signal=0;

do i=1 to 1000000000 until (signal=1); 
	        yij = j(n,1,.);
		    eij=j(n,1,.);
	        * Generating a Phase II sample;

            call randgen(eij,'NORMAL',0,1);
	        yij=eij+a0*j(n,1,1)+a1*xi;
            * Regression Analysis;
            ybar=mean(yij);
			*xbar=mean(xi);
			xbarv=xbar*j(n,1,1);
			devx=xi-xbarv;
			sxy=sum(devx#yij);
			sxx=sum(devx##2);
			a1j=sxy/sxx;
			a0j=ybar-a1j*xbar;

			var0=sigma0**2/n+xbar**2*sigma0**2/sxx;var1=sigma0**2/sxx;cov01=-xbar*sigma0**2/sxx;cov10=cov01;
			Sig0={var0 cov01,cov10 var1};
			/*b0j=a0j+a1j*xbar;b1j=a1j;
			*mub0=b0; Varb0=sigma0**2/n;mub1=b1;Varb1=sigma0**2/sxx;*/
			b0j=(a0j+delta*stdev)+(a1j+beta0*stdev)*xbar;b1j=a1j+(beta0*stdev);
			mub0=a0j+a1j*xbar; Varb0=sigma0**2/n;mub1=a1j;Varb1=sigma0**2/sxx;
			*Residuals;
			Vare=Sigma0**2/n;
			eij=yij-b0j*j(n,1,1)-b1j*devx;
			ebar=mean(eij);
			msej=sum((eij##2))/(n-2);
            *Control limits;
     		* Asymptotic;

		   a=(l1)**2;

           Vz = lambda/(2-lambda);
            Vzb0=Vz*Varb0;Vzb1=Vz*Varb1;Vze=(2/(n-2)+2/(n-2)**2+4/(3*(n-2)**3)-16/(15*(n-2)**5))*Vz;
		    lclb0=b0-Lb0*sqrt(Vzb0);uclb0=b0+Lb0*sqrt(Vzb0);
            lclb1=b1-Lb1*sqrt(Vzb1);uclb1=b1+Lb1*sqrt(Vzb1);
            ucle=Le*sqrt(Vze);
		   *Obtaining the charting statistics;
           *EWMA for the Intercept Plotting Statistic;
			zib0 = lambda*b0j + (1-lambda)*zib0_1; 
	        zib0_1 = zib0;

           *EWMA for the Slope Plotting Statistic;
			zib1 = lambda*b1j + (1-lambda)*zib1_1; 
            zib1_1 = zib1;

           *EWMA for the error Plotting Statistic;
			if msej=0 then zie = max(lambda*log(msej+0.00001) + (1-lambda)*zie_1,log(sigma0**2)); 
			else zie = max(lambda*log(msej) + (1-lambda)*zie_1,log(sigma0**2)); 
			zie_1 = zie;
	        *Comparing the plotting statistics to the control limits;           
			vecb0=vecb0//zib0;vecb1=vecb1//zib1;vece=vece//zie;

            if ((zib0>=uclb0)|(zib0<=lclb0)|(zib1>=uclb1)|(zib1<=lclb1)|(zie>=ucle)) then signal=1;
			else signal=0;

            count=count+1;
            rlvec[k,shift]=count;
			end;
   end;
end;
results=rlvec;
	name1={"0.0" "0.1" "0.2" "0.3" "0.4" "0.5" "0.6" "0.7" "0.8" "0.9" "1.0"};	
    varNames="ARL"||"SDRL"||"MRL";
    results11=mean(rlvec);
	results12=std(rlvec);
	results13=median(rlvec);
	results21=results11`;
	results22=results12`;
	results23=results13`;
	results2=results21||results22||results23;

    EARL=mean(results2[2:11,1]);
	ESDRL=mean(results2[2:11,2]);
	EMRL=mean(results2[2:11,3]);

print results2[colname=varNames rowname=name1][format=10.2];
print EARL[format=10.2] ESDRL[format=10.2] EMRL[format=10.2];
print 'EWMA3 SLOPE ASYMPTOTIC All Shifts';
print sim n lambda Lb0 Lb1 Le;
*create runlength_prec_normal from results[colname={d000 d025 d050 d075 d100 d125 d150 d175 d200 d225 d250 d275 d300}];
*append from results;
quit;

