*'ModHWMA Fixed Explanatory variable Asymptotic case-(lambda=0.05 and k=-lambda/4)';
proc iml;
******RUN_LENGTH CHARACTERISTICS OF THE EWMA3 SCHEME USING MONTE CARLO SIMULATIONS********;
*Number of simulations;
sim = 20000;
call randseed(123);
*Run-length vector;
rlvec = j(sim,11,.);

*Size of the Phase II sample;
n = 4; *(must be bigger than 2 or p in multivariate);

*lambda is the smothing constant;
lambda = 0.05; 
l1 = 1-lambda;
k=-lambda/4;

*Process parameters: Case K i.e. distribution parameters are known;
 mu0 = 0;
 sigma0 = 1;

stdev = sigma0;

*IC regression parameters from Equation (1);
a0=3;
a1=2;
*Fixed explanatory variable;
xi={2,4,6,8};
xbar=mean(xi);
xp=xi-xbar*j(n,1,1);
b0=a0+a1*xbar;*b0 and b1 according to Equation (5);
b1=a1;

* Control limits constants;
Lb0=28.2223; Lb1=28.2223; Le=20.2633;



*Shift, delta, in the intercept with no shift in the slope and error variance;

shift=0;
do gam=1 to 2 by 0.2;
shift=shift+1;
delta=0.0 ;/*Shift in the intercept*/;


beta0 = 0.00; /*Shift in the slope*/;

*gam=1.0; /*Shift in the error variance*/;
Sige_shift=gam*stdev;/*Shift in the error variance in expressed in standard deviation*/;

do t=1 to sim;/*k represents the row number- simulation No = 20000*/;

MHjb0_1 = b0;MHjb1_1 = b1;MHje_1 = log(Sigma0**2);*Starting values for Z0(b0), Z0(b1) and Z0(e) as defined below Equations (7b and 12)/ Note Z0(e) = ln(sigma^2)=0;

b0j_1=b0; b1j_1=b1; msej_1=log(Sigma0**2);

vecb0={};vecb1={};vece={};
count=0;signal=0; *Define an indicator such that signal=0 for the IC and Signal=1 for the OOC situation/ we also need a count of the number of rational subgroups;

do i=1 to 1000000000 until (signal=1); * i reprents the time- Plotting Stat are computed for i=1,2,3, .../ stopp if signal = 1;
	        yij = j(n,1,.);
		    eij=j(n,1,.);
	        * Generating a Phase II sample;

            call randgen(eij,'NORMAL',0,Sige_shift);
	        yij=eij+a0*j(n,1,1)+a1*xi;*From Equation (1) where yi represents the profile i.e. dependent variable;

            * Regression Analysis: see Equation (2);
            ybar=mean(yij);
			*xbar=mean(xi);
			xbarv=xbar*j(n,1,1);
			devx=xi-xbarv;
			sxy=sum(devx#yij);
			sxx=sum(devx##2);
			a1j=sxy/sxx;                   *Equation (2);
			a0j=ybar-a1j*xbar;             *Equation (2);

            *Below Equation (3);
			var0=sigma0**2/n+xbar**2*sigma0**2/sxx;
            var1=sigma0**2/sxx;
            cov01=-xbar*sigma0**2/sxx;cov10=cov01;

			Sig0={var0 cov01,cov10 var1}; *Covariance matrix in Equation (3);

			b0j=(a0j+delta*stdev)+(a1j+beta0*stdev)*xbar;
            b1j=a1j+(beta0*stdev);*Equation (5) with a shift in the parameters- Note the shifts are expressed in standard deviations;

            *Expected values and Cov-matrix of the regression parameters (see Equation 6);
			mub0=a0j+a1j*xbar;
            mub1=a1j; 
            Varb0=sigma0**2/n;
            Varb1=sigma0**2/sxx;

			*Residuals;
			Vare=Sigma0**2/n;
			* Computing the mean squared error;
			eij=yij-b0j*j(n,1,1)-b1j*devx;
			ebar=mean(eij);
			msej=sum((eij##2))/(n-2);
            
            
           vecb0=vecb0//b0j;vecb1=vecb1//b1j;vece=vece//msej; *Dr JC i moved this up because i kept getting an error where it was. It was an error of placement;
            
            *Control limits;
     		* Asymptotic;

		
           Vz = ((lambda+k)##2+k##2);*Asymptotic case;
            Vzb0=Vz*Varb0; *variance for the intercept;
            Vzb1=Vz*Varb1; *variance for the slope;

            Vze=(2/(n-2)+2/(n-2)**2+4/(3*(n-2)**3)-16/(15*(n-2)**5))*Vz;* see equation (13);

		    lclb0=b0-Lb0*sqrt(Vzb0);uclb0=b0+Lb0*sqrt(Vzb0); *control limits for the intercept;
            lclb1=b1-Lb1*sqrt(Vzb1);uclb1=b1+Lb1*sqrt(Vzb1);*control limits for the slope;
            ucle=Le*sqrt(Vze);*control limits for the error variance;

		   *Obtaining the charting statistics;
           *Modified HWMA for the Intercept Plotting Statistic;
			MHjb0 = lambda*b0j + k*(b0j-b0j_1) + (1-lambda)*MHjb0_1; 
			b0j_1=b0j;
	        MHjb0_1 = mean(vecb0);

           *Modified HWMA for the Slope Plotting Statistic;
			MHjb1 = lambda*b1j + k*(b1j-b1j_1) + (1-lambda)*MHjb1_1; 
			b1j_1=b1j;
	        MHjb1_1 = mean(vecb1);

           *Modified HWMA for the error Plotting Statistic: See Equation (12);
			if (msej=0 | msej_1=0) then MHje = max(lambda*log(msej+0.00001) + k*(log(msej+0.00001)-log(msej_1+0.00001)) + (1-lambda)*MHje_1,log(sigma0**2)); 
			else MHje = max(lambda*log(msej) + k*(log(msej)-log(msej_1)) + (1-lambda)*MHje_1,log(sigma0**2));
			msej_1=msej;
			MHje_1 = mean(vece);

	        *Comparing the plotting statistics to the control limits;           
			

            if ((MHjb0>=uclb0)|(MHjb0<=lclb0)|(MHjb1>=uclb1)|(MHjb1<=lclb1)|(MHje>=ucle)) then signal=1;
			else signal=0;

            count=count+1;*counting the number of rational; 
            rlvec[t,shift]=count;
			end;
  
  end;
end;
results=rlvec;
	name1={"0.00" "0.20" "0.40" "0.60" "0.80" "1.00" "1.20" "1.40" "1.60" "1.80" "2.00"};
	name2={"1.00" "1.20" "1.40" "1.60" "1.80" "2.00"};
    varNames="ARL"||"SDRL"||"MRL";
    results11=mean(rlvec);
	results12=std(rlvec);
	results13=median(rlvec);
	results21=results11`;
	results22=results12`;
	results23=results13`;
	results2=results21||results22||results23;


create runlengths from rlvec;
append from rlvec;
print 'ModHWMA Fixed Explanatory variable Time varying case-(lambda=0.05 and k=-lambda/4)';
print sim n lambda k Lb0 Lb1 Le beta0 gam;
print results2[colname=varNames rowname=name2][format=10.2];
quit;







*'ModHWMA Fixed Explanatory variable Asymptotic case-(lambda=0.05 and k=-lambda/2)';
proc iml;
******RUN_LENGTH CHARACTERISTICS OF THE EWMA3 SCHEME USING MONTE CARLO SIMULATIONS********;
*Number of simulations;
sim = 20000;
call randseed(123);
*Run-length vector;
rlvec = j(sim,11,.);

*Size of the Phase II sample;
n = 4; *(must be bigger than 2 or p in multivariate);

*lambda is the smothing constant;
lambda = 0.05; 
l1 = 1-lambda;
k=-lambda/2;

*Process parameters: Case K i.e. distribution parameters are known;
 mu0 = 0;
 sigma0 = 1;

stdev = sigma0;

*IC regression parameters from Equation (1);
a0=3;
a1=2;
*Fixed explanatory variable;
xi={2,4,6,8};
xbar=mean(xi);
xp=xi-xbar*j(n,1,1);
b0=a0+a1*xbar;*b0 and b1 according to Equation (5);
b1=a1;

* Control limits constants;
Lb0=43.4555; Lb1=43.4555; Le=22.453;



*Shift, delta, in the intercept with no shift in the slope and error variance;

shift=0;
do gam=1 to 2 by 0.2;
shift=shift+1;
delta=0.0 ;/*Shift in the intercept*/;


beta0 = 0.00; /*Shift in the slope*/;

*gam=1.0; /*Shift in the error variance*/;
Sige_shift=gam*stdev;/*Shift in the error variance in expressed in standard deviation*/;

do t=1 to sim;/*k represents the row number- simulation No = 20000*/;

MHjb0_1 = b0;MHjb1_1 = b1;MHje_1 = log(Sigma0**2);*Starting values for Z0(b0), Z0(b1) and Z0(e) as defined below Equations (7b and 12)/ Note Z0(e) = ln(sigma^2)=0;

b0j_1=b0; b1j_1=b1; msej_1=log(Sigma0**2);

vecb0={};vecb1={};vece={};
count=0;signal=0; *Define an indicator such that signal=0 for the IC and Signal=1 for the OOC situation/ we also need a count of the number of rational subgroups;

do i=1 to 1000000000 until (signal=1); * i reprents the time- Plotting Stat are computed for i=1,2,3, .../ stopp if signal = 1;
	        yij = j(n,1,.);
		    eij=j(n,1,.);
	        * Generating a Phase II sample;

            call randgen(eij,'NORMAL',0,Sige_shift);
	        yij=eij+a0*j(n,1,1)+a1*xi;*From Equation (1) where yi represents the profile i.e. dependent variable;

            * Regression Analysis: see Equation (2);
            ybar=mean(yij);
			*xbar=mean(xi);
			xbarv=xbar*j(n,1,1);
			devx=xi-xbarv;
			sxy=sum(devx#yij);
			sxx=sum(devx##2);
			a1j=sxy/sxx;                   *Equation (2);
			a0j=ybar-a1j*xbar;             *Equation (2);

            *Below Equation (3);
			var0=sigma0**2/n+xbar**2*sigma0**2/sxx;
            var1=sigma0**2/sxx;
            cov01=-xbar*sigma0**2/sxx;cov10=cov01;

			Sig0={var0 cov01,cov10 var1}; *Covariance matrix in Equation (3);

			b0j=(a0j+delta*stdev)+(a1j+beta0*stdev)*xbar;
            b1j=a1j+(beta0*stdev);*Equation (5) with a shift in the parameters- Note the shifts are expressed in standard deviations;

            *Expected values and Cov-matrix of the regression parameters (see Equation 6);
			mub0=a0j+a1j*xbar;
            mub1=a1j; 
            Varb0=sigma0**2/n;
            Varb1=sigma0**2/sxx;

			*Residuals;
			Vare=Sigma0**2/n;
			* Computing the mean squared error;
			eij=yij-b0j*j(n,1,1)-b1j*devx;
			ebar=mean(eij);
			msej=sum((eij##2))/(n-2);
            
            
           vecb0=vecb0//b0j;vecb1=vecb1//b1j;vece=vece//msej; *Dr JC i moved this up because i kept getting an error where it was. It was an error of placement;
            
            *Control limits;
     		* Asymptotic;

		
           Vz = ((lambda+k)##2+k##2);*Asymptotic case;
            Vzb0=Vz*Varb0; *variance for the intercept;
            Vzb1=Vz*Varb1; *variance for the slope;

            Vze=(2/(n-2)+2/(n-2)**2+4/(3*(n-2)**3)-16/(15*(n-2)**5))*Vz;* see equation (13);

		    lclb0=b0-Lb0*sqrt(Vzb0);uclb0=b0+Lb0*sqrt(Vzb0); *control limits for the intercept;
            lclb1=b1-Lb1*sqrt(Vzb1);uclb1=b1+Lb1*sqrt(Vzb1);*control limits for the slope;
            ucle=Le*sqrt(Vze);*control limits for the error variance;

		   *Obtaining the charting statistics;
           *Modified HWMA for the Intercept Plotting Statistic;
			MHjb0 = lambda*b0j + k*(b0j-b0j_1) + (1-lambda)*MHjb0_1; 
			b0j_1=b0j;
	        MHjb0_1 = mean(vecb0);

           *Modified HWMA for the Slope Plotting Statistic;
			MHjb1 = lambda*b1j + k*(b1j-b1j_1) + (1-lambda)*MHjb1_1; 
			b1j_1=b1j;
	        MHjb1_1 = mean(vecb1);

           *Modified HWMA for the error Plotting Statistic: See Equation (12);
			if (msej=0 | msej_1=0) then MHje = max(lambda*log(msej+0.00001) + k*(log(msej+0.00001)-log(msej_1+0.00001)) + (1-lambda)*MHje_1,log(sigma0**2)); 
			else MHje = max(lambda*log(msej) + k*(log(msej)-log(msej_1)) + (1-lambda)*MHje_1,log(sigma0**2));
			msej_1=msej;
			MHje_1 = mean(vece);

	        *Comparing the plotting statistics to the control limits;           
			

            if ((MHjb0>=uclb0)|(MHjb0<=lclb0)|(MHjb1>=uclb1)|(MHjb1<=lclb1)|(MHje>=ucle)) then signal=1;
			else signal=0;

            count=count+1;*counting the number of rational; 
            rlvec[t,shift]=count;
			end;
  
  end;
end;
results=rlvec;
	name1={"0.00" "0.20" "0.40" "0.60" "0.80" "1.00" "1.20" "1.40" "1.60" "1.80" "2.00"};
	name2={"1.00" "1.20" "1.40" "1.60" "1.80" "2.00"};
    varNames="ARL"||"SDRL"||"MRL";
    results11=mean(rlvec);
	results12=std(rlvec);
	results13=median(rlvec);
	results21=results11`;
	results22=results12`;
	results23=results13`;
	results2=results21||results22||results23;

create runlengths from rlvec;
append from rlvec;
print 'ModHWMA Fixed Explanatory variable Asymptotic case-(lambda=0.05 and k=-lambda/2)';
print sim n lambda k Lb0 Lb1 Le beta0 gam;
print results2[colname=varNames rowname=name2][format=10.2];
quit;















*'ModHWMA Fixed Explanatory variable Asymptotic case-(lambda=0.05 and k=0)';
proc iml;
******RUN_LENGTH CHARACTERISTICS OF THE EWMA3 SCHEME USING MONTE CARLO SIMULATIONS********;
*Number of simulations;
sim = 20000;
call randseed(123);
*Run-length vector;
rlvec = j(sim,11,.);

*Size of the Phase II sample;
n = 4; *(must be bigger than 2 or p in multivariate);

*lambda is the smothing constant;
lambda = 0.05; 
l1 = 1-lambda;
k=0;

*Process parameters: Case K i.e. distribution parameters are known;
 mu0 = 0;
 sigma0 = 1;

stdev = sigma0;

*IC regression parameters from Equation (1);
a0=3;
a1=2;
*Fixed explanatory variable;
xi={2,4,6,8};
xbar=mean(xi);
xp=xi-xbar*j(n,1,1);
b0=a0+a1*xbar;*b0 and b1 according to Equation (5);
b1=a1;

* Control limits constants;
Lb0=32.1853; Lb1=32.1853; Le=16.175;



*Shift, delta, in the intercept with no shift in the slope and error variance;

shift=0;
do gam=1 to 2 by 0.2;
shift=shift+1;
delta=0.0 ;/*Shift in the intercept*/;


beta0 = 0.00; /*Shift in the slope*/;

*gam=1.0; /*Shift in the error variance*/;
Sige_shift=gam*stdev;/*Shift in the error variance in expressed in standard deviation*/;

do t=1 to sim;/*k represents the row number- simulation No = 20000*/;

MHjb0_1 = b0;MHjb1_1 = b1;MHje_1 = log(Sigma0**2);*Starting values for Z0(b0), Z0(b1) and Z0(e) as defined below Equations (7b and 12)/ Note Z0(e) = ln(sigma^2)=0;

b0j_1=b0; b1j_1=b1; msej_1=log(Sigma0**2);

vecb0={};vecb1={};vece={};
count=0;signal=0; *Define an indicator such that signal=0 for the IC and Signal=1 for the OOC situation/ we also need a count of the number of rational subgroups;

do i=1 to 1000000000 until (signal=1); * i reprents the time- Plotting Stat are computed for i=1,2,3, .../ stopp if signal = 1;
	        yij = j(n,1,.);
		    eij=j(n,1,.);
	        * Generating a Phase II sample;

            call randgen(eij,'NORMAL',0,Sige_shift);
	        yij=eij+a0*j(n,1,1)+a1*xi;*From Equation (1) where yi represents the profile i.e. dependent variable;

            * Regression Analysis: see Equation (2);
            ybar=mean(yij);
			*xbar=mean(xi);
			xbarv=xbar*j(n,1,1);
			devx=xi-xbarv;
			sxy=sum(devx#yij);
			sxx=sum(devx##2);
			a1j=sxy/sxx;                   *Equation (2);
			a0j=ybar-a1j*xbar;             *Equation (2);

            *Below Equation (3);
			var0=sigma0**2/n+xbar**2*sigma0**2/sxx;
            var1=sigma0**2/sxx;
            cov01=-xbar*sigma0**2/sxx;cov10=cov01;

			Sig0={var0 cov01,cov10 var1}; *Covariance matrix in Equation (3);

			b0j=(a0j+delta*stdev)+(a1j+beta0*stdev)*xbar;
            b1j=a1j+(beta0*stdev);*Equation (5) with a shift in the parameters- Note the shifts are expressed in standard deviations;

            *Expected values and Cov-matrix of the regression parameters (see Equation 6);
			mub0=a0j+a1j*xbar;
            mub1=a1j; 
            Varb0=sigma0**2/n;
            Varb1=sigma0**2/sxx;

			*Residuals;
			Vare=Sigma0**2/n;
			* Computing the mean squared error;
			eij=yij-b0j*j(n,1,1)-b1j*devx;
			ebar=mean(eij);
			msej=sum((eij##2))/(n-2);
            
            
           vecb0=vecb0//b0j;vecb1=vecb1//b1j;vece=vece//msej; *Dr JC i moved this up because i kept getting an error where it was. It was an error of placement;
            
            *Control limits;
     		* Asymptotic;

		
           Vz = ((lambda+k)##2+k##2);*Asymptotic case;
            Vzb0=Vz*Varb0; *variance for the intercept;
            Vzb1=Vz*Varb1; *variance for the slope;

            Vze=(2/(n-2)+2/(n-2)**2+4/(3*(n-2)**3)-16/(15*(n-2)**5))*Vz;* see equation (13);

		    lclb0=b0-Lb0*sqrt(Vzb0);uclb0=b0+Lb0*sqrt(Vzb0); *control limits for the intercept;
            lclb1=b1-Lb1*sqrt(Vzb1);uclb1=b1+Lb1*sqrt(Vzb1);*control limits for the slope;
            ucle=Le*sqrt(Vze);*control limits for the error variance;

		   *Obtaining the charting statistics;
           *Modified HWMA for the Intercept Plotting Statistic;
			MHjb0 = lambda*b0j + k*(b0j-b0j_1) + (1-lambda)*MHjb0_1; 
			b0j_1=b0j;
	        MHjb0_1 = mean(vecb0);

           *Modified HWMA for the Slope Plotting Statistic;
			MHjb1 = lambda*b1j + k*(b1j-b1j_1) + (1-lambda)*MHjb1_1; 
			b1j_1=b1j;
	        MHjb1_1 = mean(vecb1);

           *Modified HWMA for the error Plotting Statistic: See Equation (12);
			if (msej=0 | msej_1=0) then MHje = max(lambda*log(msej+0.00001) + k*(log(msej+0.00001)-log(msej_1+0.00001)) + (1-lambda)*MHje_1,log(sigma0**2)); 
			else MHje = max(lambda*log(msej) + k*(log(msej)-log(msej_1)) + (1-lambda)*MHje_1,log(sigma0**2));
			msej_1=msej;
			MHje_1 = mean(vece);

	        *Comparing the plotting statistics to the control limits;           
			

            if ((MHjb0>=uclb0)|(MHjb0<=lclb0)|(MHjb1>=uclb1)|(MHjb1<=lclb1)|(MHje>=ucle)) then signal=1;
			else signal=0;

            count=count+1;*counting the number of rational; 
            rlvec[t,shift]=count;
			end;
  
  end;
end;
results=rlvec;
	name1={"0.00" "0.20" "0.40" "0.60" "0.80" "1.00" "1.20" "1.40" "1.60" "1.80" "2.00"};
	name2={"1.00" "1.20" "1.40" "1.60" "1.80" "2.00"};
    varNames="ARL"||"SDRL"||"MRL";
    results11=mean(rlvec);
	results12=std(rlvec);
	results13=median(rlvec);
	results21=results11`;
	results22=results12`;
	results23=results13`;
	results2=results21||results22||results23;

create runlengths from rlvec;
append from rlvec;
print 'ModHWMA Fixed Explanatory variable Asymptotic case-(lambda=0.05 and k=0)';
print sim n lambda k Lb0 Lb1 Le beta0 gam;
print results2[colname=varNames rowname=name2][format=10.2];
quit;











*'ModHWMA Fixed Explanatory variable Asymptotic case-(lambda=0.05 and k=lambda/2)';
proc iml;
******RUN_LENGTH CHARACTERISTICS OF THE EWMA3 SCHEME USING MONTE CARLO SIMULATIONS********;
*Number of simulations;
sim = 20000;
call randseed(123);
*Run-length vector;
rlvec = j(sim,11,.);

*Size of the Phase II sample;
n = 4; *(must be bigger than 2 or p in multivariate);

*lambda is the smothing constant;
lambda = 0.05; 
l1 = 1-lambda;
k=lambda/2;

*Process parameters: Case K i.e. distribution parameters are known;
 mu0 = 0;
 sigma0 = 1;

stdev = sigma0;

*IC regression parameters from Equation (1);
a0=3;
a1=2;
*Fixed explanatory variable;
xi={2,4,6,8};
xbar=mean(xi);
xp=xi-xbar*j(n,1,1);
b0=a0+a1*xbar;*b0 and b1 according to Equation (5);
b1=a1;

* Control limits constants;
Lb0=22.278; Lb1=22.278; Le=11.273;



*Shift, delta, in the intercept with no shift in the slope and error variance;

shift=0;
do gam=1 to 2 by 0.2;
shift=shift+1;
delta=0.0 ;/*Shift in the intercept*/;


beta0 = 0.00; /*Shift in the slope*/;

*gam=1.0; /*Shift in the error variance*/;
Sige_shift=gam*stdev;/*Shift in the error variance in expressed in standard deviation*/;

do t=1 to sim;/*k represents the row number- simulation No = 20000*/;

MHjb0_1 = b0;MHjb1_1 = b1;MHje_1 = log(Sigma0**2);*Starting values for Z0(b0), Z0(b1) and Z0(e) as defined below Equations (7b and 12)/ Note Z0(e) = ln(sigma^2)=0;

b0j_1=b0; b1j_1=b1; msej_1=log(Sigma0**2);

vecb0={};vecb1={};vece={};
count=0;signal=0; *Define an indicator such that signal=0 for the IC and Signal=1 for the OOC situation/ we also need a count of the number of rational subgroups;

do i=1 to 1000000000 until (signal=1); * i reprents the time- Plotting Stat are computed for i=1,2,3, .../ stopp if signal = 1;
	        yij = j(n,1,.);
		    eij=j(n,1,.);
	        * Generating a Phase II sample;

            call randgen(eij,'NORMAL',0,Sige_shift);
	        yij=eij+a0*j(n,1,1)+a1*xi;*From Equation (1) where yi represents the profile i.e. dependent variable;

            * Regression Analysis: see Equation (2);
            ybar=mean(yij);
			*xbar=mean(xi);
			xbarv=xbar*j(n,1,1);
			devx=xi-xbarv;
			sxy=sum(devx#yij);
			sxx=sum(devx##2);
			a1j=sxy/sxx;                   *Equation (2);
			a0j=ybar-a1j*xbar;             *Equation (2);

            *Below Equation (3);
			var0=sigma0**2/n+xbar**2*sigma0**2/sxx;
            var1=sigma0**2/sxx;
            cov01=-xbar*sigma0**2/sxx;cov10=cov01;

			Sig0={var0 cov01,cov10 var1}; *Covariance matrix in Equation (3);

			b0j=(a0j+delta*stdev)+(a1j+beta0*stdev)*xbar;
            b1j=a1j+(beta0*stdev);*Equation (5) with a shift in the parameters- Note the shifts are expressed in standard deviations;

            *Expected values and Cov-matrix of the regression parameters (see Equation 6);
			mub0=a0j+a1j*xbar;
            mub1=a1j; 
            Varb0=sigma0**2/n;
            Varb1=sigma0**2/sxx;

			*Residuals;
			Vare=Sigma0**2/n;
			* Computing the mean squared error;
			eij=yij-b0j*j(n,1,1)-b1j*devx;
			ebar=mean(eij);
			msej=sum((eij##2))/(n-2);
            
            
           vecb0=vecb0//b0j;vecb1=vecb1//b1j;vece=vece//msej; *Dr JC i moved this up because i kept getting an error where it was. It was an error of placement;
            
            *Control limits;
     		* Asymptotic;

		
           Vz = ((lambda+k)##2+k##2);*Asymptotic case;
            Vzb0=Vz*Varb0; *variance for the intercept;
            Vzb1=Vz*Varb1; *variance for the slope;

            Vze=(2/(n-2)+2/(n-2)**2+4/(3*(n-2)**3)-16/(15*(n-2)**5))*Vz;* see equation (13);

		    lclb0=b0-Lb0*sqrt(Vzb0);uclb0=b0+Lb0*sqrt(Vzb0); *control limits for the intercept;
            lclb1=b1-Lb1*sqrt(Vzb1);uclb1=b1+Lb1*sqrt(Vzb1);*control limits for the slope;
            ucle=Le*sqrt(Vze);*control limits for the error variance;

		   *Obtaining the charting statistics;
           *Modified HWMA for the Intercept Plotting Statistic;
			MHjb0 = lambda*b0j + k*(b0j-b0j_1) + (1-lambda)*MHjb0_1; 
			b0j_1=b0j;
	        MHjb0_1 = mean(vecb0);

           *Modified HWMA for the Slope Plotting Statistic;
			MHjb1 = lambda*b1j + k*(b1j-b1j_1) + (1-lambda)*MHjb1_1; 
			b1j_1=b1j;
	        MHjb1_1 = mean(vecb1);

           *Modified HWMA for the error Plotting Statistic: See Equation (12);
			if (msej=0 | msej_1=0) then MHje = max(lambda*log(msej+0.00001) + k*(log(msej+0.00001)-log(msej_1+0.00001)) + (1-lambda)*MHje_1,log(sigma0**2)); 
			else MHje = max(lambda*log(msej) + k*(log(msej)-log(msej_1)) + (1-lambda)*MHje_1,log(sigma0**2));
			msej_1=msej;
			MHje_1 = mean(vece);

	        *Comparing the plotting statistics to the control limits;           
			

            if ((MHjb0>=uclb0)|(MHjb0<=lclb0)|(MHjb1>=uclb1)|(MHjb1<=lclb1)|(MHje>=ucle)) then signal=1;
			else signal=0;

            count=count+1;*counting the number of rational; 
            rlvec[t,shift]=count;
			end;
  
  end;
end;
results=rlvec;
	name1={"0.00" "0.20" "0.40" "0.60" "0.80" "1.00" "1.20" "1.40" "1.60" "1.80" "2.00"};
	name2={"1.00" "1.20" "1.40" "1.60" "1.80" "2.00"};
    varNames="ARL"||"SDRL"||"MRL";
    results11=mean(rlvec);
	results12=std(rlvec);
	results13=median(rlvec);
	results21=results11`;
	results22=results12`;
	results23=results13`;
	results2=results21||results22||results23;

create runlengths from rlvec;
append from rlvec;
print 'ModHWMA Fixed Explanatory variable Asymptotic case-(lambda=0.05 and k=lambda/2)';
print sim n lambda k Lb0 Lb1 Le beta0 gam;
print results2[colname=varNames rowname=name2][format=10.2];
quit;














*'ModHWMA Fixed Explanatory variable Asymptotic case-(lambda=0.05 and k=lambda/4)';
proc iml;
******RUN_LENGTH CHARACTERISTICS OF THE EWMA3 SCHEME USING MONTE CARLO SIMULATIONS********;
*Number of simulations;
sim = 20000;
call randseed(123);
*Run-length vector;
rlvec = j(sim,11,.);

*Size of the Phase II sample;
n = 4; *(must be bigger than 2 or p in multivariate);

*lambda is the smothing constant;
lambda = 0.05; 
l1 = 1-lambda;
k=lambda/4;

*Process parameters: Case K i.e. distribution parameters are known;
 mu0 = 0;
 sigma0 = 1;

stdev = sigma0;

*IC regression parameters from Equation (1);
a0=3;
a1=2;
*Fixed explanatory variable;
xi={2,4,6,8};
xbar=mean(xi);
xp=xi-xbar*j(n,1,1);
b0=a0+a1*xbar;*b0 and b1 according to Equation (5);
b1=a1;

* Control limits constants;
Lb0=23.798; Lb1=23.798; Le=13.245;



*Shift, delta, in the intercept with no shift in the slope and error variance;

shift=0;
do gam=1 to 2 by 0.2;
shift=shift+1;
delta=0.0 ;/*Shift in the intercept*/;


beta0 = 0.00; /*Shift in the slope*/;

*gam=1.0; /*Shift in the error variance*/;
Sige_shift=gam*stdev;/*Shift in the error variance in expressed in standard deviation*/;

do t=1 to sim;/*k represents the row number- simulation No = 20000*/;

MHjb0_1 = b0;MHjb1_1 = b1;MHje_1 = log(Sigma0**2);*Starting values for Z0(b0), Z0(b1) and Z0(e) as defined below Equations (7b and 12)/ Note Z0(e) = ln(sigma^2)=0;

b0j_1=b0; b1j_1=b1; msej_1=log(Sigma0**2);

vecb0={};vecb1={};vece={};
count=0;signal=0; *Define an indicator such that signal=0 for the IC and Signal=1 for the OOC situation/ we also need a count of the number of rational subgroups;

do i=1 to 1000000000 until (signal=1); * i reprents the time- Plotting Stat are computed for i=1,2,3, .../ stopp if signal = 1;
	        yij = j(n,1,.);
		    eij=j(n,1,.);
	        * Generating a Phase II sample;

            call randgen(eij,'NORMAL',0,Sige_shift);
	        yij=eij+a0*j(n,1,1)+a1*xi;*From Equation (1) where yi represents the profile i.e. dependent variable;

            * Regression Analysis: see Equation (2);
            ybar=mean(yij);
			*xbar=mean(xi);
			xbarv=xbar*j(n,1,1);
			devx=xi-xbarv;
			sxy=sum(devx#yij);
			sxx=sum(devx##2);
			a1j=sxy/sxx;                   *Equation (2);
			a0j=ybar-a1j*xbar;             *Equation (2);

            *Below Equation (3);
			var0=sigma0**2/n+xbar**2*sigma0**2/sxx;
            var1=sigma0**2/sxx;
            cov01=-xbar*sigma0**2/sxx;cov10=cov01;

			Sig0={var0 cov01,cov10 var1}; *Covariance matrix in Equation (3);

			b0j=(a0j+delta*stdev)+(a1j+beta0*stdev)*xbar;
            b1j=a1j+(beta0*stdev);*Equation (5) with a shift in the parameters- Note the shifts are expressed in standard deviations;

            *Expected values and Cov-matrix of the regression parameters (see Equation 6);
			mub0=a0j+a1j*xbar;
            mub1=a1j; 
            Varb0=sigma0**2/n;
            Varb1=sigma0**2/sxx;

			*Residuals;
			Vare=Sigma0**2/n;
			* Computing the mean squared error;
			eij=yij-b0j*j(n,1,1)-b1j*devx;
			ebar=mean(eij);
			msej=sum((eij##2))/(n-2);
            
            
           vecb0=vecb0//b0j;vecb1=vecb1//b1j;vece=vece//msej; *Dr JC i moved this up because i kept getting an error where it was. It was an error of placement;
            
            *Control limits;
     		* Asymptotic;

		
           Vz = ((lambda+k)##2+k##2);*Asymptotic case;
            Vzb0=Vz*Varb0; *variance for the intercept;
            Vzb1=Vz*Varb1; *variance for the slope;

            Vze=(2/(n-2)+2/(n-2)**2+4/(3*(n-2)**3)-16/(15*(n-2)**5))*Vz;* see equation (13);

		    lclb0=b0-Lb0*sqrt(Vzb0);uclb0=b0+Lb0*sqrt(Vzb0); *control limits for the intercept;
            lclb1=b1-Lb1*sqrt(Vzb1);uclb1=b1+Lb1*sqrt(Vzb1);*control limits for the slope;
            ucle=Le*sqrt(Vze);*control limits for the error variance;

		   *Obtaining the charting statistics;
           *Modified HWMA for the Intercept Plotting Statistic;
			MHjb0 = lambda*b0j + k*(b0j-b0j_1) + (1-lambda)*MHjb0_1; 
			b0j_1=b0j;
	        MHjb0_1 = mean(vecb0);

           *Modified HWMA for the Slope Plotting Statistic;
			MHjb1 = lambda*b1j + k*(b1j-b1j_1) + (1-lambda)*MHjb1_1; 
			b1j_1=b1j;
	        MHjb1_1 = mean(vecb1);

           *Modified HWMA for the error Plotting Statistic: See Equation (12);
			if (msej=0 | msej_1=0) then MHje = max(lambda*log(msej+0.00001) + k*(log(msej+0.00001)-log(msej_1+0.00001)) + (1-lambda)*MHje_1,log(sigma0**2)); 
			else MHje = max(lambda*log(msej) + k*(log(msej)-log(msej_1)) + (1-lambda)*MHje_1,log(sigma0**2));
			msej_1=msej;
			MHje_1 = mean(vece);

	        *Comparing the plotting statistics to the control limits;           
			

            if ((MHjb0>=uclb0)|(MHjb0<=lclb0)|(MHjb1>=uclb1)|(MHjb1<=lclb1)|(MHje>=ucle)) then signal=1;
			else signal=0;

            count=count+1;*counting the number of rational; 
            rlvec[t,shift]=count;
			end;
  
  end;
end;
results=rlvec;
	name1={"0.00" "0.20" "0.40" "0.60" "0.80" "1.00" "1.20" "1.40" "1.60" "1.80" "2.00"};
	name2={"1.00" "1.20" "1.40" "1.60" "1.80" "2.00"};
    varNames="ARL"||"SDRL"||"MRL";
    results11=mean(rlvec);
	results12=std(rlvec);
	results13=median(rlvec);
	results21=results11`;
	results22=results12`;
	results23=results13`;
	results2=results21||results22||results23;

create runlengths from rlvec;
append from rlvec;
print 'ModHWMA Fixed Explanatory variable Asymptotic case-(lambda=0.05 and k=lambda/4)';
print sim n lambda k Lb0 Lb1 Le beta0 gam;
print results2[colname=varNames rowname=name2][format=10.2];
quit;














*'ModHWMA Fixed Explanatory variable Asymptotic case-(lambda=0.1 and k=-lambda/4)';
proc iml;
******RUN_LENGTH CHARACTERISTICS OF THE EWMA3 SCHEME USING MONTE CARLO SIMULATIONS********;
*Number of simulations;
sim = 20000;
call randseed(123);
*Run-length vector;
rlvec = j(sim,11,.);

*Size of the Phase II sample;
n = 4; *(must be bigger than 2 or p in multivariate);

*lambda is the smothing constant;
lambda = 0.1; 
l1 = 1-lambda;
k=-lambda/4;

*Process parameters: Case K i.e. distribution parameters are known;
 mu0 = 0;
 sigma0 = 1;

stdev = sigma0;

*IC regression parameters from Equation (1);
a0=3;
a1=2;
*Fixed explanatory variable;
xi={2,4,6,8};
xbar=mean(xi);
xp=xi-xbar*j(n,1,1);
b0=a0+a1*xbar;*b0 and b1 according to Equation (5);
b1=a1;

* Control limits constants;
Lb0=21.2458; Lb1=21.2459; Le=10.2516;



*Shift, delta, in the intercept with no shift in the slope and error variance;

shift=0;
do gam=1 to 2 by 0.2;
shift=shift+1;
delta=0.0 ;/*Shift in the intercept*/;


beta0 = 0.00; /*Shift in the slope*/;

*gam=1.0; /*Shift in the error variance*/;
Sige_shift=gam*stdev;/*Shift in the error variance in expressed in standard deviation*/;

do t=1 to sim;/*k represents the row number- simulation No = 20000*/;

MHjb0_1 = b0;MHjb1_1 = b1;MHje_1 = log(Sigma0**2);*Starting values for Z0(b0), Z0(b1) and Z0(e) as defined below Equations (7b and 12)/ Note Z0(e) = ln(sigma^2)=0;

b0j_1=b0; b1j_1=b1; msej_1=log(Sigma0**2);

vecb0={};vecb1={};vece={};
count=0;signal=0; *Define an indicator such that signal=0 for the IC and Signal=1 for the OOC situation/ we also need a count of the number of rational subgroups;

do i=1 to 1000000000 until (signal=1); * i reprents the time- Plotting Stat are computed for i=1,2,3, .../ stopp if signal = 1;
	        yij = j(n,1,.);
		    eij=j(n,1,.);
	        * Generating a Phase II sample;

            call randgen(eij,'NORMAL',0,Sige_shift);
	        yij=eij+a0*j(n,1,1)+a1*xi;*From Equation (1) where yi represents the profile i.e. dependent variable;

            * Regression Analysis: see Equation (2);
            ybar=mean(yij);
			*xbar=mean(xi);
			xbarv=xbar*j(n,1,1);
			devx=xi-xbarv;
			sxy=sum(devx#yij);
			sxx=sum(devx##2);
			a1j=sxy/sxx;                   *Equation (2);
			a0j=ybar-a1j*xbar;             *Equation (2);

            *Below Equation (3);
			var0=sigma0**2/n+xbar**2*sigma0**2/sxx;
            var1=sigma0**2/sxx;
            cov01=-xbar*sigma0**2/sxx;cov10=cov01;

			Sig0={var0 cov01,cov10 var1}; *Covariance matrix in Equation (3);

			b0j=(a0j+delta*stdev)+(a1j+beta0*stdev)*xbar;
            b1j=a1j+(beta0*stdev);*Equation (5) with a shift in the parameters- Note the shifts are expressed in standard deviations;

            *Expected values and Cov-matrix of the regression parameters (see Equation 6);
			mub0=a0j+a1j*xbar;
            mub1=a1j; 
            Varb0=sigma0**2/n;
            Varb1=sigma0**2/sxx;

			*Residuals;
			Vare=Sigma0**2/n;
			* Computing the mean squared error;
			eij=yij-b0j*j(n,1,1)-b1j*devx;
			ebar=mean(eij);
			msej=sum((eij##2))/(n-2);
            
            
           vecb0=vecb0//b0j;vecb1=vecb1//b1j;vece=vece//msej; *Dr JC i moved this up because i kept getting an error where it was. It was an error of placement;
            
            *Control limits;
     		* Asymptotic;

		
           Vz = ((lambda+k)##2+k##2);*Asymptotic case;
            Vzb0=Vz*Varb0; *variance for the intercept;
            Vzb1=Vz*Varb1; *variance for the slope;

            Vze=(2/(n-2)+2/(n-2)**2+4/(3*(n-2)**3)-16/(15*(n-2)**5))*Vz;* see equation (13);

		    lclb0=b0-Lb0*sqrt(Vzb0);uclb0=b0+Lb0*sqrt(Vzb0); *control limits for the intercept;
            lclb1=b1-Lb1*sqrt(Vzb1);uclb1=b1+Lb1*sqrt(Vzb1);*control limits for the slope;
            ucle=Le*sqrt(Vze);*control limits for the error variance;

		   *Obtaining the charting statistics;
           *Modified HWMA for the Intercept Plotting Statistic;
			MHjb0 = lambda*b0j + k*(b0j-b0j_1) + (1-lambda)*MHjb0_1; 
			b0j_1=b0j;
	        MHjb0_1 = mean(vecb0);

           *Modified HWMA for the Slope Plotting Statistic;
			MHjb1 = lambda*b1j + k*(b1j-b1j_1) + (1-lambda)*MHjb1_1; 
			b1j_1=b1j;
	        MHjb1_1 = mean(vecb1);

           *Modified HWMA for the error Plotting Statistic: See Equation (12);
			if (msej=0 | msej_1=0) then MHje = max(lambda*log(msej+0.00001) + k*(log(msej+0.00001)-log(msej_1+0.00001)) + (1-lambda)*MHje_1,log(sigma0**2)); 
			else MHje = max(lambda*log(msej) + k*(log(msej)-log(msej_1)) + (1-lambda)*MHje_1,log(sigma0**2));
			msej_1=msej;
			MHje_1 = mean(vece);

	        *Comparing the plotting statistics to the control limits;           
			

            if ((MHjb0>=uclb0)|(MHjb0<=lclb0)|(MHjb1>=uclb1)|(MHjb1<=lclb1)|(MHje>=ucle)) then signal=1;
			else signal=0;

            count=count+1;*counting the number of rational; 
            rlvec[t,shift]=count;
			end;
  
  end;
end;
results=rlvec;
	name1={"0.00" "0.20" "0.40" "0.60" "0.80" "1.00" "1.20" "1.40" "1.60" "1.80" "2.00"};
	name2={"1.00" "1.20" "1.40" "1.60" "1.80" "2.00"};
    varNames="ARL"||"SDRL"||"MRL";
    results11=mean(rlvec);
	results12=std(rlvec);
	results13=median(rlvec);
	results21=results11`;
	results22=results12`;
	results23=results13`;
	results2=results21||results22||results23;

create runlengths from rlvec;
append from rlvec;
print 'ModHWMA Fixed Explanatory variable Asymptotic case-(lambda=0.1 and k=-lambda/4)';
print sim n lambda k Lb0 Lb1 Le beta0 gam;
print results2[colname=varNames rowname=name2][format=10.2];
quit;







*'ModHWMA Fixed Explanatory variable Asymptotic case-(lambda=0.1 and k=-lambda/2)';
proc iml;
******RUN_LENGTH CHARACTERISTICS OF THE EWMA3 SCHEME USING MONTE CARLO SIMULATIONS********;
*Number of simulations;
sim = 20000;
call randseed(123);
*Run-length vector;
rlvec = j(sim,11,.);

*Size of the Phase II sample;
n = 4; *(must be bigger than 2 or p in multivariate);

*lambda is the smothing constant;
lambda = 0.1; 
l1 = 1-lambda;
k=-lambda/2;

*Process parameters: Case K i.e. distribution parameters are known;
 mu0 = 0;
 sigma0 = 1;

stdev = sigma0;

*IC regression parameters from Equation (1);
a0=3;
a1=2;
*Fixed explanatory variable;
xi={2,4,6,8};
xbar=mean(xi);
xp=xi-xbar*j(n,1,1);
b0=a0+a1*xbar;*b0 and b1 according to Equation (5);
b1=a1;

* Control limits constants;
Lb0=22.7899; Lb1=22.7899; Le=11.408;



*Shift, delta, in the intercept with no shift in the slope and error variance;

shift=0;
do gam=1 to 2 by 0.2;
shift=shift+1;
delta=0.0 ;/*Shift in the intercept*/;


beta0 = 0.00; /*Shift in the slope*/;

*gam=1.0; /*Shift in the error variance*/;
Sige_shift=gam*stdev;/*Shift in the error variance in expressed in standard deviation*/;

do t=1 to sim;/*k represents the row number- simulation No = 20000*/;

MHjb0_1 = b0;MHjb1_1 = b1;MHje_1 = log(Sigma0**2);*Starting values for Z0(b0), Z0(b1) and Z0(e) as defined below Equations (7b and 12)/ Note Z0(e) = ln(sigma^2)=0;

b0j_1=b0; b1j_1=b1; msej_1=log(Sigma0**2);

vecb0={};vecb1={};vece={};
count=0;signal=0; *Define an indicator such that signal=0 for the IC and Signal=1 for the OOC situation/ we also need a count of the number of rational subgroups;

do i=1 to 1000000000 until (signal=1); * i reprents the time- Plotting Stat are computed for i=1,2,3, .../ stopp if signal = 1;
	        yij = j(n,1,.);
		    eij=j(n,1,.);
	        * Generating a Phase II sample;

            call randgen(eij,'NORMAL',0,Sige_shift);
	        yij=eij+a0*j(n,1,1)+a1*xi;*From Equation (1) where yi represents the profile i.e. dependent variable;

            * Regression Analysis: see Equation (2);
            ybar=mean(yij);
			*xbar=mean(xi);
			xbarv=xbar*j(n,1,1);
			devx=xi-xbarv;
			sxy=sum(devx#yij);
			sxx=sum(devx##2);
			a1j=sxy/sxx;                   *Equation (2);
			a0j=ybar-a1j*xbar;             *Equation (2);

            *Below Equation (3);
			var0=sigma0**2/n+xbar**2*sigma0**2/sxx;
            var1=sigma0**2/sxx;
            cov01=-xbar*sigma0**2/sxx;cov10=cov01;

			Sig0={var0 cov01,cov10 var1}; *Covariance matrix in Equation (3);

			b0j=(a0j+delta*stdev)+(a1j+beta0*stdev)*xbar;
            b1j=a1j+(beta0*stdev);*Equation (5) with a shift in the parameters- Note the shifts are expressed in standard deviations;

            *Expected values and Cov-matrix of the regression parameters (see Equation 6);
			mub0=a0j+a1j*xbar;
            mub1=a1j; 
            Varb0=sigma0**2/n;
            Varb1=sigma0**2/sxx;

			*Residuals;
			Vare=Sigma0**2/n;
			* Computing the mean squared error;
			eij=yij-b0j*j(n,1,1)-b1j*devx;
			ebar=mean(eij);
			msej=sum((eij##2))/(n-2);
            
            
           vecb0=vecb0//b0j;vecb1=vecb1//b1j;vece=vece//msej; *Dr JC i moved this up because i kept getting an error where it was. It was an error of placement;
            
            *Control limits;
     		* Asymptotic;

		
           Vz = ((lambda+k)##2+k##2);*Asymptotic case;
            Vzb0=Vz*Varb0; *variance for the intercept;
            Vzb1=Vz*Varb1; *variance for the slope;

            Vze=(2/(n-2)+2/(n-2)**2+4/(3*(n-2)**3)-16/(15*(n-2)**5))*Vz;* see equation (13);

		    lclb0=b0-Lb0*sqrt(Vzb0);uclb0=b0+Lb0*sqrt(Vzb0); *control limits for the intercept;
            lclb1=b1-Lb1*sqrt(Vzb1);uclb1=b1+Lb1*sqrt(Vzb1);*control limits for the slope;
            ucle=Le*sqrt(Vze);*control limits for the error variance;

		   *Obtaining the charting statistics;
           *Modified HWMA for the Intercept Plotting Statistic;
			MHjb0 = lambda*b0j + k*(b0j-b0j_1) + (1-lambda)*MHjb0_1; 
			b0j_1=b0j;
	        MHjb0_1 = mean(vecb0);

           *Modified HWMA for the Slope Plotting Statistic;
			MHjb1 = lambda*b1j + k*(b1j-b1j_1) + (1-lambda)*MHjb1_1; 
			b1j_1=b1j;
	        MHjb1_1 = mean(vecb1);

           *Modified HWMA for the error Plotting Statistic: See Equation (12);
			if (msej=0 | msej_1=0) then MHje = max(lambda*log(msej+0.00001) + k*(log(msej+0.00001)-log(msej_1+0.00001)) + (1-lambda)*MHje_1,log(sigma0**2)); 
			else MHje = max(lambda*log(msej) + k*(log(msej)-log(msej_1)) + (1-lambda)*MHje_1,log(sigma0**2));
			msej_1=msej;
			MHje_1 = mean(vece);

	        *Comparing the plotting statistics to the control limits;           
			

            if ((MHjb0>=uclb0)|(MHjb0<=lclb0)|(MHjb1>=uclb1)|(MHjb1<=lclb1)|(MHje>=ucle)) then signal=1;
			else signal=0;

            count=count+1;*counting the number of rational; 
            rlvec[t,shift]=count;
			end;
  
  end;
end;
results=rlvec;
	name1={"0.00" "0.20" "0.40" "0.60" "0.80" "1.00" "1.20" "1.40" "1.60" "1.80" "2.00"};
	name2={"1.00" "1.20" "1.40" "1.60" "1.80" "2.00"};
    varNames="ARL"||"SDRL"||"MRL";
    results11=mean(rlvec);
	results12=std(rlvec);
	results13=median(rlvec);
	results21=results11`;
	results22=results12`;
	results23=results13`;
	results2=results21||results22||results23;

create runlengths from rlvec;
append from rlvec;
print 'ModHWMA Fixed Explanatory variable Asymptotic case-(lambda=0.1 and k=-lambda/2)';
print sim n lambda k Lb0 Lb1 Le beta0 gam;
print results2[colname=varNames rowname=name2][format=10.2];
quit;















*'ModHWMA Fixed Explanatory variable Asymptotic case-(lambda=0.1 and k=0)';
proc iml;
******RUN_LENGTH CHARACTERISTICS OF THE EWMA3 SCHEME USING MONTE CARLO SIMULATIONS********;
*Number of simulations;
sim = 20000;
call randseed(123);
*Run-length vector;
rlvec = j(sim,11,.);

*Size of the Phase II sample;
n = 4; *(must be bigger than 2 or p in multivariate);

*lambda is the smothing constant;
lambda = 0.1; 
l1 = 1-lambda;
k=0;

*Process parameters: Case K i.e. distribution parameters are known;
 mu0 = 0;
 sigma0 = 1;

stdev = sigma0;

*IC regression parameters from Equation (1);
a0=3;
a1=2;
*Fixed explanatory variable;
xi={2,4,6,8};
xbar=mean(xi);
xp=xi-xbar*j(n,1,1);
b0=a0+a1*xbar;*b0 and b1 according to Equation (5);
b1=a1;

* Control limits constants;
Lb0=16.6329; Lb1=16.6329; Le=8.3644;



*Shift, delta, in the intercept with no shift in the slope and error variance;

shift=0;
do gam=1 to 2 by 0.2;
shift=shift+1;
delta=0.0 ;/*Shift in the intercept*/;


beta0 = 0.00; /*Shift in the slope*/;

*gam=1.0; /*Shift in the error variance*/;
Sige_shift=gam*stdev;/*Shift in the error variance in expressed in standard deviation*/;

do t=1 to sim;/*k represents the row number- simulation No = 20000*/;

MHjb0_1 = b0;MHjb1_1 = b1;MHje_1 = log(Sigma0**2);*Starting values for Z0(b0), Z0(b1) and Z0(e) as defined below Equations (7b and 12)/ Note Z0(e) = ln(sigma^2)=0;

b0j_1=b0; b1j_1=b1; msej_1=log(Sigma0**2);

vecb0={};vecb1={};vece={};
count=0;signal=0; *Define an indicator such that signal=0 for the IC and Signal=1 for the OOC situation/ we also need a count of the number of rational subgroups;

do i=1 to 1000000000 until (signal=1); * i reprents the time- Plotting Stat are computed for i=1,2,3, .../ stopp if signal = 1;
	        yij = j(n,1,.);
		    eij=j(n,1,.);
	        * Generating a Phase II sample;

            call randgen(eij,'NORMAL',0,Sige_shift);
	        yij=eij+a0*j(n,1,1)+a1*xi;*From Equation (1) where yi represents the profile i.e. dependent variable;

            * Regression Analysis: see Equation (2);
            ybar=mean(yij);
			*xbar=mean(xi);
			xbarv=xbar*j(n,1,1);
			devx=xi-xbarv;
			sxy=sum(devx#yij);
			sxx=sum(devx##2);
			a1j=sxy/sxx;                   *Equation (2);
			a0j=ybar-a1j*xbar;             *Equation (2);

            *Below Equation (3);
			var0=sigma0**2/n+xbar**2*sigma0**2/sxx;
            var1=sigma0**2/sxx;
            cov01=-xbar*sigma0**2/sxx;cov10=cov01;

			Sig0={var0 cov01,cov10 var1}; *Covariance matrix in Equation (3);

			b0j=(a0j+delta*stdev)+(a1j+beta0*stdev)*xbar;
            b1j=a1j+(beta0*stdev);*Equation (5) with a shift in the parameters- Note the shifts are expressed in standard deviations;

            *Expected values and Cov-matrix of the regression parameters (see Equation 6);
			mub0=a0j+a1j*xbar;
            mub1=a1j; 
            Varb0=sigma0**2/n;
            Varb1=sigma0**2/sxx;

			*Residuals;
			Vare=Sigma0**2/n;
			* Computing the mean squared error;
			eij=yij-b0j*j(n,1,1)-b1j*devx;
			ebar=mean(eij);
			msej=sum((eij##2))/(n-2);
            
            
           vecb0=vecb0//b0j;vecb1=vecb1//b1j;vece=vece//msej; *Dr JC i moved this up because i kept getting an error where it was. It was an error of placement;
            
            *Control limits;
     		* Asymptotic;

		
           Vz = ((lambda+k)##2+k##2);*Asymptotic case;
            Vzb0=Vz*Varb0; *variance for the intercept;
            Vzb1=Vz*Varb1; *variance for the slope;

            Vze=(2/(n-2)+2/(n-2)**2+4/(3*(n-2)**3)-16/(15*(n-2)**5))*Vz;* see equation (13);

		    lclb0=b0-Lb0*sqrt(Vzb0);uclb0=b0+Lb0*sqrt(Vzb0); *control limits for the intercept;
            lclb1=b1-Lb1*sqrt(Vzb1);uclb1=b1+Lb1*sqrt(Vzb1);*control limits for the slope;
            ucle=Le*sqrt(Vze);*control limits for the error variance;

		   *Obtaining the charting statistics;
           *Modified HWMA for the Intercept Plotting Statistic;
			MHjb0 = lambda*b0j + k*(b0j-b0j_1) + (1-lambda)*MHjb0_1; 
			b0j_1=b0j;
	        MHjb0_1 = mean(vecb0);

           *Modified HWMA for the Slope Plotting Statistic;
			MHjb1 = lambda*b1j + k*(b1j-b1j_1) + (1-lambda)*MHjb1_1; 
			b1j_1=b1j;
	        MHjb1_1 = mean(vecb1);

           *Modified HWMA for the error Plotting Statistic: See Equation (12);
			if (msej=0 | msej_1=0) then MHje = max(lambda*log(msej+0.00001) + k*(log(msej+0.00001)-log(msej_1+0.00001)) + (1-lambda)*MHje_1,log(sigma0**2)); 
			else MHje = max(lambda*log(msej) + k*(log(msej)-log(msej_1)) + (1-lambda)*MHje_1,log(sigma0**2));
			msej_1=msej;
			MHje_1 = mean(vece);

	        *Comparing the plotting statistics to the control limits;           
			

            if ((MHjb0>=uclb0)|(MHjb0<=lclb0)|(MHjb1>=uclb1)|(MHjb1<=lclb1)|(MHje>=ucle)) then signal=1;
			else signal=0;

            count=count+1;*counting the number of rational; 
            rlvec[t,shift]=count;
			end;
  
  end;
end;
results=rlvec;
	name1={"0.00" "0.20" "0.40" "0.60" "0.80" "1.00" "1.20" "1.40" "1.60" "1.80" "2.00"};
	name2={"1.00" "1.20" "1.40" "1.60" "1.80" "2.00"};
    varNames="ARL"||"SDRL"||"MRL";
    results11=mean(rlvec);
	results12=std(rlvec);
	results13=median(rlvec);
	results21=results11`;
	results22=results12`;
	results23=results13`;
	results2=results21||results22||results23;

create runlengths from rlvec;
append from rlvec;
print 'ModHWMA Fixed Explanatory variable Asymptotic case-(lambda=0.1 and k=0)';
print sim n lambda k Lb0 Lb1 Le beta0 gam;
print results2[colname=varNames rowname=name2][format=10.2];
quit;











*'ModHWMA Fixed Explanatory variable Asymptotic case-(lambda=0.1 and k=lambda/2)';
proc iml;
******RUN_LENGTH CHARACTERISTICS OF THE EWMA3 SCHEME USING MONTE CARLO SIMULATIONS********;
*Number of simulations;
sim = 20000;
call randseed(123);
*Run-length vector;
rlvec = j(sim,11,.);

*Size of the Phase II sample;
n = 4; *(must be bigger than 2 or p in multivariate);

*lambda is the smothing constant;
lambda = 0.1; 
l1 = 1-lambda;
k=lambda/2;

*Process parameters: Case K i.e. distribution parameters are known;
 mu0 = 0;
 sigma0 = 1;

stdev = sigma0;

*IC regression parameters from Equation (1);
a0=3;
a1=2;
*Fixed explanatory variable;
xi={2,4,6,8};
xbar=mean(xi);
xp=xi-xbar*j(n,1,1);
b0=a0+a1*xbar;*b0 and b1 according to Equation (5);
b1=a1;

* Control limits constants;
Lb0=12.55; Lb1=12.55; Le=6.264;



*Shift, delta, in the intercept with no shift in the slope and error variance;

shift=0;
do gam=1 to 2 by 0.2;
shift=shift+1;
delta=0.0 ;/*Shift in the intercept*/;


beta0 = 0.00; /*Shift in the slope*/;

*gam=1.0; /*Shift in the error variance*/;
Sige_shift=gam*stdev;/*Shift in the error variance in expressed in standard deviation*/;

do t=1 to sim;/*k represents the row number- simulation No = 20000*/;

MHjb0_1 = b0;MHjb1_1 = b1;MHje_1 = log(Sigma0**2);*Starting values for Z0(b0), Z0(b1) and Z0(e) as defined below Equations (7b and 12)/ Note Z0(e) = ln(sigma^2)=0;

b0j_1=b0; b1j_1=b1; msej_1=log(Sigma0**2);

vecb0={};vecb1={};vece={};
count=0;signal=0; *Define an indicator such that signal=0 for the IC and Signal=1 for the OOC situation/ we also need a count of the number of rational subgroups;

do i=1 to 1000000000 until (signal=1); * i reprents the time- Plotting Stat are computed for i=1,2,3, .../ stopp if signal = 1;
	        yij = j(n,1,.);
		    eij=j(n,1,.);
	        * Generating a Phase II sample;

            call randgen(eij,'NORMAL',0,Sige_shift);
	        yij=eij+a0*j(n,1,1)+a1*xi;*From Equation (1) where yi represents the profile i.e. dependent variable;

            * Regression Analysis: see Equation (2);
            ybar=mean(yij);
			*xbar=mean(xi);
			xbarv=xbar*j(n,1,1);
			devx=xi-xbarv;
			sxy=sum(devx#yij);
			sxx=sum(devx##2);
			a1j=sxy/sxx;                   *Equation (2);
			a0j=ybar-a1j*xbar;             *Equation (2);

            *Below Equation (3);
			var0=sigma0**2/n+xbar**2*sigma0**2/sxx;
            var1=sigma0**2/sxx;
            cov01=-xbar*sigma0**2/sxx;cov10=cov01;

			Sig0={var0 cov01,cov10 var1}; *Covariance matrix in Equation (3);

			b0j=(a0j+delta*stdev)+(a1j+beta0*stdev)*xbar;
            b1j=a1j+(beta0*stdev);*Equation (5) with a shift in the parameters- Note the shifts are expressed in standard deviations;

            *Expected values and Cov-matrix of the regression parameters (see Equation 6);
			mub0=a0j+a1j*xbar;
            mub1=a1j; 
            Varb0=sigma0**2/n;
            Varb1=sigma0**2/sxx;

			*Residuals;
			Vare=Sigma0**2/n;
			* Computing the mean squared error;
			eij=yij-b0j*j(n,1,1)-b1j*devx;
			ebar=mean(eij);
			msej=sum((eij##2))/(n-2);
            
            
           vecb0=vecb0//b0j;vecb1=vecb1//b1j;vece=vece//msej; *Dr JC i moved this up because i kept getting an error where it was. It was an error of placement;
            
            *Control limits;
     		* Asymptotic;

		
           Vz = ((lambda+k)##2+k##2);*Asymptotic case;
            Vzb0=Vz*Varb0; *variance for the intercept;
            Vzb1=Vz*Varb1; *variance for the slope;

            Vze=(2/(n-2)+2/(n-2)**2+4/(3*(n-2)**3)-16/(15*(n-2)**5))*Vz;* see equation (13);

		    lclb0=b0-Lb0*sqrt(Vzb0);uclb0=b0+Lb0*sqrt(Vzb0); *control limits for the intercept;
            lclb1=b1-Lb1*sqrt(Vzb1);uclb1=b1+Lb1*sqrt(Vzb1);*control limits for the slope;
            ucle=Le*sqrt(Vze);*control limits for the error variance;

		   *Obtaining the charting statistics;
           *Modified HWMA for the Intercept Plotting Statistic;
			MHjb0 = lambda*b0j + k*(b0j-b0j_1) + (1-lambda)*MHjb0_1; 
			b0j_1=b0j;
	        MHjb0_1 = mean(vecb0);

           *Modified HWMA for the Slope Plotting Statistic;
			MHjb1 = lambda*b1j + k*(b1j-b1j_1) + (1-lambda)*MHjb1_1; 
			b1j_1=b1j;
	        MHjb1_1 = mean(vecb1);

           *Modified HWMA for the error Plotting Statistic: See Equation (12);
			if (msej=0 | msej_1=0) then MHje = max(lambda*log(msej+0.00001) + k*(log(msej+0.00001)-log(msej_1+0.00001)) + (1-lambda)*MHje_1,log(sigma0**2)); 
			else MHje = max(lambda*log(msej) + k*(log(msej)-log(msej_1)) + (1-lambda)*MHje_1,log(sigma0**2));
			msej_1=msej;
			MHje_1 = mean(vece);

	        *Comparing the plotting statistics to the control limits;           
			

            if ((MHjb0>=uclb0)|(MHjb0<=lclb0)|(MHjb1>=uclb1)|(MHjb1<=lclb1)|(MHje>=ucle)) then signal=1;
			else signal=0;

            count=count+1;*counting the number of rational; 
            rlvec[t,shift]=count;
			end;
  
  end;
end;
results=rlvec;
	name1={"0.00" "0.20" "0.40" "0.60" "0.80" "1.00" "1.20" "1.40" "1.60" "1.80" "2.00"};
	name2={"1.00" "1.20" "1.40" "1.60" "1.80" "2.00"};
    varNames="ARL"||"SDRL"||"MRL";
    results11=mean(rlvec);
	results12=std(rlvec);
	results13=median(rlvec);
	results21=results11`;
	results22=results12`;
	results23=results13`;
	results2=results21||results22||results23;

create runlengths from rlvec;
append from rlvec;
print 'ModHWMA Fixed Explanatory variable Asymptotic case-(lambda=0.1 and k=lambda/2)';
print sim n lambda k Lb0 Lb1 Le beta0 gam;
print results2[colname=varNames rowname=name2][format=10.2];
quit;














*'ModHWMA Fixed Explanatory variable Asymptotic case-(lambda=0.1 and k=lambda/4)';
proc iml;
******RUN_LENGTH CHARACTERISTICS OF THE EWMA3 SCHEME USING MONTE CARLO SIMULATIONS********;
*Number of simulations;
sim = 20000;
call randseed(123);
*Run-length vector;
rlvec = j(sim,11,.);

*Size of the Phase II sample;
n = 4; *(must be bigger than 2 or p in multivariate);

*lambda is the smothing constant;
lambda = 0.1; 
l1 = 1-lambda;
k=lambda/4;

*Process parameters: Case K i.e. distribution parameters are known;
 mu0 = 0;
 sigma0 = 1;

stdev = sigma0;

*IC regression parameters from Equation (1);
a0=3;
a1=2;
*Fixed explanatory variable;
xi={2,4,6,8};
xbar=mean(xi);
xp=xi-xbar*j(n,1,1);
b0=a0+a1*xbar;*b0 and b1 according to Equation (5);
b1=a1;

* Control limits constants;
Lb0=14.0500; Lb1=14.0500; Le=7.0628;



*Shift, delta, in the intercept with no shift in the slope and error variance;

shift=0;
do gam=1 to 2 by 0.2;
shift=shift+1;
delta=0.0 ;/*Shift in the intercept*/;


beta0 = 0.00; /*Shift in the slope*/;

*gam=1.0; /*Shift in the error variance*/;
Sige_shift=gam*stdev;/*Shift in the error variance in expressed in standard deviation*/;

do t=1 to sim;/*k represents the row number- simulation No = 20000*/;

MHjb0_1 = b0;MHjb1_1 = b1;MHje_1 = log(Sigma0**2);*Starting values for Z0(b0), Z0(b1) and Z0(e) as defined below Equations (7b and 12)/ Note Z0(e) = ln(sigma^2)=0;

b0j_1=b0; b1j_1=b1; msej_1=log(Sigma0**2);

vecb0={};vecb1={};vece={};
count=0;signal=0; *Define an indicator such that signal=0 for the IC and Signal=1 for the OOC situation/ we also need a count of the number of rational subgroups;

do i=1 to 1000000000 until (signal=1); * i reprents the time- Plotting Stat are computed for i=1,2,3, .../ stopp if signal = 1;
	        yij = j(n,1,.);
		    eij=j(n,1,.);
	        * Generating a Phase II sample;

            call randgen(eij,'NORMAL',0,Sige_shift);
	        yij=eij+a0*j(n,1,1)+a1*xi;*From Equation (1) where yi represents the profile i.e. dependent variable;

            * Regression Analysis: see Equation (2);
            ybar=mean(yij);
			*xbar=mean(xi);
			xbarv=xbar*j(n,1,1);
			devx=xi-xbarv;
			sxy=sum(devx#yij);
			sxx=sum(devx##2);
			a1j=sxy/sxx;                   *Equation (2);
			a0j=ybar-a1j*xbar;             *Equation (2);

            *Below Equation (3);
			var0=sigma0**2/n+xbar**2*sigma0**2/sxx;
            var1=sigma0**2/sxx;
            cov01=-xbar*sigma0**2/sxx;cov10=cov01;

			Sig0={var0 cov01,cov10 var1}; *Covariance matrix in Equation (3);

			b0j=(a0j+delta*stdev)+(a1j+beta0*stdev)*xbar;
            b1j=a1j+(beta0*stdev);*Equation (5) with a shift in the parameters- Note the shifts are expressed in standard deviations;

            *Expected values and Cov-matrix of the regression parameters (see Equation 6);
			mub0=a0j+a1j*xbar;
            mub1=a1j; 
            Varb0=sigma0**2/n;
            Varb1=sigma0**2/sxx;

			*Residuals;
			Vare=Sigma0**2/n;
			* Computing the mean squared error;
			eij=yij-b0j*j(n,1,1)-b1j*devx;
			ebar=mean(eij);
			msej=sum((eij##2))/(n-2);
            
            
           vecb0=vecb0//b0j;vecb1=vecb1//b1j;vece=vece//msej; *Dr JC i moved this up because i kept getting an error where it was. It was an error of placement;
            
            *Control limits;
     		* Asymptotic;

		
           Vz = ((lambda+k)##2+k##2);*Asymptotic case;
            Vzb0=Vz*Varb0; *variance for the intercept;
            Vzb1=Vz*Varb1; *variance for the slope;

            Vze=(2/(n-2)+2/(n-2)**2+4/(3*(n-2)**3)-16/(15*(n-2)**5))*Vz;* see equation (13);

		    lclb0=b0-Lb0*sqrt(Vzb0);uclb0=b0+Lb0*sqrt(Vzb0); *control limits for the intercept;
            lclb1=b1-Lb1*sqrt(Vzb1);uclb1=b1+Lb1*sqrt(Vzb1);*control limits for the slope;
            ucle=Le*sqrt(Vze);*control limits for the error variance;

		   *Obtaining the charting statistics;
           *Modified HWMA for the Intercept Plotting Statistic;
			MHjb0 = lambda*b0j + k*(b0j-b0j_1) + (1-lambda)*MHjb0_1; 
			b0j_1=b0j;
	        MHjb0_1 = mean(vecb0);

           *Modified HWMA for the Slope Plotting Statistic;
			MHjb1 = lambda*b1j + k*(b1j-b1j_1) + (1-lambda)*MHjb1_1; 
			b1j_1=b1j;
	        MHjb1_1 = mean(vecb1);

           *Modified HWMA for the error Plotting Statistic: See Equation (12);
			if (msej=0 | msej_1=0) then MHje = max(lambda*log(msej+0.00001) + k*(log(msej+0.00001)-log(msej_1+0.00001)) + (1-lambda)*MHje_1,log(sigma0**2)); 
			else MHje = max(lambda*log(msej) + k*(log(msej)-log(msej_1)) + (1-lambda)*MHje_1,log(sigma0**2));
			msej_1=msej;
			MHje_1 = mean(vece);

	        *Comparing the plotting statistics to the control limits;           
			

            if ((MHjb0>=uclb0)|(MHjb0<=lclb0)|(MHjb1>=uclb1)|(MHjb1<=lclb1)|(MHje>=ucle)) then signal=1;
			else signal=0;

            count=count+1;*counting the number of rational; 
            rlvec[t,shift]=count;
			end;
  
  end;
end;
results=rlvec;
	name1={"0.00" "0.20" "0.40" "0.60" "0.80" "1.00" "1.20" "1.40" "1.60" "1.80" "2.00"};
	name2={"1.00" "1.20" "1.40" "1.60" "1.80" "2.00"};
    varNames="ARL"||"SDRL"||"MRL";
    results11=mean(rlvec);
	results12=std(rlvec);
	results13=median(rlvec);
	results21=results11`;
	results22=results12`;
	results23=results13`;
	results2=results21||results22||results23;

create runlengths from rlvec;
append from rlvec;
print 'ModHWMA Fixed Explanatory variable Asymptotic case-(lambda=0.1 and k=lambda/4)';
print sim n lambda k Lb0 Lb1 Le beta0 gam;
print results2[colname=varNames rowname=name2][format=10.2];
quit;











*'ModHWMA Fixed Explanatory variable Asymptotic case-(lambda=0.2 and k=-lambda/4)';
proc iml;
******RUN_LENGTH CHARACTERISTICS OF THE EWMA3 SCHEME USING MONTE CARLO SIMULATIONS********;
*Number of simulations;
sim = 20000;
call randseed(123);
*Run-length vector;
rlvec = j(sim,11,.);

*Size of the Phase II sample;
n = 4; *(must be bigger than 2 or p in multivariate);

*lambda is the smothing constant;
lambda = 0.2; 
l1 = 1-lambda;
k=-lambda/4;

*Process parameters: Case K i.e. distribution parameters are known;
 mu0 = 0;
 sigma0 = 1;

stdev = sigma0;

*IC regression parameters from Equation (1);
a0=3;
a1=2;
*Fixed explanatory variable;
xi={2,4,6,8};
xbar=mean(xi);
xp=xi-xbar*j(n,1,1);
b0=a0+a1*xbar;*b0 and b1 according to Equation (5);
b1=a1;

* Control limits constants;
Lb0=7.3569; Lb1=7.3569; Le=5.3615;



*Shift, delta, in the intercept with no shift in the slope and error variance;

shift=0;
do gam=1 to 2 by 0.2;
shift=shift+1;
delta=0.0 ;/*Shift in the intercept*/;


beta0 = 0.00; /*Shift in the slope*/;

*gam=1.0; /*Shift in the error variance*/;
Sige_shift=gam*stdev;/*Shift in the error variance in expressed in standard deviation*/;

do t=1 to sim;/*k represents the row number- simulation No = 20000*/;

MHjb0_1 = b0;MHjb1_1 = b1;MHje_1 = log(Sigma0**2);*Starting values for Z0(b0), Z0(b1) and Z0(e) as defined below Equations (7b and 12)/ Note Z0(e) = ln(sigma^2)=0;

b0j_1=b0; b1j_1=b1; msej_1=log(Sigma0**2);

vecb0={};vecb1={};vece={};
count=0;signal=0; *Define an indicator such that signal=0 for the IC and Signal=1 for the OOC situation/ we also need a count of the number of rational subgroups;

do i=1 to 1000000000 until (signal=1); * i reprents the time- Plotting Stat are computed for i=1,2,3, .../ stopp if signal = 1;
	        yij = j(n,1,.);
		    eij=j(n,1,.);
	        * Generating a Phase II sample;

            call randgen(eij,'NORMAL',0,Sige_shift);
	        yij=eij+a0*j(n,1,1)+a1*xi;*From Equation (1) where yi represents the profile i.e. dependent variable;

            * Regression Analysis: see Equation (2);
            ybar=mean(yij);
			*xbar=mean(xi);
			xbarv=xbar*j(n,1,1);
			devx=xi-xbarv;
			sxy=sum(devx#yij);
			sxx=sum(devx##2);
			a1j=sxy/sxx;                   *Equation (2);
			a0j=ybar-a1j*xbar;             *Equation (2);

            *Below Equation (3);
			var0=sigma0**2/n+xbar**2*sigma0**2/sxx;
            var1=sigma0**2/sxx;
            cov01=-xbar*sigma0**2/sxx;cov10=cov01;

			Sig0={var0 cov01,cov10 var1}; *Covariance matrix in Equation (3);

			b0j=(a0j+delta*stdev)+(a1j+beta0*stdev)*xbar;
            b1j=a1j+(beta0*stdev);*Equation (5) with a shift in the parameters- Note the shifts are expressed in standard deviations;

            *Expected values and Cov-matrix of the regression parameters (see Equation 6);
			mub0=a0j+a1j*xbar;
            mub1=a1j; 
            Varb0=sigma0**2/n;
            Varb1=sigma0**2/sxx;

			*Residuals;
			Vare=Sigma0**2/n;
			* Computing the mean squared error;
			eij=yij-b0j*j(n,1,1)-b1j*devx;
			ebar=mean(eij);
			msej=sum((eij##2))/(n-2);
            
            
           vecb0=vecb0//b0j;vecb1=vecb1//b1j;vece=vece//msej; *Dr JC i moved this up because i kept getting an error where it was. It was an error of placement;
            
            *Control limits;
     		* Asymptotic;

		
           Vz = ((lambda+k)##2+k##2);*Asymptotic case;
            Vzb0=Vz*Varb0; *variance for the intercept;
            Vzb1=Vz*Varb1; *variance for the slope;

            Vze=(2/(n-2)+2/(n-2)**2+4/(3*(n-2)**3)-16/(15*(n-2)**5))*Vz;* see equation (13);

		    lclb0=b0-Lb0*sqrt(Vzb0);uclb0=b0+Lb0*sqrt(Vzb0); *control limits for the intercept;
            lclb1=b1-Lb1*sqrt(Vzb1);uclb1=b1+Lb1*sqrt(Vzb1);*control limits for the slope;
            ucle=Le*sqrt(Vze);*control limits for the error variance;

		   *Obtaining the charting statistics;
           *Modified HWMA for the Intercept Plotting Statistic;
			MHjb0 = lambda*b0j + k*(b0j-b0j_1) + (1-lambda)*MHjb0_1; 
			b0j_1=b0j;
	        MHjb0_1 = mean(vecb0);

           *Modified HWMA for the Slope Plotting Statistic;
			MHjb1 = lambda*b1j + k*(b1j-b1j_1) + (1-lambda)*MHjb1_1; 
			b1j_1=b1j;
	        MHjb1_1 = mean(vecb1);

           *Modified HWMA for the error Plotting Statistic: See Equation (12);
			if (msej=0 | msej_1=0) then MHje = max(lambda*log(msej+0.00001) + k*(log(msej+0.00001)-log(msej_1+0.00001)) + (1-lambda)*MHje_1,log(sigma0**2)); 
			else MHje = max(lambda*log(msej) + k*(log(msej)-log(msej_1)) + (1-lambda)*MHje_1,log(sigma0**2));
			msej_1=msej;
			MHje_1 = mean(vece);

	        *Comparing the plotting statistics to the control limits;           
			

            if ((MHjb0>=uclb0)|(MHjb0<=lclb0)|(MHjb1>=uclb1)|(MHjb1<=lclb1)|(MHje>=ucle)) then signal=1;
			else signal=0;

            count=count+1;*counting the number of rational; 
            rlvec[t,shift]=count;
			end;
  
  end;
end;
results=rlvec;
	name1={"0.00" "0.20" "0.40" "0.60" "0.80" "1.00" "1.20" "1.40" "1.60" "1.80" "2.00"};
	name2={"1.00" "1.20" "1.40" "1.60" "1.80" "2.00"};
    varNames="ARL"||"SDRL"||"MRL";
    results11=mean(rlvec);
	results12=std(rlvec);
	results13=median(rlvec);
	results21=results11`;
	results22=results12`;
	results23=results13`;
	results2=results21||results22||results23;

create runlengths from rlvec;
append from rlvec;
print 'ModHWMA Fixed Explanatory variable Asymptotic case-(lambda=0.2 and k=-lambda/4)';
print sim n lambda k Lb0 Lb1 Le beta0 gam;
print results2[colname=varNames rowname=name2][format=10.2];
quit;







*'ModHWMA Fixed Explanatory variable Asymptotic case-(lambda=0.2 and k=-lambda/2)';
proc iml;
******RUN_LENGTH CHARACTERISTICS OF THE EWMA3 SCHEME USING MONTE CARLO SIMULATIONS********;
*Number of simulations;
sim = 20000;
call randseed(123);
*Run-length vector;
rlvec = j(sim,11,.);

*Size of the Phase II sample;
n = 4; *(must be bigger than 2 or p in multivariate);

*lambda is the smothing constant;
lambda = 0.2; 
l1 = 1-lambda;
k=-lambda/2;

*Process parameters: Case K i.e. distribution parameters are known;
 mu0 = 0;
 sigma0 = 1;

stdev = sigma0;

*IC regression parameters from Equation (1);
a0=3;
a1=2;
*Fixed explanatory variable;
xi={2,4,6,8};
xbar=mean(xi);
xp=xi-xbar*j(n,1,1);
b0=a0+a1*xbar;*b0 and b1 according to Equation (5);
b1=a1;

* Control limits constants;
Lb0=9.3568; Lb1=9.3568; Le=5.891;



*Shift, delta, in the intercept with no shift in the slope and error variance;

shift=0;
do gam=1 to 2 by 0.2;
shift=shift+1;
delta=0.0 ;/*Shift in the intercept*/;


beta0 = 0.00; /*Shift in the slope*/;

*gam=1.0; /*Shift in the error variance*/;
Sige_shift=gam*stdev;/*Shift in the error variance in expressed in standard deviation*/;

do t=1 to sim;/*k represents the row number- simulation No = 20000*/;

MHjb0_1 = b0;MHjb1_1 = b1;MHje_1 = log(Sigma0**2);*Starting values for Z0(b0), Z0(b1) and Z0(e) as defined below Equations (7b and 12)/ Note Z0(e) = ln(sigma^2)=0;

b0j_1=b0; b1j_1=b1; msej_1=log(Sigma0**2);

vecb0={};vecb1={};vece={};
count=0;signal=0; *Define an indicator such that signal=0 for the IC and Signal=1 for the OOC situation/ we also need a count of the number of rational subgroups;

do i=1 to 1000000000 until (signal=1); * i reprents the time- Plotting Stat are computed for i=1,2,3, .../ stopp if signal = 1;
	        yij = j(n,1,.);
		    eij=j(n,1,.);
	        * Generating a Phase II sample;

            call randgen(eij,'NORMAL',0,Sige_shift);
	        yij=eij+a0*j(n,1,1)+a1*xi;*From Equation (1) where yi represents the profile i.e. dependent variable;

            * Regression Analysis: see Equation (2);
            ybar=mean(yij);
			*xbar=mean(xi);
			xbarv=xbar*j(n,1,1);
			devx=xi-xbarv;
			sxy=sum(devx#yij);
			sxx=sum(devx##2);
			a1j=sxy/sxx;                   *Equation (2);
			a0j=ybar-a1j*xbar;             *Equation (2);

            *Below Equation (3);
			var0=sigma0**2/n+xbar**2*sigma0**2/sxx;
            var1=sigma0**2/sxx;
            cov01=-xbar*sigma0**2/sxx;cov10=cov01;

			Sig0={var0 cov01,cov10 var1}; *Covariance matrix in Equation (3);

			b0j=(a0j+delta*stdev)+(a1j+beta0*stdev)*xbar;
            b1j=a1j+(beta0*stdev);*Equation (5) with a shift in the parameters- Note the shifts are expressed in standard deviations;

            *Expected values and Cov-matrix of the regression parameters (see Equation 6);
			mub0=a0j+a1j*xbar;
            mub1=a1j; 
            Varb0=sigma0**2/n;
            Varb1=sigma0**2/sxx;

			*Residuals;
			Vare=Sigma0**2/n;
			* Computing the mean squared error;
			eij=yij-b0j*j(n,1,1)-b1j*devx;
			ebar=mean(eij);
			msej=sum((eij##2))/(n-2);
            
            
           vecb0=vecb0//b0j;vecb1=vecb1//b1j;vece=vece//msej; *Dr JC i moved this up because i kept getting an error where it was. It was an error of placement;
            
            *Control limits;
     		* Asymptotic;

		
           Vz = ((lambda+k)##2+k##2);*Asymptotic case;
            Vzb0=Vz*Varb0; *variance for the intercept;
            Vzb1=Vz*Varb1; *variance for the slope;

            Vze=(2/(n-2)+2/(n-2)**2+4/(3*(n-2)**3)-16/(15*(n-2)**5))*Vz;* see equation (13);

		    lclb0=b0-Lb0*sqrt(Vzb0);uclb0=b0+Lb0*sqrt(Vzb0); *control limits for the intercept;
            lclb1=b1-Lb1*sqrt(Vzb1);uclb1=b1+Lb1*sqrt(Vzb1);*control limits for the slope;
            ucle=Le*sqrt(Vze);*control limits for the error variance;

		   *Obtaining the charting statistics;
           *Modified HWMA for the Intercept Plotting Statistic;
			MHjb0 = lambda*b0j + k*(b0j-b0j_1) + (1-lambda)*MHjb0_1; 
			b0j_1=b0j;
	        MHjb0_1 = mean(vecb0);

           *Modified HWMA for the Slope Plotting Statistic;
			MHjb1 = lambda*b1j + k*(b1j-b1j_1) + (1-lambda)*MHjb1_1; 
			b1j_1=b1j;
	        MHjb1_1 = mean(vecb1);

           *Modified HWMA for the error Plotting Statistic: See Equation (12);
			if (msej=0 | msej_1=0) then MHje = max(lambda*log(msej+0.00001) + k*(log(msej+0.00001)-log(msej_1+0.00001)) + (1-lambda)*MHje_1,log(sigma0**2)); 
			else MHje = max(lambda*log(msej) + k*(log(msej)-log(msej_1)) + (1-lambda)*MHje_1,log(sigma0**2));
			msej_1=msej;
			MHje_1 = mean(vece);

	        *Comparing the plotting statistics to the control limits;           
			

            if ((MHjb0>=uclb0)|(MHjb0<=lclb0)|(MHjb1>=uclb1)|(MHjb1<=lclb1)|(MHje>=ucle)) then signal=1;
			else signal=0;

            count=count+1;*counting the number of rational; 
            rlvec[t,shift]=count;
			end;
  
  end;
end;
results=rlvec;
	name1={"0.00" "0.20" "0.40" "0.60" "0.80" "1.00" "1.20" "1.40" "1.60" "1.80" "2.00"};
	name2={"1.00" "1.20" "1.40" "1.60" "1.80" "2.00"};
    varNames="ARL"||"SDRL"||"MRL";
    results11=mean(rlvec);
	results12=std(rlvec);
	results13=median(rlvec);
	results21=results11`;
	results22=results12`;
	results23=results13`;
	results2=results21||results22||results23;

create runlengths from rlvec;
append from rlvec;
print 'ModHWMA Fixed Explanatory variable Asymptotic case-(lambda=0.2 and k=-lambda/2)';
print sim n lambda k Lb0 Lb1 Le beta0 gam;
print results2[colname=varNames rowname=name2][format=10.2];
quit;















*'ModHWMA Fixed Explanatory variable Asymptotic case-(lambda=0.2 and k=0)';
proc iml;
******RUN_LENGTH CHARACTERISTICS OF THE EWMA3 SCHEME USING MONTE CARLO SIMULATIONS********;
*Number of simulations;
sim = 20000;
call randseed(123);
*Run-length vector;
rlvec = j(sim,11,.);

*Size of the Phase II sample;
n = 4; *(must be bigger than 2 or p in multivariate);

*lambda is the smothing constant;
lambda = 0.2; 
l1 = 1-lambda;
k=0;

*Process parameters: Case K i.e. distribution parameters are known;
 mu0 = 0;
 sigma0 = 1;

stdev = sigma0;

*IC regression parameters from Equation (1);
a0=3;
a1=2;
*Fixed explanatory variable;
xi={2,4,6,8};
xbar=mean(xi);
xp=xi-xbar*j(n,1,1);
b0=a0+a1*xbar;*b0 and b1 according to Equation (5);
b1=a1;

* Control limits constants;
Lb0=8.8050; Lb1=8.8050; Le=4.4430;



*Shift, delta, in the intercept with no shift in the slope and error variance;

shift=0;
do gam=1 to 2 by 0.2;
shift=shift+1;
delta=0.0 ;/*Shift in the intercept*/;


beta0 = 0.00; /*Shift in the slope*/;

*gam=1.0; /*Shift in the error variance*/;
Sige_shift=gam*stdev;/*Shift in the error variance in expressed in standard deviation*/;

do t=1 to sim;/*k represents the row number- simulation No = 20000*/;

MHjb0_1 = b0;MHjb1_1 = b1;MHje_1 = log(Sigma0**2);*Starting values for Z0(b0), Z0(b1) and Z0(e) as defined below Equations (7b and 12)/ Note Z0(e) = ln(sigma^2)=0;

b0j_1=b0; b1j_1=b1; msej_1=log(Sigma0**2);

vecb0={};vecb1={};vece={};
count=0;signal=0; *Define an indicator such that signal=0 for the IC and Signal=1 for the OOC situation/ we also need a count of the number of rational subgroups;

do i=1 to 1000000000 until (signal=1); * i reprents the time- Plotting Stat are computed for i=1,2,3, .../ stopp if signal = 1;
	        yij = j(n,1,.);
		    eij=j(n,1,.);
	        * Generating a Phase II sample;

            call randgen(eij,'NORMAL',0,Sige_shift);
	        yij=eij+a0*j(n,1,1)+a1*xi;*From Equation (1) where yi represents the profile i.e. dependent variable;

            * Regression Analysis: see Equation (2);
            ybar=mean(yij);
			*xbar=mean(xi);
			xbarv=xbar*j(n,1,1);
			devx=xi-xbarv;
			sxy=sum(devx#yij);
			sxx=sum(devx##2);
			a1j=sxy/sxx;                   *Equation (2);
			a0j=ybar-a1j*xbar;             *Equation (2);

            *Below Equation (3);
			var0=sigma0**2/n+xbar**2*sigma0**2/sxx;
            var1=sigma0**2/sxx;
            cov01=-xbar*sigma0**2/sxx;cov10=cov01;

			Sig0={var0 cov01,cov10 var1}; *Covariance matrix in Equation (3);

			b0j=(a0j+delta*stdev)+(a1j+beta0*stdev)*xbar;
            b1j=a1j+(beta0*stdev);*Equation (5) with a shift in the parameters- Note the shifts are expressed in standard deviations;

            *Expected values and Cov-matrix of the regression parameters (see Equation 6);
			mub0=a0j+a1j*xbar;
            mub1=a1j; 
            Varb0=sigma0**2/n;
            Varb1=sigma0**2/sxx;

			*Residuals;
			Vare=Sigma0**2/n;
			* Computing the mean squared error;
			eij=yij-b0j*j(n,1,1)-b1j*devx;
			ebar=mean(eij);
			msej=sum((eij##2))/(n-2);
            
            
           vecb0=vecb0//b0j;vecb1=vecb1//b1j;vece=vece//msej; *Dr JC i moved this up because i kept getting an error where it was. It was an error of placement;
            
            *Control limits;
     		* Asymptotic;

		
           Vz = ((lambda+k)##2+k##2);*Asymptotic case;
            Vzb0=Vz*Varb0; *variance for the intercept;
            Vzb1=Vz*Varb1; *variance for the slope;

            Vze=(2/(n-2)+2/(n-2)**2+4/(3*(n-2)**3)-16/(15*(n-2)**5))*Vz;* see equation (13);

		    lclb0=b0-Lb0*sqrt(Vzb0);uclb0=b0+Lb0*sqrt(Vzb0); *control limits for the intercept;
            lclb1=b1-Lb1*sqrt(Vzb1);uclb1=b1+Lb1*sqrt(Vzb1);*control limits for the slope;
            ucle=Le*sqrt(Vze);*control limits for the error variance;

		   *Obtaining the charting statistics;
           *Modified HWMA for the Intercept Plotting Statistic;
			MHjb0 = lambda*b0j + k*(b0j-b0j_1) + (1-lambda)*MHjb0_1; 
			b0j_1=b0j;
	        MHjb0_1 = mean(vecb0);

           *Modified HWMA for the Slope Plotting Statistic;
			MHjb1 = lambda*b1j + k*(b1j-b1j_1) + (1-lambda)*MHjb1_1; 
			b1j_1=b1j;
	        MHjb1_1 = mean(vecb1);

           *Modified HWMA for the error Plotting Statistic: See Equation (12);
			if (msej=0 | msej_1=0) then MHje = max(lambda*log(msej+0.00001) + k*(log(msej+0.00001)-log(msej_1+0.00001)) + (1-lambda)*MHje_1,log(sigma0**2)); 
			else MHje = max(lambda*log(msej) + k*(log(msej)-log(msej_1)) + (1-lambda)*MHje_1,log(sigma0**2));
			msej_1=msej;
			MHje_1 = mean(vece);

	        *Comparing the plotting statistics to the control limits;           
			

            if ((MHjb0>=uclb0)|(MHjb0<=lclb0)|(MHjb1>=uclb1)|(MHjb1<=lclb1)|(MHje>=ucle)) then signal=1;
			else signal=0;

            count=count+1;*counting the number of rational; 
            rlvec[t,shift]=count;
			end;
  
  end;
end;
results=rlvec;
	name1={"0.00" "0.20" "0.40" "0.60" "0.80" "1.00" "1.20" "1.40" "1.60" "1.80" "2.00"};	
	name2={"1.00" "1.20" "1.40" "1.60" "1.80" "2.00"};
    varNames="ARL"||"SDRL"||"MRL";
    results11=mean(rlvec);
	results12=std(rlvec);
	results13=median(rlvec);
	results21=results11`;
	results22=results12`;
	results23=results13`;
	results2=results21||results22||results23;

create runlengths from rlvec;
append from rlvec;
print 'ModHWMA Fixed Explanatory variable Asymptotic case-(lambda=0.2 and k=0)';
print sim n lambda k Lb0 Lb1 Le beta0 gam;
print results2[colname=varNames rowname=name2][format=10.2];
quit;











*'ModHWMA Fixed Explanatory variable Asymptotic case-(lambda=0.2 and k=lambda/2)';
proc iml;
******RUN_LENGTH CHARACTERISTICS OF THE EWMA3 SCHEME USING MONTE CARLO SIMULATIONS********;
*Number of simulations;
sim = 20000;
call randseed(123);
*Run-length vector;
rlvec = j(sim,11,.);

*Size of the Phase II sample;
n = 4; *(must be bigger than 2 or p in multivariate);

*lambda is the smothing constant;
lambda = 0.2; 
l1 = 1-lambda;
k=lambda/2;

*Process parameters: Case K i.e. distribution parameters are known;
 mu0 = 0;
 sigma0 = 1;

stdev = sigma0;

*IC regression parameters from Equation (1);
a0=3;
a1=2;
*Fixed explanatory variable;
xi={2,4,6,8};
xbar=mean(xi);
xp=xi-xbar*j(n,1,1);
b0=a0+a1*xbar;*b0 and b1 according to Equation (5);
b1=a1;

* Control limits constants;
Lb0=6.78; Lb1=6.78; Le=3.7403;



*Shift, delta, in the intercept with no shift in the slope and error variance;

shift=0;
do gam=1 to 2 by 0.2;
shift=shift+1;
delta=0.0 ;/*Shift in the intercept*/;


beta0 = 0.00; /*Shift in the slope*/;

*gam=1.0; /*Shift in the error variance*/;
Sige_shift=gam*stdev;/*Shift in the error variance in expressed in standard deviation*/;

do t=1 to sim;/*k represents the row number- simulation No = 20000*/;

MHjb0_1 = b0;MHjb1_1 = b1;MHje_1 = log(Sigma0**2);*Starting values for Z0(b0), Z0(b1) and Z0(e) as defined below Equations (7b and 12)/ Note Z0(e) = ln(sigma^2)=0;

b0j_1=b0; b1j_1=b1; msej_1=log(Sigma0**2);

vecb0={};vecb1={};vece={};
count=0;signal=0; *Define an indicator such that signal=0 for the IC and Signal=1 for the OOC situation/ we also need a count of the number of rational subgroups;

do i=1 to 1000000000 until (signal=1); * i reprents the time- Plotting Stat are computed for i=1,2,3, .../ stopp if signal = 1;
	        yij = j(n,1,.);
		    eij=j(n,1,.);
	        * Generating a Phase II sample;

            call randgen(eij,'NORMAL',0,Sige_shift);
	        yij=eij+a0*j(n,1,1)+a1*xi;*From Equation (1) where yi represents the profile i.e. dependent variable;

            * Regression Analysis: see Equation (2);
            ybar=mean(yij);
			*xbar=mean(xi);
			xbarv=xbar*j(n,1,1);
			devx=xi-xbarv;
			sxy=sum(devx#yij);
			sxx=sum(devx##2);
			a1j=sxy/sxx;                   *Equation (2);
			a0j=ybar-a1j*xbar;             *Equation (2);

            *Below Equation (3);
			var0=sigma0**2/n+xbar**2*sigma0**2/sxx;
            var1=sigma0**2/sxx;
            cov01=-xbar*sigma0**2/sxx;cov10=cov01;

			Sig0={var0 cov01,cov10 var1}; *Covariance matrix in Equation (3);

			b0j=(a0j+delta*stdev)+(a1j+beta0*stdev)*xbar;
            b1j=a1j+(beta0*stdev);*Equation (5) with a shift in the parameters- Note the shifts are expressed in standard deviations;

            *Expected values and Cov-matrix of the regression parameters (see Equation 6);
			mub0=a0j+a1j*xbar;
            mub1=a1j; 
            Varb0=sigma0**2/n;
            Varb1=sigma0**2/sxx;

			*Residuals;
			Vare=Sigma0**2/n;
			* Computing the mean squared error;
			eij=yij-b0j*j(n,1,1)-b1j*devx;
			ebar=mean(eij);
			msej=sum((eij##2))/(n-2);
            
            
           vecb0=vecb0//b0j;vecb1=vecb1//b1j;vece=vece//msej; *Dr JC i moved this up because i kept getting an error where it was. It was an error of placement;
            
            *Control limits;
     		* Asymptotic;

		
           Vz = ((lambda+k)##2+k##2);*Asymptotic case;
            Vzb0=Vz*Varb0; *variance for the intercept;
            Vzb1=Vz*Varb1; *variance for the slope;

            Vze=(2/(n-2)+2/(n-2)**2+4/(3*(n-2)**3)-16/(15*(n-2)**5))*Vz;* see equation (13);

		    lclb0=b0-Lb0*sqrt(Vzb0);uclb0=b0+Lb0*sqrt(Vzb0); *control limits for the intercept;
            lclb1=b1-Lb1*sqrt(Vzb1);uclb1=b1+Lb1*sqrt(Vzb1);*control limits for the slope;
            ucle=Le*sqrt(Vze);*control limits for the error variance;

		   *Obtaining the charting statistics;
           *Modified HWMA for the Intercept Plotting Statistic;
			MHjb0 = lambda*b0j + k*(b0j-b0j_1) + (1-lambda)*MHjb0_1; 
			b0j_1=b0j;
	        MHjb0_1 = mean(vecb0);

           *Modified HWMA for the Slope Plotting Statistic;
			MHjb1 = lambda*b1j + k*(b1j-b1j_1) + (1-lambda)*MHjb1_1; 
			b1j_1=b1j;
	        MHjb1_1 = mean(vecb1);

           *Modified HWMA for the error Plotting Statistic: See Equation (12);
			if (msej=0 | msej_1=0) then MHje = max(lambda*log(msej+0.00001) + k*(log(msej+0.00001)-log(msej_1+0.00001)) + (1-lambda)*MHje_1,log(sigma0**2)); 
			else MHje = max(lambda*log(msej) + k*(log(msej)-log(msej_1)) + (1-lambda)*MHje_1,log(sigma0**2));
			msej_1=msej;
			MHje_1 = mean(vece);

	        *Comparing the plotting statistics to the control limits;           
			

            if ((MHjb0>=uclb0)|(MHjb0<=lclb0)|(MHjb1>=uclb1)|(MHjb1<=lclb1)|(MHje>=ucle)) then signal=1;
			else signal=0;

            count=count+1;*counting the number of rational; 
            rlvec[t,shift]=count;
			end;
  
  end;
end;
results=rlvec;
	name1={"0.00" "0.20" "0.40" "0.60" "0.80" "1.00" "1.20" "1.40" "1.60" "1.80" "2.00"};
	name2={"1.00" "1.20" "1.40" "1.60" "1.80" "2.00"};
    varNames="ARL"||"SDRL"||"MRL";
    results11=mean(rlvec);
	results12=std(rlvec);
	results13=median(rlvec);
	results21=results11`;
	results22=results12`;
	results23=results13`;
	results2=results21||results22||results23;

create runlengths from rlvec;
append from rlvec;
print 'ModHWMA Fixed Explanatory variable Asymptotic case-(lambda=0.2 and k=lambda/2)';
print sim n lambda k Lb0 Lb1 Le beta0 gam;
print results2[colname=varNames rowname=name2][format=10.2];
quit;














*'ModHWMA Fixed Explanatory variable Asymptotic case-(lambda=0.2 and k=lambda/4)';
proc iml;
******RUN_LENGTH CHARACTERISTICS OF THE EWMA3 SCHEME USING MONTE CARLO SIMULATIONS********;
*Number of simulations;
sim = 20000;
call randseed(123);
*Run-length vector;
rlvec = j(sim,11,.);

*Size of the Phase II sample;
n = 4; *(must be bigger than 2 or p in multivariate);

*lambda is the smothing constant;
lambda = 0.2; 
l1 = 1-lambda;
k=lambda/4;

*Process parameters: Case K i.e. distribution parameters are known;
 mu0 = 0;
 sigma0 = 1;

stdev = sigma0;

*IC regression parameters from Equation (1);
a0=3;
a1=2;
*Fixed explanatory variable;
xi={2,4,6,8};
xbar=mean(xi);
xp=xi-xbar*j(n,1,1);
b0=a0+a1*xbar;*b0 and b1 according to Equation (5);
b1=a1;

* Control limits constants;
Lb0=6.78; Lb1=6.78; Le=3.9767;



*Shift, delta, in the intercept with no shift in the slope and error variance;

shift=0;
do gam=1 to 2 by 0.2;
shift=shift+1;
delta=0.0 ;/*Shift in the intercept*/;


beta0 = 0.00; /*Shift in the slope*/;

*gam=1.0; /*Shift in the error variance*/;
Sige_shift=gam*stdev;/*Shift in the error variance in expressed in standard deviation*/;

do t=1 to sim;/*k represents the row number- simulation No = 20000*/;

MHjb0_1 = b0;MHjb1_1 = b1;MHje_1 = log(Sigma0**2);*Starting values for Z0(b0), Z0(b1) and Z0(e) as defined below Equations (7b and 12)/ Note Z0(e) = ln(sigma^2)=0;

b0j_1=b0; b1j_1=b1; msej_1=log(Sigma0**2);

vecb0={};vecb1={};vece={};
count=0;signal=0; *Define an indicator such that signal=0 for the IC and Signal=1 for the OOC situation/ we also need a count of the number of rational subgroups;

do i=1 to 1000000000 until (signal=1); * i reprents the time- Plotting Stat are computed for i=1,2,3, .../ stopp if signal = 1;
	        yij = j(n,1,.);
		    eij=j(n,1,.);
	        * Generating a Phase II sample;

            call randgen(eij,'NORMAL',0,Sige_shift);
	        yij=eij+a0*j(n,1,1)+a1*xi;*From Equation (1) where yi represents the profile i.e. dependent variable;

            * Regression Analysis: see Equation (2);
            ybar=mean(yij);
			*xbar=mean(xi);
			xbarv=xbar*j(n,1,1);
			devx=xi-xbarv;
			sxy=sum(devx#yij);
			sxx=sum(devx##2);
			a1j=sxy/sxx;                   *Equation (2);
			a0j=ybar-a1j*xbar;             *Equation (2);

            *Below Equation (3);
			var0=sigma0**2/n+xbar**2*sigma0**2/sxx;
            var1=sigma0**2/sxx;
            cov01=-xbar*sigma0**2/sxx;cov10=cov01;

			Sig0={var0 cov01,cov10 var1}; *Covariance matrix in Equation (3);

			b0j=(a0j+delta*stdev)+(a1j+beta0*stdev)*xbar;
            b1j=a1j+(beta0*stdev);*Equation (5) with a shift in the parameters- Note the shifts are expressed in standard deviations;

            *Expected values and Cov-matrix of the regression parameters (see Equation 6);
			mub0=a0j+a1j*xbar;
            mub1=a1j; 
            Varb0=sigma0**2/n;
            Varb1=sigma0**2/sxx;

			*Residuals;
			Vare=Sigma0**2/n;
			* Computing the mean squared error;
			eij=yij-b0j*j(n,1,1)-b1j*devx;
			ebar=mean(eij);
			msej=sum((eij##2))/(n-2);
            
            
           vecb0=vecb0//b0j;vecb1=vecb1//b1j;vece=vece//msej; *Dr JC i moved this up because i kept getting an error where it was. It was an error of placement;
            
            *Control limits;
     		* Asymptotic;

		
           Vz = ((lambda+k)##2+k##2);*Asymptotic case;
            Vzb0=Vz*Varb0; *variance for the intercept;
            Vzb1=Vz*Varb1; *variance for the slope;

            Vze=(2/(n-2)+2/(n-2)**2+4/(3*(n-2)**3)-16/(15*(n-2)**5))*Vz;* see equation (13);

		    lclb0=b0-Lb0*sqrt(Vzb0);uclb0=b0+Lb0*sqrt(Vzb0); *control limits for the intercept;
            lclb1=b1-Lb1*sqrt(Vzb1);uclb1=b1+Lb1*sqrt(Vzb1);*control limits for the slope;
            ucle=Le*sqrt(Vze);*control limits for the error variance;

		   *Obtaining the charting statistics;
           *Modified HWMA for the Intercept Plotting Statistic;
			MHjb0 = lambda*b0j + k*(b0j-b0j_1) + (1-lambda)*MHjb0_1; 
			b0j_1=b0j;
	        MHjb0_1 = mean(vecb0);

           *Modified HWMA for the Slope Plotting Statistic;
			MHjb1 = lambda*b1j + k*(b1j-b1j_1) + (1-lambda)*MHjb1_1; 
			b1j_1=b1j;
	        MHjb1_1 = mean(vecb1);

           *Modified HWMA for the error Plotting Statistic: See Equation (12);
			if (msej=0 | msej_1=0) then MHje = max(lambda*log(msej+0.00001) + k*(log(msej+0.00001)-log(msej_1+0.00001)) + (1-lambda)*MHje_1,log(sigma0**2)); 
			else MHje = max(lambda*log(msej) + k*(log(msej)-log(msej_1)) + (1-lambda)*MHje_1,log(sigma0**2));
			msej_1=msej;
			MHje_1 = mean(vece);

	        *Comparing the plotting statistics to the control limits;           
			

            if ((MHjb0>=uclb0)|(MHjb0<=lclb0)|(MHjb1>=uclb1)|(MHjb1<=lclb1)|(MHje>=ucle)) then signal=1;
			else signal=0;

            count=count+1;*counting the number of rational; 
            rlvec[t,shift]=count;
			end;
  
  end;
end;
results=rlvec;
	name1={"0.00" "0.20" "0.40" "0.60" "0.80" "1.00" "1.20" "1.40" "1.60" "1.80" "2.00"};
	name2={"1.00" "1.20" "1.40" "1.60" "1.80" "2.00"};
    varNames="ARL"||"SDRL"||"MRL";
    results11=mean(rlvec);
	results12=std(rlvec);
	results13=median(rlvec);
	results21=results11`;
	results22=results12`;
	results23=results13`;
	results2=results21||results22||results23;

create runlengths from rlvec;
append from rlvec;
print 'ModHWMA Fixed Explanatory variable Asymptotic case-(lambda=0.2 and k=lambda/4)';
print sim n lambda k Lb0 Lb1 Le beta0 gam;
print results2[colname=varNames rowname=name2][format=10.2];
quit;



















































































































































































