proc iml;
**********************************************************************;
*Number of simulations;
call randseed(123);
sim = 10000;
rlvec = j(sim,11,.);

*Size of the Phase II sample;
n = 4;

*lambda is the smothing constant;
lambda = 0.2; 
l1 = 1-lambda;
k=lambda/2;

*Process parameters: Case K i.e. distribution parameters are known;
 p=2;
 mu0 = {0,0};
 sigma0 = {1 0 , 0 1};
sig2=1;
*Shift: for the IC case, delta=0-for the OOC case, vary delta from 0.25 to 3 with an increment of 0.25;
stdev = 1;
*Model Parameters (a0 and a1) of the in-control process;
a0=3;
a1=2;

*Explanatory variable;
xi={2, 4, 6, 8};                    /*Fixed explanatory variable*/;

xbar=mean(xi);
x=j(n,1,1)||xi;

Lz=13.149;                           /*Control limit constant*/;

shift=0;
do beta0=0.00 to 2 by 0.2 ;  /*Shift in the slope*/;
shift=shift+1;                     /*Indicate the column of the run-length matrix*/;

delta = 0.00;
b0j=(a0+delta*stdev);
b1j=a1+(beta0*stdev);
betashift=b0j//b1j;
df=n-p;                     /*Degrees of freedom*/;
do t=1 to sim;

MHj_1 = {0,0,0};            /*Starting point of the EWMA statistic for profile*/;

vecz={};vecu={};
count=0;signal=0;

do i=1 to 1000000000 until (signal=1); /*When the processm is IC "signal=0" and when it is OOC, "signal=1"*/;
	        yij = j(n,1,.);
	        *y2ij = j(n,1,.);
		    eij=j(n,1,.);
		    *e2ij=j(n,1,.);
	        * Generating a Phase II sample;

            call randgen(eij,'NORMAL',0,1);          /*Error vector*/;
            *call randgen(e2ij,'NORMAL',0,1);
	        yij=eij+a0*j(n,1,1)+a1*xi;               /* Equation (15)*/;
	       
            * Regression Analysis;                  /*Below Equestion (16b)*/;
         inv_x=inv(x`*x);
		 xyj=x`*yij;
		 beta_hatj=inv_x*xyj;             
         sig2j=(1/(n-p))*(yij-x*beta_hatj)`*(yij-x*beta_hatj);/*Below equation (16b)*/;
         tq=sig2j/sig2;
		 npsig=df*tq;
         f=probchi(npsig,df);
         zjsig=quantile('normal',f);             /*Equation (16b)*/;
         zjb=(beta_hatj-betashift)/sqrt(sig2);  /*Equation (16a)*/;
         zj=zjb//zjsig; 
         zj_1=zj;
         /*Below Equation (16b)*/;
         *Getting the covariance matrix below Equation (16b)-its dimension will be p+1 since we added the component of the error variance to be monitored;
		 zerocol={0, 0};
		 zerorow=zerocol`;
		 unit1={1};

		 row1=inv_x||zerocol;
         row2=zerorow||unit1;
		 InvSig0=inv(row1//row2);

            * Asymptotic;

           Vz = ((lambda+k)##2+k##2);
            ucl=Lz*vz;;
 
		   *Obtaining the charting statistics;
           *EWMA for the Intercept Plotting Statistic;
           MHj = lambda*zj + k*(zj-zj_1) + (1-lambda)*MHj_1; /*Equation (17*/;
           MHj_1=MHj;
           uj = MHj`*InvSig0*MHj;			 /* Equation (18)- LHS*/;
			
			                

	        *Comparing the plotting statistics to the control limits;           
			vecz=vecz//zj;vecu=vecu//uj;

            if uj>=ucl then signal=1;         /*Equation (18) - RHS*/;
			else signal=0;

            count=count+1;
            rlvec[t,shift]=count;
			end;
   end;
end;
results=rlvec;
	name1={"0.0" "0.2" "0.4" "0.6" "0.8" "1.0" "1.2" "1.4" "1.6" "1.8" "2.0" };	
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

print 'shift in beta';
print results2[colname=varNames rowname=name1][format=10.2];
print EARL[format=10.2] ESDRL[format=10.2] EMRL[format=10.2];
print 'MEWMA3 SLOPE Asymptotic All Shifts';
print sim n lambda Lz delta;

quit;




proc iml;
**********************************************************************;
*Number of simulations;
call randseed(123);
sim = 10000;
rlvec = j(sim,11,.);

*Size of the Phase II sample;
n = 4;

*lambda is the smothing constant;
lambda = 0.2; 
l1 = 1-lambda;
k=lambda/4;

*Process parameters: Case K i.e. distribution parameters are known;
 p=2;
 mu0 = {0,0};
 sigma0 = {1 0 , 0 1};
sig2=1;
*Shift: for the IC case, delta=0-for the OOC case, vary delta from 0.25 to 3 with an increment of 0.25;
stdev = 1;
*Model Parameters (a0 and a1) of the in-control process;
a0=3;
a1=2;

*Explanatory variable;
xi={2, 4, 6, 8};                    /*Fixed explanatory variable*/;

xbar=mean(xi);
x=j(n,1,1)||xi;

Lz=20.24;                          /*Control limit constant*/;

shift=0;
do beta0=0.00 to 2 by 0.2 ;  /*Shift in the slope*/;
shift=shift+1;                     /*Indicate the column of the run-length matrix*/;

delta = 0.00;
b0j=(a0+delta*stdev);
b1j=a1+(beta0*stdev);
betashift=b0j//b1j;
df=n-p;                     /*Degrees of freedom*/;
do t=1 to sim;

MHj_1 = {0,0,0};            /*Starting point of the EWMA statistic for profile*/;

vecz={};vecu={};
count=0;signal=0;

do i=1 to 1000000000 until (signal=1); /*When the processm is IC "signal=0" and when it is OOC, "signal=1"*/;
	        yij = j(n,1,.);
	        *y2ij = j(n,1,.);
		    eij=j(n,1,.);
		    *e2ij=j(n,1,.);
	        * Generating a Phase II sample;

            call randgen(eij,'NORMAL',0,1);          /*Error vector*/;
            *call randgen(e2ij,'NORMAL',0,1);
	        yij=eij+a0*j(n,1,1)+a1*xi;               /* Equation (15)*/;
	       
            * Regression Analysis;                  /*Below Equestion (16b)*/;
         inv_x=inv(x`*x);
		 xyj=x`*yij;
		 beta_hatj=inv_x*xyj;             
         sig2j=(1/(n-p))*(yij-x*beta_hatj)`*(yij-x*beta_hatj);/*Below equation (16b)*/;
         tq=sig2j/sig2;
		 npsig=df*tq;
         f=probchi(npsig,df);
         zjsig=quantile('normal',f);             /*Equation (16b)*/;
         zjb=(beta_hatj-betashift)/sqrt(sig2);  /*Equation (16a)*/;
         zj=zjb//zjsig; 
         zj_1=zj;
         /*Below Equation (16b)*/;
         *Getting the covariance matrix below Equation (16b)-its dimension will be p+1 since we added the component of the error variance to be monitored;
		 zerocol={0, 0};
		 zerorow=zerocol`;
		 unit1={1};

		 row1=inv_x||zerocol;
         row2=zerorow||unit1;
		 InvSig0=inv(row1//row2);

            * Asymptotic;

           Vz = ((lambda+k)##2+k##2);
            ucl=Lz*vz;;
 
		   *Obtaining the charting statistics;
           *EWMA for the Intercept Plotting Statistic;
           MHj = lambda*zj + k*(zj-zj_1) + (1-lambda)*MHj_1; /*Equation (17*/;
           MHj_1=MHj;
           uj = MHj`*InvSig0*MHj;			 /* Equation (18)- LHS*/;
			
			                

	        *Comparing the plotting statistics to the control limits;           
			vecz=vecz//zj;vecu=vecu//uj;

            if uj>=ucl then signal=1;         /*Equation (18) - RHS*/;
			else signal=0;

            count=count+1;
            rlvec[t,shift]=count;
			end;
   end;
end;
results=rlvec;
	name1={"0.0" "0.2" "0.4" "0.6" "0.8" "1.0" "1.2" "1.4" "1.6" "1.8" "2.0" };	
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

print 'shift in beta';
print results2[colname=varNames rowname=name1][format=10.2];
print EARL[format=10.2] ESDRL[format=10.2] EMRL[format=10.2];
print 'MEWMA3 SLOPE Asymptotic All Shifts';
print sim n lambda Lz delta;

quit;





proc iml;
**********************************************************************;
*Number of simulations;
call randseed(123);
sim = 10000;
rlvec = j(sim,11,.);

*Size of the Phase II sample;
n = 4;

*lambda is the smothing constant;
lambda = 0.2; 
l1 = 1-lambda;
k=0;

*Process parameters: Case K i.e. distribution parameters are known;
 p=2;
 mu0 = {0,0};
 sigma0 = {1 0 , 0 1};
sig2=1;
*Shift: for the IC case, delta=0-for the OOC case, vary delta from 0.25 to 3 with an increment of 0.25;
stdev = 1;
*Model Parameters (a0 and a1) of the in-control process;
a0=3;
a1=2;

*Explanatory variable;
xi={2, 4, 6, 8};                    /*Fixed explanatory variable*/;

xbar=mean(xi);
x=j(n,1,1)||xi;

Lz=32.873;                          /*Control limit constant*/;

shift=0;
do beta0=0.00 to 2 by 0.2 ;  /*Shift in the slope*/;
shift=shift+1;                     /*Indicate the column of the run-length matrix*/;

delta = 0.00;
b0j=(a0+delta*stdev);
b1j=a1+(beta0*stdev);
betashift=b0j//b1j;
df=n-p;                     /*Degrees of freedom*/;
do t=1 to sim;

MHj_1 = {0,0,0};            /*Starting point of the EWMA statistic for profile*/;

vecz={};vecu={};
count=0;signal=0;

do i=1 to 1000000000 until (signal=1); /*When the processm is IC "signal=0" and when it is OOC, "signal=1"*/;
	        yij = j(n,1,.);
	        *y2ij = j(n,1,.);
		    eij=j(n,1,.);
		    *e2ij=j(n,1,.);
	        * Generating a Phase II sample;

            call randgen(eij,'NORMAL',0,1);          /*Error vector*/;
            *call randgen(e2ij,'NORMAL',0,1);
	        yij=eij+a0*j(n,1,1)+a1*xi;               /* Equation (15)*/;
	       
            * Regression Analysis;                  /*Below Equestion (16b)*/;
         inv_x=inv(x`*x);
		 xyj=x`*yij;
		 beta_hatj=inv_x*xyj;             
         sig2j=(1/(n-p))*(yij-x*beta_hatj)`*(yij-x*beta_hatj);/*Below equation (16b)*/;
         tq=sig2j/sig2;
		 npsig=df*tq;
         f=probchi(npsig,df);
         zjsig=quantile('normal',f);             /*Equation (16b)*/;
         zjb=(beta_hatj-betashift)/sqrt(sig2);  /*Equation (16a)*/;
         zj=zjb//zjsig; 
         zj_1=zj;
         /*Below Equation (16b)*/;
         *Getting the covariance matrix below Equation (16b)-its dimension will be p+1 since we added the component of the error variance to be monitored;
		 zerocol={0, 0};
		 zerorow=zerocol`;
		 unit1={1};

		 row1=inv_x||zerocol;
         row2=zerorow||unit1;
		 InvSig0=inv(row1//row2);

            * Asymptotic;

           Vz = ((lambda+k)##2+k##2);
            ucl=Lz*vz;;
 
		   *Obtaining the charting statistics;
           *EWMA for the Intercept Plotting Statistic;
           MHj = lambda*zj + k*(zj-zj_1) + (1-lambda)*MHj_1; /*Equation (17*/;
           MHj_1=MHj;
           uj = MHj`*InvSig0*MHj;			 /* Equation (18)- LHS*/;
			
			                

	        *Comparing the plotting statistics to the control limits;           
			vecz=vecz//zj;vecu=vecu//uj;

            if uj>=ucl then signal=1;         /*Equation (18) - RHS*/;
			else signal=0;

            count=count+1;
            rlvec[t,shift]=count;
			end;
   end;
end;
results=rlvec;
	name1={"0.0" "0.2" "0.4" "0.6" "0.8" "1.0" "1.2" "1.4" "1.6" "1.8" "2.0" };	
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

print 'shift in beta';
print results2[colname=varNames rowname=name1][format=10.2];
print EARL[format=10.2] ESDRL[format=10.2] EMRL[format=10.2];
print 'MEWMA3 SLOPE Asymptotic All Shifts';
print sim n lambda Lz delta;

quit;





proc iml;
**********************************************************************;
*Number of simulations;
call randseed(123);
sim = 10000;
rlvec = j(sim,11,.);

*Size of the Phase II sample;
n = 4;

*lambda is the smothing constant;
lambda = 0.2; 
l1 = 1-lambda;
k=-lambda/4;

*Process parameters: Case K i.e. distribution parameters are known;
 p=2;
 mu0 = {0,0};
 sigma0 = {1 0 , 0 1};
sig2=1;
*Shift: for the IC case, delta=0-for the OOC case, vary delta from 0.25 to 3 with an increment of 0.25;
stdev = 1;
*Model Parameters (a0 and a1) of the in-control process;
a0=3;
a1=2;

*Explanatory variable;
xi={2, 4, 6, 8};                    /*Fixed explanatory variable*/;

xbar=mean(xi);
x=j(n,1,1)||xi;

Lz=52.6;                           /*Control limit constant*/;

shift=0;
do beta0=0.00 to 2 by 0.2 ;  /*Shift in the slope*/;
shift=shift+1;                     /*Indicate the column of the run-length matrix*/;

delta = 0.00;
b0j=(a0+delta*stdev);
b1j=a1+(beta0*stdev);
betashift=b0j//b1j;
df=n-p;                     /*Degrees of freedom*/;
do t=1 to sim;

MHj_1 = {0,0,0};            /*Starting point of the EWMA statistic for profile*/;

vecz={};vecu={};
count=0;signal=0;

do i=1 to 1000000000 until (signal=1); /*When the processm is IC "signal=0" and when it is OOC, "signal=1"*/;
	        yij = j(n,1,.);
	        *y2ij = j(n,1,.);
		    eij=j(n,1,.);
		    *e2ij=j(n,1,.);
	        * Generating a Phase II sample;

            call randgen(eij,'NORMAL',0,1);          /*Error vector*/;
            *call randgen(e2ij,'NORMAL',0,1);
	        yij=eij+a0*j(n,1,1)+a1*xi;               /* Equation (15)*/;
	       
            * Regression Analysis;                  /*Below Equestion (16b)*/;
         inv_x=inv(x`*x);
		 xyj=x`*yij;
		 beta_hatj=inv_x*xyj;             
         sig2j=(1/(n-p))*(yij-x*beta_hatj)`*(yij-x*beta_hatj);/*Below equation (16b)*/;
         tq=sig2j/sig2;
		 npsig=df*tq;
         f=probchi(npsig,df);
         zjsig=quantile('normal',f);             /*Equation (16b)*/;
         zjb=(beta_hatj-betashift)/sqrt(sig2);  /*Equation (16a)*/;
         zj=zjb//zjsig; 
         zj_1=zj;
         /*Below Equation (16b)*/;
         *Getting the covariance matrix below Equation (16b)-its dimension will be p+1 since we added the component of the error variance to be monitored;
		 zerocol={0, 0};
		 zerorow=zerocol`;
		 unit1={1};

		 row1=inv_x||zerocol;
         row2=zerorow||unit1;
		 InvSig0=inv(row1//row2);

            * Asymptotic;

           Vz = ((lambda+k)##2+k##2);
            ucl=Lz*vz;;
 
		   *Obtaining the charting statistics;
           *EWMA for the Intercept Plotting Statistic;
           MHj = lambda*zj + k*(zj-zj_1) + (1-lambda)*MHj_1; /*Equation (17*/;
           MHj_1=MHj;
           uj = MHj`*InvSig0*MHj;			 /* Equation (18)- LHS*/;
			
			                

	        *Comparing the plotting statistics to the control limits;           
			vecz=vecz//zj;vecu=vecu//uj;

            if uj>=ucl then signal=1;         /*Equation (18) - RHS*/;
			else signal=0;

            count=count+1;
            rlvec[t,shift]=count;
			end;
   end;
end;
results=rlvec;
	name1={"0.0" "0.2" "0.4" "0.6" "0.8" "1.0" "1.2" "1.4" "1.6" "1.8" "2.0" };	
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

print 'shift in beta';
print results2[colname=varNames rowname=name1][format=10.2];
print EARL[format=10.2] ESDRL[format=10.2] EMRL[format=10.2];
print 'MEWMA3 SLOPE Asymptotic All Shifts';
print sim n lambda Lz delta;

quit;





proc iml;
**********************************************************************;
*Number of simulations;
call randseed(123);
sim = 10000;
rlvec = j(sim,11,.);

*Size of the Phase II sample;
n = 4;

*lambda is the smothing constant;
lambda = 0.2; 
l1 = 1-lambda;
k=-lambda/2;

*Process parameters: Case K i.e. distribution parameters are known;
 p=2;
 mu0 = {0,0};
 sigma0 = {1 0 , 0 1};
sig2=1;
*Shift: for the IC case, delta=0-for the OOC case, vary delta from 0.25 to 3 with an increment of 0.25;
stdev = 1;
*Model Parameters (a0 and a1) of the in-control process;
a0=3;
a1=2;

*Explanatory variable;
xi={2, 4, 6, 8};                    /*Fixed explanatory variable*/;

xbar=mean(xi);
x=j(n,1,1)||xi;

Lz=65.765;                           /*Control limit constant*/;

shift=0;
do beta0=0.00 to 2 by 0.2 ;  /*Shift in the slope*/;
shift=shift+1;                     /*Indicate the column of the run-length matrix*/;

delta = 0.00;
b0j=(a0+delta*stdev);
b1j=a1+(beta0*stdev);
betashift=b0j//b1j;
df=n-p;                     /*Degrees of freedom*/;
do t=1 to sim;

MHj_1 = {0,0,0};            /*Starting point of the EWMA statistic for profile*/;

vecz={};vecu={};
count=0;signal=0;

do i=1 to 1000000000 until (signal=1); /*When the processm is IC "signal=0" and when it is OOC, "signal=1"*/;
	        yij = j(n,1,.);
	        *y2ij = j(n,1,.);
		    eij=j(n,1,.);
		    *e2ij=j(n,1,.);
	        * Generating a Phase II sample;

            call randgen(eij,'NORMAL',0,1);          /*Error vector*/;
            *call randgen(e2ij,'NORMAL',0,1);
	        yij=eij+a0*j(n,1,1)+a1*xi;               /* Equation (15)*/;
	       
            * Regression Analysis;                  /*Below Equestion (16b)*/;
         inv_x=inv(x`*x);
		 xyj=x`*yij;
		 beta_hatj=inv_x*xyj;             
         sig2j=(1/(n-p))*(yij-x*beta_hatj)`*(yij-x*beta_hatj);/*Below equation (16b)*/;
         tq=sig2j/sig2;
		 npsig=df*tq;
         f=probchi(npsig,df);
         zjsig=quantile('normal',f);             /*Equation (16b)*/;
         zjb=(beta_hatj-betashift)/sqrt(sig2);  /*Equation (16a)*/;
         zj=zjb//zjsig; 
         zj_1=zj;
         /*Below Equation (16b)*/;
         *Getting the covariance matrix below Equation (16b)-its dimension will be p+1 since we added the component of the error variance to be monitored;
		 zerocol={0, 0};
		 zerorow=zerocol`;
		 unit1={1};

		 row1=inv_x||zerocol;
         row2=zerorow||unit1;
		 InvSig0=inv(row1//row2);

            * Asymptotic;

           Vz = ((lambda+k)##2+k##2);
            ucl=Lz*vz;;
 
		   *Obtaining the charting statistics;
           *EWMA for the Intercept Plotting Statistic;
           MHj = lambda*zj + k*(zj-zj_1) + (1-lambda)*MHj_1; /*Equation (17*/;
           MHj_1=MHj;
           uj = MHj`*InvSig0*MHj;			 /* Equation (18)- LHS*/;
			
			                

	        *Comparing the plotting statistics to the control limits;           
			vecz=vecz//zj;vecu=vecu//uj;

            if uj>=ucl then signal=1;         /*Equation (18) - RHS*/;
			else signal=0;

            count=count+1;
            rlvec[t,shift]=count;
			end;
   end;
end;
results=rlvec;
	name1={"0.0" "0.2" "0.4" "0.6" "0.8" "1.0" "1.2" "1.4" "1.6" "1.8" "2.0" };	
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

print 'shift in beta';
print results2[colname=varNames rowname=name1][format=10.2];
print EARL[format=10.2] ESDRL[format=10.2] EMRL[format=10.2];
print 'MEWMA3 SLOPE Asymptotic All Shifts';
print sim n lambda Lz delta;

quit;
