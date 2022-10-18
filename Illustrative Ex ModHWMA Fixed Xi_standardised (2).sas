proc iml;
**********************************************************************;

*Size of the Phase II sample;
n = 11;

*lambda is the smothing constant;
lambda = 0.2; 
l1 = 1-lambda;
k=lambda/2;

*Process parameters: Case K i.e. distribution parameters are known;
 mu0 = 0;
 sigma0 = 1;

*Shift: for the IC case;
stdev = 1;
df=5;
v=1;
wp=2;
a0=-0.000004167;
a1=-0.5084333;
delta=0;
beta=0;
xi={-1.507556723,-1.206045378,-0.904534034,-0.603022689,-0.301511345,0,	0.301511344577764,0.603022689155527,0.904534033733291,1.20604537831105,1.50755672288882
};
y={1.708891952	1.479229139	1.446932806	1.285451142	1.127557958	0.844067924	0.661055371	0.445746484	0.356034448	0.115606192	0.029482637,
0.043836563	-0.131999028	-0.279126767	-0.508789579	-0.681036688	-0.594913133	-0.738452391	-0.96093824	-0.907111019	-0.93581887	-1.18342409,
1.963674134	1.762719173	1.533056361	1.242389364	1.339278363	1.242389364	0.987607182	1.080907699	0.829713999	0.582108779	0.521104595,
1.450521288	1.482817621	1.350043807	1.15626581	1.027080478	0.865598813	0.83330248	0.610816631	0.302207227	0.147902525	0.025894156,
1.080907699	0.793829184	0.607228149	0.542635483	0.173021895	0.054602007	0.083309859	-0.185826249	-0.304246137	-0.393958173	-0.627209466,
-0.968115203	-1.136773831	-1.305432459	-1.387967532	-1.578157048	-1.793465935	-1.872412526	-1.98006697	-2.216906745	-2.371211447	-2.475277408,
0.524693076	0.277087857	0.158667969	-0.053052436	-0.203768656	-0.408312098	-0.584147689	-0.78510265	-0.957349759	-0.982469129	-1.262370681,
1.396694066	1.045022885	0.951722368	0.765121333	0.639524482	0.33450356	0.244791524	0.126371636	0.0366596	-0.081760287	-0.372427284,
-0.042286992	-0.33654247	-0.566205282	-0.677448207	-0.874814686	-1.043473313	-0.98605761	-1.20854346	-1.291078533	-1.588922492	-1.614041863,
1.428990399	1.364397733	1.20650455	0.923014516	0.643112964	0.675409297	0.312972671	0.133548599	-0.074583324	0.015128712	-0.246830434,
-0.828164427	-1.050650276	-1.161893201	-1.298255496	-1.488445012	-1.596099455	-1.675046047	-1.80781986	-2.044659635	-2.094898376	-2.198964337,
0.668232334	0.452923447	0.269910894	-0.203768656	-0.232476508	-0.494435653	-0.663094281	-0.637974911	-0.946584315	-1.194189534	-1.35925968,
1.080907699	1.145500365	0.987607182	0.829713999	0.829713999	0.70770563	0.524693076	0.463688892	0.370388374	0.251968487	0.072544415,
1.13114644	0.976841738	0.790240703	0.668232334	0.452923447	0.33450356	0.212495191	0.0366596	-0.196591693	-0.340130951	-0.397546654,
1.436167362	1.30698203	1.177796698	1.116792514	0.804594628	0.926602997	0.743590444	0.589285742	0.43498104	0.198141265	0.076132896,
-1.255193718	-1.391556013	-1.527918308	-1.664280603	-1.746815676	-1.958536081	-2.062602043	-2.249203078	-2.274322448	-2.42862715	-2.42862715,
0.539047002	0.115606192	0.025894156	-0.092525732	-0.268361322	-0.41190058	-0.526731986	-0.756394798	-0.896345574	-1.007588499	-1.036296351,
0.467277373	0.194552784	0.076132896	-0.070994843	-0.279126767	-0.566205282	-0.70974454	-0.763571761	-0.903522537	-1.126008387	-1.126008387,
0.887129702	1.045022885	0.793829184	0.639524482	0.549812446	0.33450356	0.305795708	0.0366596	-0.053052436	-0.167883842	-0.401135136,
1.066553774	0.78306374	0.693351704	0.567754853	0.234026079	-0.035110029	-0.210945619	-0.297069174	-0.268361322	-0.383192728	-0.498024134,
2.00314743	1.80578095	1.608414471	1.38234014	1.350043807	1.188562143	0.969664775	0.750767407	0.628759038	0.478042817	0.266322412,
1.626356878	1.360809252	1.167031254	0.940956923	0.90866059	0.750767407	0.625170556	0.28426482	0.101252266	0.01154023	-0.13917599,
1.486406102	1.454109769	1.102438588	0.822537036	0.729236518	0.452923447	0.363211411	0.183787339	0.033071119	0.004363267	-0.081760287,
0.35962293	0.119194673	0.029482637	-0.149941435	-0.236064989	-0.383192728	-0.612855541	-0.838929871	-0.978880648	-1.061415721	-1.197778015	

};
xbar=mean(xi);
xp=xi-xbar*j(n,1,1);
b0=a0+a1*xbar;
b1=a1;
Lb0=13.83; Lb1=14.23; Le=7.268;

MHjb0_1 = b0;MHjb1_1 = b1;MHje_1 = log(Sigma0**2);*Starting values for Z0(b0), Z0(b1) and Z0(e) as defined below Equations (7b and 12)/ Note Z0(e) = ln(sigma^2)=0;

b0j_1=b0; b1j_1=b1; msej_1=log(Sigma0**2);



vecb0={};vecb1={};vece={};
veclimb01={};veclimb02={};
veclimb11={};veclimb12={};
veclime={};
no=0;vecno={};

do i=1 to 24 ; 
            no=no+1;
		    eij=j(n,1,.);
	        * Generating a Phase II sample;

            call randgen(eij,'NORMAL',0,1);
	        y1=y[i,];
			yij=y1`;
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
			b0j=a0j+a1j*xbar;b1j=a1j;
			mub0=b0; Varb0=sigma0**2/n;mub1=b1;Varb1=sigma0**2/sxx;
			b0j=(a0j+delta*stdev)+a1j*xbar;b1j=a1j+(beta*stdev); b0j=a0j;b1j=a1j;
			mub0=a0j+a1j*xbar; Varb0=sigma0**2/n;mub1=a1j;Varb1=sigma0**2/sxx;
			*Residuals;
			Vare=stdev**2/n;
			eij=yij-b0j*j(n,1,1)-b1j*devx;
			ebar=mean(eij);
			msej=sum((eij##2))/(n-2);
            
            
            vecb0=vecb0//b0j;vecb1=vecb1//b1j;vece=vece//msej;
            *Control limits;
            * Time-varying;
		   a=(l1)**2;
			
            if i=1 then vz=(lambda+k)**2;
		   else Vz = ((lambda+k)**2) + (((l1/(i-1))+k)**2) + ((i-2)*(l1/(i-1))**2);
            Vzb0=Vz*Varb0;Vzb1=Vz*Varb1;Vze=(2/(n-2)+2/(n-2)**2+4/(3*(n-2)**3)-16/(15*(n-2)**5))*Vz;
		    lclb0=b0-Lb0*sqrt(Vzb0);uclb0=b0+Lb0*sqrt(Vzb0);
            lclb1=b1-Lb1*sqrt(Vzb1);uclb1=b1+Lb1*sqrt(Vzb1);
            ucle=Le*sqrt(Vze);
			veclimb01=veclimb01//lclb0;
			veclimb02=veclimb02//uclb0;
			veclimb11=veclimb11//lclb1;
			veclimb12=veclimb12//uclb1;
			veclime=veclime//ucle;
			vecno=vecno//no;
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
			vecb0=vecb0//MHjb0;vecb1=vecb1//MHjb1;vece=vece//MHje;

		end;
		
		_vecb0={};_vecb1={};_vece={};
do i=2 to 48 by 2;
_vecb0=_vecb0//vecb0[i,];
_vecb1=_vecb1//vecb1[i,];
_vece=_vece//vece[i,];
end;






print 'EWMA3 scheme for profile';
print n lambda Lb0 Lb1 Le;
print vecno _vecb0 veclimb01 veclimb02 _vecb1 veclimb11 veclimb12 _vece veclime;
dat=vecno|| _vecb0|| veclimb01|| veclimb02|| _vecb1|| veclimb11|| veclimb12|| _vece|| veclime;
cols={'vecno' '_vecb0' 'veclimb01' 'veclimb02' '_vecb1' 'veclimb11' 'veclimb12' '_vece' 'veclime'};

create chart from dat[colname=cols];
append from dat;
quit;





proc sgplot data=chart;
series x=vecno y=veclimb01;
series x=vecno y=veclimb02;
scatter x=vecno y=_vecb0;
run;


proc sgplot data=chart;
series x=vecno y=veclimb11;
series x=vecno y=veclimb12;
scatter x=vecno y=_vecb1;
run;

proc sgplot data=chart;
series x=vecno y=veclime;
scatter x=vecno y=_vece;
run;
