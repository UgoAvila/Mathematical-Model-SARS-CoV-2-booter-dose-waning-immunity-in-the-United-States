load('BestFitParameters')

%% Setting initial values for Dec-13-2020:
IniVal = [S1(57); E1(57); I1(57); A1(57); H1_exp(57); 0; 0; 0; 0; 0; 0; ...
    R(57); W(57); 0; 0; D_exp(57); Icumul_exp(57); Acumul_exp(57); 0; 0; D_exp(57); 0];

dateslong = dates(327:714);  % Dec-13-2020 to Jan-04-2022
tsim=linspace(0,length(dateslong)-1,length(dateslong));
Icumul_exp_vacc=C_exp_US(57+initdate-1:714);
D_exp_vacc=D_exp_US(57+initdate-1:714);
H1_exp_vacc=Hosp_exp_US(57+initdate-diffDates-1:end);

global alpha1 chi omega beta1 beta2 Cm epsilonm kappa epsa1 epsa2 epsa3 ...
    SumEpsLi SumEpsL2i Sumthetai Kv Kv1 Kv2 eps3i p2 gamma2 alphaB chiB omegaB Cr

%% Parameter settings

% Parameter values obtained from data fitting:
beta1orig = x(1)/( 1 - 0.5*0.5); % Transmission rate for Wuhan variant (symptomatic)
beta2orig = x(2)/( 1 - 0.5*0.5); % Transmission rate for Wuhan variant (asymptomatic)
alpha1 = x(3);       % Hospitalization rate
chi = x(4);          % Death probability for H
omega = x(5);        % 1/(average hospitalization period)

% Time-varying transmission rates:
beta1 = @(X) beta1orig*(X<=125) ...                 % Second wave
          + 2.8*beta1orig.*( (X>125)&(X<=364) ) ... % Third wave
          + 5.0*beta1orig.*( (X>364) ) ;            % Fourth wave (Omicron)

beta2 = @(X) beta2orig*(X<=125) ...                 % Second wave
          + 2.8*beta2orig.*( (X>125)&(X<=364) ) ... % Third wave
          + 5.0*beta2orig.*( (X>364) ) ;            % Fourth wave (Omicron)

% Other parameter values:

Cm = 0.5;   % Face mask compliance in the community
Cr=0;
epsilonm = 0.5;     % Efficacy of face masks
kappa = 0.52;   % Relative infectivity of breakthrough infections
epsa1 = 7.23e-4;    % All-or-nothing (Moderna)
epsa2 = 3.77e-4;    % All-or-nothing (Pfizer)
epsa3 = 0.0091;     % All-or-nothing (Janssen J&J)
Kv1 = 1/60;     % Waning rate of immunity after 1 dose
Kv2 = 1/180;    % Waning rate of immunity after 2 doses
Kv = 0.0037;    % Waning rate of acquired immunity
Kn=1/270;     % Waning rate of natural immunity
eps3i = 1-0.37;   % 1 - (effectiveness of booster dose)
EpsL = 0.26;    % Effectiveness for 1 dose
EpsL2 = 0.19;   % Effectiveness for 2 doses
SumEpsLi = 1 - EpsL;    % Sum of 1-epsilon_L,i (i=1,2,3)
SumEpsL2i = 1 - EpsL2;   % Sum of 1-epsilon_L2,i (i=1,2,3)
Sumthetai = 2/180;      % Sum of theta1 and theta2

chiB = 0.000172;          % Death probability for HB (breakthrough infections)
omegaB = 2; %0.0514;       % 1/(average hospitalization period) (breakthrough infections)
alphaB = 0.0118;        % Hospitalization rate (breakthrough infections)
p2 = 1/9;   % Proportion of symptomatic infections (breakthrough infections)
gamma2 = 1/14;  % Recovery rate (breakthrough infections)

%% Defining vaccination rates
rho1 = @(x) (2.1268e+04)*x.*exp(-0.015479*x);
rho2 = @(x) (1.8017e+04)*x.*exp(-0.010932*x);
rho3 = @(x) (8.6062e+03)*(x-81).*exp(-0.023304*(x-81)).*( x>=81 );
rho4 = @(x) (1e-03).*( x>=242 );    % Booster doses began on Aug-13-3021

%% Solving the system
[tsim,S1a,E1a,I1a,A1a,H1a,V1a,V2a,EBa,IBa,ABa,HBa,Ra,Wa,Wva,V3a,Da, ...
    I1cumula,A1cumula,IBcumula,ABcumula,D1a,DBa]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, tsim, IniVal);

%% Graphs
figure
subplot(3,1,1)
plot(dateslong,I1cumula+IBcumula,'-r','LineWidth',3)
hold on
plot(dateslong,Icumul_exp_vacc,'--b')
legend({'Cumul. symptomatic infected','Reported data'},'Location','northwest')
grid on
subplot(3,1,2)
plot(dateslong,H1a+HBa,'-r','LineWidth',3)
hold on
plot(dateslong(1:380),H1_exp_vacc,'--b')
legend({'Hospitalizations (H1+HB)','Reported data'},'Location','northwest')
grid on
subplot(3,1,3)
plot(dateslong,Da,'-r','LineWidth',3)
hold on
plot(dateslong,D_exp_vacc,'--b')
legend({'Deaths (D)','Reported data'},'Location','northwest')
grid on

figure
subplot(2,2,1)
plot(dateslong,I1a,'-r','LineWidth',3)
hold on
plot(dateslong,IBa,'--b','LineWidth',3)
legend({'I1','IB'},'Location','northwest')
grid on

subplot(2,2,2)
plot(dateslong,A1a,'-r','LineWidth',3)
hold on
plot(dateslong,ABa,'--b','LineWidth',3)
legend({'A1','AB'},'Location','northwest')
grid on

subplot(2,2,3)
plot(dateslong,H1a,'-r','LineWidth',3)
hold on
plot(dateslong,HBa,'--b','LineWidth',3)
legend({'H1','HB'},'Location','northwest')
grid on

subplot(2,2,4)
plot(dateslong,V1a,'-r','LineWidth',3)
hold on
plot(dateslong,V2a,'--b','LineWidth',3)
plot(dateslong,V3a,'-.k','LineWidth',3)
legend({'V1','V2','V3'},'Location','northwest')
grid on

figure
plot(dateslong,D1a,'-r','LineWidth',3)
hold on
plot(dateslong,DBa,'--b','LineWidth',3)
legend({'Deaths from COVID-19 Unvaccinated','Deaths from COVID-19 Vaccinated'},'Location','northwest')
grid on

%% Seting initial conditions for simulations

S10 = S1a(end);
E10 = E1a(end);
I10 = I1a(end);
A10 = A1a(end);
H10 = H1a(end);
V10 = V1a(end);
V20 = V2a(end);
EB0 = EBa(end);
IB0 = IBa(end);
AB0 = ABa(end);
HB0 = HBa(end);
R0 = Ra(end);
W0 = Wa(end);
WV0= Wva(end);
V30= V3a(end);
D0 = 832180;
I1cumula0 = I1cumula(end);
A1cumula0 = A1cumula(end);
IBcumula0 = IBcumula(end);
ABcumula0 = ABcumula(end);
D10=D0*D1a(end)/( D1a(end)+DBa(end) );
DB0=D0*DBa(end)/( D1a(end)+DBa(end) );

IniVal = [S10;E10;I10;A10;H10;V10;V20;EB0;IB0;AB0;HB0;R0;...
    W0;WV0;V30;D0;I1cumula0;A1cumula0;IBcumula0;ABcumula0;D10;DB0];

newinitdate=714;   % 4 Jan 2022

numberofdays = 300;
t=linspace(tsim(end), tsim(end)+numberofdays,numberofdays+1);
newdata_long=dateshift(dateslong(388),'start','day',0:numberofdays);

%% Defining vaccination rates
rho1 = @(x) (2.1268e+04)*x.*exp(-0.015479*x);
rho2 = @(x) (1.8017e+04)*x.*exp(-0.010932*x);
rho3 = @(x) (8.6062e+03)*(x-81).*exp(-0.023304*(x-81)).*( x>=81 );
rho4 = @(x) (1e-03).*( x>=242 );    % Booster doses began on Aug-13-3021

%% Solving the system
[t,S10,E10,I10,A10,H10,V10,V20,EB0,IB0,AB0,HB0,R0,W0,Wv0,V30,D0, ...
    I1cumula0,A1cumula0,IBcumula0,ABcumula0,D10,DB0]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);

figure
subplot(2,2,1)
plot(newdata_long,I10,'-r','LineWidth',3)
hold on
plot(newdata_long,IB0,'--b','LineWidth',3)
title('Symptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'Infected Unvaccinated','Infected Vaccinated'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,A10,'-r','LineWidth',3)
hold on
plot(newdata_long,AB0,'--b','LineWidth',3)
title('Asymptomatic Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
ylabel('Individuals')
legend({'Asymptomatic Unvaccinated','Asymptomatic Vaccinated'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,H10,'-r','LineWidth',3)
hold on
plot(newdata_long,HB0,'--b','LineWidth',3)
title('Hospitalized Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
ylabel('Individuals')
legend({'Hospitalized Unvaccinated','Hospitalized Vaccinated'},'Location','northwest')
grid on

subplot(2,2,4)
plot(newdata_long,V10,'-r','LineWidth',3)
hold on
plot(newdata_long,V20,'--b','LineWidth',3)
plot(newdata_long,V30,'-.k','LineWidth',3)
title('Vaccinated Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfD','FontSize',12)
legend({'One Dose','Two Dose','Booster'},'Location','northwest')
grid on

figure
plot(newdata_long,H10./E10,'-r','LineWidth',3)
hold on
plot(newdata_long,HB0./EB0,'--b','LineWidth',3)
ylabel('Hospitalization rate')
legend({'Hospitalized/Exposed Unvaccinated','Hospitalized/Exposed Vaccinated'},'Location','northwest')
grid on

%%%%%%%%%%%%%%

figure
subplot(2,2,1)
plot(newdata_long,R0,'-g','LineWidth',3)
hold on
text(0.5,0.5,'\bfA','FontSize',12)
title('Recuperation Individuals')
ylabel('Individuals')
legend({'Recuperated regardless of vaccination status'},'Location','northwest')
grid on

subplot(2,2,2)
plot(newdata_long,W0,'-r','LineWidth',3)
hold on
plot(newdata_long,Wv0,'--b','LineWidth',3)
title('Waning Immunity')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfA','FontSize',12)
legend({'Waning Natural Immunity','Waning Vaccine Immunity'},'Location','northwest')
grid on

subplot(2,2,3)
plot(newdata_long,D0,'-r','LineWidth',3)
hold on
title('Death Toll')
ylabel('Individuals')
legend({'Deaths from COVID-19'},'Location','northeast')
text(0.5,0.5,'\bfB','FontSize',12)
grid on

subplot(2,2,4)
plot(newdata_long,D10,'-r','LineWidth',3)
hold on
plot(newdata_long,DB0,'--b','LineWidth',3)
title('Death Toll Vaccination Status')
ylabel('Individuals')
legend({'Deaths from COVID-19 Unvaccinated','Deaths from COVID-19 Vaccinated'},'Location','northeast')
text(0.5,0.5,'\bfB','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
grid on


%% Lets vary the usage of face masks where they provide a 30% of protection
Cm=0; 
epsilonm=0.3;
[t,S1o,E1o,I1o,A1o,H1o,V1o,V2o,EBo,IBo,ABo,HBo,Ro,Wo,Wvo,V3o,Do, ...
    I1cumulao,A1cumulao,IBcumulao,ABcumulao,D1o,DBo]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.25; 
epsilonm=0.3;
[t,S1a,E1a,I1a,A1a,H1a,V1a,V2a,EBa,IBa,ABa,HBa,Ra,Wa,Wva,V3a,Da, ...
    I1cumulaa,A1cumulaa,IBcumulaa,ABcumulaa,D1a,DBa]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.50; 
epsilonm=0.3;
[t,S1b,E1b,I1b,A1b,H1b,V1b,V2b,EBb,IBb,ABb,HBb,Rb,Wb,Wvb,V3b,Db, ...
    I1cumulab,A1cumulab,IBcumulab,ABcumulab,D1b,DBb]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.75; 
epsilonm=0.3;
[t,S1c,E1c,I1c,A1c,H1c,V1c,V2c,EBc,IBc,ABc,HBc,Rc,Wc,Wvc,V3c,Dc, ...
    I1cumulac,A1cumulac,IBcumulac,ABcumulac,D1c,DBc]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=1.0; 
epsilonm=0.3;
[t,S1d,E1d,I1d,A1d,H1d,V1d,V2d,EBd,IBd,ABd,HBd,Rd,Wd,Wvd,V3d,Dd, ...
    I1cumulad,A1cumulad,IBcumulad,ABcumulad,D1d,DBd]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);


figure
subplot(2,2,1)
plot(newdata_long,I1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,I1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,I1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,I1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,I1d,'-.k','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,A1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,A1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,A1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,A1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,A1d,'-.k','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,H1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,H1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,H1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,H1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,H1d,'-.k','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,4)
plot(newdata_long,D1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,D1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,D1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,D1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,D1d,'-.k','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

% Plot for Vaccinated Individuals 
figure
subplot(2,2,1)
plot(newdata_long,IBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,IBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,IBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,IBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,IBd,'-.k','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,ABo,'-.r','LineWidth',2)
hold on
plot(newdata_long,ABa,'-.m','LineWidth',2)
hold on
plot(newdata_long,ABb,'-.b','LineWidth',2)
hold on
plot(newdata_long,ABc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,ABd,'-.k','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,HBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,HBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,HBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,HBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,HBd,'-.k','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,4)
plot(newdata_long,DBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,DBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,DBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,DBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,DBd,'-.k','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

%% Lets vary the usage of face masks where they provide a 50% of protection
Cm=0; 
epsilonm=0.5;
[t,S1o,E1o,I1o,A1o,H1o,V1o,V2o,EBo,IBo,ABo,HBo,Ro,Wo,Wvo,V3o,Do, ...
    I1cumulao,A1cumulao,IBcumulao,ABcumulao,D1o,DBo]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.25; 
epsilonm=0.5;
[t,S1a,E1a,I1a,A1a,H1a,V1a,V2a,EBa,IBa,ABa,HBa,Ra,Wa,Wva,V3a,Da, ...
    I1cumulaa,A1cumulaa,IBcumulaa,ABcumulaa,D1a,DBa]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.50; 
epsilonm=0.5;
[t,S1b,E1b,I1b,A1b,H1b,V1b,V2b,EBb,IBb,ABb,HBb,Rb,Wb,Wvb,V3b,Db, ...
    I1cumulab,A1cumulab,IBcumulab,ABcumulab,D1b,DBb]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.75; 
epsilonm=0.5;
[t,S1c,E1c,I1c,A1c,H1c,V1c,V2c,EBc,IBc,ABc,HBc,Rc,Wc,Wvc,V3c,Dc, ...
    I1cumulac,A1cumulac,IBcumulac,ABcumulac,D1c,DBc]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=1.0; 
epsilonm=0.5;
[t,S1d,E1d,I1d,A1d,H1d,V1d,V2d,EBd,IBd,ABd,HBd,Rd,Wd,Wvd,V3d,Dd, ...
    I1cumulad,A1cumulad,IBcumulad,ABcumulad,D1d,DBd]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);


figure
subplot(2,2,1)
plot(newdata_long,I1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,I1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,I1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,I1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,I1d,'-.k','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,A1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,A1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,A1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,A1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,A1d,'-.k','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,H1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,H1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,H1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,H1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,H1d,'-.k','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,4)
plot(newdata_long,D1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,D1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,D1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,D1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,D1d,'-.k','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

% Plot for Vaccinated Individuals 
figure
subplot(2,2,1)
plot(newdata_long,IBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,IBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,IBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,IBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,IBd,'-.k','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,ABo,'-.r','LineWidth',2)
hold on
plot(newdata_long,ABa,'-.m','LineWidth',2)
hold on
plot(newdata_long,ABb,'-.b','LineWidth',2)
hold on
plot(newdata_long,ABc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,ABd,'-.k','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,HBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,HBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,HBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,HBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,HBd,'-.k','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,4)
plot(newdata_long,DBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,DBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,DBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,DBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,DBd,'-.k','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

%% Lets vary the usage of face masks where they provide a 70% of protection
Cm=0; 
epsilonm=0.7;
[t,S1o,E1o,I1o,A1o,H1o,V1o,V2o,EBo,IBo,ABo,HBo,Ro,Wo,Wvo,V3o,Do, ...
    I1cumulao,A1cumulao,IBcumulao,ABcumulao,D1o,DBo]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.25; 
epsilonm=0.7;
[t,S1a,E1a,I1a,A1a,H1a,V1a,V2a,EBa,IBa,ABa,HBa,Ra,Wa,Wva,V3a,Da, ...
    I1cumulaa,A1cumulaa,IBcumulaa,ABcumulaa,D1a,DBa]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.50; 
epsilonm=0.7;
[t,S1b,E1b,I1b,A1b,H1b,V1b,V2b,EBb,IBb,ABb,HBb,Rb,Wb,Wvb,V3b,Db, ...
    I1cumulab,A1cumulab,IBcumulab,ABcumulab,D1b,DBb]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.75; 
epsilonm=0.7;
[t,S1c,E1c,I1c,A1c,H1c,V1c,V2c,EBc,IBc,ABc,HBc,Rc,Wc,Wvc,V3c,Dc, ...
    I1cumulac,A1cumulac,IBcumulac,ABcumulac,D1c,DBc]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=1.0; 
epsilonm=0.7;
[t,S1d,E1d,I1d,A1d,H1d,V1d,V2d,EBd,IBd,ABd,HBd,Rd,Wd,Wvd,V3d,Dd, ...
    I1cumulad,A1cumulad,IBcumulad,ABcumulad,D1d,DBd]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);


figure
subplot(2,2,1)
plot(newdata_long,I1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,I1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,I1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,I1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,I1d,'-.k','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,A1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,A1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,A1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,A1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,A1d,'-.k','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,H1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,H1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,H1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,H1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,H1d,'-.k','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,4)
plot(newdata_long,D1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,D1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,D1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,D1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,D1d,'-.k','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

% Plot for Vaccinated Individuals 
figure
subplot(2,2,1)
plot(newdata_long,IBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,IBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,IBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,IBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,IBd,'-.k','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,ABo,'-.r','LineWidth',2)
hold on
plot(newdata_long,ABa,'-.m','LineWidth',2)
hold on
plot(newdata_long,ABb,'-.b','LineWidth',2)
hold on
plot(newdata_long,ABc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,ABd,'-.k','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,HBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,HBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,HBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,HBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,HBd,'-.k','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,4)
plot(newdata_long,DBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,DBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,DBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,DBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,DBd,'-.k','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

%% Lets vary the usage of face masks where they provide a 95% of protection
Cm=0; 
epsilonm=0.95;
[t,S1o,E1o,I1o,A1o,H1o,V1o,V2o,EBo,IBo,ABo,HBo,Ro,Wo,Wvo,V3o,Do, ...
    I1cumulao,A1cumulao,IBcumulao,ABcumulao,D1o,DBo]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.25; 
epsilonm=0.95;
[t,S1a,E1a,I1a,A1a,H1a,V1a,V2a,EBa,IBa,ABa,HBa,Ra,Wa,Wva,V3a,Da, ...
    I1cumulaa,A1cumulaa,IBcumulaa,ABcumulaa,D1a,DBa]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.50; 
epsilonm=0.95;
[t,S1b,E1b,I1b,A1b,H1b,V1b,V2b,EBb,IBb,ABb,HBb,Rb,Wb,Wvb,V3b,Db, ...
    I1cumulab,A1cumulab,IBcumulab,ABcumulab,D1b,DBb]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.75; 
epsilonm=0.95;
[t,S1c,E1c,I1c,A1c,H1c,V1c,V2c,EBc,IBc,ABc,HBc,Rc,Wc,Wvc,V3c,Dc, ...
    I1cumulac,A1cumulac,IBcumulac,ABcumulac,D1c,DBc]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=1.0; 
epsilonm=0.95;
[t,S1d,E1d,I1d,A1d,H1d,V1d,V2d,EBd,IBd,ABd,HBd,Rd,Wd,Wvd,V3d,Dd, ...
    I1cumulad,A1cumulad,IBcumulad,ABcumulad,D1d,DBd]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);


figure
subplot(2,2,1)
plot(newdata_long,I1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,I1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,I1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,I1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,I1d,'-.k','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,A1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,A1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,A1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,A1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,A1d,'-.k','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,H1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,H1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,H1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,H1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,H1d,'-.k','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,4)
plot(newdata_long,D1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,D1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,D1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,D1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,D1d,'-.k','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

% Plot for Vaccinated Individuals 
figure
subplot(2,2,1)
plot(newdata_long,IBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,IBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,IBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,IBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,IBd,'-.k','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,ABo,'-.r','LineWidth',2)
hold on
plot(newdata_long,ABa,'-.m','LineWidth',2)
hold on
plot(newdata_long,ABb,'-.b','LineWidth',2)
hold on
plot(newdata_long,ABc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,ABd,'-.k','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,HBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,HBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,HBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,HBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,HBd,'-.k','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,4)
plot(newdata_long,DBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,DBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,DBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,DBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,DBd,'-.k','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on




%% Modify the vaccination rate, first enhancing the one doses/two doses and then boosters.
rho1a = @(x) 0.5*rho1(x);
rho2a = @(x) 0.5*rho2(x);
rho3a = @(x) 0.5*rho3(x);
rho4 = @(x) (1e-03).*( x>=242 );    % Booster doses began on Aug-13-3021
[t,S1o,E1o,I1o,A1o,H1o,V1o,V2o,EBo,IBo,ABo,HBo,Ro,Wo,Wvo,V3o,Do, ...
    I1cumulao,A1cumulao,IBcumulao,ABcumulao,D1o,DBo]= ...
   USA_covid_vaccine_solver(rho1a, rho2a, rho3a, rho4, t, IniVal);

rho1 = @(x) (2.1268e+04)*x.*exp(-0.015479*x);
rho2 = @(x) (1.8017e+04)*x.*exp(-0.010932*x);
rho3 = @(x) (8.6062e+03)*(x-81).*exp(-0.023304*(x-81)).*( x>=81 );
rho4 = @(x) (1e-03).*( x>=242 );    % Booster doses began on Aug-13-3021 
[t,S1a,E1a,I1a,A1a,H1a,V1a,V2a,EBa,IBa,ABa,HBa,Ra,Wa,Wva,V3a,Da, ...
    I1cumulaa,A1cumulaa,IBcumulaa,ABcumulaa,D1a,DBa]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);

rho1b = @(x) 2*rho1(x);
rho2b = @(x) 2*rho2(x);
rho3b = @(x) 2*rho3(x);
rho4 = @(x) (1e-03).*( x>=242 );    % Booster doses began on Aug-13-3021
[t,S1b,E1b,I1b,A1b,H1b,V1b,V2b,EBb,IBb,ABb,HBb,Rb,Wb,Wvb,V3b,Db, ...
    I1cumulab,A1cumulab,IBcumulab,ABcumulab,D1b,DBb]= ...
   USA_covid_vaccine_solver(rho1b, rho2b, rho3b, rho4, t, IniVal);

% Behaviour of vaccination 
tStart=datetime('2022-01-04'); tEnd=datetime('2022-08-30');
figure
subplot(1,3,1)
plot(newdata_long,V1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,V1a,'-.b','LineWidth',2)
hold on
plot(newdata_long,V1b,'-.k','LineWidth',2)
hold on
 title('One dose vaccination')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
text(0.5,0.5,'\bfB','FontSize',12)
text(0.5,0.5,'\bfC','FontSize',12)
legend({'50%','Baseline', '200%'},'Location','northeast')
xlim([tStart tEnd])
grid on

subplot(1,3,2)
plot(newdata_long,V2o,'-.r','LineWidth',2)
hold on
plot(newdata_long,V2a,'-.b','LineWidth',2)
hold on
plot(newdata_long,V2b,'-.k','LineWidth',2)
hold on
 title('Two dose vaccination')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'50%','Baseline', '200%'},'Location','northeast')
xlim([tStart tEnd])
grid on

subplot(1,3,3)
plot(newdata_long,V3o,'-.r','LineWidth',2)
hold on
plot(newdata_long,V3a,'-.b','LineWidth',2)
hold on
plot(newdata_long,V3b,'-.k','LineWidth',2)
hold on
 title('Booster dose vaccination')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
legend({'50%','Baseline', '200%'},'Location','northeast')
xlim([tStart tEnd])
grid on


%Unvaccinated
tStart=datetime('2022-01-04'); tEnd=datetime('2022-06-30');
figure
subplot(2,2,1)
plot(newdata_long,I1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,I1a,'-.b','LineWidth',2)
hold on
plot(newdata_long,I1b,'-.k','LineWidth',2)
hold on
 title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'50%','Baseline', '200%'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,A1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,A1a,'-.b','LineWidth',2)
hold on
plot(newdata_long,A1b,'-.k','LineWidth',2)
hold on
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'50%','Baseline', '200%'},'Location','northeast')
xlim([tStart tEnd])
grid on

subplot(2,2,3)
plot(newdata_long,H1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,H1a,'-.b','LineWidth',2)
hold on
plot(newdata_long,H1b,'-.k','LineWidth',2)
hold on
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'50%','Baseline', '200%'},'Location','northeast')
xlim([tStart tEnd])
grid on

subplot(2,2,4)
plot(newdata_long,D1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,D1a,'-.b','LineWidth',2)
hold on
plot(newdata_long,D1b,'-.k','LineWidth',2)
hold on
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'50%','Baseline', '200%'},'Location','northeast')
xlim([tStart tEnd])
grid on

%Vaccinated
tStart=datetime('2022-01-04'); tEnd=datetime('2022-06-30');
figure
subplot(2,2,1)
plot(newdata_long,IBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,IBa,'-.b','LineWidth',2)
hold on
plot(newdata_long,IBb,'-.k','LineWidth',2)
hold on
 title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'200%','Baseline', '50%'},'Location','northeast')
xlim([tStart tEnd])
grid on

subplot(2,2,2)
plot(newdata_long,ABo,'-.r','LineWidth',2)
hold on
plot(newdata_long,ABa,'-.b','LineWidth',2)
hold on
plot(newdata_long,ABb,'-.k','LineWidth',2)
hold on
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'200%','Baseline', '50%'},'Location','northeast')
xlim([tStart tEnd])
grid on

subplot(2,2,3)
plot(newdata_long,HBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,HBa,'-.b','LineWidth',2)
hold on
plot(newdata_long,HBb,'-.k','LineWidth',2)
hold on
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'200%','Baseline', '50%'},'Location','northeast')
xlim([tStart tEnd])
grid on

subplot(2,2,4)
plot(newdata_long,DBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,DBa,'-.b','LineWidth',2)
hold on
plot(newdata_long,DBb,'-.k','LineWidth',2)
hold on
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'200%','Baseline', '50%'},'Location','northeast')
xlim([tStart tEnd])
grid on


%% Lets try by enhancing booster vaccine
rho1 = @(x) (2.1268e+04)*x.*exp(-0.015479*x);
rho2 = @(x) (1.8017e+04)*x.*exp(-0.010932*x);
rho3 = @(x) (8.6062e+03)*(x-81).*exp(-0.023304*(x-81)).*( x>=81 );
rho4a = @(x) (1e-03).*( (x>=242)&(x<387) )...
    + (5.0000e-04).*(x>=387);% Booster doses began on Aug-13-3021
[t,S1o,E1o,I1o,A1o,H1o,V1o,V2o,EBo,IBo,ABo,HBo,Ro,Wo,Wvo,V3o,Do, ...
    I1cumulao,A1cumulao,IBcumulao,ABcumulao,D1o,DBo]= ...
   USA_covid_vaccine_solver(rho1, rho2, rho3, rho4a, t, IniVal);

rho1 = @(x) (2.1268e+04)*x.*exp(-0.015479*x);
rho2 = @(x) (1.8017e+04)*x.*exp(-0.010932*x);
rho3 = @(x) (8.6062e+03)*(x-81).*exp(-0.023304*(x-81)).*( x>=81 );
rho4 = @(x) (1e-03).*( x>=242 );    % Booster doses began on Aug-13-3021 
[t,S1a,E1a,I1a,A1a,H1a,V1a,V2a,EBa,IBa,ABa,HBa,Ra,Wa,Wva,V3a,Da, ...
    I1cumulaa,A1cumulaa,IBcumulaa,ABcumulaa,D1a,DBa]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);

rho1 = @(x) (2.1268e+04)*x.*exp(-0.015479*x);
rho2 = @(x) (1.8017e+04)*x.*exp(-0.010932*x);
rho3 = @(x) (8.6062e+03)*(x-81).*exp(-0.023304*(x-81)).*( x>=81 );
rho4b = @(x) (1e-03).*( (x>=242)&(x<387) )...
    + (0.0020).*(x>=387);% Booster doses began on Aug-13-3021
[t,S1b,E1b,I1b,A1b,H1b,V1b,V2b,EBb,IBb,ABb,HBb,Rb,Wb,Wvb,V3b,Db, ...
    I1cumulab,A1cumulab,IBcumulab,ABcumulab,D1b,DBb]= ...
   USA_covid_vaccine_solver(rho1, rho2, rho3, rho4b, t, IniVal);

%Behaviour of vaccination rates 
figure
subplot(1,3,1)
plot(newdata_long,V1a,'-.b','LineWidth',2)
title('One dose vaccination')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'Baseline'},'Location','northeast')
xlim([tStart tEnd])
grid on

subplot(1,3,2)
plot(newdata_long,V2a,'-.b','LineWidth',2)
title('Two dose vaccination')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'Baseline'},'Location','northeast')
xlim([tStart tEnd])
grid on

subplot(1,3,3)
plot(newdata_long,V3o,'-.r','LineWidth',2)
hold on
plot(newdata_long,V3a,'-.b','LineWidth',2)
hold on
plot(newdata_long,V3b,'-.k','LineWidth',2)
hold on
 title('Booster dose vaccination')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
legend({'50%','Baseline', '200%'},'Location','northeast')
xlim([tStart tEnd])
grid on

%Vaccinated
figure
subplot(2,2,1)
plot(newdata_long,IBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,IBa,'-.b','LineWidth',2)
hold on
plot(newdata_long,IBb,'-.k','LineWidth',2)
hold on
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'200%','Baseline', '50%'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,ABo,'-.r','LineWidth',2)
hold on
plot(newdata_long,ABa,'-.b','LineWidth',2)
hold on
plot(newdata_long,ABb,'-.k','LineWidth',2)
hold on
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'200%','Baseline', '50%'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,HBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,HBa,'-.b','LineWidth',2)
hold on
plot(newdata_long,HBb,'-.k','LineWidth',2)
hold on
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'200%','Baseline', '50%'},'Location','northeast')
grid on

subplot(2,2,4)
plot(newdata_long,DBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,DBa,'-.b','LineWidth',2)
hold on
plot(newdata_long,DBb,'-.k','LineWidth',2)
hold on
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'200%','Baseline', '50%'},'Location','northeast')
grid on

%% Varying vaccine efficiency (Baseline Vaccination)
%Baseline Vaccination
rho1 = @(x) (2.1268e+04)*x.*exp(-0.015479*x);
rho2 = @(x) (1.8017e+04)*x.*exp(-0.010932*x);
rho3 = @(x) (8.6062e+03)*(x-81).*exp(-0.023304*(x-81)).*( x>=81 );
rho4 = @(x) (1e-03).*( x>=242 );    % Booster doses began on Aug-13-3021 
EpsL=0.;
EpsL2=0.1;
eps3i=0.55;
% Baseline Efficiency 
[t,S1a,E1a,I1a,A1a,H1a,V1a,V2a,EBa,IBa,ABa,HBa,Ra,Wa,Wva,V3a,Da, ...
    I1cumulaa,A1cumulaa,IBcumulaa,ABcumulaa,D1a,DBa]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
%Low Efficiency 
EpsL=0.26;
EpsL2=0.19;
eps3i=0.63;
[t,S1o,E1o,I1o,A1o,H1o,V1o,V2o,EBo,IBo,ABo,HBo,Ro,Wo,Wvo,V3o,Do, ...
    I1cumulao,A1cumulao,IBcumulao,ABcumulao,D1o,DBo]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
%High Efficiency 
EpsL=0.3;
EpsL2=0.35;
eps3i=0.75;
[t,S1b,E1b,I1b,A1b,H1b,V1b,V2b,EBb,IBb,ABb,HBb,Rb,Wb,Wvb,V3b,Db, ...
    I1cumulab,A1cumulab,IBcumulab,ABcumulab,D1b,DBb]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);

figure
subplot(2,2,1)
plot(newdata_long,I1o,'-.r','LineWidth',1)
hold on
plot(newdata_long,I1a,'-.b','LineWidth',1)
hold on
plot(newdata_long,I1b,'-.k','LineWidth',2)
hold on
 title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'\epsilon_L=0 , \epsilon_L2=0.1,\epsilon_3=0.55',...
    '\epsilon_L=0.26 , \epsilon_L2=0.19,\epsilon_3=0.63',...
    '\epsilon_L=0.3 , \epsilon_L2=0.35 ,\epsilon_3=0.75'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,A1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,A1a,'-.b','LineWidth',2)
hold on
plot(newdata_long,A1b,'-.k','LineWidth',2)
hold on
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'\epsilon_L=0 , \epsilon_L2=0.1,\epsilon_3=0.55',...
    '\epsilon_L=0.26 , \epsilon_L2=0.19,\epsilon_3=0.63',...
    '\epsilon_L=0.3 , \epsilon_L2=0.35 ,\epsilon_3=0.75'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,H1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,H1a,'-.b','LineWidth',2)
hold on
plot(newdata_long,H1b,'-.k','LineWidth',2)
hold on
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'\epsilon_L=0 , \epsilon_L2=0.1,\epsilon_3=0.55',...
    '\epsilon_L=0.26 , \epsilon_L2=0.19,\epsilon_3=0.63',...
    '\epsilon_L=0.3 , \epsilon_L2=0.35 ,\epsilon_3=0.75'},'Location','northeast')
grid on

subplot(2,2,4)
plot(newdata_long,D1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,D1a,'-.b','LineWidth',2)
hold on
plot(newdata_long,D1b,'-.k','LineWidth',2)
hold on
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'\epsilon_L=0 , \epsilon_L2=0.1,\epsilon_3=0.55',...
    '\epsilon_L=0.26 , \epsilon_L2=0.19,\epsilon_3=0.63',...
    '\epsilon_L=0.3 , \epsilon_L2=0.35 ,\epsilon_3=0.75'},'Location','northeast')
grid on

%Plot Vaccinated Individuals 
figure
subplot(2,2,1)
plot(newdata_long,IBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,IBa,'-.b','LineWidth',2)
hold on
plot(newdata_long,IBb,'-.k','LineWidth',2)
hold on
 title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'\epsilon_L=0 , \epsilon_L2=0.1,\epsilon_3=0.55',...
    '\epsilon_L=0.26 , \epsilon_L2=0.19,\epsilon_3=0.63',...
    '\epsilon_L=0.3 , \epsilon_L2=0.35 ,\epsilon_3=0.75'},'Location','northeast')
xlim([tStart tEnd])
grid on

subplot(2,2,2)
plot(newdata_long,ABo,'-.r','LineWidth',2)
hold on
plot(newdata_long,ABa,'-.b','LineWidth',2)
hold on
plot(newdata_long,ABb,'-.k','LineWidth',2)
hold on
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'\epsilon_L=0 , \epsilon_L2=0.1,\epsilon_3=0.55',...
    '\epsilon_L=0.26 , \epsilon_L2=0.19,\epsilon_3=0.63',...
    '\epsilon_L=0.3 , \epsilon_L2=0.35 ,\epsilon_3=0.75'},'Location','northeast')
xlim([tStart tEnd])
grid on

subplot(2,2,3)
plot(newdata_long,HBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,HBa,'-.b','LineWidth',2)
hold on
plot(newdata_long,HBb,'-.k','LineWidth',2)
hold on
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'\epsilon_L=0 , \epsilon_L2=0.1,\epsilon_3=0.55',...
    '\epsilon_L=0.26 , \epsilon_L2=0.19,\epsilon_3=0.63',...
    '\epsilon_L=0.3 , \epsilon_L2=0.35 ,\epsilon_3=0.75'},'Location','northeast')
xlim([tStart tEnd])
grid on

subplot(2,2,4)
plot(newdata_long,DBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,DBa,'-.b','LineWidth',2)
hold on
plot(newdata_long,DBb,'-.k','LineWidth',2)
hold on
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'\epsilon_L=0 , \epsilon_L2=0.1,\epsilon_3=0.55',...
    '\epsilon_L=0.26 , \epsilon_L2=0.19,\epsilon_3=0.63',...
    '\epsilon_L=0.3 , \epsilon_L2=0.35 ,\epsilon_3=0.75'},'Location','northeast')
xlim([tStart tEnd])
grid on

%% Varying vaccine efficiency (50% Vaccination)
%Baseline Vaccination
rho1a = @(x) 0.5*rho1(x);
rho2a = @(x) 0.5*rho2(x);
rho3a = @(x) 0.5*rho3(x);
rho4 = @(x) (1e-03).*( x>=242 );   % Booster doses began on Aug-13-3021 
%Low Efficiency 
EpsL=0.;
EpsL2=0.1;
eps3i=0.55;
% Baseline Efficiency 
[t,S1a,E1a,I1a,A1a,H1a,V1a,V2a,EBa,IBa,ABa,HBa,Ra,Wa,Wva,V3a,Da, ...
    I1cumulaa,A1cumulaa,IBcumulaa,ABcumulaa,D1a,DBa]= ...
    USA_covid_vaccine_solver(rho1a,rho2a,rho3a,rho4, t, IniVal);
%Low Efficiency 
EpsL=0.26;
EpsL2=0.19;
eps3i=0.63;
[t,S1o,E1o,I1o,A1o,H1o,V1o,V2o,EBo,IBo,ABo,HBo,Ro,Wo,Wvo,V3o,Do, ...
    I1cumulao,A1cumulao,IBcumulao,ABcumulao,D1o,DBo]= ...
    USA_covid_vaccine_solver(rho1a,rho2a,rho3a,rho4, t, IniVal);
%High Efficiency 
EpsL=03.;
EpsL2=0.35;
eps3i=0.75;
[t,S1b,E1b,I1b,A1b,H1b,V1b,V2b,EBb,IBb,ABb,HBb,Rb,Wb,Wvb,V3b,Db, ...
    I1cumulab,A1cumulab,IBcumulab,ABcumulab,D1b,DBb]= ...
    USA_covid_vaccine_solver(rho1a,rho2a,rho3a,rho4, t, IniVal);

figure
subplot(2,2,1)
plot(newdata_long,I1o,'-.r','LineWidth',1)
hold on
plot(newdata_long,I1a,'-.b','LineWidth',1)
hold on
plot(newdata_long,I1b,'-.k','LineWidth',2)
hold on
 title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'\epsilon_L=0 , \epsilon_L2=0.1,\epsilon_3=0.55',...
    '\epsilon_L=0.26 , \epsilon_L2=0.19,\epsilon_3=0.63',...
    '\epsilon_L=0.3 , \epsilon_L2=0.35 ,\epsilon_3=0.75'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,A1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,A1a,'-.b','LineWidth',2)
hold on
plot(newdata_long,A1b,'-.k','LineWidth',2)
hold on
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'\epsilon_L=0 , \epsilon_L2=0.1,\epsilon_3=0.55',...
    '\epsilon_L=0.26 , \epsilon_L2=0.19,\epsilon_3=0.63',...
    '\epsilon_L=0.3 , \epsilon_L2=0.35 ,\epsilon_3=0.75'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,H1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,H1a,'-.b','LineWidth',2)
hold on
plot(newdata_long,H1b,'-.k','LineWidth',2)
hold on
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'\epsilon_L=0 , \epsilon_L2=0.1,\epsilon_3=0.55',...
    '\epsilon_L=0.26 , \epsilon_L2=0.19,\epsilon_3=0.63',...
    '\epsilon_L=0.3 , \epsilon_L2=0.35 ,\epsilon_3=0.75'},'Location','northeast')
grid on

subplot(2,2,4)
plot(newdata_long,D1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,D1a,'-.b','LineWidth',2)
hold on
plot(newdata_long,D1b,'-.k','LineWidth',2)
hold on
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'\epsilon_L=0 , \epsilon_L2=0.1,\epsilon_3=0.55',...
    '\epsilon_L=0.26 , \epsilon_L2=0.19,\epsilon_3=0.63',...
    '\epsilon_L=0.3 , \epsilon_L2=0.35 ,\epsilon_3=0.75'},'Location','northeast')
grid on

%Plot Vaccinated Individuals 
figure
subplot(2,2,1)
plot(newdata_long,IBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,IBa,'-.b','LineWidth',1)
hold on
plot(newdata_long,IBb,'-.k','LineWidth',2)
hold on
 title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'\epsilon_L=0 , \epsilon_L2=0.1,\epsilon_3=0.55',...
    '\epsilon_L=0.26 , \epsilon_L2=0.19,\epsilon_3=0.63',...
    '\epsilon_L=0.3 , \epsilon_L2=0.35 ,\epsilon_3=0.75'},'Location','northeast')
xlim([tStart tEnd])
grid on

subplot(2,2,2)
plot(newdata_long,ABo,'-.r','LineWidth',2)
hold on
plot(newdata_long,ABa,'-.b','LineWidth',1)
hold on
plot(newdata_long,ABb,'-.k','LineWidth',2)
hold on
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'\epsilon_L=0 , \epsilon_L2=0.1,\epsilon_3=0.55',...
    '\epsilon_L=0.26 , \epsilon_L2=0.19,\epsilon_3=0.63',...
    '\epsilon_L=0.3 , \epsilon_L2=0.35 ,\epsilon_3=0.75'},'Location','northeast')
xlim([tStart tEnd])
grid on

subplot(2,2,3)
plot(newdata_long,HBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,HBa,'-.b','LineWidth',1)
hold on
plot(newdata_long,HBb,'-.k','LineWidth',2)
hold on
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'\epsilon_L=0 , \epsilon_L2=0.1,\epsilon_3=0.55',...
    '\epsilon_L=0.26 , \epsilon_L2=0.19,\epsilon_3=0.63',...
    '\epsilon_L=0.3 , \epsilon_L2=0.35 ,\epsilon_3=0.75'},'Location','northeast')
xlim([tStart tEnd])
grid on

subplot(2,2,4)
plot(newdata_long,DBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,DBa,'-.b','LineWidth',1)
hold on
plot(newdata_long,DBb,'-.k','LineWidth',2)
hold on
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'\epsilon_L=-0.05 , \epsilon_L2=-0.01,\epsilon_3=0.68',...
    '\epsilon_L=0 , \epsilon_L2=0.07,\epsilon_3=0.66',...
    '\epsilon_L=0.05 , \epsilon_L2=0.13 ,\epsilon_3=0.64'},'Location','northeast')
xlim([tStart tEnd])
grid on

%% Varying vaccine efficiency (200% Vaccination)
%Baseline Vaccination
rho1b = @(x) 2*rho1(x);
rho2b = @(x) 2*rho2(x);
rho3b = @(x) 2*rho3(x);
rho4 = @(x) (1e-03).*( x>=242 );   % Booster doses began on Aug-13-3021 
EpsL=0.;
EpsL2=0.1;
eps3i=0.55;
% Baseline Efficiency 
[t,S1a,E1a,I1a,A1a,H1a,V1a,V2a,EBa,IBa,ABa,HBa,Ra,Wa,Wva,V3a,Da, ...
    I1cumulaa,A1cumulaa,IBcumulaa,ABcumulaa,D1a,DBa]= ...
    USA_covid_vaccine_solver(rho1b,rho2b,rho3b,rho4, t, IniVal);
%Low Efficiency 
EpsL=0.26;
EpsL2=0.19;
eps3i=0.63;
[t,S1o,E1o,I1o,A1o,H1o,V1o,V2o,EBo,IBo,ABo,HBo,Ro,Wo,Wvo,V3o,Do, ...
    I1cumulao,A1cumulao,IBcumulao,ABcumulao,D1o,DBo]= ...
    USA_covid_vaccine_solver(rho1b,rho2b,rho3b,rho4, t, IniVal);
%High Efficiency 
EpsL=03.;
EpsL2=0.35;
eps3i=0.75;
[t,S1b,E1b,I1b,A1b,H1b,V1b,V2b,EBb,IBb,ABb,HBb,Rb,Wb,Wvb,V3b,Db, ...
    I1cumulab,A1cumulab,IBcumulab,ABcumulab,D1b,DBb]= ...
    USA_covid_vaccine_solver(rho1b,rho2b,rho3b,rho4, t, IniVal);
%Plot Vaccinated Individuals 
figure
subplot(2,2,1)
plot(newdata_long,IBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,IBa,'-.b','LineWidth',1)
hold on
plot(newdata_long,IBb,'-.k','LineWidth',2)
hold on
 title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'\epsilon_L=0.3 , \epsilon_L2=0.35 ,\epsilon_3=0.75',...
    '\epsilon_L=0.26 , \epsilon_L2=0.19,\epsilon_3=0.63',...
    '\epsilon_L=0 , \epsilon_L2=0.1,\epsilon_3=0.55'},'Location','northeast')
xlim([tStart tEnd])
grid on

subplot(2,2,2)
plot(newdata_long,ABo,'-.r','LineWidth',2)
hold on
plot(newdata_long,ABa,'-.b','LineWidth',1)
hold on
plot(newdata_long,ABb,'-.k','LineWidth',2)
hold on
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'\epsilon_L=0.3 , \epsilon_L2=0.35 ,\epsilon_3=0.75',...
    '\epsilon_L=0.26 , \epsilon_L2=0.19,\epsilon_3=0.63',...
    '\epsilon_L=0 , \epsilon_L2=0.1,\epsilon_3=0.55'},'Location','northeast')
xlim([tStart tEnd])
grid on

subplot(2,2,3)
plot(newdata_long,HBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,HBa,'-.b','LineWidth',1)
hold on
plot(newdata_long,HBb,'-.k','LineWidth',2)
hold on
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'\epsilon_L=0.3 , \epsilon_L2=0.35 ,\epsilon_3=0.75',...
    '\epsilon_L=0.26 , \epsilon_L2=0.19,\epsilon_3=0.63',...
    '\epsilon_L=0 , \epsilon_L2=0.1,\epsilon_3=0.55'},'Location','northeast')
xlim([tStart tEnd])
grid on

subplot(2,2,4)
plot(newdata_long,DBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,DBa,'-.b','LineWidth',1)
hold on
plot(newdata_long,DBb,'-.k','LineWidth',2)
hold on
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'\epsilon_L=0.3 , \epsilon_L2=0.35 ,\epsilon_3=0.75',...
    '\epsilon_L=0.26 , \epsilon_L2=0.19,\epsilon_3=0.63',...
    '\epsilon_L=0 , \epsilon_L2=0.1,\epsilon_3=0.55'},'Location','northeast')
xlim([tStart tEnd])
grid on

%% Varying Booster Vaccine (Baseline Booster)
rho1 = @(x) (2.1268e+04)*x.*exp(-0.015479*x);
rho2 = @(x) (1.8017e+04)*x.*exp(-0.010932*x);
rho3 = @(x) (8.6062e+03)*(x-81).*exp(-0.023304*(x-81)).*( x>=81 );
rho4 = @(x) (1e-03).*( x>=242 );    % Booster doses began on Aug-13-3021 

EpsL=0.;
EpsL2=0.1;
eps3i=0.55;
% Baseline Efficiency 
[t,S1a,E1a,I1a,A1a,H1a,V1a,V2a,EBa,IBa,ABa,HBa,Ra,Wa,Wva,V3a,Da, ...
    I1cumulaa,A1cumulaa,IBcumulaa,ABcumulaa,D1a,DBa]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
%Low Efficiency 
EpsL=0.26;
EpsL2=0.19;
eps3i=0.63;
[t,S1o,E1o,I1o,A1o,H1o,V1o,V2o,EBo,IBo,ABo,HBo,Ro,Wo,Wvo,V3o,Do, ...
    I1cumulao,A1cumulao,IBcumulao,ABcumulao,D1o,DBo]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
%High Efficiency 
EpsL=03.;
EpsL2=0.35;
eps3i=0.75;
[t,S1b,E1b,I1b,A1b,H1b,V1b,V2b,EBb,IBb,ABb,HBb,Rb,Wb,Wvb,V3b,Db, ...
    I1cumulab,A1cumulab,IBcumulab,ABcumulab,D1b,DBb]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
%Plot Vaccinated Individuals 
figure
subplot(2,2,1)
plot(newdata_long,IBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,IBa,'-.b','LineWidth',1)
hold on
plot(newdata_long,IBb,'-.k','LineWidth',2)
hold on
 title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'\epsilon_L=0.3 , \epsilon_L2=0.35 ,\epsilon_3=0.75',...
    '\epsilon_L=0.26 , \epsilon_L2=0.19,\epsilon_3=0.63',...
    '\epsilon_L=0 , \epsilon_L2=0.1,\epsilon_3=0.55'},'Location','northeast')
xlim([tStart tEnd])
grid on

subplot(2,2,2)
plot(newdata_long,ABo,'-.r','LineWidth',2)
hold on
plot(newdata_long,ABa,'-.b','LineWidth',1)
hold on
plot(newdata_long,ABb,'-.k','LineWidth',2)
hold on
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'\epsilon_L=0.3 , \epsilon_L2=0.35 ,\epsilon_3=0.75',...
    '\epsilon_L=0.26 , \epsilon_L2=0.19,\epsilon_3=0.63',...
    '\epsilon_L=0 , \epsilon_L2=0.1,\epsilon_3=0.55'},'Location','northeast')
xlim([tStart tEnd])
grid on

subplot(2,2,3)
plot(newdata_long,HBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,HBa,'-.b','LineWidth',1)
hold on
plot(newdata_long,HBb,'-.k','LineWidth',2)
hold on
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'\epsilon_L=0.3 , \epsilon_L2=0.35 ,\epsilon_3=0.75',...
    '\epsilon_L=0.26 , \epsilon_L2=0.19,\epsilon_3=0.63',...
    '\epsilon_L=0 , \epsilon_L2=0.1,\epsilon_3=0.55'},'Location','northeast')
xlim([tStart tEnd])
grid on

subplot(2,2,4)
plot(newdata_long,DBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,DBa,'-.b','LineWidth',1)
hold on
plot(newdata_long,DBb,'-.k','LineWidth',2)
hold on
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'\epsilon_L=0.3 , \epsilon_L2=0.35 ,\epsilon_3=0.75',...
    '\epsilon_L=0.26 , \epsilon_L2=0.19,\epsilon_3=0.63',...
    '\epsilon_L=0 , \epsilon_L2=0.1,\epsilon_3=0.55'},'Location','northeast')
xlim([tStart tEnd])
grid on


%% Varying Booster Vaccine (50% Baseline Booster)
rho1 = @(x) (2.1268e+04)*x.*exp(-0.015479*x);
rho2 = @(x) (1.8017e+04)*x.*exp(-0.010932*x);
rho3 = @(x) (8.6062e+03)*(x-81).*exp(-0.023304*(x-81)).*( x>=81 );
rho4a = @(x) (1e-03).*( (x>=242)&(x<388) )...
      + (5.0000e-04).*(x>=388);% Booster doses began on Aug-13-3021


% Baseline Efficiency 
EpsL=0.26;
EpsL2=0.19;
eps3i=0.63;
[t,S1a,E1a,I1a,A1a,H1a,V1a,V2a,EBa,IBa,ABa,HBa,Ra,Wa,Wva,V3a,Da, ...
    I1cumulaa,A1cumulaa,IBcumulaa,ABcumulaa,D1a,DBa]= ...
    USA_covid_vaccine_solver(rho1b,rho2b,rho3b,rho4a, t, IniVal);
%Low Efficiency 
EpsL=0.;
EpsL2=0.1;
eps3i=0.55;
[t,S1o,E1o,I1o,A1o,H1o,V1o,V2o,EBo,IBo,ABo,HBo,Ro,Wo,Wvo,V3o,Do, ...
    I1cumulao,A1cumulao,IBcumulao,ABcumulao,D1o,DBo]= ...
    USA_covid_vaccine_solver(rho1b,rho2b,rho3b,rho4a, t, IniVal);
%High Efficiency 
EpsL=0.3;
EpsL2=0.35;
eps3i=0.75;
[t,S1b,E1b,I1b,A1b,H1b,V1b,V2b,EBb,IBb,ABb,HBb,Rb,Wb,Wvb,V3b,Db, ...
    I1cumulab,A1cumulab,IBcumulab,ABcumulab,D1b,DBb]= ...
    USA_covid_vaccine_solver(rho1b,rho2b,rho3b,rho4a, t, IniVal);
%Plot Vaccinated Individuals 
figure
subplot(2,2,1)
plot(newdata_long,IBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,IBa,'-.b','LineWidth',2)
hold on
plot(newdata_long,IBb,'-.k','LineWidth',2)
hold on
 title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'\epsilon_L=0.3 , \epsilon_L2=0.35 ,\epsilon_3=0.75',...
    '\epsilon_L=0.26 , \epsilon_L2=0.19,\epsilon_3=0.63',...
    '\epsilon_L=0 , \epsilon_L2=0.1,\epsilon_3=0.55'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,ABo,'-.r','LineWidth',2)
hold on
plot(newdata_long,ABa,'-.b','LineWidth',2)
hold on
plot(newdata_long,ABb,'-.k','LineWidth',2)
hold on
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'\epsilon_L=0.3 , \epsilon_L2=0.35 ,\epsilon_3=0.75',...
    '\epsilon_L=0.26 , \epsilon_L2=0.19,\epsilon_3=0.63',...
    '\epsilon_L=0 , \epsilon_L2=0.1,\epsilon_3=0.55'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,HBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,HBa,'-.b','LineWidth',2)
hold on
plot(newdata_long,HBb,'-.k','LineWidth',2)
hold on
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'\epsilon_L=0.3 , \epsilon_L2=0.35 ,\epsilon_3=0.75',...
    '\epsilon_L=0.26 , \epsilon_L2=0.19,\epsilon_3=0.63',...
    '\epsilon_L=0 , \epsilon_L2=0.1,\epsilon_3=0.55'},'Location','northeast')
grid on

subplot(2,2,4)
plot(newdata_long,DBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,DBa,'-.b','LineWidth',2)
hold on
plot(newdata_long,DBb,'-.k','LineWidth',2)
hold on
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'\epsilon_L=0.3 , \epsilon_L2=0.35 ,\epsilon_3=0.75',...
    '\epsilon_L=0.26 , \epsilon_L2=0.19,\epsilon_3=0.63',...
    '\epsilon_L=0 , \epsilon_L2=0.1,\epsilon_3=0.55'},'Location','northeast')
grid on

%% Varying Booster Vaccine (200% Baseline Booster)
rho1 = @(x) (2.1268e+04)*x.*exp(-0.015479*x);
rho2 = @(x) (1.8017e+04)*x.*exp(-0.010932*x);
rho3 = @(x) (8.6062e+03)*(x-81).*exp(-0.023304*(x-81)).*( x>=81 );
rho4b = @(x) (1e-03).*( (x>=242)&(x<388) )...
      + (0.0020).*(x>=388);% Booster doses began on Aug-13-3021


% Baseline Efficiency 
EpsL=0.26;
EpsL2=0.19;
eps3i=0.63;
[t,S1a,E1a,I1a,A1a,H1a,V1a,V2a,EBa,IBa,ABa,HBa,Ra,Wa,Wva,V3a,Da, ...
    I1cumulaa,A1cumulaa,IBcumulaa,ABcumulaa,D1a,DBa]= ...
    USA_covid_vaccine_solver(rho1b,rho2b,rho3b,rho4b, t, IniVal);
%Low Efficiency 
EpsL=0.;
EpsL2=0.1;
eps3i=0.55;
[t,S1o,E1o,I1o,A1o,H1o,V1o,V2o,EBo,IBo,ABo,HBo,Ro,Wo,Wvo,V3o,Do, ...
    I1cumulao,A1cumulao,IBcumulao,ABcumulao,D1o,DBo]= ...
    USA_covid_vaccine_solver(rho1b,rho2b,rho3b,rho4b, t, IniVal);
%High Efficiency 
EpsL=0.3;
EpsL2=0.35;
eps3i=0.75;
[t,S1b,E1b,I1b,A1b,H1b,V1b,V2b,EBb,IBb,ABb,HBb,Rb,Wb,Wvb,V3b,Db, ...
    I1cumulab,A1cumulab,IBcumulab,ABcumulab,D1b,DBb]= ...
    USA_covid_vaccine_solver(rho1b,rho2b,rho3b,rho4b, t, IniVal);
%Plot Vaccinated Individuals 
figure
subplot(2,2,1)
plot(newdata_long,IBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,IBa,'-.b','LineWidth',2)
hold on
plot(newdata_long,IBb,'-.k','LineWidth',2)
hold on
 title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'\epsilon_L=0.3 , \epsilon_L2=0.35 ,\epsilon_3=0.75',...
    '\epsilon_L=0.26 , \epsilon_L2=0.19,\epsilon_3=0.63',...
    '\epsilon_L=0 , \epsilon_L2=0.1,\epsilon_3=0.55'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,ABo,'-.r','LineWidth',2)
hold on
plot(newdata_long,ABa,'-.b','LineWidth',2)
hold on
plot(newdata_long,ABb,'-.k','LineWidth',2)
hold on
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'\epsilon_L=0.3 , \epsilon_L2=0.35 ,\epsilon_3=0.75',...
    '\epsilon_L=0.26 , \epsilon_L2=0.19,\epsilon_3=0.63',...
    '\epsilon_L=0 , \epsilon_L2=0.1,\epsilon_3=0.55'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,HBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,HBa,'-.b','LineWidth',2)
hold on
plot(newdata_long,HBb,'-.k','LineWidth',2)
hold on
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'\epsilon_L=0.3 , \epsilon_L2=0.35 ,\epsilon_3=0.75',...
    '\epsilon_L=0.26 , \epsilon_L2=0.19,\epsilon_3=0.63',...
    '\epsilon_L=0 , \epsilon_L2=0.1,\epsilon_3=0.55'},'Location','northeast')
grid on

subplot(2,2,4)
plot(newdata_long,DBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,DBa,'-.b','LineWidth',2)
hold on
plot(newdata_long,DBb,'-.k','LineWidth',2)
hold on
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'\epsilon_L=0.3 , \epsilon_L2=0.35 ,\epsilon_3=0.75',...
    '\epsilon_L=0.26 , \epsilon_L2=0.19,\epsilon_3=0.63',...
    '\epsilon_L=0 , \epsilon_L2=0.1,\epsilon_3=0.55'},'Location','northeast')
grid on

%%
load('BestFitParameters')

%% Setting initial values for Dec-13-2020:
IniVal = [S1(57); E1(57); I1(57); A1(57); H1_exp(57); 0; 0; 0; 0; 0; 0; ...
    R(57); W(57); 0; 0; D_exp(57); Icumul_exp(57); Acumul_exp(57); 0; 0; D_exp(57); 0];

dateslong = dates(327:714);  % Dec-13-2020 to Jan-04-2022
tsim=linspace(0,length(dateslong)-1,length(dateslong));
Icumul_exp_vacc=C_exp_US(57+initdate-1:714);
D_exp_vacc=D_exp_US(57+initdate-1:714);
H1_exp_vacc=Hosp_exp_US(57+initdate-diffDates-1:end);

global alpha1 chi omega beta1 beta2 Cm epsilonm kappa epsa1 epsa2 epsa3 ...
    SumEpsLi SumEpsL2i Sumthetai Kv Kn Kv1 Kv2 eps3i p2 gamma2 alphaB chiB omegaB 

%% Parameter settings

% Parameter values obtained from data fitting:
beta1orig = x(1)/( 1 - 0.5*0.5); % Transmission rate for Wuhan variant (symptomatic)
beta2orig = x(2)/( 1 - 0.5*0.5); % Transmission rate for Wuhan variant (asymptomatic)
alpha1 = x(3);       % Hospitalization rate
chi = x(4);          % Death probability for H
omega = x(5);        % 1/(average hospitalization period)

% Time-varying transmission rates:
beta1 = @(X) beta1orig*(X<=125) ...                 % Second wave
          + 2.8*beta1orig.*( (X>125)&(X<=364) ) ... % Third wave
          + 5.0*beta1orig.*( (X>364) ) ;            % Fourth wave (Omicron)

beta2 = @(X) beta2orig*(X<=125) ...                 % Second wave
          + 2.8*beta2orig.*( (X>125)&(X<=364) ) ... % Third wave
          + 5.0*beta2orig.*( (X>364) ) ;            % Fourth wave (Omicron)

% Other parameter values:

Cm = 0.5;   % Face mask compliance in the community
epsilonm = 0.5;     % Efficacy of face masks
kappa = 0.52;   % Relative infectivity of breakthrough infections
epsa1 = 7.23e-4;    % All-or-nothing (Moderna)
epsa2 = 3.77e-4;    % All-or-nothing (Pfizer)
epsa3 = 0.0091;     % All-or-nothing (Janssen J&J)
Kv1 = 1/60;     % Waning rate of immunity after 1 dose
Kv2 = 1/180;    % Waning rate of immunity after 2 doses
Kv = 1;    % Waning rate of acquired immunity
Kn=1/270;     % Waning rate of natural immunity
eps3i = 1-0.37;   % 1 - (effectiveness of booster dose)
EpsL = 0.26;    % Effectiveness for 1 dose
EpsL2 = 0.19;   % Effectiveness for 2 doses
SumEpsLi = 1 - EpsL;    % Sum of 1-epsilon_L,i (i=1,2,3)
SumEpsL2i = 1 - EpsL2;   % Sum of 1-epsilon_L2,i (i=1,2,3)
Sumthetai = 2/180;      % Sum of theta1 and theta2

chiB = 0.000172;          % Death probability for HB (breakthrough infections)
omegaB = 2; %0.0514;       % 1/(average hospitalization period) (breakthrough infections)
alphaB = 0.0118;        % Hospitalization rate (breakthrough infections)
p2 = 1/9;   % Proportion of symptomatic infections (breakthrough infections)
gamma2 = 1/14;  % Recovery rate (breakthrough infections)

%% Defining vaccination rates
rho1 = @(x) (2.1268e+04)*x.*exp(-0.015479*x);
rho2 = @(x) (1.8017e+04)*x.*exp(-0.010932*x);
rho3 = @(x) (8.6062e+03)*(x-81).*exp(-0.023304*(x-81)).*( x>=81 );
rho4 = @(x) (1e-03).*( x>=242 );    % Booster doses began on Aug-13-3021

%% Solving the system
[tsim,S1a,E1a,I1a,A1a,H1a,V1a,V2a,EBa,IBa,ABa,HBa,Ra,Wa,Wva,V3a,Da, ...
    I1cumula,A1cumula,IBcumula,ABcumula,D1a,DBa]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, tsim, IniVal);

%% Seting initial conditions for simulations

S10 = S1a(end);
E10 = E1a(end);
I10 = I1a(end);
A10 = A1a(end);
H10 = H1a(end);
V10 = V1a(end);
V20 = V2a(end);
EB0 = EBa(end);
IB0 = IBa(end);
AB0 = ABa(end);
HB0 = HBa(end);
R0 = Ra(end);
W0 = Wa(end);
WV0= Wva(end);
V30= V3a(end);
D0 = 832180;
I1cumula0 = I1cumula(end);
A1cumula0 = A1cumula(end);
IBcumula0 = IBcumula(end);
ABcumula0 = ABcumula(end);
D10=D0*D1a(end)/( D1a(end)+DBa(end) );
DB0=D0*DBa(end)/( D1a(end)+DBa(end) );

IniVal = [S10;E10;I10;A10;H10;V10;V20;EB0;IB0;AB0;HB0;R0;...
    W0;WV0;V30;D0;I1cumula0;A1cumula0;IBcumula0;ABcumula0;D10;DB0];

newinitdate=714;   % 4 Jan 2022

numberofdays = 300;
t=linspace(tsim(end), tsim(end)+numberofdays,numberofdays+1);
newdata_long=dateshift(dateslong(388),'start','day',0:numberofdays);

%% Defining vaccination rates
rho1 = @(x) (2.1268e+04)*x.*exp(-0.015479*x);
rho2 = @(x) (1.8017e+04)*x.*exp(-0.010932*x);
rho3 = @(x) (8.6062e+03)*(x-81).*exp(-0.023304*(x-81)).*( x>=81 );
rho4 = @(x) (1e-03).*( x>=242 );    % Booster doses began on Aug-13-3021

%% Solving the system
[t,S10,E10,I10,A10,H10,V10,V20,EB0,IB0,AB0,HB0,R0,W0,Wv0,V30,D0, ...
    I1cumula0,A1cumula0,IBcumula0,ABcumula0,D10,DB0]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);


%% Varying Natural Waning Immunity with face mask usage of 50& and 50& efficiency
Cm=0.5;
epsilonm=0.5;
Kn=0.0048;
[t,S1o,E1o,I1o,A1o,H1o,V1o,V2o,EBo,IBo,ABo,HBo,Ro,Wo,Wvo,V3o,Do, ...
    I1cumulao,A1cumulao,IBcumulao,ABcumulao,D1o,DBo]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.5;
epsilonm=0.5;
Kn=0.0037;
[t,S1a,E1a,I1a,A1a,H1a,V1a,V2a,EBa,IBa,ABa,HBa,Ra,Wa,Wva,V3a,Da, ...
    I1cumulaa,A1cumulaa,IBcumulaa,ABcumulaa,D1a,DBa]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.5;
epsilonm=0.5;
Kn=0.0030;
[t,S1b,E1b,I1b,A1b,H1b,V1b,V2b,EBb,IBb,ABb,HBb,Rb,Wb,Wvb,V3b,Db, ...
    I1cumulab,A1cumulab,IBcumulab,ABcumulab,D1b,DBb]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.5;
epsilonm=0.5;
Kv=0.0;
Kn=0.0;
[t,S1c,E1c,I1c,A1c,H1c,V1c,V2c,EBc,IBc,ABc,HBc,Rc,Wc,Wvc,V3c,Dc, ...
    I1cumulac,A1cumulac,IBcumulac,ABcumulac,D1c,DBc]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
%Plot Unvaccinated Individuals 
figure
subplot(2,2,1)
plot(newdata_long,I1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,I1a,'-.b','LineWidth',2)
hold on
plot(newdata_long,I1b,'-.k','LineWidth',2)
hold on
plot(newdata_long,I1c,'-.m','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'K_v=0.0048','K_v=0.0037','K_v=0.0030','K_v=0'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,A1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,A1a,'-.b','LineWidth',2)
hold on
plot(newdata_long,A1b,'-.k','LineWidth',2)
hold on
plot(newdata_long,A1c,'-.m','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'K_v=0.0048','K_v=0.0037','K_v=0.0030','K_v=0'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,H1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,H1a,'-.b','LineWidth',2)
hold on
plot(newdata_long,H1b,'-.k','LineWidth',2)
hold on
plot(newdata_long,H1c,'-.m','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'K_v=0.0048','K_v=0.0037','K_v=0.0030','K_v=0'},'Location','northeast')
grid on

subplot(2,2,4)
plot(newdata_long,D1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,D1a,'-.b','LineWidth',2)
hold on
plot(newdata_long,D1b,'-.k','LineWidth',2)
hold on
plot(newdata_long,D1c,'-.m','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'K_v=0.0048','K_v=0.0037','K_v=0.0030','K_v=0'},'Location','northeast')
grid on

%% Varying Acquired Waning Immunity with face mask usage of 50& and 50& efficiency

Cm=0.5;
epsilonm=0.5;
Kv1=1/30;
Kv2=1/150;
[t,S1o,E1o,I1o,A1o,H1o,V1o,V2o,EBo,IBo,ABo,HBo,Ro,Wo,Wvo,V3o,Do, ...
    I1cumulao,A1cumulao,IBcumulao,ABcumulao,D1o,DBo]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.5;
epsilonm=0.5;
Kv1=1/60;
Kv2=1/180;
[t,S1a,E1a,I1a,A1a,H1a,V1a,V2a,EBa,IBa,ABa,HBa,Ra,Wa,Wva,V3a,Da, ...
    I1cumulaa,A1cumulaa,IBcumulaa,ABcumulaa,D1a,DBa]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.5;
epsilonm=0.5;
Kv1=1/90;
Kv2=1/210;
[t,S1b,E1b,I1b,A1b,H1b,V1b,V2b,EBb,IBb,ABb,HBb,Rb,Wb,Wvb,V3b,Db, ...
    I1cumulab,A1cumulab,IBcumulab,ABcumulab,D1b,DBb]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.5;
epsilonm=0.5;
Kv1=0;
Kv2=0;
[t,S1c,E1c,I1c,A1c,H1c,V1c,V2c,EBc,IBc,ABc,HBc,Rc,Wc,Wvc,V3c,Dc, ...
    I1cumulac,A1cumulac,IBcumulac,ABcumulac,D1c,DBc]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
%Plot Unvaccinated Individuals 
figure
subplot(2,2,1)
plot(newdata_long,IBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,IBa,'-.b','LineWidth',2)
hold on
plot(newdata_long,IBb,'-.k','LineWidth',2)
hold on
plot(newdata_long,IBc,'-.m','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'K_v_1=0.0333,K_v_2=0.0067',...
    'K_v_1=0.0167,K_v_2=0.0056',...
    'K_v_1= 0.0111,K_v_2=0.0048',...
    'K_v_1=0.,K_v_2=0'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,ABo,'-.r','LineWidth',2)
hold on
plot(newdata_long,ABa,'-.b','LineWidth',2)
hold on
plot(newdata_long,ABb,'-.k','LineWidth',2)
hold on
plot(newdata_long,ABc,'-.m','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'K_v_1=0.0333,K_v_2=0.0067',...
    'K_v_1=0.0167,K_v_2=0.0056',...
    'K_v_1= 0.0111,K_v_2=0.0048',...
    'K_v_1=0.,K_v_2=0'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,HBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,HBa,'-.b','LineWidth',2)
hold on
plot(newdata_long,HBb,'-.k','LineWidth',2)
hold on
plot(newdata_long,HBc,'-.m','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'K_v_1=0.0333,K_v_2=0.0067',...
    'K_v_1=0.0167,K_v_2=0.0056',...
    'K_v_1= 0.0111,K_v_2=0.0048',...
    'K_v_1=0.,K_v_2=0'},'Location','northeast')
grid on

subplot(2,2,4)
plot(newdata_long,DBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,DBa,'-.b','LineWidth',2)
hold on
plot(newdata_long,DBb,'-.k','LineWidth',2)
hold on
plot(newdata_long,DBc,'-.m','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'K_v_1=0.0333,K_v_2=0.0067',...
    'K_v_1=0.0167,K_v_2=0.0056',...
    'K_v_1= 0.0111,K_v_2=0.0048',...
    'K_v_1=0.,K_v_2=0'},'Location','northeast')
grid on

%% Varying Acquired Waning Immunity with face mask usage of 50& and 50& efficiency
Kn=1/270;
Cm=0.5;
epsilonm=0.5;
Kv=0.0666666;
[t,S1o,E1o,I1o,A1o,H1o,V1o,V2o,EBo,IBo,ABo,HBo,Ro,Wo,Wvo,V3o,Do, ...
    I1cumulao,A1cumulao,IBcumulao,ABcumulao,D1o,DBo]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.5;
epsilonm=0.5;
Kv=0.033333;
[t,S1a,E1a,I1a,A1a,H1a,V1a,V2a,EBa,IBa,ABa,HBa,Ra,Wa,Wva,V3a,Da, ...
    I1cumulaa,A1cumulaa,IBcumulaa,ABcumulaa,D1a,DBa]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.5;
epsilonm=0.5;
Kv=0.022222;
[t,S1b,E1b,I1b,A1b,H1b,V1b,V2b,EBb,IBb,ABb,HBb,Rb,Wb,Wvb,V3b,Db, ...
    I1cumulab,A1cumulab,IBcumulab,ABcumulab,D1b,DBb]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.5;
epsilonm=0.5;
Kv=0.0;
[t,S1c,E1c,I1c,A1c,H1c,V1c,V2c,EBc,IBc,ABc,HBc,Rc,Wc,Wvc,V3c,Dc, ...
    I1cumulac,A1cumulac,IBcumulac,ABcumulac,D1c,DBc]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
%Plot Unvaccinated Individuals 
figure
subplot(2,2,1)
plot(newdata_long,IBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,IBa,'-.b','LineWidth',2)
hold on
plot(newdata_long,IBb,'-.k','LineWidth',2)
hold on
plot(newdata_long,IBc,'-.m','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'K_v=0.0048','K_v=0.0037','K_v=0.0030','K_v=0'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,A1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,A1a,'-.b','LineWidth',2)
hold on
plot(newdata_long,A1b,'-.k','LineWidth',2)
hold on
plot(newdata_long,A1c,'-.m','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'K_v=0.0048','K_v=0.0037','K_v=0.0030','K_v=0'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,H1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,H1a,'-.b','LineWidth',2)
hold on
plot(newdata_long,H1b,'-.k','LineWidth',2)
hold on
plot(newdata_long,H1c,'-.m','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'K_v=0.0048','K_v=0.0037','K_v=0.0030','K_v=0'},'Location','northeast')
grid on

subplot(2,2,4)
plot(newdata_long,D1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,D1a,'-.b','LineWidth',2)
hold on
plot(newdata_long,D1b,'-.k','LineWidth',2)
hold on
plot(newdata_long,D1c,'-.m','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'K_v=0.0048','K_v=0.0037','K_v=0.0030','K_v=0'},'Location','northeast')
grid on


%% Natural Immunity at baseline with 30% protection and coverage
Cm=0; 
epsilonm=0.3;
Kn=0.0037;
[t,S1o,E1o,I1o,A1o,H1o,V1o,V2o,EBo,IBo,ABo,HBo,Ro,Wo,Wvo,V3o,Do, ...
    I1cumulao,A1cumulao,IBcumulao,ABcumulao,D1o,DBo]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.25; 
epsilonm=0.3;
Kn=0.0037;
[t,S1a,E1a,I1a,A1a,H1a,V1a,V2a,EBa,IBa,ABa,HBa,Ra,Wa,Wva,V3a,Da, ...
    I1cumulaa,A1cumulaa,IBcumulaa,ABcumulaa,D1a,DBa]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.50; 
epsilonm=0.3;
Kn=0.0037;
[t,S1b,E1b,I1b,A1b,H1b,V1b,V2b,EBb,IBb,ABb,HBb,Rb,Wb,Wvb,V3b,Db, ...
    I1cumulab,A1cumulab,IBcumulab,ABcumulab,D1b,DBb]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.75; 
epsilonm=0.3;
Kn=0.0037;
[t,S1c,E1c,I1c,A1c,H1c,V1c,V2c,EBc,IBc,ABc,HBc,Rc,Wc,Wvc,V3c,Dc, ...
    I1cumulac,A1cumulac,IBcumulac,ABcumulac,D1c,DBc]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=1.0; 
epsilonm=0.3;
Kn=0.0037;
[t,S1d,E1d,I1d,A1d,H1d,V1d,V2d,EBd,IBd,ABd,HBd,Rd,Wd,Wvd,V3d,Dd, ...
    I1cumulad,A1cumulad,IBcumulad,ABcumulad,D1d,DBd]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);


figure
subplot(2,2,1)
plot(newdata_long,I1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,I1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,I1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,I1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,I1d,'-.k','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,A1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,A1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,A1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,A1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,A1d,'-.k','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,H1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,H1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,H1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,H1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,H1d,'-.k','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,4)
plot(newdata_long,D1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,D1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,D1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,D1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,D1d,'-.k','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

% Plot for Vaccinated Individuals 
figure
subplot(2,2,1)
plot(newdata_long,IBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,IBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,IBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,IBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,IBd,'-.k','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,ABo,'-.r','LineWidth',2)
hold on
plot(newdata_long,ABa,'-.m','LineWidth',2)
hold on
plot(newdata_long,ABb,'-.b','LineWidth',2)
hold on
plot(newdata_long,ABc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,ABd,'-.k','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,HBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,HBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,HBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,HBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,HBd,'-.k','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,4)
plot(newdata_long,DBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,DBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,DBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,DBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,DBd,'-.k','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

%% Natural Immunity at baseline with 70% protection and coverage
Cm=0; 
epsilonm=0.7;
Kn=0.0037;
[t,S1o,E1o,I1o,A1o,H1o,V1o,V2o,EBo,IBo,ABo,HBo,Ro,Wo,Wvo,V3o,Do, ...
    I1cumulao,A1cumulao,IBcumulao,ABcumulao,D1o,DBo]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.25; 
epsilonm=0.7;
Kn=0.0037;
[t,S1a,E1a,I1a,A1a,H1a,V1a,V2a,EBa,IBa,ABa,HBa,Ra,Wa,Wva,V3a,Da, ...
    I1cumulaa,A1cumulaa,IBcumulaa,ABcumulaa,D1a,DBa]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.50; 
epsilonm=0.7;
Kn=0.0037;
[t,S1b,E1b,I1b,A1b,H1b,V1b,V2b,EBb,IBb,ABb,HBb,Rb,Wb,Wvb,V3b,Db, ...
    I1cumulab,A1cumulab,IBcumulab,ABcumulab,D1b,DBb]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.75; 
epsilonm=0.7;
Kn=0.0037;
[t,S1c,E1c,I1c,A1c,H1c,V1c,V2c,EBc,IBc,ABc,HBc,Rc,Wc,Wvc,V3c,Dc, ...
    I1cumulac,A1cumulac,IBcumulac,ABcumulac,D1c,DBc]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=1.0; 
epsilonm=0.7;
Kn=0.0037;
[t,S1d,E1d,I1d,A1d,H1d,V1d,V2d,EBd,IBd,ABd,HBd,Rd,Wd,Wvd,V3d,Dd, ...
    I1cumulad,A1cumulad,IBcumulad,ABcumulad,D1d,DBd]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);


figure
subplot(2,2,1)
plot(newdata_long,I1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,I1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,I1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,I1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,I1d,'-.k','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,A1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,A1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,A1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,A1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,A1d,'-.k','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,H1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,H1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,H1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,H1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,H1d,'-.k','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,4)
plot(newdata_long,D1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,D1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,D1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,D1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,D1d,'-.k','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

% Plot for Vaccinated Individuals 
figure
subplot(2,2,1)
plot(newdata_long,IBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,IBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,IBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,IBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,IBd,'-.k','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,ABo,'-.r','LineWidth',2)
hold on
plot(newdata_long,ABa,'-.m','LineWidth',2)
hold on
plot(newdata_long,ABb,'-.b','LineWidth',2)
hold on
plot(newdata_long,ABc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,ABd,'-.k','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,HBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,HBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,HBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,HBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,HBd,'-.k','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,4)
plot(newdata_long,DBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,DBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,DBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,DBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,DBd,'-.k','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

%% Natural Immunity at baseline with 95% protection and coverage
Cm=0; 
epsilonm=0.95;
Kn=0.0037;
[t,S1o,E1o,I1o,A1o,H1o,V1o,V2o,EBo,IBo,ABo,HBo,Ro,Wo,Wvo,V3o,Do, ...
    I1cumulao,A1cumulao,IBcumulao,ABcumulao,D1o,DBo]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.25; 
epsilonm=0.95;
Kn=0.0037;
[t,S1a,E1a,I1a,A1a,H1a,V1a,V2a,EBa,IBa,ABa,HBa,Ra,Wa,Wva,V3a,Da, ...
    I1cumulaa,A1cumulaa,IBcumulaa,ABcumulaa,D1a,DBa]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.50; 
epsilonm=0.95;
Kn=0.0037;
[t,S1b,E1b,I1b,A1b,H1b,V1b,V2b,EBb,IBb,ABb,HBb,Rb,Wb,Wvb,V3b,Db, ...
    I1cumulab,A1cumulab,IBcumulab,ABcumulab,D1b,DBb]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.75; 
epsilonm=0.95;
Kn=0.0037;
[t,S1c,E1c,I1c,A1c,H1c,V1c,V2c,EBc,IBc,ABc,HBc,Rc,Wc,Wvc,V3c,Dc, ...
    I1cumulac,A1cumulac,IBcumulac,ABcumulac,D1c,DBc]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=1.0; 
epsilonm=0.95;
Kn=0.0037;
[t,S1d,E1d,I1d,A1d,H1d,V1d,V2d,EBd,IBd,ABd,HBd,Rd,Wd,Wvd,V3d,Dd, ...
    I1cumulad,A1cumulad,IBcumulad,ABcumulad,D1d,DBd]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);


figure
subplot(2,2,1)
plot(newdata_long,I1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,I1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,I1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,I1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,I1d,'-.k','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,A1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,A1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,A1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,A1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,A1d,'-.k','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,H1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,H1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,H1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,H1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,H1d,'-.k','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,4)
plot(newdata_long,D1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,D1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,D1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,D1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,D1d,'-.k','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

% Plot for Vaccinated Individuals 
figure
subplot(2,2,1)
plot(newdata_long,IBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,IBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,IBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,IBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,IBd,'-.k','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,ABo,'-.r','LineWidth',2)
hold on
plot(newdata_long,ABa,'-.m','LineWidth',2)
hold on
plot(newdata_long,ABb,'-.b','LineWidth',2)
hold on
plot(newdata_long,ABc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,ABd,'-.k','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,HBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,HBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,HBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,HBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,HBd,'-.k','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on

subplot(2,2,4)
plot(newdata_long,DBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,DBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,DBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,DBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,DBd,'-.k','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_m=0','C_m=0.25','C_m=0.50','C_m=0.75','C_m=1.0'},'Location','northeast')
grid on


%% Control Reproduction Number
lambda = 3747540/365; % Birth rate (people/day)     Obtained from https://www.cdc.gov/nchs/fastats/births.htm
d = (869.7/100000)/365; % Natural death rate (1/day)    Obtained from https://www.cdc.gov/nchs/fastats/deaths.htm
w = 1/4.;     % Transfer from exposed to infectious (1/day)
     % Waning rate of natural immunity (1/day)
p = 1/9;   % Proportion of symptomatic infections
gamma = 1/14;
beta1=2.0989e-09;
beta2=2.2965e-09;
alpha=0.1180;
% x: proportion of face masks usage
% y: face mask efficiencyt
[X,Y]=meshgrid(0:0.01:1.0, 0:0.01:1.0);
epsilon=0.1;

RC = @(x,y) (lambda.*(1-epsilon).*beta1.*(1-x.*y).*p.*w)/(d.*(w+d).*(alpha+gamma+d)) ...
     +(lambda.*(1-epsilon).*beta2.*(1-x.*y).*(1-p).*w)/(d.*(w+d).*(gamma+d));

    
z = RC(X,Y);

figure
subplot(1,3,1)
[C,h]=contourf(X,Y,z, [0.0 : 0.5 : 20.0]);
clabel(C,h,'FontSize',7,'Color','k','FontWeight','bold')
xlabel('Usage of Face Mask')
ylabel('Face Mask Efficiency')
title('Low Vaccine Efficiency, \epsilon_L_2=0.1')
text(0.5,0.5,'\bfA','FontSize',12)

[X,Y]=meshgrid(0:0.01:1.0, 0:0.01:1.0);
epsilon=0.2;

RC = @(x,y) (lambda.*(1-epsilon).*beta1.*(1-x.*y).*p.*w)/(d.*(w+d).*(alpha+gamma+d)) ...
     +(lambda.*(1-epsilon).*beta2.*(1-x.*y).*(1-p).*w)/(d.*(w+d).*(gamma+d));

    
z = RC(X,Y);
subplot(1,3,2)
[C,h]=contourf(X,Y,z, [0.0 : 0.5 : 20.0]);
clabel(C,h,'FontSize',7,'Color','k','FontWeight','bold')
xlabel('Usage of Face Mask')
ylabel('Face Mask Efficiency')
title('Baseline Vaccine Efficiency, \epsilon_L_2=0.2')
text(0.5,0.5,'\bfB','FontSize',12)

[X,Y]=meshgrid(0:0.01:1.0, 0:0.01:1.0);
epsilon=0.35;

RC = @(x,y) (lambda.*(1-epsilon).*beta1.*(1-x.*y).*p.*w)/(d.*(w+d).*(alpha+gamma+d)) ...
     +(lambda.*(1-epsilon).*beta2.*(1-x.*y).*(1-p).*w)/(d.*(w+d).*(gamma+d));

    
z = RC(X,Y);
subplot(1,3,3)
[C,h]=contourf(X,Y,z, [0.0 : 0.5 : 20.0]);
clabel(C,h,'FontSize',7,'Color','k','FontWeight','bold')
xlabel('Usage of Face Mask')
ylabel('Face Mask Efficiency')
title('High Vaccine Efficiency, \epsilon_L_2=0.35')
text(0.5,0.5,'\bfC','FontSize',12)

%% Control Reproduction Number with booster dose efficiency
lambda = 3747540/365; % Birth rate (people/day)     Obtained from https://www.cdc.gov/nchs/fastats/births.htm
d = (869.7/100000)/365; % Natural death rate (1/day)    Obtained from https://www.cdc.gov/nchs/fastats/deaths.htm
w = 1/4.;     % Transfer from exposed to infectious (1/day)
     % Waning rate of natural immunity (1/day)
p = 1/9;   % Proportion of symptomatic infections
gamma = 1/14;
beta1=2.0989e-09;
beta2=2.2965e-09;
alpha=0.1180;
% x: proportion of face masks usage
% y: face mask efficiencyt
[X,Y]=meshgrid(0:0.01:1.0, 0:0.01:1.0);
epsilon=0.55;

RC = @(x,y) (lambda.*(1-epsilon).*beta1.*(1-x.*y).*p.*w)/(d.*(w+d).*(alpha+gamma+d)) ...
     +(lambda.*(1-epsilon).*beta2.*(1-x.*y).*(1-p).*w)/(d.*(w+d).*(gamma+d));

    
z = RC(X,Y);

figure
subplot(1,3,1)
[C,h]=contourf(X,Y,z, [0.0 : 0.5 : 20.0]);
clabel(C,h,'FontSize',7,'Color','k','FontWeight','bold')
xlabel('Usage of Face Mask')
ylabel('Face Mask Efficiency')
title('Booster Low Efficency, \epsilon_L_2=0.55')
text(0.5,0.5,'\bfA','FontSize',12)

[X,Y]=meshgrid(0:0.01:1.0, 0:0.01:1.0);
epsilon=0.63;

RC = @(x,y) (lambda.*(1-epsilon).*beta1.*(1-x.*y).*p.*w)/(d.*(w+d).*(alpha+gamma+d)) ...
     +(lambda.*(1-epsilon).*beta2.*(1-x.*y).*(1-p).*w)/(d.*(w+d).*(gamma+d));

    
z = RC(X,Y);
subplot(1,3,2)
[C,h]=contourf(X,Y,z, [0.0 : 0.5 : 20.0]);
clabel(C,h,'FontSize',7,'Color','k','FontWeight','bold')
xlabel('Usage of Face Mask')
ylabel('Face Mask Efficiency')
title('Booster Baseline Efficiency, \epsilon_L_2=0.63')
text(0.5,0.5,'\bfB','FontSize',12)

[X,Y]=meshgrid(0:0.01:1.0, 0:0.01:1.0);
epsilon=0.75;

RC = @(x,y) (lambda.*(1-epsilon).*beta1.*(1-x.*y).*p.*w)/(d.*(w+d).*(alpha+gamma+d)) ...
     +(lambda.*(1-epsilon).*beta2.*(1-x.*y).*(1-p).*w)/(d.*(w+d).*(gamma+d));

    
z = RC(X,Y);
subplot(1,3,3)
[C,h]=contourf(X,Y,z, [0.0 : 0.5 : 20.0]);
clabel(C,h,'FontSize',7,'Color','k','FontWeight','bold')
xlabel('Usage of Face Mask')
ylabel('Face Mask Efficiency')
title('Booster High Efficiency, \epsilon_L_2=0.75')
text(0.5,0.5,'\bfC','FontSize',12)
%% What will happen if other non pharmaceutical strategies were applied?? 
% Lets vary the usage of face masks where they provide a 30% of protection
Cm=0.5; 
epsilonm=0.3;
Cr=0.0;
[t,S1o,E1o,I1o,A1o,H1o,V1o,V2o,EBo,IBo,ABo,HBo,Ro,Wo,Wvo,V3o,Do, ...
    I1cumulao,A1cumulao,IBcumulao,ABcumulao,D1o,DBo]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.5; 
epsilonm=0.3;
Cr=0.05;
[t,S1a,E1a,I1a,A1a,H1a,V1a,V2a,EBa,IBa,ABa,HBa,Ra,Wa,Wva,V3a,Da, ...
    I1cumulaa,A1cumulaa,IBcumulaa,ABcumulaa,D1a,DBa]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.50; 
epsilonm=0.3;
Cr=0.1;
[t,S1b,E1b,I1b,A1b,H1b,V1b,V2b,EBb,IBb,ABb,HBb,Rb,Wb,Wvb,V3b,Db, ...
    I1cumulab,A1cumulab,IBcumulab,ABcumulab,D1b,DBb]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.5; 
epsilonm=0.3;
Cr=0.15;
[t,S1c,E1c,I1c,A1c,H1c,V1c,V2c,EBc,IBc,ABc,HBc,Rc,Wc,Wvc,V3c,Dc, ...
    I1cumulac,A1cumulac,IBcumulac,ABcumulac,D1c,DBc]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.5; 
epsilonm=0.3;
Cr=0.2;
[t,S1d,E1d,I1d,A1d,H1d,V1d,V2d,EBd,IBd,ABd,HBd,Rd,Wd,Wvd,V3d,Dd, ...
    I1cumulad,A1cumulad,IBcumulad,ABcumulad,D1d,DBd]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);


figure
subplot(2,2,1)
plot(newdata_long,I1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,I1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,I1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,I1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,I1d,'-.k','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,A1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,A1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,A1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,A1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,A1d,'-.k','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,H1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,H1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,H1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,H1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,H1d,'-.k','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

subplot(2,2,4)
plot(newdata_long,D1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,D1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,D1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,D1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,D1d,'-.k','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

%Vaccinated Individuals
figure
subplot(2,2,1)
plot(newdata_long,IBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,IBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,IBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,IBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,IBd,'-.k','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,ABo,'-.r','LineWidth',2)
hold on
plot(newdata_long,ABa,'-.m','LineWidth',2)
hold on
plot(newdata_long,ABb,'-.b','LineWidth',2)
hold on
plot(newdata_long,ABc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,ABd,'-.k','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,HBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,HBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,HBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,HBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,HBd,'-.k','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

subplot(2,2,4)
plot(newdata_long,DBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,DBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,DBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,DBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,DBd,'-.k','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

%Waning Immunity Vaccinated and Unvaccinted 
subplot(1,2,1)
plot(newdata_long,Wo,'-.r','LineWidth',2)
hold on
plot(newdata_long,Wa,'-.m','LineWidth',2)
hold on
plot(newdata_long,Wb,'-.b','LineWidth',2)
hold on
plot(newdata_long,Wc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,Wd,'-.k','LineWidth',2)
title('Waning Natural Immunity Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

subplot(1,2,2)
plot(newdata_long,Wvo,'-.r','LineWidth',2)
hold on
plot(newdata_long,Wva,'-.m','LineWidth',2)
hold on
plot(newdata_long,Wvb,'-.b','LineWidth',2)
hold on
plot(newdata_long,Wvc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,Wvd,'-.k','LineWidth',2)
title('Waning Acquired Immunity  Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
text(0.5,0.5,'\bfB','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

%% What will happen if other non pharmaceutical strategies were applied?? 
% Lets vary the usage of face masks where they provide a 50% of protection
Cm=0.5; 
epsilonm=0.5;
Cr=0.0;
[t,S1o,E1o,I1o,A1o,H1o,V1o,V2o,EBo,IBo,ABo,HBo,Ro,Wo,Wvo,V3o,Do, ...
    I1cumulao,A1cumulao,IBcumulao,ABcumulao,D1o,DBo]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.5; 
epsilonm=0.5;
Cr=0.05;
[t,S1a,E1a,I1a,A1a,H1a,V1a,V2a,EBa,IBa,ABa,HBa,Ra,Wa,Wva,V3a,Da, ...
    I1cumulaa,A1cumulaa,IBcumulaa,ABcumulaa,D1a,DBa]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.50; 
epsilonm=0.5;
Cr=0.1;
[t,S1b,E1b,I1b,A1b,H1b,V1b,V2b,EBb,IBb,ABb,HBb,Rb,Wb,Wvb,V3b,Db, ...
    I1cumulab,A1cumulab,IBcumulab,ABcumulab,D1b,DBb]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.5; 
epsilonm=0.5;
Cr=0.15;
[t,S1c,E1c,I1c,A1c,H1c,V1c,V2c,EBc,IBc,ABc,HBc,Rc,Wc,Wvc,V3c,Dc, ...
    I1cumulac,A1cumulac,IBcumulac,ABcumulac,D1c,DBc]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.5; 
epsilonm=0.5;
Cr=0.2;
[t,S1d,E1d,I1d,A1d,H1d,V1d,V2d,EBd,IBd,ABd,HBd,Rd,Wd,Wvd,V3d,Dd, ...
    I1cumulad,A1cumulad,IBcumulad,ABcumulad,D1d,DBd]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);


figure
subplot(2,2,1)
plot(newdata_long,I1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,I1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,I1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,I1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,I1d,'-.k','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,A1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,A1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,A1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,A1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,A1d,'-.k','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,H1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,H1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,H1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,H1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,H1d,'-.k','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

subplot(2,2,4)
plot(newdata_long,D1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,D1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,D1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,D1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,D1d,'-.k','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

%Vaccinated Individuals
figure
subplot(2,2,1)
plot(newdata_long,IBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,IBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,IBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,IBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,IBd,'-.k','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,ABo,'-.r','LineWidth',2)
hold on
plot(newdata_long,ABa,'-.m','LineWidth',2)
hold on
plot(newdata_long,ABb,'-.b','LineWidth',2)
hold on
plot(newdata_long,ABc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,ABd,'-.k','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,HBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,HBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,HBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,HBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,HBd,'-.k','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

subplot(2,2,4)
plot(newdata_long,DBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,DBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,DBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,DBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,DBd,'-.k','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

%Waning Immunity Vaccinated and Unvaccinted 
subplot(1,2,1)
plot(newdata_long,Wo,'-.r','LineWidth',2)
hold on
plot(newdata_long,Wa,'-.m','LineWidth',2)
hold on
plot(newdata_long,Wb,'-.b','LineWidth',2)
hold on
plot(newdata_long,Wc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,Wd,'-.k','LineWidth',2)
title('Waning Natural Immunity Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

subplot(1,2,2)
plot(newdata_long,Wvo,'-.r','LineWidth',2)
hold on
plot(newdata_long,Wva,'-.m','LineWidth',2)
hold on
plot(newdata_long,Wvb,'-.b','LineWidth',2)
hold on
plot(newdata_long,Wvc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,Wvd,'-.k','LineWidth',2)
title('Waning Acquired Immunity  Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
text(0.5,0.5,'\bfB','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

%% What will happen if other non pharmaceutical strategies were applied?? 
% Lets vary the usage of face masks where they provide a 70% of protection
Cm=0.5; 
epsilonm=0.7;
Cr=0.0;
[t,S1o,E1o,I1o,A1o,H1o,V1o,V2o,EBo,IBo,ABo,HBo,Ro,Wo,Wvo,V3o,Do, ...
    I1cumulao,A1cumulao,IBcumulao,ABcumulao,D1o,DBo]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.5; 
epsilonm=0.7;
Cr=0.05;
[t,S1a,E1a,I1a,A1a,H1a,V1a,V2a,EBa,IBa,ABa,HBa,Ra,Wa,Wva,V3a,Da, ...
    I1cumulaa,A1cumulaa,IBcumulaa,ABcumulaa,D1a,DBa]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.50; 
epsilonm=0.7;
Cr=0.1;
[t,S1b,E1b,I1b,A1b,H1b,V1b,V2b,EBb,IBb,ABb,HBb,Rb,Wb,Wvb,V3b,Db, ...
    I1cumulab,A1cumulab,IBcumulab,ABcumulab,D1b,DBb]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.5; 
epsilonm=0.7;
Cr=0.15;
[t,S1c,E1c,I1c,A1c,H1c,V1c,V2c,EBc,IBc,ABc,HBc,Rc,Wc,Wvc,V3c,Dc, ...
    I1cumulac,A1cumulac,IBcumulac,ABcumulac,D1c,DBc]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.5; 
epsilonm=0.7;
Cr=0.2;
[t,S1d,E1d,I1d,A1d,H1d,V1d,V2d,EBd,IBd,ABd,HBd,Rd,Wd,Wvd,V3d,Dd, ...
    I1cumulad,A1cumulad,IBcumulad,ABcumulad,D1d,DBd]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);


figure
subplot(2,2,1)
plot(newdata_long,I1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,I1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,I1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,I1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,I1d,'-.k','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,A1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,A1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,A1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,A1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,A1d,'-.k','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,H1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,H1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,H1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,H1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,H1d,'-.k','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

subplot(2,2,4)
plot(newdata_long,D1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,D1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,D1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,D1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,D1d,'-.k','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

%Vaccinated Individuals
figure
subplot(2,2,1)
plot(newdata_long,IBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,IBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,IBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,IBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,IBd,'-.k','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,ABo,'-.r','LineWidth',2)
hold on
plot(newdata_long,ABa,'-.m','LineWidth',2)
hold on
plot(newdata_long,ABb,'-.b','LineWidth',2)
hold on
plot(newdata_long,ABc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,ABd,'-.k','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,HBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,HBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,HBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,HBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,HBd,'-.k','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

subplot(2,2,4)
plot(newdata_long,DBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,DBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,DBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,DBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,DBd,'-.k','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

%Waning Immunity Vaccinated and Unvaccinted 
subplot(1,2,1)
plot(newdata_long,Wo,'-.r','LineWidth',2)
hold on
plot(newdata_long,Wa,'-.m','LineWidth',2)
hold on
plot(newdata_long,Wb,'-.b','LineWidth',2)
hold on
plot(newdata_long,Wc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,Wd,'-.k','LineWidth',2)
title('Waning Natural Immunity Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

subplot(1,2,2)
plot(newdata_long,Wvo,'-.r','LineWidth',2)
hold on
plot(newdata_long,Wva,'-.m','LineWidth',2)
hold on
plot(newdata_long,Wvb,'-.b','LineWidth',2)
hold on
plot(newdata_long,Wvc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,Wvd,'-.k','LineWidth',2)
title('Waning Acquired Immunity  Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
text(0.5,0.5,'\bfB','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

%% What will happen if other non pharmaceutical strategies were applied?? 
% Lets vary the usage of face masks where they provide a 70% of protection
Cm=0.5; 
epsilonm=0.95;
Cr=0.0;
[t,S1o,E1o,I1o,A1o,H1o,V1o,V2o,EBo,IBo,ABo,HBo,Ro,Wo,Wvo,V3o,Do, ...
    I1cumulao,A1cumulao,IBcumulao,ABcumulao,D1o,DBo]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.5; 
epsilonm=0.95;
Cr=0.05;
[t,S1a,E1a,I1a,A1a,H1a,V1a,V2a,EBa,IBa,ABa,HBa,Ra,Wa,Wva,V3a,Da, ...
    I1cumulaa,A1cumulaa,IBcumulaa,ABcumulaa,D1a,DBa]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.50; 
epsilonm=0.95;
Cr=0.1;
[t,S1b,E1b,I1b,A1b,H1b,V1b,V2b,EBb,IBb,ABb,HBb,Rb,Wb,Wvb,V3b,Db, ...
    I1cumulab,A1cumulab,IBcumulab,ABcumulab,D1b,DBb]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.5; 
epsilonm=0.95;
Cr=0.15;
[t,S1c,E1c,I1c,A1c,H1c,V1c,V2c,EBc,IBc,ABc,HBc,Rc,Wc,Wvc,V3c,Dc, ...
    I1cumulac,A1cumulac,IBcumulac,ABcumulac,D1c,DBc]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.5; 
epsilonm=0.95;
Cr=0.2;
[t,S1d,E1d,I1d,A1d,H1d,V1d,V2d,EBd,IBd,ABd,HBd,Rd,Wd,Wvd,V3d,Dd, ...
    I1cumulad,A1cumulad,IBcumulad,ABcumulad,D1d,DBd]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);


figure
subplot(2,2,1)
plot(newdata_long,I1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,I1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,I1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,I1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,I1d,'-.k','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,A1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,A1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,A1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,A1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,A1d,'-.k','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,H1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,H1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,H1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,H1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,H1d,'-.k','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

subplot(2,2,4)
plot(newdata_long,D1o,'-.r','LineWidth',2)
hold on
plot(newdata_long,D1a,'-.m','LineWidth',2)
hold on
plot(newdata_long,D1b,'-.b','LineWidth',2)
hold on
plot(newdata_long,D1c,'-.c','LineWidth',2)
hold on 
plot(newdata_long,D1d,'-.k','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

%Vaccinated Individuals
figure
subplot(2,2,1)
plot(newdata_long,IBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,IBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,IBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,IBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,IBd,'-.k','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,ABo,'-.r','LineWidth',2)
hold on
plot(newdata_long,ABa,'-.m','LineWidth',2)
hold on
plot(newdata_long,ABb,'-.b','LineWidth',2)
hold on
plot(newdata_long,ABc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,ABd,'-.k','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,HBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,HBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,HBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,HBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,HBd,'-.k','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

subplot(2,2,4)
plot(newdata_long,DBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,DBa,'-.m','LineWidth',2)
hold on
plot(newdata_long,DBb,'-.b','LineWidth',2)
hold on
plot(newdata_long,DBc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,DBd,'-.k','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

%Waning Immunity Vaccinated and Unvaccinted 
subplot(1,2,1)
plot(newdata_long,Wo,'-.r','LineWidth',2)
hold on
plot(newdata_long,Wa,'-.m','LineWidth',2)
hold on
plot(newdata_long,Wb,'-.b','LineWidth',2)
hold on
plot(newdata_long,Wc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,Wd,'-.k','LineWidth',2)
title('Waning Natural Immunity Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

subplot(1,2,2)
plot(newdata_long,Wvo,'-.r','LineWidth',2)
hold on
plot(newdata_long,Wva,'-.m','LineWidth',2)
hold on
plot(newdata_long,Wvb,'-.b','LineWidth',2)
hold on
plot(newdata_long,Wvc,'-.c','LineWidth',2)
hold on 
plot(newdata_long,Wvd,'-.k','LineWidth',2)
title('Waning Acquired Immunity  Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
text(0.5,0.5,'\bfB','FontSize',12)
legend({'C_r=0','C_r=0.05','C_r=0.1','C_r=0.15','C_r=0.2'},'Location','northeast')
grid on

%% Seting initial conditions for simulations

S10 = S1a(end);
E10 = E1a(end);
I10 = I1a(end);
A10 = A1a(end);
H10 = H1a(end);
V10 = V1a(end);
V20 = V2a(end);
EB0 = EBa(end);
IB0 = IBa(end);
AB0 = ABa(end);
HB0 = HBa(end);
R0 = Ra(end);
W0 = Wa(end);
WV0= Wva(end);
V30= V3a(end);
D0 = 832180;
I1cumula0 = I1cumula(end);
A1cumula0 = A1cumula(end);
IBcumula0 = IBcumula(end);
ABcumula0 = ABcumula(end);
D10=D0*D1a(end)/( D1a(end)+DBa(end) );
DB0=D0*DBa(end)/( D1a(end)+DBa(end) );

IniVal = [S10;E10;I10;A10;H10;V10;V20;EB0;IB0;AB0;HB0;R0;...
    W0;WV0;V30;D0;I1cumula0;A1cumula0;IBcumula0;ABcumula0;D10;DB0];

newinitdate=714;   % 4 Jan 2022

numberofdays = 300;
t=linspace(tsim(end), tsim(end)+numberofdays,numberofdays+1);
newdata_long=dateshift(dateslong(388),'start','day',0:numberofdays);

%% Defining vaccination rates
rho1 = @(x) (2.1268e+04)*x.*exp(-0.015479*x);
rho2 = @(x) (1.8017e+04)*x.*exp(-0.010932*x);
rho3 = @(x) (8.6062e+03)*(x-81).*exp(-0.023304*(x-81)).*( x>=81 );
rho4 = @(x) (1e-03).*( x>=242 );    % Booster doses began on Aug-13-3021
date1 = datenum('2022/04/01');
date2 = datenum('2022/10/01');
Cm = @(x) 0.5.*(x < date1-datenum(dateslong(1))) ...
       + 0.25.*( (x >= date1-datenum(dateslong(1)))&(x < date2-datenum(dateslong(1))) ) ...
       + 0.*( x >= date2-datenum(dateslong(1)) );

%% Solving the system
[t,S10,E10,I10,A10,H10,V10,V20,EB0,IB0,AB0,HB0,R0,W0,Wv0,V30,D0, ...
    I1cumula0,A1cumula0,IBcumula0,ABcumula0,D10,DB0]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);

%%Solving the system without varying the use of face masks
Cm=0.5;
epsilonm=0.3;
[t,S1d,E1d,I1d,A1d,H1d,V1d,V2d,EBd,IBd,ABd,HBd,Rd,Wd,Wvd,V3d,Dd, ...
    I1cumulad,A1cumulad,IBcumulad,ABcumulad,D1d,DBd]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);



%% Lets vary the usage of face masks with 30 % efficiency
epsilonm=0.3;
date1 = datenum('2022/04/01');
date2 = datenum('2022/10/01');
Cm = @(x) 0.5.*(x < date1-datenum(dateslong(1))) ...
       + 0.25.*( (x >= date1-datenum(dateslong(1)))&(x < date2-datenum(dateslong(1))) ) ...
       + 0.*( x >= date2-datenum(dateslong(1)) );
[t,S1o,E1o,I1o,A1o,H1o,V1o,V2o,EBo,IBo,ABo,HBo,Ro,Wo,Wvo,V3o,Do, ...
    I1cumulao,A1cumulao,IBcumulao,ABcumulao,D1o,DBo]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.5;
epsilonm=0.3;
[t,S1d,E1d,I1d,A1d,H1d,V1d,V2d,EBd,IBd,ABd,HBd,Rd,Wd,Wvd,V3d,Dd, ...
    I1cumulad,A1cumulad,IBcumulad,ABcumulad,D1d,DBd]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
figure
subplot(2,2,1)
plot(newdata_long,I1o,'-.r','LineWidth',2)
hold on 
plot(newdata_long,I1d,'-.k','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'\epsilon_m=0.3','No Variation'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,A1o,'-.k','LineWidth',2)
hold on 
plot(newdata_long,A1d,'-.r','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'\epsilon_m=0.3','No Variation'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,H1o,'-.k','LineWidth',2)
hold on
plot(newdata_long,H1d,'-.r','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'\epsilon_m=0.3','No Variation'},'Location','northeast')
grid on

tStart=datetime('2022-01-04'); tEnd=datetime('2022-04-30');
subplot(2,2,4)
plot(newdata_long,D1o,'-.k','LineWidth',2)
hold on
plot(newdata_long,D1d,'-.r','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfD','FontSize',12)
xlim([tStart tEnd])
legend({'\epsilon_m=0.3','No Variation'},'Location','northeast')
grid on

% Plot for Vaccinated Individuals 
figure
subplot(2,2,1)
plot(newdata_long,IBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,IBd,'-.k','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'\epsilon_m=0.3','No Variation'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,ABo,'-.r','LineWidth',2)
hold on
plot(newdata_long,ABd,'-.k','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'\epsilon_m=0.3','No Variation'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,HBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,HBd,'-.k','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'\epsilon_m=0.3','No Variation'},'Location','northeast')
grid on

tStart=datetime('2022-01-04'); tEnd=datetime('2022-04-30');
subplot(2,2,4)
plot(newdata_long,DBo,'-.r','LineWidth',2)
hold on
plot(newdata_long,DBd,'-.k','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
xlim([tStart tEnd])
text(0.5,0.5,'\bfD','FontSize',12)
legend({'\epsilon_m=0.3','No Variation'},'Location','northeast')
grid on

%Waning Immunity Vaccinated and Unvaccinted 
figure
subplot(1,2,1)
plot(newdata_long,Wo,'-.r','LineWidth',2)
hold on
plot(newdata_long,Wd,'-.k','LineWidth',2)
title('Waning Natural Immunity Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'\epsilon_m=0.3','No Variation'},'Location','northeast')
grid on

subplot(1,2,2)
plot(newdata_long,Wvo,'-.r','LineWidth',2)
hold on
plot(newdata_long,Wvd,'-.k','LineWidth',2)
title('Waning Acquired Immunity  Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
text(0.5,0.5,'\bfB','FontSize',12)
legend({'\epsilon_m=0.3','No Variation'},'Location','northeast')
grid on

%Deaths a partir del primer de abril
tStart=datetime('2022-04-01'); tEnd=datetime('2022-08-30');
figure
subplot(1,2,1)
plot(newdata_long,D1o,'-.r','LineWidth',2)
hold on 
plot(newdata_long,D1d,'-.k','LineWidth',2)
title('Death of the Unvaccinated Individuals')
ylabel('Individuals')
xlim([tStart tEnd])
text(1500,1500,'\bfA','FontSize',12)
legend({'\epsilon_m=0.3','No Variation'},'Location','northwest')
grid on

subplot(1,2,2)
plot(newdata_long,DBo,'-.r','LineWidth',2)
hold on 
plot(newdata_long,DBd,'-.k','LineWidth',2)
title('Death of the Vaccinated Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
text(0.5,0.5,'\bfB','FontSize',12)
legend({'\epsilon_m=0.3','No Variation'},'Location','northwest')
grid on

%% Lets vary the usage of face masks with 50% efficiency
epsilonm=0.5;
date1 = datenum('2022/04/01');
date2 = datenum('2022/10/01');
Cm = @(x) 0.5.*(x < date1-datenum(dateslong(1))) ...
       + 0.25.*( (x >= date1-datenum(dateslong(1)))&(x < date2-datenum(dateslong(1))) ) ...
       + 0.*( x >= date2-datenum(dateslong(1)) );
[t,S1o,E1o,I1o,A1o,H1o,V1o,V2o,EBo,IBo,ABo,HBo,Ro,Wo,Wvo,V3o,Do, ...
    I1cumulao,A1cumulao,IBcumulao,ABcumulao,D1o,DBo]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.5;
epsilonm=0.5;
[t,S1d,E1d,I1d,A1d,H1d,V1d,V2d,EBd,IBd,ABd,HBd,Rd,Wd,Wvd,V3d,Dd, ...
    I1cumulad,A1cumulad,IBcumulad,ABcumulad,D1d,DBd]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
figure
subplot(2,2,1)
plot(newdata_long,I1o,'-.m','LineWidth',2)
hold on 
plot(newdata_long,I1d,'-.k','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'\epsilon_m=0.5','No Variation'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,A1o,'-.m','LineWidth',2)
hold on 
plot(newdata_long,A1d,'-.k','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'\epsilon_m=0.5','No Variation'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,H1o,'-.m','LineWidth',2)
hold on
plot(newdata_long,H1d,'-.k','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'\epsilon_m=0.5','No Variation'},'Location','northeast')
grid on

tStart=datetime('2022-01-04'); tEnd=datetime('2022-04-30');
subplot(2,2,4)
plot(newdata_long,D1o,'-.m','LineWidth',2)
hold on
plot(newdata_long,D1d,'-.k','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfD','FontSize',12)
xlim([tStart tEnd])
legend({'\epsilon_m=0.5','No Variation'},'Location','northeast')
grid on

% Plot for Vaccinated Individuals 
figure
subplot(2,2,1)
plot(newdata_long,IBo,'-.m','LineWidth',2)
hold on
plot(newdata_long,IBd,'-.k','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'\epsilon_m=0.5','No Variation'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,ABo,'-.m','LineWidth',2)
hold on
plot(newdata_long,ABd,'-.k','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'\epsilon_m=0.5','No Variation'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,HBo,'-.m','LineWidth',2)
hold on
plot(newdata_long,HBd,'-.k','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'\epsilon_m=0.5','No Variation'},'Location','northeast')
grid on

tStart=datetime('2022-01-04'); tEnd=datetime('2022-04-30');
subplot(2,2,4)
plot(newdata_long,DBo,'-.m','LineWidth',2)
hold on
plot(newdata_long,DBd,'-.k','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
xlim([tStart tEnd])
text(0.5,0.5,'\bfD','FontSize',12)
legend({'\epsilon_m=0.5','No Variation'},'Location','northeast')
grid on

%Waning Immunity Vaccinated and Unvaccinted 
figure
subplot(1,2,1)
plot(newdata_long,Wo,'-.m','LineWidth',2)
hold on
plot(newdata_long,Wd,'-.k','LineWidth',2)
title('Waning Natural Immunity Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'\epsilon_m=0.5','No Variation'},'Location','northeast')
grid on

subplot(1,2,2)
plot(newdata_long,Wvo,'-.m','LineWidth',2)
hold on
plot(newdata_long,Wvd,'-.k','LineWidth',2)
title('Waning Acquired Immunity  Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
text(0.5,0.5,'\bfB','FontSize',12)
legend({'\epsilon_m=0.5','No Variation'},'Location','northeast')
grid on

%Deaths a partir del primer de abril
tStart=datetime('2022-04-01'); tEnd=datetime('2022-08-30');
figure
subplot(1,2,1)
plot(newdata_long,D1o,'-.m','LineWidth',2)
hold on 
plot(newdata_long,D1d,'-.k','LineWidth',2)
title('Death of the Unvaccinated Individuals')
ylabel('Individuals')
xlim([tStart tEnd])
text(1500,1500,'\bfA','FontSize',12)
legend({'\epsilon_m=0.5','No Variation'},'Location','northwest')
grid on

subplot(1,2,2)
plot(newdata_long,DBo,'-.m','LineWidth',2)
hold on 
plot(newdata_long,DBd,'-.k','LineWidth',2)
title('Death of the Vaccinated Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
text(0.5,0.5,'\bfB','FontSize',12)
legend({'\epsilon_m=0.5','No Variation'},'Location','northwest')
grid on

%% Lets vary the usage of face masks with 70% efficiency
epsilonm=0.7;
date1 = datenum('2022/04/01');
date2 = datenum('2022/10/01');
Cm = @(x) 0.5.*(x < date1-datenum(dateslong(1))) ...
       + 0.25.*( (x >= date1-datenum(dateslong(1)))&(x < date2-datenum(dateslong(1))) ) ...
       + 0.*( x >= date2-datenum(dateslong(1)) );
[t,S1o,E1o,I1o,A1o,H1o,V1o,V2o,EBo,IBo,ABo,HBo,Ro,Wo,Wvo,V3o,Do, ...
    I1cumulao,A1cumulao,IBcumulao,ABcumulao,D1o,DBo]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.5;
epsilonm=0.7;
[t,S1d,E1d,I1d,A1d,H1d,V1d,V2d,EBd,IBd,ABd,HBd,Rd,Wd,Wvd,V3d,Dd, ...
    I1cumulad,A1cumulad,IBcumulad,ABcumulad,D1d,DBd]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
figure
subplot(2,2,1)
plot(newdata_long,I1o,'-.b','LineWidth',2)
hold on 
plot(newdata_long,I1d,'-.k','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'\epsilon_m=0.7','No Variation'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,A1o,'-.b','LineWidth',2)
hold on 
plot(newdata_long,A1d,'-.k','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'\epsilon_m=0.7','No Variation'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,H1o,'-.b','LineWidth',2)
hold on
plot(newdata_long,H1d,'-.k','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'\epsilon_m=0.7','No Variation'},'Location','northeast')
grid on

tStart=datetime('2022-01-04'); tEnd=datetime('2022-04-30');
subplot(2,2,4)
plot(newdata_long,D1o,'-.b','LineWidth',2)
hold on
plot(newdata_long,D1d,'-.k','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfD','FontSize',12)
xlim([tStart tEnd])
legend({'\epsilon_m=0.7','No Variation'},'Location','northeast')
grid on

% Plot for Vaccinated Individuals 
figure
subplot(2,2,1)
plot(newdata_long,IBo,'-.b','LineWidth',2)
hold on
plot(newdata_long,IBd,'-.k','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'\epsilon_m=0.7','No Variation'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,ABo,'-.b','LineWidth',2)
hold on
plot(newdata_long,ABd,'-.k','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'\epsilon_m=0.7','No Variation'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,HBo,'-.b','LineWidth',2)
hold on
plot(newdata_long,HBd,'-.k','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'\epsilon_m=0.7','No Variation'},'Location','northeast')
grid on

tStart=datetime('2022-01-04'); tEnd=datetime('2022-04-30');
subplot(2,2,4)
plot(newdata_long,DBo,'-.b','LineWidth',2)
hold on
plot(newdata_long,DBd,'-.k','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
xlim([tStart tEnd])
text(0.5,0.5,'\bfD','FontSize',12)
legend({'\epsilon_m=0.7','No Variation'},'Location','northeast')
grid on

%Waning Immunity Vaccinated and Unvaccinted 
figure
subplot(1,2,1)
plot(newdata_long,Wo,'-.b','LineWidth',2)
hold on
plot(newdata_long,Wd,'-.k','LineWidth',2)
title('Waning Natural Immunity Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'\epsilon_m=0.7','No Variation'},'Location','northeast')
grid on

subplot(1,2,2)
plot(newdata_long,Wvo,'-.b','LineWidth',2)
hold on
plot(newdata_long,Wvd,'-.k','LineWidth',2)
title('Waning Acquired Immunity  Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
text(0.5,0.5,'\bfB','FontSize',12)
legend({'\epsilon_m=0.7','No Variation'},'Location','northeast')
grid on

%Deaths a partir del primer de abril

tStart=datetime('2022-04-01'); tEnd=datetime('2022-08-30');
figure
subplot(1,2,1)
plot(newdata_long,D1o,'-.b','LineWidth',2)
hold on 
plot(newdata_long,D1d,'-.k','LineWidth',2)
title('Death of the Unvaccinated Individuals')
ylabel('Individuals')
xlim([tStart tEnd])
text(1500,1500,'\bfA','FontSize',12)
legend({'\epsilon_m=0.7','No Variation'},'Location','northwest')
grid on

subplot(1,2,2)
plot(newdata_long,DBo,'-.b','LineWidth',2)
hold on 
plot(newdata_long,DBd,'-.k','LineWidth',2)
title('Death of the Vaccinated Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
text(0.5,0.5,'\bfB','FontSize',12)
legend({'\epsilon_m=0.7','No Variation'},'Location','northwest')
grid on

%% Lets vary the usage of face masks with 50% efficiency
epsilonm=0.95;
date1 = datenum('2022/04/01');
date2 = datenum('2022/10/01');
Cm = @(x) 0.5.*(x < date1-datenum(dateslong(1))) ...
       + 0.25.*( (x >= date1-datenum(dateslong(1)))&(x < date2-datenum(dateslong(1))) ) ...
       + 0.*( x >= date2-datenum(dateslong(1)) );
[t,S1o,E1o,I1o,A1o,H1o,V1o,V2o,EBo,IBo,ABo,HBo,Ro,Wo,Wvo,V3o,Do, ...
    I1cumulao,A1cumulao,IBcumulao,ABcumulao,D1o,DBo]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
Cm=0.5;
epsilonm=0.95;
[t,S1d,E1d,I1d,A1d,H1d,V1d,V2d,EBd,IBd,ABd,HBd,Rd,Wd,Wvd,V3d,Dd, ...
    I1cumulad,A1cumulad,IBcumulad,ABcumulad,D1d,DBd]= ...
    USA_covid_vaccine_solver(rho1,rho2,rho3,rho4, t, IniVal);
figure
subplot(2,2,1)
plot(newdata_long,I1o,'-.c','LineWidth',2)
hold on 
plot(newdata_long,I1d,'-.k','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'\epsilon_m=0.95','No Variation'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,A1o,'-.c','LineWidth',2)
hold on 
plot(newdata_long,A1d,'-.k','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'\epsilon_m=0.95','No Variation'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,H1o,'-.c','LineWidth',2)
hold on
plot(newdata_long,H1d,'-.k','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'\epsilon_m=0.95','No Variation'},'Location','northeast')
grid on

tStart=datetime('2022-01-04'); tEnd=datetime('2022-04-30');
subplot(2,2,4)
plot(newdata_long,D1o,'-.c','LineWidth',2)
hold on
plot(newdata_long,D1d,'-.k','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfD','FontSize',12)
xlim([tStart tEnd])
legend({'\epsilon_m=0.95','No Variation'},'Location','northeast')
grid on

% Plot for Vaccinated Individuals 
figure
subplot(2,2,1)
plot(newdata_long,IBo,'-.c','LineWidth',2)
hold on
plot(newdata_long,IBd,'-.k','LineWidth',2)
title('Infected Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
legend({'\epsilon_m=0.95','No Variation'},'Location','northeast')
grid on

subplot(2,2,2)
plot(newdata_long,ABo,'-.c','LineWidth',2)
hold on
plot(newdata_long,ABd,'-.k','LineWidth',2)
title('Asymptomatic Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfB','FontSize',12)
legend({'\epsilon_m=0.95','No Variation'},'Location','northeast')
grid on

subplot(2,2,3)
plot(newdata_long,HBo,'-.c','LineWidth',2)
hold on
plot(newdata_long,HBd,'-.k','LineWidth',2)
title('Hospitalized Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'\epsilon_m=0.95','No Variation'},'Location','northeast')
grid on

subplot(2,2,4)
plot(newdata_long,DBo,'-.c','LineWidth',2)
hold on
plot(newdata_long,DBd,'-.k','LineWidth',2)
title('Death Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfD','FontSize',12)
legend({'\epsilon_m=0.95','No Variation'},'Location','northeast')
grid on

%Waning Immunity Vaccinated and Unvaccinted 
figure
subplot(1,2,1)
plot(newdata_long,Wo,'-.c','LineWidth',2)
hold on
plot(newdata_long,Wd,'-.k','LineWidth',2)
title('Waning Natural Immunity Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfC','FontSize',12)
text(0.5,0.5,'\bfD','FontSize',12)
legend({'\epsilon_m=0.95','No Variation'},'Location','northeast')
grid on

subplot(1,2,2)
plot(newdata_long,Wvo,'-.c','LineWidth',2)
hold on
plot(newdata_long,Wvd,'-.k','LineWidth',2)
title('Waning Acquired Immunity  Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
text(0.5,0.5,'\bfB','FontSize',12)
legend({'\epsilon_m=0.95','No Variation'},'Location','northeast')
grid on

%Deaths a partir del primer de abril
tStart=datetime('2022-04-01'); tEnd=datetime('2022-08-30');
figure
subplot(1,2,1)
plot(newdata_long,D1o,'-.c','LineWidth',2)
hold on 
plot(newdata_long,D1d,'-.k','LineWidth',2)
title('Death of the Unvaccinated Individuals')
ylabel('Individuals')
xlim([tStart tEnd])
text(1500,1500,'\bfA','FontSize',12)
legend({'\epsilon_m=0.95','No Variation'},'Location','northwest')
grid on

subplot(1,2,2)
plot(newdata_long,DBo,'-.c','LineWidth',2)
hold on 
plot(newdata_long,DBd,'-.k','LineWidth',2)
title('Death of the Vaccinated Individuals')
ylabel('Individuals')
text(0.5,0.5,'\bfA','FontSize',12)
text(0.5,0.5,'\bfB','FontSize',12)
legend({'\epsilon_m=0.95','No Variation'},'Location','northwest')
grid on