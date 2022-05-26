function [t,S1,E1,I1,A1,H1,V1,V2,EB,IB,AB,HB,R,W,Wv,V3,D,I1cumul,A1cumul,IBcumul,ABcumul,D1,DB]...
    =USA_covid_vaccine_solver(rho1,rho2,rho3,rho4,t,IV)
% This function computes the solutions of the model with vaccination rates
% (rho1,rho2,rho3,rho4) and initial values IV over the time range t

%% Solving the model
options=odeset('NonNegative',(1:20));
[t,y] = ode45(@(t,y)USA_covid_vaccine_odes2(t,y,rho1,rho2,rho3,rho4), ...
    t,IV,options);

S1=y(:,1);
E1=y(:,2);
I1=y(:,3);
A1=y(:,4);
H1=y(:,5);
V1=y(:,6);
V2=y(:,7);
EB=y(:,8);
IB=y(:,9);
AB=y(:,10);
HB=y(:,11);
R =y(:,12);
W =y(:,13);
Wv=y(:,14);
V3=y(:,15);
D =y(:,16);
I1cumul=y(:,17);
A1cumul=y(:,18);
IBcumul=y(:,19);
ABcumul=y(:,20);
D1=y(:,21);
DB=y(:,22);


end