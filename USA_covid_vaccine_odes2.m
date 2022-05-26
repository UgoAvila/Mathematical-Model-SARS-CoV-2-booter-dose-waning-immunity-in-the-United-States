function dydt = USA_covid_vaccine_odes(t,y,rho1,rho2,rho3,rho4)

% We define fixed parameter values:

Lambda = 3747540/365; % Birth rate (people/day)     Obtained from https://www.cdc.gov/nchs/fastats/births.htm
d = (869.7/100000)/365; % Natural death rate (1/day)    Obtained from https://www.cdc.gov/nchs/fastats/deaths.htm
w = 1/4.;     % Transfer from exposed to infectious (1/day)
Kn = 1/270;     % Waning rate of natural immunity (1/day)
p1 = 1/9;   % Proportion of symptomatic infections
gamma1 = 1/14;  % Recovery rate (1/day)

global alpha1 chi omega beta1 beta2 Cm epsilonm kappa epsa1 epsa2 epsa3 ...
    SumEpsLi SumEpsL2i Sumthetai Kv Kv1 Kv2 eps3i p2 gamma2 alphaB chiB omegaB Kn Cr

S1=y(1);
E1=y(2);
I1=y(3);
A1=y(4);
H1=y(5);
V1=y(6);
V2=y(7);
EB=y(8);
IB=y(9);
AB=y(10);
HB=y(11);
R =y(12);
W =y(13);
Wv=y(14);
V3=y(15);
D =y(16);
I1cumul=y(17);
A1cumul=y(18);
IBcumul=y(19);
ABcumul=y(20);
D1=y(21);
DB=y(22);

lambd = @(x) beta1(x)*(I1 + kappa*IB) + beta2(x)*(A1 + kappa*AB);   % Transmission rate at time x
SumEpsai = @(x) (1-epsa1)*rho1(x) + (1-epsa2)*rho2(x) + (1-epsa3)*rho3(x);  %Sum of vaccin. rates at time x
% SumEpsLi = (1-epsL1) + (1-epsL2) + (1-epsL3);
% SumEpsL2i = (1-epsL21) + (1-epsL22) + (1-epsL23);

dS1_dt = Lambda + 2*Kv*Wv + 2*Kn*W - lambd(t).*(1-Cm*epsilonm).*(1-Cr).*S1 ...
    - SumEpsai(t) - d*S1;
dE1_dt = lambd(t).*(1-Cm*epsilonm).*(1-Cr).*S1 - w*E1 - d*E1;
dI1_dt = p1*w*E1 - (alpha1+gamma1)*I1 - d*I1;
dA1_dt = (1-p1)*w*E1 - gamma1*A1 - d*A1;
dH1_dt = alpha1*I1 - omega*H1 - d*H1;
dV1_dt = SumEpsai(t) - SumEpsLi.*lambd(t).*(1-Cm*epsilonm).*(1-Cr).*V1 ...
    - Sumthetai.*V1 - 2*Kv1*V1 - d*V1;
dV2_dt = Sumthetai.*V1 - SumEpsL2i.*lambd(t).*(1-Cm*epsilonm).*(1-Cr).*V2 ...
    - 2*Kv2*V2 - d*V2;
dEB_dt = SumEpsLi.*lambd(t).*(1-Cm*epsilonm).*(1-Cr).*V1 + SumEpsL2i.*lambd(t).*(1-Cm*epsilonm).*(1-Cr).*V2 ...
    + lambd(t).*(1-Cm*epsilonm).*(1-Cr).*eps3i.*V3 - w*EB - d*EB;
dIB_dt = p2*w*EB - (alphaB+gamma2)*IB - d*IB;
dAB_dt = (1-p2)*w*EB - gamma2*AB - d*AB;
dHB_dt = alphaB*IB - omegaB*HB - d*HB;
dR_dt = gamma1*(I1 + A1) + gamma2*(IB + AB) + omega*(1-chi)*H1...
    + omegaB*(1-chiB)*HB - 2*Kn*R - d*R;
dW_dt = 2*Kn*R - 2*Kn*W - d*W;
dWv_dt = 2*Kv1*V1 + 2*Kv2*V2 - 2*Kv*Wv - rho4(t)*Wv - d*Wv;
dV3_dt = rho4(t)*Wv - lambd(t).*(1-Cm*epsilonm).*(1-Cr)*eps3i.*V3 - d*V3;
dD_dt = omega*chi*H1 + omegaB*chiB*HB;

dI1cumul_dt = p1*w*E1;  % Cumulative symptomatic infections (unvaccinated)
dA1cumul_dt = (1-p1)*w*E1;  % Cumulative asymptomatic infections (unvaccinated)
dIBcumul_dt = p2*w*EB;  % Cumulative symptomatic infections (breakthrough)
dABcumul_dt = (1-p2)*w*EB;  % Cumulative asymptomatic infections (breakthrough)
dD1_dt=omega*chi*H1;
dDB_dt=omegaB*chiB*HB;


dydt = [dS1_dt; dE1_dt; dI1_dt; dA1_dt; dH1_dt; dV1_dt; dV2_dt; dEB_dt; ...
    dIB_dt; dAB_dt; dHB_dt; dR_dt; dW_dt; dWv_dt; dV3_dt; dD_dt; ...
    dI1cumul_dt; dA1cumul_dt; dIBcumul_dt; dABcumul_dt; dD1_dt; dDB_dt];
end