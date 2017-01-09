function [y]=iksode(t,States,v)
Para = [0.001, 0.01, 0.001, -0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001, 0.01, 1, 1, 0.01, -0.001];
F = 96487;%Faraday constant (C/mol)
R = 8314;%Gas constant (J/kmol/K)
T = 310;%Temperature (K)

alphaA = Para(1);
alphaB = Para(2);
betaA = Para(3);
betaB = Para(4);
gammaA = Para(5);
gammaB = Para(6);
thetaA = Para(7);
etaA = Para(8);
etaB = Para(9);
deltaA = Para(10);
deltaB = Para(11);
omegaA = Para(12);
omgeaB = Para(13);
psiA = Para(14);
psiB = Para(15);
    
% alpha=alphaA*exp(v*F/(R*T)*alphaB);
% beta=betaA*exp(v*F/(R*T)*betaB);
% gamma=gammaA*exp(v*F/(R*T)*gammaB);
% delta=deltaA*exp(v*F/(R*T)*deltaB);
% theta=thetaA;
% eta=etaA*exp(v*F/(R*T)*etaB);
% psi=psiA*exp(v*F/(R*T)*psiB);
% omega=omegaA*exp(v*F/(R*T)*omgeaB);

alpha=3.72e-003*exp(v*F/(R*T)*2.10e-001);
beta=2.35e-004*exp(v*F/(R*T)*-2.42e-001);
gamma=7.25e-003*exp(v*F/(R*T)*2.43e+000);
delta=1.53e-003*exp(v*F/(R*T)*-6.26e-001);
theta=1.96e-003;
eta=1.67e-002*exp(v*F/(R*T)*-1.34e+000);
psi=3.41e-004*exp(v*F/(R*T)*1.24e+000);
omega=2.62e-004*exp(v*F/(R*T)*-8.07e-001);

c1 = States(1);
c2 = States(2);
c3 = States(3);
c4 = States(4);
c5 = States(5);
c6 = States(6);
c7 = States(7);
c8 = States(8);
c9 = States(9);
c10 = States(10);
c11 = States(11);
c12 = States(12);
c13 = States(13);
c14 = States(14);
c15 = States(15);
o1 = States(16);
o2 = States(17);

do2dt = psi*o1-omega*o2;
do1dt = theta*c15+omega*o2-(psi+eta)*o1;
dc15dt = gamma*c14+eta*o1-(4*delta+theta)*c15;
dc14dt = alpha*c13+4*delta*c15+2*gamma*c12-(beta+3*delta+gamma)*c14;
dc13dt = beta*c14+gamma*c11-(alpha+3*delta)*c13;
dc12dt = alpha*c11+3*delta*c14+3*gamma*c9-(2*beta+2*delta+2*gamma)*c12;
dc11dt = 2*alpha*c10+2*beta*c12+2*gamma*c8+3*delta*c13-(beta+alpha+gamma+2*delta)*c11;
dc10dt = beta*c11+gamma*c7-(2*alpha+2*delta)*c10;
dc9dt = alpha*c8+2*delta*c12+4*gamma*c5-(3*beta+delta+3*gamma)*c9;
dc8dt = 2*alpha*c7+3*beta*c9+3*gamma*c4+2*delta*c11-(2*beta+alpha+2*gamma+delta)*c8;
dc7dt = 3*alpha*c6+2*beta*c8+2*gamma*c3+2*delta*c10-(beta+2*alpha+delta+gamma)*c7;
dc6dt = beta*c7+gamma*c2-(3*alpha+delta)*c6;
dc5dt = alpha*c4+delta*c9-(4*beta+4*gamma)*c5;
dc4dt = 2*alpha*c3+4*beta*c5+delta*c8-(3*beta+alpha+3*gamma)*c4;
dc3dt = 3*alpha*c2+3*beta*c4+delta*c7-(2*beta+2*alpha+2*gamma)*c3;
dc2dt = 4*alpha*c1+2*beta*c3+delta*c6-(beta+3*alpha+gamma)*c2;
dc1dt = beta*c2-4*alpha*c1;

y(1,1)=dc1dt;
y(2,1)=dc2dt;
y(3,1)=dc3dt;
y(4,1)=dc4dt;
y(5,1)=dc5dt;
y(6,1)=dc6dt;
y(7,1)=dc7dt;
y(8,1)=dc8dt;
y(9,1)=dc9dt;
y(10,1)=dc10dt;
y(11,1)=dc11dt;
y(12,1)=dc12dt;
y(13,1)=dc13dt;
y(14,1)=dc14dt;
y(15,1)=dc15dt;
y(16,1)=do1dt;
y(17,1)=do2dt;

return

if abs(sum(y))>1e-6
    error('sum of derivatives greater than 0')
end