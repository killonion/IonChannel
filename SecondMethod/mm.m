function [I, o] = mm(Vm, Para, States, Time)
    %parameter vector contains 15 elements and States contains 17
    %Define constant
    F = 96487;%Faraday constant (C/mol)
    R = 8314;%Gas constant (J/kmol/K)
    T = 310;%Temperature (K)
    Ca_in = 1.3394*0.01;%calcium concentration in myoplasm (mmol/L)
    K_out = 5.4;%Extracellular K+ concentration (mmol/L)
    K_in = 1.4562*100;%K+ concentration in myoplasm (mmol/L)
    Na_out = 140;%Extracellular Na+ concentration (mmol/L)
    Na_in = 6.8909;%Na+ concentration in myoplasm (mmol/L)
    P_Na_K = 0.01833;%permeability
    Gbar = 0.19561*(1+0.6/(1+(3.8*0.000001/Ca_in)));
    Eks = R*T*log((K_out+P_Na_K*Na_out)/(K_in+P_Na_K*Na_in))/F;
    t = length(Vm);
    %parameters
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
    %States
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
    %calculate transition states
    alpha = alphaA*exp(alphaB*Vm*F/(R*T));
    beta = betaA*exp(betaB*Vm*F/(R*T));
    gamma = gammaA*exp(gammaB*Vm*F/(R*T));
    theta = thetaA;
    eta = etaA*exp(etaB*Vm*F/(R*T));
    delta = deltaA*exp(deltaB*Vm*F/(R*T));
    omega = omegaA*exp(omgeaB*Vm*F/(R*T));
    psi = psiA*exp(psiB*Vm*F/(R*T));
    %loop through all the parameters and compute the corresponding simulated
    %current
        TR = [-4*alpha, beta, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
            4*alpha, -(beta+3*alpha+gamma), 2*beta, 0, 0, delta, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
            0, 3*alpha, -(2*beta+2*alpha+2*gamma), 3*beta, 0, 0, delta, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 2*alpha, -(3*beta+alpha+3*gamma), 4*beta, 0, 0, delta, 0, 0, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 0, alpha, -(4*beta+4*gamma), 0, 0, 0, delta, 0, 0, 0, 0, 0, 0, 0, 0;
            0, gamma, 0, 0, 0, -(delta+3*alpha), beta, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 2*gamma, 0, 0, 3*alpha, -(delta+beta+2*alpha+gamma), 2*beta, 0, 2*delta, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 0, 3*gamma, 0, 0, 2*alpha, -(delta+2*beta+alpha+2*gamma), 3*beta, 0, 2*delta, 0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 4*gamma, 0, 0, alpha, -(delta+3*beta+3*gamma), 0, 0, 2*delta, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0, gamma, 0, 0, -(2*delta+2*alpha), beta, 0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0, 0, 2*gamma, 0, 2*alpha, -(2*delta+beta+alpha+gamma), 2*beta, 3*delta, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0, 0, 0, 3*gamma, 0, alpha, 0, 0, 3*delta, 0, 0, 0;
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, gamma, 0, -(3*delta+alpha), beta, 0, 0, 0;
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2*gamma, alpha, -(3*delta+alpha+gamma), 4*delta, 0, 0;
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, gamma, -(4*delta+theta), eta, 0;
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, theta, -(eta+psi), omega;
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, psi, -omega];
        %solve the system
        for i = 1:length(Time);
            [P,D]=eig(TR);
            etD=exp(D*i);
            etTR =P*etD*inv(P);
            StatesX = etTR*transpose(States);
            o1 = StatesX(16);
            o2 = StatesX(17);
            %to compute the simulated current
            Iks(i) = Gbar*(o1+o2)*(Vm-Eks);
            o(i) = o1+o2;
        end
    I = Iks;
end