    %parameter vector contains 15 elements and States contains 17
    %Define constant
    Vm = -87.491;%membrane potential (mV)
    F = 96487;%Faraday constant (C/mol)
    R = 8314;%Gas constant (J/kmol/K)
    T = 310;%Temperature (K)
    Ca_in = 1.3394*0.01;%calcium concentration in myoplasm (mmol/L)
    K_out = 5.4;%Extracellular K+ concentration (mmol/L)
    K_in = 1.4562*100;%K+ concentration in myoplasm (mmol/L)
    Na_out = 140;%Extracellular Na+ concentration (mmol/L)
    Na_in = 6.8909;%Na+ concentration in myoplasm (mmol/L)
    P_Na_K = 0.01833;%permeability
    Gbar = 0.19561*(1+0.6/(1+(3.8*0.00001/Ca_in)));
    Eks = R*T*log((K_out+P_Na_K*Na_out)/(K_in+P_Na_K*Na_in))/F;
    %define initial values
    States = [1 zeros(1,16)];
    %Para = [1.4864*0.01, 2.9877*0.01, 8.3986*0.01, -5.5461*0.01, 1.4601*0.01, 2.4465*0.1, 8.9538*0.02, 7.7320*0.02, -6.4726*0.01, 3.1173*0.001, -4.2625*0.1, 7.9405*0.1, -8.0174*0.01, 5.8638*0.1, 2.8206*0.1];
    Para = [0.1, 0.01, 0.001, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001, 0.01, 0.001, 0.0001, 0.1, 0.1];
    %define protocol
    v = [-70:10:90]';
    vh = -60;
    
    %[s1, s2] = ode23s(iksrates,[0 10000], States, [], vh);
    [s1, s2] = ode23s(@iksode,[0 100000], States, [], vh);
    sp=s2(length(s2),:);
    
    for i=1:length(v)
        [s1, s2]=ode23s(@iksode, [0:3000], sp, [], v(i));
        Iact(:,i)=Gbar*sum(s2(:,16:17),2)*(v(i)-Eks);
        [s3, s4]=ode23s(@iksode, [0:1000], s2(length(s2),:),[],vh);
        Itail(:,i)=Gbar*sum(s4(:,16:17),2)*(vh-Eks);
    end
    
Isim(:,1)=[s1;s3+max(s1)];
Isim(:,2:length(v)+1)=[Iact;Itail]/max(max([Iact;Itail]));

fileID = fopen('average_A341V_het_cAMP.txt','r');
formatSpec = '%f %f';
RD = textscan(fileID,formatSpec);
Time = RD{1}*1000;
RC = RD{2};
RC = RC/max(RC);

%plot
figure                                 
plot(Isim(:,1),Isim(:,2:length(v)+1))
title('IKs IVcurve')
xlabel('Time (ms)')
ylabel('Current')
% plot(RC)
% title('Normalized Reference Current')
% xlabel('Time (ms)')
% ylabel('Current')