%% Program Heat exhanger 
% program to test Epsilon - NTU method.
% assuming a counterflow heat exhanger with known properties, known inlet
% temperatures and mass flowrates.
clc;
clear all;
close all;
format longg;
%% Input Variable
value1=input("Enter Specific Heat Capacity Of Hot Fluid:"); %1.2KJ/kg-K
c_p_hot=input(value1); % KJ/kg-K
value2=input("Enter Mass Flow Rate Of The Hot Fluid:"); %2.222kg/s
m_rate_hot=input(value2); % kg/s
value3=input("Enter Temperature Of Hot Fluid:"); %378K
T_hot_in=input(value3); % K
value4=input("Enter Specific Heat Capacity Of Cold Fluid:"); %4.18KJ/kg-K
c_p_cold=input(value4); % KJ/kg-K 
value5=input("Enter Mass Flow Rate Of Cold Fluid:"); %2.08333kg/s
m_rate_cold=input(value5); %kg/s
value6=input("Enter Temperature Of Cold Fluid:"); %288K
T_cold_in=input(value6);
value7=input("Enter Overall Heat Transfer Coefficent:"); %0.145kW/m2-K
U=input(value7); %kW/m2.K (Overall heat transfer coefficient)
value8=input("Enter Area Of Cross Section:"); %20m2
A=input(value8);  %m2 
value9=input("Enter Flow Type:","s"); %Counter Flow
HE_Type=input(value9);
%% calculation related to Heat Exchanger
C_hot = m_rate_hot*c_p_hot; %{heat capacity of hot fluid}
C_cold = m_rate_cold*c_p_cold; %{heat capacity of cold fluid}
C_min = min(C_hot,C_cold); % finds the flow with lower heat capacity and higher temperature change.
C_max = max(C_hot,C_cold); % finds the flow with higher heat capacity and lower temperature change. 
C_r=C_min/C_max; %[capacity rate ratio of the heat exchanger]
NTU = U*A/C_min;
%% Effectiveness calculation
% Special case of boiling or condensing:
if C_r==0
    epsilon=1-exp(-NTU);
    return;
end
%general cases please refer the TABLE from NTU method 
switch HE_Type
    case 'Parallel Flow'
        epsilon=(1-exp(-NTU*(1+C_r)))/(1+C_r);
    case 'Counter Flow'
        if C_r==1
            epsilon=NTU/(1+NTU);
        else
            epsilon=(1-exp(-NTU*(1-C_r)))/(1-C_r*exp(-NTU*(1-C_r)));
        end
    case 'One Shell Pass'
        epsilon=2/(1+C_r+sqrt(1+C_r^2)*(1+exp(-NTU*sqrt(1+C_r^2)))/(1+exp(-NTU*sqrt(1+C_r^2))));
    case 'N Shell Pass'
        NTUN=NTU/N;
        epsilon1=2/(1+C_r+sqrt(1+C_r^2)*(1+exp(-NTUN*sqrt(1+C_r^2)))/(1+exp(-NTUN*sqrt(1+C_r^2))));
        epsilon=(((1-epsilon1*C_r)/(1-epsilon1))^N-1) / (((1-epsilon1*C_r)/(1-epsilon1))^N-C_r);
    case 'Cross Both Unmixed'
        epsilon=1-exp(1/C_r * NTU^0.22 * (exp(-C_r*NTU^0.78)-1));
    case 'Cross Cmax Mixed'
        epsilon=1/C_r*(1-exp(-C_r*(1-exp(-NTU))));
    case 'Cross Cmin Mixed'
        epsilon=1-epx(-1/C_r*(1-exp(-C_r*NTU)));
    otherwise % the type is not in the list, therefore we assume there's no heat exchanger.
        epsilon=0; 
end
%% Final output of the results 
Q_max=C_min*(T_hot_in-T_cold_in);
Q=epsilon*Q_max ;
T_hot_out=T_hot_in-Q/C_hot; 
T_cold_out=T_cold_in+Q/C_cold; 
Q_hot=m_rate_hot*c_p_hot*(T_hot_in-T_hot_out);
Q_cold=m_rate_cold*c_p_cold*(T_cold_in-T_cold_out);
LMTD=((T_hot_in-T_cold_out)-(T_hot_out-T_cold_in))/log((T_hot_in-T_cold_out)/(T_hot_out-T_cold_in));
Q_exchange=U*A*LMTD;
Effectiveness_Heat_exchanger=epsilon;