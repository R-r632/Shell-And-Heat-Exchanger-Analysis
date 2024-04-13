%% Program Heat exhanger 
% program to test Epsilon - NTU method.
% assuming a counterflow heat exhanger with known properties, known inlet
% temperatures and mass flowrates.
clc;
clear all;
close all;
format longg;
%% Input Variable
c_p_hot = 1.2; % KJ/kg-K
m_rate_hot =2.222; % kg/s
T_hot_in = 378; % K
c_p_cold = 4.18; % KJ/kg-K 
m_rate_cold = 2.08333; %kg/s
T_cold_in = 288;
U=0.145; %kW/m2.K (Overall heat transfer coefficient)
A=20;  %m2 
HE_Type = 'Counter Flow'
%% calculation related to Heat Exchanger
C_hot = m_rate_hot*c_p_hot; %{heat capacity of hot fluid}
C_cold = m_rate_cold*c_p_cold; %{heat capacity of cold fluid}
C_min = min(C_hot,C_cold); % finds the flow with lower heat capacity and higher temperature change.
C_max = max(C_hot,C_cold); % finds the flow with higher heat capacity and lower temperature change. 
C_r=C_min/C_max; %[capacity rate ratio of the heat exchanger]
NTU = U*A/C_min
%% Effectiveness calculation
% Special case of boiling or condensing:
if C_r == 0
    epsilon = 1-exp(-NTU);
    return;
end
%general cases please refer the TABLE from NTU method 
switch HE_Type
    case 'Parallel Flow'
        epsilon = (1-exp(-NTU*(1+C_r)))/(1+C_r);
    case 'Counter Flow'
        if C_r==1
            epsilon = NTU/(1+NTU);
        else
            epsilon = (1-exp(-NTU*(1-C_r)))/(1-C_r*exp(-NTU*(1-C_r)));
        end
    case 'One Shell Pass'
        epsilon = 2/(1+C_r+sqrt(1+C_r^2)*(1+exp(-NTU*sqrt(1+C_r^2)))/(1+exp(-NTU*sqrt(1+C_r^2))));
    case 'N Shell Pass'
        NTUN = NTU/N;
        epsilon1 = 2/(1+C_r+sqrt(1+C_r^2)*(1+exp(-NTUN*sqrt(1+C_r^2)))/(1+exp(-NTUN*sqrt(1+C_r^2))));
        epsilon = (((1-epsilon1*C_r)/(1-epsilon1))^N-1) / (((1-epsilon1*C_r)/(1-epsilon1))^N-C_r);
    case 'Cross Both Unmixed'
        epsilon = 1-exp(1/C_r * NTU^0.22 * (exp(-C_r*NTU^0.78)-1));
    case 'Cross Cmax Mixed'
        epsilon = 1/C_r*(1-exp(-C_r*(1-exp(-NTU))));
    case 'Cross Cmin Mixed'
        epsilon = 1 - epx(-1/C_r*(1-exp(-C_r*NTU)));
    otherwise % the type is not in the list, therefore we assume there's no heat exchanger.
        epsilon = 0; 
end
%% Final output of the results 
Q_max = C_min*(T_hot_in-T_cold_in);
Q = epsilon * Q_max ;
T_hot_out = T_hot_in - Q/C_hot 
T_cold_out = T_cold_in + Q/C_cold 
Q_hot=m_rate_hot*c_p_hot*(T_hot_in-T_hot_out)
Q_cold=m_rate_cold*c_p_cold*(T_cold_in-T_cold_out)
LMTD = ((T_hot_in-T_cold_out)-(T_hot_out-T_cold_in))/log((T_hot_in-T_cold_out)/(T_hot_out-T_cold_in))
Q_exchange = U*A*LMTD
Effectiveness_Heat_exchanger = epsilon