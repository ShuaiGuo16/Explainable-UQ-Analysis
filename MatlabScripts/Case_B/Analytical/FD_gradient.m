% FD way to calculate gradient
clear
clc

% Basic parameters
l_Plen = 0.17;   c_up = 343.1143;  area_ratio_1 = 29.76;
l_Swir = 0.18;   area_ratio_2 = 0.13;  l_com = 0.6;
c_down = 880.61;  temp_jump = 5.59;   imp_ratio = 2.57;
r_out = -0.6351;

syms omega F
syms freq growth real

%% Construct upstream acoustic matrix

% Plenum & 1st area jump
T_upstream(1,:) = 0.5*[(1+area_ratio_1)*exp(-omega*1i*l_Plen/c_up),(1-area_ratio_1)*exp(omega*1i*l_Plen/c_up)];
T_upstream(2,:) = 0.5*[(1-area_ratio_1)*exp(-omega*1i*l_Plen/c_up),(1+area_ratio_1)*exp(omega*1i*l_Plen/c_up)];

% Swirler tube
T_swirler = [exp(-omega*1i*l_Swir/c_up),0;0,exp(omega*1i*l_Swir/c_up)];

% Flame & 2st area jump
T_flame(1,:) = [0.5*(imp_ratio+area_ratio_2+area_ratio_2*temp_jump*F),0.5*(imp_ratio-area_ratio_2-area_ratio_2*temp_jump*F)];
T_flame(2,:) = [0.5*(imp_ratio-area_ratio_2-area_ratio_2*temp_jump*F),0.5*(imp_ratio+area_ratio_2+area_ratio_2*temp_jump*F)];

% Combustor
T_comb = [exp(-omega*1i*l_com/c_down),0;0,exp(omega*1i*l_com/c_down)];

% Full system without B.C.
T = T_comb*T_flame*T_swirler*T_upstream;

%% Construct characteristic equation
CharEqn = T(2,2)-r_out*T(1,2)+T(2,1)-r_out*T(1,1);

H = solve(CharEqn==0,F);
H_subs = subs(H,omega,freq-growth*1i);
a_ref = double(subs(subs(real(H_subs),freq,97.46*2*pi),growth,-4));
b_ref = double(subs(subs(imag(H_subs),freq,97.46*2*pi),growth,-4));

a_increase = double(subs(subs(real(H_subs),freq,98.46*2*pi),growth,-4));
b_increase = double(subs(subs(imag(H_subs),freq,98.46*2*pi),growth,-4));
a_omega = (a_increase-a_ref)/(2*pi)
b_omega = (b_increase-b_ref)/(2*pi)