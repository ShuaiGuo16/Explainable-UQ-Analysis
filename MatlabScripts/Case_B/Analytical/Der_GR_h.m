clear
clc

load 'model_Luis_30kW_9.5%_4_21_350ms.mat'
load 'eigenvector_70.mat'

temp_mean = model.B(4:19);
temp_cov = getcov(model);   
mean_translation = temp_mean;   
variance_translation = sqrt(diag(temp_cov)); 
% Results from active subspace
a = V(:,end)./variance_translation;
ratio_AS = a(5)/a(4)

delta_Ts = 6.25e-4;
FR_ref = 97.46*2*pi;
[ a_omega_value, b_omega_value ] = Direction_Cal_B( FR_ref, -4 );

A = 0;   B = 0;
for i = 1:16
    A = A + mean_translation(i)*(i+2)*delta_Ts*exp(-(i+2)*delta_Ts*(-4))*cos((i+2)*delta_Ts*FR_ref);
    B = B + mean_translation(i)*(i+2)*delta_Ts*exp(-(i+2)*delta_Ts*(-4))*sin((i+2)*delta_Ts*FR_ref);
end

GR_h4 = cos(6*FR_ref*delta_Ts)*(b_omega_value+A)+sin(6*FR_ref*delta_Ts)*(a_omega_value+B);
GR_h5 = cos(7*FR_ref*delta_Ts)*(b_omega_value+A)+sin(7*FR_ref*delta_Ts)*(a_omega_value+B);
ratio_ana = GR_h5/GR_h4



%% Approximate through FD
h4_per = model.B(4:19);
h4_per(4) = h4_per(4)*1.01;
GR_4 = solver_GR(h4_per);
diff_h4 = (GR_4+4)/(0.01*temp_mean(4));

h5_per = model.B(4:19);
h5_per(5) = h5_per(5)*1.01;
GR_5 = solver_GR(h5_per);
diff_h5 = (GR_5+4)/(0.01*temp_mean(5));

ratio_FD = diff_h5/diff_h4