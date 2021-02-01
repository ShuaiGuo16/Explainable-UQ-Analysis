clear
clc

load 'model_Luis_30kW_9.5%_4_21_350ms.mat'   
temp_mean = model.B(4:19); 
options = optimoptions('fsolve','Display','off');
initial_value = [107.6*2*pi,-4]; 

EigenFun = @(omega) Eigenmode_solver(omega,temp_mean);
Eigen = fsolve(EigenFun, initial_value(1)-initial_value(2)*1i,options);    % solving characteristic equation
f_FR = real(Eigen);
f_GR = -imag(Eigen);


omega = f_FR - f_GR*1i; 
delta_t = 6.25e-4;
F = 0;
for k = 1:16
    F = F + temp_mean(k)*exp(-omega*1i*(k+2)*delta_t);
end

A = Acoustic_term( f_FR, f_GR );
result = double(A*F)
