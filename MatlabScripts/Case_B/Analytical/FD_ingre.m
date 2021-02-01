function [ freq_diff, growth_diff, F, f_FR, f_GR] = FD_ingre( h )

initial_value = [612.4085, -4.0032];
delta_t = 6.25e-4;
options = optimoptions('fsolve','Display','off');
EigenFun = @(omega) Eigenmode_solver(omega,h);
Eigen = fsolve(EigenFun, initial_value(1)-initial_value(2)*1i,options);    % solving characteristic equation
f_FR = real(Eigen);
f_GR = -imag(Eigen);

freq_diff = f_FR-initial_value(1);
growth_diff = f_GR-initial_value(2);

omega = f_FR-f_GR*1i;
F = 0;
for k = 1:16
    F = F + h(k)*exp(-omega*1i*(k+2)*delta_t);
end



end
