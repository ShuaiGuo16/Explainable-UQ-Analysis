function [ GR ] = solver_GR( h )

options = optimoptions('fsolve','Display','off');
initial_value = [107.6*2*pi,-4];

EigenFun = @(omega) Eigenmode_solver(omega,h);
Eigen = fsolve(EigenFun, initial_value(1)-initial_value(2)*1i,options);    % solving characteristic equation
GR = -imag(Eigen);


end

