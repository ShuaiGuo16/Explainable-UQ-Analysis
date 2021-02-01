function [ f ] = TargetGrowthRate( x,vector )

initial_value = [107.6*2*pi,-4;];
options = optimoptions('fsolve','Display','off');

target_GrowthRate = -4;

h(1:15) = vector;
h(16) = x;

        EigenFun = @(omega) Eigenmode_solver(omega,h);
        Eigen = fsolve(EigenFun, initial_value(1)-initial_value(2)*1i,options);    % solving characteristic equation
        f_GR = -imag(Eigen);

f = (f_GR-target_GrowthRate)^2;

end

