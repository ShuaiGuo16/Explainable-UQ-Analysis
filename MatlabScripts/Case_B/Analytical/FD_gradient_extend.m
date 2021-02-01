% Strict FD way to calculate gradient
clear
clc

delta_t = 6.25e-4;
ref = [612.4085, -4.0032];
load 'model_Luis_30kW_9.5%_4_21_350ms.mat'
h_ref = model.B(4:19);
h_per_1 = model.B(4:19);
h_per_1(2) = h_per_1(2)*1.04;

h_per_2 = model.B(4:19);
h_per_2(10) = h_per_2(10)*1.04;

omega = ref(1)-ref(2)*1i;
F_ref = 0;
for k = 1:16
    F_ref = F_ref + h_ref(k)*exp(-omega*1i*(k+2)*delta_t);
end

[ lin_A(1,1), lin_A(1,2), F_per_1, f_FR, f_GR] = FD_ingre( h_per_1 );
[ lin_A(2,1), lin_A(2,2), F_per_2, f_FR, f_GR] = FD_ingre( h_per_2 );
lin_b = [real(F_per_1-F_ref);real(F_per_2-F_ref)];
a_diff = lin_A\lin_b