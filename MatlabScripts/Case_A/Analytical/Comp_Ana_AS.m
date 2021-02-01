clear
clc

load 'eigenvector_70.mat'
load 'model_Luis_30kW_9.5%_4_21_350ms.mat'

delta_Ts = 6.25e-4;
FR_ref = 434.2*2*pi;
temp_mean = model.B(4:19);
temp_cov = getcov(model);   

mean_translation = temp_mean;   
variance_translation = sqrt(diag(temp_cov)); 

% Results from active subspace
a = V(:,end)./variance_translation;

% Results from analytical solution
% 1-determine the direction
 [ a_omega_value, b_omega_value ] = Direction_Cal_A( FR_ref, -4 );
 
  A = 0;   B = 0;
for i = 1:16
    A = A + mean_translation(i)*(i+2)*delta_Ts*cos((i+2)*delta_Ts*FR_ref);
    B = B + mean_translation(i)*(i+2)*delta_Ts*sin((i+2)*delta_Ts*FR_ref);
end
 
 direction = [b_omega_value+A,-a_omega_value-B];
% 2-projection
FTF_individual = zeros(16,2);
for i = 1:16
    FTF_individual(i,1) = cos((i+2)*delta_Ts*FR_ref);
    FTF_individual(i,2) = -sin((i+2)*delta_Ts*FR_ref);
end

FTF_projection = zeros(16,1);
for i = 1:16
    vector_a = direction;
    vector_b = FTF_individual(i,:);
    FTF_projection(i) = sum(vector_a.*vector_b)/sqrt(vector_a(1)^2+vector_a(2)^2);
end
save 'Analytical_Coeff.mat' FTF_projection
% Compare
figure(1)
plot(FTF_projection,a,'ok','MarkerSize',6,'LineWidth',1.2)
hold on

new_X = ones(16,2);
new_X(:,2) = FTF_projection;
new_beta = regress(a,new_X);

new_coor_x = -1:0.1:1;
new_coor_y = new_beta(2)*new_coor_x;
plot(new_coor_x,new_coor_y,'--r','LineWidth',1.2)
axis([-1 1 -12 12])

% Plot axis
plot([-1,1],[0,0],'--k','LineWidth',1.2)
plot([0,0],[-12,12],'--k','LineWidth',1.2)

% xlabel('Real','FontSize',14)
% ylabel('Imag','FontSize',14)
set(gca,'FontSize',12)

hold off

% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print('Analytical_measure_compare','-dtiff','-r600')