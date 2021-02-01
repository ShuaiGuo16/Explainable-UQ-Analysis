% Perform UQ by using analytical derived coefficients
clear
clc

load 'Analytical_Coeff.mat'
load 'f_GR.mat'
load 'Samples_N.mat'
load 'model_Luis_30kW_9.5%_4_21_350ms.mat'

mean_translation = model.B(4:19); 
temp_cov = getcov(model);
variance_translation = sqrt(diag(temp_cov));

% Results from analytical solution
coeff = FTF_projection;

n = 1;
m = size(coeff,1);

for i = 1:size(Samples_N,1)
    
    delta_h = Samples_N(i,:) - mean_translation;
    AV(i,1) =  coeff'*delta_h';
    
end

figure(1)
plot(AV, f_GR,'ko')
hold on

index = [5,12,19,22,31,37,91];    % Ready for fitting

X = ones(7,n+2);
X(:,2) = AV(index);
X(:,3) = AV(index).^2;

gamma = regress(f_GR(index),X); 
save 'gamma.mat' gamma
%-----------------------------------------------------


coor_x = -0.3:0.05:0.3;
coor_y = gamma(1)+gamma(2)*coor_x+gamma(3)*coor_x.^2;
plot(coor_x,coor_y)
hold off




