% To calculate the projection

clear
clc
% Load data
load 'Datasample_004.mat'
load 'eigenvector_70.mat'
load 'model_Luis_30kW_9.5%_4_21_350ms.mat'
delta_Ts = 6.25e-4;
FR_ref = 434.2*2*pi;
temp_mean = model.B(4:19);
temp_cov = getcov(model);   
variance_translation = sqrt(diag(temp_cov)); 

% measurement - Active Subspace
a = V(:,end)./variance_translation; 

% Calculate representations in FTF-plane
sample_num = size(Datasample_004,1);
P_total_ref = zeros(sample_num,2);
for i = 1:size(Datasample_004,1)
    
    for j = 1:16
        P_total_ref(i,1) = P_total_ref(i,1)+Datasample_004(i,j)*cos((j+2)*delta_Ts*FR_ref);
        P_total_ref(i,2) = P_total_ref(i,2)-Datasample_004(i,j)*sin((j+2)*delta_Ts*FR_ref);
    end
    
end

% Regression
X = ones(sample_num,2);
X(:,2) = P_total_ref(:,1);
Y = P_total_ref(:,2);
beta = regress(Y,X); 

% Orthogonal direction
direction = [-1,1/beta(2)];

% Individual vector
FTF_individual = zeros(16,2);
for i = 1:16
    FTF_individual(i,1) = cos((i+2)*delta_Ts*FR_ref);
    FTF_individual(i,2) = -sin((i+2)*delta_Ts*FR_ref);
end

% Calculate projection
FTF_projection = zeros(16,1);
for i = 1:16
    vector_a = direction;
    vector_b = FTF_individual(i,:);
    FTF_projection(i) = sum(vector_a.*vector_b)/sqrt(vector_a(1)^2+vector_a(2)^2);
end

% Compare results
entry_x = a;
entry_y = FTF_individual(:,2);

figure(1)
plot(entry_x,entry_y,'ok','MarkerSize',6,'LineWidth',1.2)
hold on

new_X = ones(16,2);
new_X(:,2) = entry_x;
new_beta = regress(entry_y,new_X);

new_coor_x = -12:1:12;
new_coor_y = new_beta(2)*new_coor_x;
plot(new_coor_x,new_coor_y,'--r','LineWidth',1.2)
axis([-12 12 -1 1])

% Plot axis
plot([-12,12],[0,0],'--k','LineWidth',1.2)
plot([0,0],[-1,1],'--k','LineWidth',1.2)

% xlabel('Real','FontSize',14)
% ylabel('Imag','FontSize',14)
set(gca,'FontSize',12)

hold off

fig = gcf;
fig.PaperPositionMode = 'auto';
print('Sensitivity_measure_compare','-dtiff','-r600')