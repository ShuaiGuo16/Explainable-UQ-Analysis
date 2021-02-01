clear
clc

%%%%%%
%%% PDF comparison between MC & AS & Ana
%%%%%%

load 'beta.mat'
load 'gamma.mat'
load 'eigenvector_70.mat'
load 'f_MC_GR.mat'
load 'Analytical_Coeff.mat'
load 'model_Luis_30kW_9.5%_4_21_350ms.mat'

m = 16;
n = 1;
eg_vector = V(:,end);
mean_translation = model.B(4:19);
temp_cov = getcov(model);
variance_translation = sqrt(diag(temp_cov));
coeff = FTF_projection;

%% AS low-order model - data generation
sampling_num=40000;
Samples = lhsnorm(zeros(m,1),diag(ones(m,1)),sampling_num); 
X_input = ones(sampling_num,n+2);
X_input(:,2) = (eg_vector'*Samples')';
X_input(:,3) = X_input(:,2).^2;

f_ROM = X_input * beta;

%% Analytical solution
for i = 1:m
    Samples_scale(:,i) = Samples(:,i)*variance_translation(i) + mean_translation(i);
end

X_input_ana = ones(sampling_num,n+2);
for i = 1:sampling_num
    delta_h = Samples_scale(i,:) - mean_translation;
    X_input_ana(i,2) = coeff'*delta_h';
    X_input_ana(i,3) = X_input_ana(i,2)^2;
end
f_ANA = X_input_ana*gamma;

pts = -30:0.1:20;
figure(1)
[f_Pre,xi_Pre] = ksdensity(f_ROM,pts);
[f_Ana,xi_Ana] = ksdensity(f_ANA,pts);
plot(xi_Ana,f_Ana,'r-','LineWidth',2)
hold on
plot(xi_Pre,f_Pre,'k--','LineWidth',2)
% [f_Full,xi_Full] = ksdensity(f_MC_GR(1:40000,3),pts);
% plot(xi_Full,f_Full,'k--','LineWidth',2)
H = histogram(f_MC_GR(1:30000,3),'Normalization','pdf');
set(H,'FaceColor','black');

h = gca;
h.FontSize = 14;
h.XLim = [-30 20]
h.YLim = [0 0.08]
xlabel('Growth Rate (rad/s)')
ylabel('PDF')
xticks(-30:10:20)
yticks(0:0.02:0.08)
legend('Analytical','Active Subspace','Direct Monte Carlo','Location','NorthWest')

% plot([0 0],[0 0.08],'k--','LineWidth',2)

hold off

fig = gcf;
fig.PaperPositionMode = 'auto';
print('PDF_Comp_all','-dtiff','-r800')