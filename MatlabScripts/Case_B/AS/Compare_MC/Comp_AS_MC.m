clear
clc

%%%%%%
%%% PDF comparison between MC & AS
%%%%%%

load 'beta.mat'
load 'eigenvector_70.mat'
load 'f_MC_GR.mat'

m = 16;
n = 1;
eg_vector = V(:,end);

%% AS low-order model - data generation
sampling_num=40000;
Samples = lhsnorm(zeros(m,1),diag(ones(m,1)),sampling_num); 
X_input = ones(sampling_num,n+2);
X_input(:,2) = (eg_vector'*Samples')';
X_input(:,3) = X_input(:,2).^2;

f_ROM = X_input * beta;

pts = -50:0.1:30;
figure(1)
[f_Pre,xi_Pre] = ksdensity(f_ROM,pts);
plot(xi_Pre,f_Pre,'k-','LineWidth',2)
hold on
% [f_Full,xi_Full] = ksdensity(f_MC_GR,pts);
% plot(xi_Full,f_Full,'k--','LineWidth',2)
H = histogram(f_MC_GR(1:40000),'Normalization','pdf');
set(H,'FaceColor','black');

h = gca;
h.FontSize = 14;
h.XLim = [-50 30]
h.YLim = [0 0.05]
xlabel('Growth Rate (rad/s)')
ylabel('PDF')
xticks(-50:10:30)
yticks(0:0.01:0.05)
legend('Active Subspace','Monte Carlo','Location','NorthWest')

plot([0 0],[0 0.05],'k--','LineWidth',2)

hold off

% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print('PDF_Comp_B','-dtiff','-r800')