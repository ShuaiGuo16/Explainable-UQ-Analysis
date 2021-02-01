% Test: what will happen if focusing on sig=10

clear
clc

load 'Datasample_10.mat'
delta_Ts = 6.25e-4;
FR_ref = 97.46*2*pi;
A = Acoustic_term( FR_ref, 0 );
options = optimoptions('fsolve','Display','off');
initial_value = [107.6*2*pi,-4]; 

%% Compute right frequency
f_FR = zeros(size(Datasample_10,1),1);
f_GR = zeros(size(Datasample_10,1),1);
for index = 1:size(Datasample_10,1)
   
        h = Datasample_10(index,:);    % current sample
        EigenFun = @(omega) Eigenmode_solver(omega,h);
        Eigen = fsolve(EigenFun, initial_value(1)-initial_value(2)*1i,options);    % solving characteristic equation
        f_FR(index) = real(Eigen);
        f_GR(index) = -imag(Eigen);
    
end

%% Compute FTF - true
FTF_true = zeros(size(Datasample_10,1),2);
for i = 1:size(Datasample_10,1)
    
    for j = 1:16
        FTF_true(i,1) = FTF_true(i,1)+Datasample_10(i,j)*cos((j+2)*delta_Ts*f_FR(i));
        FTF_true(i,2) = FTF_true(i,2)-Datasample_10(i,j)*sin((j+2)*delta_Ts*f_FR(i));
    end
    
end

OLTF_real_true = FTF_true(:,1)*real(A)-FTF_true(:,2)*imag(A);
OLTF_imag_true = FTF_true(:,1)*imag(A)+FTF_true(:,2)*real(A);

plot(OLTF_real_true,OLTF_imag_true,'sk','LineWidth',1.5)
hold on

%% Compute FTF - wrong
FTF = zeros(size(Datasample_10,1),2);
for i = 1:size(Datasample_10,1)
    
    for j = 1:16
        FTF(i,1) = FTF(i,1)+Datasample_10(i,j)*cos((j+2)*delta_Ts*FR_ref);
        FTF(i,2) = FTF(i,2)-Datasample_10(i,j)*sin((j+2)*delta_Ts*FR_ref);
    end
    
end

OLTF_real = FTF(:,1)*real(A)-FTF(:,2)*imag(A);
OLTF_imag = FTF(:,1)*imag(A)+FTF(:,2)*real(A);

plot(OLTF_real,OLTF_imag,'or','LineWidth',1.5)

%% Plot Axis
plot([-1.8,0.2],[0,0],'--k','LineWidth',1.2)
plot([0,0],[-1.8,0.2],'--k','LineWidth',1.2)

%% Plot critical point
plot(-1,0,'g>','MarkerSize',8,'MarkerFaceColor','g')
hold off

axis([-1.4 -0.6 -0.2 0.2])
% axis([-1.3 -0.7 -1.6 -0.85])
title('BRS Burner+shrink')
xlabel('Real','FontSize',14)
ylabel('Imag','FontSize',14)
set(gca,'FontSize',12)

% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print('BRS Burner-shrink','-dtiff','-r600')