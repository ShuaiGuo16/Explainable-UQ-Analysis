% OLTF plotting
clear
clc

load 'Datasample_00.mat'
load 'Datasample_10.mat'
load 'Datasample_010.mat'
load 'Datasample_020.mat'
load 'model_Luis_30kW_9.5%_4_21_350ms.mat'


delta_Ts = 6.25e-4;
FR_ref = 97.46*2*pi;
temp_mean = model.B(4:19);
 


%% Sigma-isoline check - part2 - Full Accuracy Results
A = Acoustic_term( FR_ref, 0 );
P_total_ref = zeros(size(Datasample_00,1),2);
for i = 1:size(Datasample_00,1)
    
    for j = 1:16
        P_total_ref(i,1) = P_total_ref(i,1)+Datasample_00(i,j)*cos((j+2)*delta_Ts*FR_ref);
        P_total_ref(i,2) = P_total_ref(i,2)-Datasample_00(i,j)*sin((j+2)*delta_Ts*FR_ref);
    end
    
end
plot(P_total_ref(:,1),P_total_ref(:,2),'sk','LineWidth',1.5)
hold on

plot(P_total_ref(:,1)*real(A)-P_total_ref(:,2)*imag(A),P_total_ref(:,1)*imag(A)+P_total_ref(:,2)*real(A),'ok','LineWidth',1.5)


P_total_ref = zeros(size(Datasample_10,1),2);
for i = 1:size(Datasample_10,1)
    
    for j = 1:16
        P_total_ref(i,1) = P_total_ref(i,1)+Datasample_10(i,j)*cos((j+2)*delta_Ts*FR_ref);
        P_total_ref(i,2) = P_total_ref(i,2)-Datasample_10(i,j)*sin((j+2)*delta_Ts*FR_ref);
    end
    
end
plot(P_total_ref(:,1),P_total_ref(:,2),'sk','LineWidth',1.5)
plot(P_total_ref(:,1)*real(A)-P_total_ref(:,2)*imag(A),P_total_ref(:,1)*imag(A)+P_total_ref(:,2)*real(A),'ok','LineWidth',1.5)

P_total_ref = zeros(size(Datasample_010,1),2);
for i = 1:size(Datasample_010,1)
    
    for j = 1:16
        P_total_ref(i,1) = P_total_ref(i,1)+Datasample_010(i,j)*cos((j+2)*delta_Ts*FR_ref);
        P_total_ref(i,2) = P_total_ref(i,2)-Datasample_010(i,j)*sin((j+2)*delta_Ts*FR_ref);
    end
    
end
plot(P_total_ref(:,1),P_total_ref(:,2),'sk','LineWidth',1.5)
plot(P_total_ref(:,1)*real(A)-P_total_ref(:,2)*imag(A),P_total_ref(:,1)*imag(A)+P_total_ref(:,2)*real(A),'ok','LineWidth',1.5)

P_total_ref = zeros(size(Datasample_020,1),2);
for i = 1:size(Datasample_020,1)
    
    for j = 1:16
        P_total_ref(i,1) = P_total_ref(i,1)+Datasample_020(i,j)*cos((j+2)*delta_Ts*FR_ref);
        P_total_ref(i,2) = P_total_ref(i,2)-Datasample_020(i,j)*sin((j+2)*delta_Ts*FR_ref);
    end
    
end
plot(P_total_ref(:,1),P_total_ref(:,2),'sk','LineWidth',1.5)
plot(P_total_ref(:,1)*real(A)-P_total_ref(:,2)*imag(A),P_total_ref(:,1)*imag(A)+P_total_ref(:,2)*real(A),'ok','LineWidth',1.5)

%% FTF mini-vector
FTF_total_ref = [0,0];
for i = 1:16
    FTF_total_ref(1) = FTF_total_ref(1)+temp_mean(i)*cos((i+2)*delta_Ts*FR_ref);
    FTF_total_ref(2) = FTF_total_ref(2)-temp_mean(i)*sin((i+2)*delta_Ts*FR_ref);
end
plot([0,FTF_total_ref(1)],[0,FTF_total_ref(2)],'-k','LineWidth',0.8)

for i = 1:16
    FTF_individual(1) = temp_mean(i)*cos((i+2)*delta_Ts*FR_ref);
    FTF_individual(2) = -temp_mean(i)*sin((i+2)*delta_Ts*FR_ref);
    plot([0,FTF_individual(1)],[0,FTF_individual(2)],'-k','LineWidth',0.8)
end


%% Plot Axis
plot([-1.8,0.2],[0,0],'--k','LineWidth',1.2)
plot([0,0],[-1.8,0.2],'--k','LineWidth',1.2)

%% Plot critical point
plot(-1,0,'r>','MarkerSize',10,'MarkerFaceColor','r')

hold off

axis([-1.8 0.2 -1.8 0.2])
% axis([-1.3 -0.7 -1.6 -0.85])
title('BRS Burner+Intrinsic Mode+Full Transformation')
xlabel('Real','FontSize',14)
ylabel('Imag','FontSize',14)
set(gca,'FontSize',12)

fig = gcf;
fig.PaperPositionMode = 'auto';
print('BRS Burner-1','-dtiff','-r600')