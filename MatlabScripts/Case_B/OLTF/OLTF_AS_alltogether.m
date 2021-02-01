% OLTF plotting
clear
clc

load 'model_Luis_30kW_9.5%_4_21_350ms.mat'
load 'beta.mat'
load 'eigenvector_70.mat'

delta_Ts = 6.25e-4;
FR_ref = 97.46*2*pi;
temp_mean = model.B(4:19);
temp_cov = getcov(model);   

mean_translation = temp_mean;   
variance_translation = sqrt(diag(temp_cov));

%% Basic Setup
N = [200,120,120,200];
m = 16;
A = Acoustic_term( FR_ref, -4 );

%% Sigma-isoline check - part1 - AS Results - transformation
figure(1)
for i = 1:4
    NumSamples = N(i);
    target_GrowthRate = (i-1)*10-20;
    AV = Reverse_AV( target_GrowthRate, beta);
    
    h1_15 = lhsnorm(zeros(m-1,1),diag(ones(m-1,1)),NumSamples); 
    h1_16_initial = zeros(NumSamples,16);
    
    index = 0;
for AV_i = 1:NumSamples
    h16_st_temp = (AV - sum(V(1:15,end)'.*h1_15(AV_i,:)))/V(16,end);
    
    if h16_st_temp>-3 && h16_st_temp<3
        
        index = index+1;
        for k = 1:m-1
            h1_16_initial(index,k) = h1_15(AV_i,k)*variance_translation(k)+mean_translation(k);
        end
        h1_16_initial(index,m) = h16_st_temp*variance_translation(m)+mean_translation(m);
        
    else
        
        continue
        
    end
    
end

P_total_AS = zeros(index,2);
for plot_i = 1:index
    
    for j = 1:16
        P_total_AS(plot_i,1) = P_total_AS(plot_i,1)+h1_16_initial(plot_i,j)*cos((j+2)*delta_Ts*FR_ref);
        P_total_AS(plot_i,2) = P_total_AS(plot_i,2)-h1_16_initial(plot_i,j)*sin((j+2)*delta_Ts*FR_ref);
    end

end
    [ Plot_P_total_AS ] = sortFTF( P_total_AS );
    
    Plot_OLTF_real = Plot_P_total_AS(:,1)*real(A)-Plot_P_total_AS(:,2)*imag(A);
    Plot_OLTF_imag = Plot_P_total_AS(:,1)*imag(A)+Plot_P_total_AS(:,2)*real(A);
    
    h1 = plot(Plot_P_total_AS(:,1),Plot_P_total_AS(:,2),'sk','LineWidth',1.5)
    hold on
    
    h2 = plot(Plot_OLTF_real,Plot_OLTF_imag,'ok','LineWidth',1.5)
    
end

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
title('BRS Burner+Intrinsic Mode+AS Transformation')
xlabel('Real','FontSize',14)
ylabel('Imag','FontSize',14)
set(gca,'FontSize',12)

fig = gcf;
fig.PaperPositionMode = 'auto';
print('BRS Burner','-dtiff','-r600')