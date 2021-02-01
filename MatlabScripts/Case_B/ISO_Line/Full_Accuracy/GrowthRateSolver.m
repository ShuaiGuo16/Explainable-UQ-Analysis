clear
clc


% Target value
target_GrowthRate = -4;

m = 16;  

load 'eigenvector_70.mat'
load 'beta.mat'
load 'model_Luis_30kW_9.5%_4_21_350ms.mat'
temp_mean = model.B(4:19);
temp_cov = getcov(model);   

mean_translation = temp_mean;   
variance_translation = sqrt(diag(temp_cov)); 

% Number of points
N = 80;

% active variable
AV = Reverse_AV( target_GrowthRate, beta);

% Initial guess
h1_15 = lhsnorm(zeros(m-1,1),diag(ones(m-1,1)),N); 
% Calculate corresponding tau values
index = 0;
for i = 1:N
    h16_st_temp = (AV - sum(V(1:15,end)'.*h1_15(i,:)))/V(16,end);
    
    if h16_st_temp>-3 && h16_st_temp<3
        
        index = index+1;
        for k = 1:m-1
            h1_16_initial(index,k) = h1_15(i,k)*variance_translation(k)+mean_translation(k);
        end
        h1_16_initial(index,m) = h16_st_temp*variance_translation(m)+mean_translation(m);
        
    else
        
        continue
        
    end
    
end

Datasample_004 = zeros(size(h1_16_initial,1),16);
Datasample_004(:,1:15) = h1_16_initial(:,1:15);
counter = size(h1_16_initial,1)
for i = 1:size(h1_16_initial,1)
    x0 = h1_16_initial(i,end);
    Datasample_004(i,16) = fminsearch(@(x) TargetGrowthRate(x,h1_16_initial(i,1:15)),x0);
    counter = counter - 1
end

% Check procedure
for i = 1:size(Datasample_004,1)
    
    options = optimoptions('fsolve','Display','off');
    initial_value = [107.6*2*pi,-4;];  
    
    EigenFun = @(omega) Eigenmode_solver(omega,Datasample_004(i,:));
    Eigen = fsolve(EigenFun, initial_value(1)-initial_value(2)*1i,options); 
    predict_GrowthRate(i) = -imag(Eigen);
    
end

save 'Datasample_004.mat' Datasample_004
