% Test: Plot paths for sig=10
clear
clc

load 'Datasample_10.mat'
delta_Ts = 6.25e-4;
FR_ref = 97.46*2*pi;
frequency_seq = linspace(0,1.2*FR_ref,30);
count_choose = [1,3,5,8,10,12,20];

for count = 1:size(count_choose,2)
    
    for i = 1:size(frequency_seq,2)
        cur_freq = frequency_seq(i);

        FTF = zeros(1,2);


            for j = 1:16
                FTF(1,1) = FTF(1,1)+Datasample_10(count_choose(count),j)*cos((j+2)*delta_Ts*cur_freq);
                FTF(1,2) = FTF(1,2)-Datasample_10(count_choose(count),j)*sin((j+2)*delta_Ts*cur_freq);
            end


        A = Acoustic_term( cur_freq, 0 );
        OLTF_real(i,1) = FTF(1,1)*real(A)-FTF(1,2)*imag(A);
        OLTF_imag(i,1) = FTF(1,1)*imag(A)+FTF(1,2)*real(A);
    end

        plot(OLTF_real,OLTF_imag,'-k','LineWidth',1.2)
        hold on

        % Plot critical point
        plot(-1,0,'>r','MarkerSize',10)
        % Plot final solution
        FTF_sol = zeros(1,2);
        for j = 1:16
                FTF_sol(1,1) = FTF_sol(1,1)+Datasample_10(count_choose(count),j)*cos((j+2)*delta_Ts*FR_ref);
                FTF_sol(1,2) = FTF_sol(1,2)-Datasample_10(count_choose(count),j)*sin((j+2)*delta_Ts*FR_ref);
        end
        A_sol = Acoustic_term( FR_ref, 0 );
        plot(FTF_sol(1,1)*real(A_sol)-FTF_sol(1,2)*imag(A_sol),FTF_sol(1,1)*imag(A_sol)+FTF_sol(1,2)*real(A_sol),'gd','MarkerSize',8,'MarkerFaceColor','g')


%         axis([-1.3 -0.5 -1.3 0.6])

end
hold off

fig = gcf;
fig.PaperPositionMode = 'auto';
print('BRS Burner-OLTF','-dtiff','-r600')