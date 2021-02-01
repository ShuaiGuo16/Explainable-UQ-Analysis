% Developing process - frequency
clear
clc

load 'Datasample_010.mat'

delta_Ts = 6.25e-4;
% FR_ref = 97.46*2*pi;
FR_ref = 610.2504;
growthrate_seq = linspace(0,-10,11);
count_choose = [1,5,10];

for count = 1:3
    
    for i = 1:size(growthrate_seq,2)
        cur_grow = growthrate_seq(i);

        FTF = zeros(1,2);


            for j = 1:16
                FTF(1,1) = FTF(1,1)+exp(-cur_grow*(j+2)*delta_Ts)*Datasample_010(count_choose(count),j)*cos((j+2)*delta_Ts*FR_ref);
                FTF(1,2) = FTF(1,2)-exp(-cur_grow*(j+2)*delta_Ts)*Datasample_010(count_choose(count),j)*sin((j+2)*delta_Ts*FR_ref);
            end


        A = Acoustic_term( FR_ref, cur_grow );
        OLTF_real = FTF(1,1)*real(A)-FTF(1,2)*imag(A);
        OLTF_imag = FTF(1,1)*imag(A)+FTF(1,2)*real(A);

        plot(OLTF_real,OLTF_imag,'sk','MarkerSize',6,'LineWidth',1.5)
        hold on

        % Plot critical point
        plot(-1,0,'g*','MarkerSize',8)


        axis([-1.3 -0.5 -1.3 0.6])
        
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,8);
        if i == 1
            imwrite(imind,cm,'ex.gif','gif','Loopcount',inf)
        else
            imwrite(imind,cm,'ex.gif','gif','WriteMode','append','Delaytime',0.8)
        end


    end
    
end
hold off