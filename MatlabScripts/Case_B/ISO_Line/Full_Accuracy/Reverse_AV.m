function [ AV ] = Reverse_AV( GrowthRate, beta )

AV_1 = (-beta(2)+sqrt(beta(2)^2-4*beta(3)*(beta(1)-GrowthRate)))/(2*beta(3));
AV_2 = (-beta(2)-sqrt(beta(2)^2-4*beta(3)*(beta(1)-GrowthRate)))/(2*beta(3));

if AV_1>-3 && AV_1<3
    AV = AV_1;
else
    AV = AV_2;
end

end

