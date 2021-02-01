function [ new_coor ] = sortFTF( coor )

[Y,index] = sort(coor(:,1));

for i = 1:size(coor,1)
    new_coor(i,:) = coor(index(i),:);
end


end

