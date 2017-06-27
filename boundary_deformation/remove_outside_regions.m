%% --copyRight @ Shuhua Li <sue142857@gmail.com> 05/17/2017--%%
%--function: remove outside regions via scale the corresponding peices
%iterately 

function controlPs = remove_outside_regions(controlPs,CT)

n = length(controlPs);
xv = [CT(:,1);CT(1,1)];
yv = [CT(:,2);CT(1,2)];
scale = 0.99;
for i = 1:n
    C = controlPs{i};
    % check whether piece i has regions outside the trunk
    [in,on] = inpolygon(C(2:end-1,1),C(2:end-1,2),xv,yv);
    if min(in-on)
        continue;
    end
    % scale piece i until it has no regions outside
    newC = C;
    while ~min(in-on)
        newC = scale_poly(newC,scale);
        [in,on] = inpolygon(newC(2:end-1,1),newC(2:end-1,2),xv,yv);     
    end
    C(2:end-1,:) = newC(2:end-1,:);
    controlPs{i} = C;
end