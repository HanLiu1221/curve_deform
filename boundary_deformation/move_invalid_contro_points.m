%% --copyRight @ Shuhua Li <sue142857@gmail.com> 05/12/2017--%%

function [controlPs,oriControlPs] = move_invalid_contro_points(controlPs,oriControlPs)

% move control points with extra-large curvature
n = length(controlPs);
for i = 1:n
    Ci = controlPs{i};
    j = 2;
    while j< size(Ci,1)
        v1 = Ci(j-1,:);
        v2 = Ci(j,:);
        v3 = Ci(j+1,:);
        [angle,~] = compute_angle(v1,v2,v3);
        if abs(angle) < 10 && angle>0
            Ci(j,:) = [];
            oriControlPs{i}(j,:) = [];
            j = j-1;
            if j <2
                j = 2;
            end
            continue;
        end
        j = j+1;
    end
    controlPs{i} = Ci;
end