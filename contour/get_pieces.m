%-----copyRight(c) Shuhua Li<sue142857@gmail.com> 04.04.2017-----%
% function: 

function pieces = get_pieces(T,m,direction)

n = length(T);
pieces = cell(n,1);
if strcmp(direction,'counter-clockwise')
    for i = 1:n
        v1 = T(i);
        if i == n
            v2 = T(1);
        else
            v2 = T(i+1);
        end
        if v1 > v2
            pieces{i} = [v1:m 1:v2];
        else
            pieces{i} = v1:v2;
        end
    end
end

if strcmp(direction,'clockwise')
    for i = 1:n
        v1 = T(i);
        if i == n
            v2 = T(1);
        else
            v2 = T(i+1);
        end
        if v1 < v2
            pieces{i} = [v1:-1:1 m:-1:v2];
        else
            pieces{i} = v1:-1:v2;
        end
    end 
end


