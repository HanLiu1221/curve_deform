%-----copyRight(c) Shuhua Li<sue142857@gmail.com> 03.23.2017-----%
%--function: compute edge length of a polygon

function L = compute_edgeLength(P)

n = size(P,1);
L = zeros(n,1);

for i = 1:n
    v1 = P(i,:);
    if i == n
        v2 = P(1,:);
    else
        v2 = P(i+1,:);
    end
    L(i) = norm(v1-v2,2);
end