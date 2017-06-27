%% --copyRight @ Shuhua Li <sue142857@gmail.com> 05/17/2017--%%

function  fea = compute_localFea(C)

n = size(C,1);
fea = zeros(n,3); % [e1 e2 theta12]
E = compute_edgeLength(C);
[A,~] = compute_curvature_poly(C); 
fea(:,2) = E;
fea(:,1) = [E(end);E(1:end-1)];
fea(:,3) = A;
