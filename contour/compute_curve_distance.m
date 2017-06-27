%% --copyRight @ Shuhua Li <sue142857@gmail.com> 05/10/2017--%%

function d = compute_curve_distance(oldC,newC)

n = size(oldC,1);
Lines = zeros(n,2);
Lines(:,1) = (1:n)';
Lines(:,2) = [2:n 1]';
k1=LineCurvature2D(oldC,Lines);
k2=LineCurvature2D(newC,Lines);
d_ = zeros(n,1);
for i = 1:n
    edi = norm(oldC(i,:)-newC(i,:),2);
    if edi>0
        d_(i) = abs(k1(i)-k2(i))*edi/sum(abs(k1-k2));
    else
        d_(i) = 0;
    end
end
d = sum(d_);