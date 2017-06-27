%% --copyRight @ Shuhua Li <sue142857@gmail.com> 06/15/2017--%%
%-function: get the fixed candidate trunk

function CT = get_candidate_trunk(curves)

n = length(curves);
CT = zeros(n,2);
for i = 1:n
    CT(i,:) = curves{i}(1,:);
end