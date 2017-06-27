%-----copyRight(c) Shuhua Li<sue142857@gmail.com> 02.28.2017-----%
% -input:
%        -poly:
% -output:
%         -points:

function s = sample_points_via_midpoint(points,index_fea,d_sample)

% sample midpoints for each edge
n_edge = length(index_fea)-1; 
n = size(points,1);
space = zeros(0,1);
index_edge = cell(0,1);
loc = zeros(0,1);

for i = 1:n_edge
    edge = [index_fea(i);index_fea(i+1)];
    d = norm(points(index_fea(i),:) - points(index_fea(i+1),:),2);
    while d > d_sample
        n_mid = size(edge,1)-1;
        index_mid = zeros(n_mid,1);
        new_edge = zeros(n_mid+n_mid+1,1);
        for j = 1:n_mid
            if edge(j) < edge(j+1)
                interval = edge(j+1) - edge(j);
            else
                interval = n - edge(j) + edge(j+1);
            end
            index_mid(j) = edge(j) + floor(interval/2);
            if index_mid(j) > n
                index_mid(j) = index_mid(j) - n;
            end
        end
        new_edge(1:2:end,1) = edge;
        new_edge(2:2:end,1) = index_mid;
        edge = new_edge;
        d = norm(points(edge(1),:) - points(edge(2),:),2);
    end
    n_sample = size(edge,1)-1;
    space(end+1:end+n_sample,1) = edge(1:end-1,1);
    if i == 1
        loc1 = {[n_edge 1]}; 
    else
        loc1 = {[i-1 i]};
    end
    index_edge(end+1,1)  = loc1;
    index_edge(end+1:end+n_sample-1,1)  = repmat({i},n_sample-1,1);
    loc(end+1,1) = 0;
    loc(end+1:end+n_sample-1,1)  = i;
end
if edge(end)~=space(1)
   space = [space;edge(end)];
end
s.index = space;
s.index_edge = index_edge;
s.loc = loc;