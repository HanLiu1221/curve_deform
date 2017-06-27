function curve = deformCurve(curve, Trunk)

% n pieces of curves
nPs = length(curve);

n = 0; % total number of points
for i = 1:nPs
    n = n + length(curve{i});
end
n2 = n * 2;

% labeling points to indicate which piece it belongs to
label = zeros(1, n);
points = zeros(n2, 1);
T = zeros(n2, n2); 

% weights
w_lap = 0.1; % laplacian smoothing
w_ar = 0.9; % attraction-repulsion

%% Laplacian force & labeling
neig = 1;
id = 1;
boundIdx = zeros(1, length(Trunk) * 2);
for i = 1:nPs
    boundIdx(i * 2 - 1) = id;
    for j = 1: length(curve{i})
        label(id) = i;
        points(id * 2 - 1) = curve{i}(j, 1);
        points(id * 2) = curve{i}(j, 2);
        if j > neig && j <= length(curve{i}) - neig
            pre = id - neig;
            post = id + neig;
            T(2 * id - 1, pre * 2 - 1) = 0.5 * w_lap;
            T(2 * id, pre * 2) = 0.5 * w_lap;
            T(2 * id - 1, post * 2 - 1) = 0.5 * w_lap;
            T(2 * id, post * 2) = 0.5 * w_lap;
            T(2 * id - 1, 2 * id - 1) = -1 * w_lap;
            T(2 * id, 2 * id) = -1 * w_lap;
        end
        id = id + 1;
    end
    boundIdx(i * 2) = id - 1;
end

% average distance
avgd = 0;
nd = 0;
maxd = 0;
for i = 1:n
    pi = points(2 * i - 1 : 2 * i);
    for j = 1:n
        if i == j
            continue;
        end
        pj = points(2 * j - 1 : 2 * j);
        d = norm(pj - pi);
        avgd = avgd + d;
        nd = nd + 1;
        if d > maxd
            maxd = d;
        end
    end
end
avgd = avgd / nd;

k0 = 0.3;
k1 = 2.0 * k0;
R0 = k0 * avgd;
sigma = R0;
R1 = k1 * avgd;

% % the maixmum distance between any two points on the curves, for
% % normalization, using the vertices on the truck for efficiency
% maxd = 0;
% for i = 1:length(Trunk)
%     for j = 1:length(Trunk)
%         d = norm(Trunk(i) - Trunk(j));
%         if d > maxd 
%             maxd = d;
%         end
%     end
% end

%% Attraction-Repulsion force
% Note the end points on the curve should be remained unchanged
k = 0.01;
% consider nearby points within radius R0
thr = 1e-4;
ncount = 0;
for i = 1:n
    pi = points(2 * i - 1 : 2 * i);
    % end points cannot be moved
    isEndPnt = find(boundIdx == i);
    if ~isempty(isEndPnt)
        continue;
    end
    % proposed move
    sum_w = 0;
    for j = 1:n
        % the force between points on the same curve is 0
        if label(i) == label(j)
            continue;
        end
        pj = points(2 * j - 1 : 2 * j);
        % distance
        d = norm(pj - pi);
%         if d / avgd < R0 && d > thr
        if d < R0 && d > thr
            d = d / maxd;
            w = k * (1 - d) * w_ar;
%             w = LJ(sigma, d / avgd);
%             w = w / d * w_ar;
            T(2 * i - 1, 2 * j - 1) = T(2 * i - 1, 2 * j - 1) + w;
            T(2 * i, 2 * j) = T(2 * i, 2 * j) + w;
            sum_w = sum_w - w;
            ncount = ncount + 1;
        end
    end % each i
    T(2 * i - 1, 2 * i - 1) = sum_w + T(2 * i - 1, 2 * i - 1);
    T(2 * i, 2 * i) = sum_w + T(2 * i, 2 * i);
end
ncount
points = points + T * points;

% re-assign
id = 1;
for i = 1:nPs
    for j= 1: length(curve{i})
        label(id) = i;
        curve{i}(j, 1) = points(id * 2 - 1);
        curve{i}(j, 2) = points(id * 2);
        id = id + 1;
    end
end

end

function w = LJ(sigma, r)
%% Lennard-Jones potential
w = power(sigma / r, 12) - power(sigma / r, 6);
end
