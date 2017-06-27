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
w_lap = 0.3; % laplacian smoothing
w_ar = 0.8; % attraction-repulsion

%% Laplacian force & labeling
neig = 1;
id = 1;
for i = 1:nPs
    for j = 1: length(curve{i})
        label(id) = i;
        points(id * 2 - 1) = curve{i}(j, 1);
        points(id * 2) = curve{i}(j, 2);
        if j > neig && j <= length(curve{i}) - neig
            pre = id - neig;
            post = id + neig;
            T(2 * id - 1, pre * 2 - 1) = 0.5 * w_lap;
            T(2 * id, pre * 2) = 0.5;
            T(2 * id - 1, post * 2 - 1) = 0.5* w_lap;
            T(2 * id, post * 2) = 0.5;
            T(2 * id - 1, 2 * id - 1) = -1 * w_lap;
            T(2 * id, 2 * id) = -1;
        end
        id = id + 1;
    end
end

% the maixmum distance between any two points on the curves, for
% normalization, using the vertices on the truck for efficiency
maxd = 0;
for i = 1:length(Trunk)
    for j = 1:length(Trunk)
        d = norm(Trunk(i) - Trunk(j));
        if d > maxd 
            maxd = d;
        end
    end
end

%% Attraction-Repulsion force
% Note the end points on the curve should be remained unchanged
k = 0.01;
% consider nearby points within radius R0
R0 = maxd / 6;
thr = 1e-3;
ncount = 0;
for i = 1:n
    % proposed move
    sum_w = 0;
    pi = points(2 * i - 1 : 2 * i);
    shouldMove = true;
    stop_id = 1;
    for j = 1:n
        % the force between points on the same curve is 0
        if label(i) == label(j)
            continue;
        end
        pj = points(2 * j - 1 : 2 * j);
        % distance
        d = norm(pj - pi);
%         if d < thr
%             % boundary point
%             shouldMove = false;
%             sum_w = 0;
%             stop_id = j;
%             break;
%         end
        if d < R0 && d > thr
            d = d / maxd;
            w = k * (1 - d) * w_ar;
            T(2 * i - 1, 2 * j - 1) = T(2 * i - 1, 2 * j - 1) + w;
            T(2 * i, 2 * j) = T(2 * i, 2 * j) + w;
            sum_w = sum_w - w;
            ncount = ncount + 1;
        end
    end % each i
    if ~shouldMove
        for j = 1:stop_id
            T(2 * i - 1, 2 * j - 1) = 0;
            T(2 * i, 2 * j) = 0;
        end
    end
    T(2 * i - 1, 2 * i - 1) = 1 - sum_w + T(2 * i - 1, 2 * i - 1);
    T(2 * i, 2 * i) = 1- sum_w + T(2 * i, 2 * i);
end
ncount
points = T * points;

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