%% split a curve into a few curve segments for deformation
function [splitPointIds, handleIds] = getSplitCurvePointIds(curve, shouldReverse)
%% OUTPUT 
% [s1, handle, e2; ...]
% simplify the polygon
n = length(curve);
d = 0.01;
simp_poly = dpsimplify(curve, d);
[angles, types] = compute_curvature_poly(simp_poly);
if shouldReverse == 1
    types = 1 - types;
end
sids = findPointIndices(curve, simp_poly);
% keep only those with curvature changes
curvature = LineCurvature2D(curve);
dt = 3;
pids = zeros(1, length(sids));
j = 1;
nthr = 15;
for i = 1:length(sids)
    id = sids(i);
    if id <= dt || id >= length(curve) - dt
        continue;
    end
    if j > 1 && sids(i) - pids(j - 1) < nthr
        continue;
    end
    if types(i) == 1
        pids(j) = sids(i);
        j = j + 1;
    end
end
% curvature(pids(1, 1 :j - 1))
splitPointIds = zeros(j, 3);
handleIds = zeros(1, 0);
for i = 1 : j
    if i == 1
        splitPointIds(i, 1) = 1;
    else
        splitPointIds(i, 1) = pids(i - 1);
    end
    if i == j
        splitPointIds(i, 3) = length(curve);
    else
        splitPointIds(i, 3) = pids(i);
    end
    % handle point
%     curvs = curvature(splitPointIds(i, 1) : splitPointIds(i, 3), 1);
%     curvs = abs(curvs);
%     [m, id] = max(curvs);
    id = floor((splitPointIds(i, 1) + splitPointIds(i, 3)) / 2);
    splitPointIds(i, 2) = splitPointIds(i, 1) + id - 1;
end
if size(splitPointIds, 1) >= 2
    return;
end
% check if it can be split into two lines
% find the contex point, i.e., the highest point
% check the projection length from the point to the line formed by two
% end points
prev = 1;
post = n;
hIds = zeros(1, length(simp_poly));
hid = 0;
thr = 0.05;
for i = 1:length(simp_poly)
    pid = i;
    cid = sids(i);
    if types(pid) == 1 || cid < nthr || cid > n - nthr
        continue;
    end
    % convex
    if i + 1 < length(simp_poly)
        post = sids(i + 1);
    else 
        post = n;
    end
    % fit two lines
    [line1, error1] = polyfit(curve(prev : cid, 1), curve(prev : cid, 2), 1);
    [line2, error2] = polyfit(curve(cid : post, 1), curve(cid : post, 2), 1);
    prev = cid;
    if error1.normr < thr && error2.normr < thr
        hid = hid + 1;
        hIds(hid) = cid;
    end
    disp(strcat('line 1 fit error: ', num2str(error1.normr)));
    disp(strcat('line 2 fit error: ', num2str(error2.normr)));
end
if hid < 3
    handleIds = hIds(1, 1 : hid);
end
end

function isLM = isLocalMaxMin(curvature, id, d)
isLM = 0;
curvature = abs(curvature);
nthr = 2 * d;
n = 0;
for i = id - d : id + d
    if curvature(id) > curvature(i)
        n = n + 1;
    end
end
if n >= nthr
    isLM = 1;
end
end

function isInfl = isInflectionPoint(curvature, id, d)
isInfl = 1;
nthr = floor(d * d / 2);
nr = 0;
for i = id - d : id - 1
    for j = id + 1 : id + d
        if curvature(i) * curvature(j) < 0
            nr = nr + 1;
        end
    end
end
if nr < nthr
    isInfl = 0;
    return;
end
nr = 0;
for i = 1 : d
    for j = 1 : d
        if curvature(id - i) * curvature(id - j) > 0
            nr = nr + 1;
        end
        if curvature(id + i) * curvature(id + j) > 0
            nr = nr + 1;
        end
    end
end
if nr - d * 2 < nthr
    isInfl = 0;
end
end

%% identify the indices of feature points after polygon simplification
function ids = findPointIndices(curve, poly) 
ids = ones(size(poly, 1), 1);
for i = 1 : size(poly, 1)
    mind = 100000;
    for j = 1 : size(curve, 1)
        d = norm(curve(j, :) - poly(i, :));
        if d < mind
            mind = d;
            ids(i) = j;
        end
    end
end
end

%% In any case, merge if there exists dense feature points.
function mergedIds = mergeSegs(pIds, n_thr)
n = 2;
ids = pIds;
for i = 2: length(pIds)
    if pIds(i) - ids(n - 1) > n_thr
        ids(n) = pIds(i);
        n = n + 1;
    end
end
ids(n - 1) = pIds(length(pIds));
mergedIds = ids(1 : n - 1);
end