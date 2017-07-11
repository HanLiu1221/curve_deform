%% split a curve into a few curve segments for deformation
function splitPointIds = getSplitCurvePointIds(curve, shouldReverse)
%% OUTPUT 
% [s1, handle, e2; ...]
% simplify the polygon
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