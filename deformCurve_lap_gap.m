function [deformed, area_gap, area_overlap] = deformCurve_lap_gap(curves, Trunk, feaIds)
%%

% n pieces of curves
nCurves = length(curves);

%% compute offset
nthr = 6;
n2static = 2;
% analyze each curve
deformed = cell(1, nCurves);
nh = 1;
handle_anchor = zeros(1, nh);
offsets = zeros(nh, 2);
deformed = curves;
%% Iterate over each curve segment
[~, area_gap, area_overlap] = compute_gap_overlap_area(curves, Trunk);
overlap_thr = 1e-6;
iter = 0;
max_iter = 1;
ratio = 1.0;
while iter < max_iter % change offset length
    for i = 1 : nCurves
        nseg = size(feaIds{i}, 1);
        deformed{i} = curves{i};
        % for each curve segment
        for s = 1 : nseg
            activePids = feaIds{i}(s, 1) : feaIds{i}(s, 3);
            if length(activePids) < nthr
                continue;
            end
            % jump end points
            s1 = activePids(1) + n2static;
            s2 = activePids(length(activePids) - n2static);
            static_anchor = zeros(1, n2static * 2);
            for k = 1 : n2static
                static_anchor(k) = activePids(k);
                static_anchor(k + n2static) = activePids(length(activePids) - n2static + k);
            end
            % points to compute
            ipoints = curves{i}(s1 : s2, :);
            max_off_s = 0;
            max_offset_s = [0, 0];
            max_handleId = floor((s1 + s2)/2);
            for ip = 1 : size(ipoints, 1)
                pi = ipoints(ip, :);
                % find nearest neighboring point
                min_dij = intmax('int64');
                min_offset_ij = [0, 0];
                for j = 1 : nCurves
                    if i == j
                        continue;
                    end
                    jpoints = curves{j};
                    for jp = 1 : size(jpoints, 1)
                        pj = jpoints(jp, :);
                        d = norm(pj - pi);
                        dir = (pj - pi) / d;
                        if d < min_dij
                            min_dij = d;
                            min_offset_ij = d / 2 * dir;
                        end
                    end
                end % j curve
                if min_dij ~= intmax('int64') && min_dij > max_off_s
                    max_off_s = min_dij;
                    max_offset_s = min_offset_ij;
                    max_handleId = feaIds{i}(s, 1) + ip - 1;
                end
            end % point ip at i curve s segment
            % deform     
            if max_off_s ~= 0
                handle_anchor(1) = max_handleId;   
                offsets(1, :) = max_offset_s * ratio;
                prevCurve = curves{i};
                curves{i} = lap2D(curves{i}, static_anchor, handle_anchor, offsets);
                [~, cur_gap, cur_overlap] = compute_gap_overlap_area(curves, Trunk);
                if cur_overlap > overlap_thr || cur_gap >= area_gap
                    curves{i} = prevCurve;
                    curves{i} = translateCurve(curves{i}, max_offset_s);
                    [~, cur_gap, cur_overlap] = compute_gap_overlap_area(curves, Trunk);
                end
                if cur_overlap > overlap_thr || cur_gap >= area_gap
                    curves{i} = prevCurve; %unchanged
                else
                    deformed{i}(s1 : s2, :) = curves{i}(s1 : s2, :);
                    area_overlap = cur_overlap;
                    area_gap = cur_gap;
                    return;
                end
            end
        end % curve seg
    end
    iter = iter + 1;
    ratio = ratio * 0.8;
end
% if the program proceeds here, meaning the gap is not changed,
% now try move only segments instead of taking the whole curve for editing
% 1. find curve intersection point to split the curve

iter = 0;
ratio = 1.0;
while iter < max_iter % change offset length
    for i = 1 : nCurves
        nseg = size(feaIds{i}, 1);
        deformed{i} = curves{i};
        % for each curve segment
        for s = 1 : nseg
            activePids = feaIds{i}(s, 1) : feaIds{i}(s, 3);
            if length(activePids) < nthr
                continue;
            end
            % jump end points
            s1 = activePids(1) + n2static;
            s2 = activePids(length(activePids) - n2static);
            static_anchor = zeros(1, n2static * 2);
            for k = 1 : n2static
                static_anchor(k) = k;
                static_anchor(k + n2static) = length(activePids) - n2static + k;
            end
            %% only this curve segment
            sCurve = curves{i}(activePids, :);
            % points to compute
            ipoints = curves{i}(s1 : s2, :);
            max_off_s = 0;
            max_offset_s = [0, 0];
            max_handleId = floor((s1 + s2)/2);
            for ip = 1 : size(ipoints, 1)
                pi = ipoints(ip, :);
                % find nearest neighboring point
                min_dij = intmax('int64');
                min_offset_ij = [0, 0];
                for j = 1 : nCurves
                    if i == j
                        continue;
                    end
                    jpoints = curves{j};
                    for jp = 1 : size(jpoints, 1)
                        pj = jpoints(jp, :);
                        d = norm(pj - pi);
                        dir = (pj - pi) / d;
                        if d < min_dij
                            min_dij = d;
                            min_offset_ij = d / 2 * dir;
                        end
                    end
                end % j curve
                if min_dij ~= intmax('int64') && min_dij > max_off_s
                    max_off_s = min_dij;
                    max_offset_s = min_offset_ij;
                    max_handleId = ip;
                end
            end % point ip at i curve s segment
            % deform     
            if max_off_s ~= 0
                handle_anchor(1) = max_handleId;   
                offsets(1, :) = max_offset_s * ratio;
                prevCurve = curves{i};
                curves{i}(activePids, :) = lap2D(sCurve, static_anchor, handle_anchor, offsets);
                [~, cur_gap, cur_overlap] = compute_gap_overlap_area(curves, Trunk);
                if cur_overlap > overlap_thr || cur_gap >= area_gap
                    curves{i} = prevCurve;
                    curves{i}(activePids, :) = translateCurve(sCurve, max_offset_s);
                    [~, cur_gap, cur_overlap] = compute_gap_overlap_area(curves, Trunk);
                end
                if cur_overlap > overlap_thr || cur_gap >= area_gap
                    curves{i} = prevCurve; %unchanged
                else
                    deformed{i}(s1 : s2, :) = curves{i}(s1 : s2, :);
                    area_overlap = cur_overlap;
                    area_gap = cur_gap;
                    return;
                end
            end
        end % curve seg
    end
    iter = iter + 1;
    ratio = ratio * 0.8;
end
end

function segIds = findCurveSegmentIds(curves)
n = length(curves);
segIds = cell(1, n);
thr = 1e-3;
nthr = 5;
ids = zeros(1, n);
for i = 1 : n
    segIds{i} = zeros(1, length(curves{i}));
    ids(i) = 0;
end
for i = 1 : n - 1
    iCurve = curves{i};
    for j = i + 1 : n
        jCurve = curves{j};
        % find the contact
        for p = 1 : length(iCurve)
            for q = 1 : length(jCurve)
                d = norm(iCurve(p, :) - jCurve(q, :));
                if d < thr
                    ids(i) = ids(i) + 1;
                    segIds{i}(ids(i)) = p;
                    ids(j) = ids(j) + 1;
                    segIds{j}(ids(j)) = q;
                    break;
                end
            end
        end
    end
end
for i = 1 : n
    segs = segIds{i}(1, 1 : ids(i));
    sort(segs);
    tmp = zeros(1, length(segs));
    tmp(1) = segs(1);
    id = 1;
    for p = 1 : length(segs)
        if segs(p) - tmp(id) > nthr
            id = id + 1;
            tmp(id) = segs(p);
        end
    end
    segIds{i} = tmp(1, 1 : id);
end
end

function deformed = translateCurve(curve, offset)
deformed = zeros(size(curve));
for i = 1 : length(curve)
    deformed(i, :) = curve(i, :) + offset;
end
end