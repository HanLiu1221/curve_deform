function [deformed, area_gap, area_overlap] =  eliminate_gaps(curves, Trunk)

deformed = cell(1, nCurves);
iter = 0;
ratio = 1.0;
ori_gap = area_gap;
max_iter = 5;
dist_threshold = 0.01;
min_num_points_in_seg = 10;
while iter < max_iter % change offset length
    % 1. find curve intersection point to split the curve
    segIds = findCurveSegmentIds(curves, dist_threshold, min_num_points_in_seg);
    for i = 1 : nCurves
        nseg = length(segIds{i}) - 1;
        deformed{i} = curves{i};
        % for each curve segment
        for s = 1 : nseg
            activePids = segIds{i}(s) : segIds{i}(s + 1);
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
                end
            end
        end % curve seg
    end
    if cur_overlap <= overlap_thr && cur_gap < ori_gap
        break;
    end
    iter = iter + 1;
    dist_threshold = dist_threshold / 2;
    min_num_points_in_seg = min_num_points_in_seg * 0.8;
end
end

function segIds = findCurveSegmentIds(curves, dist_threshold, min_num_points_in_seg)
n = length(curves);
segIds = cell(1, n);
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
                if d < dist_threshold
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
    segs = sort(segs);
    tmp = zeros(1, length(segs) + 2);
    tmp(1) = segs(1);
    id = 1;
    for p = 2 : length(segs)
        if segs(p) - tmp(id) > min_num_points_in_seg
            id = id + 1;
            tmp(id) = segs(p);
        end
    end
    if length(curves{i}) - tmp(id) > min_num_points_in_seg
        id = id + 1;
        tmp(id) = length(curves{i}) - min_num_points_in_seg;
    end
    segIds{i} = tmp(1, 1 : id);
end
end