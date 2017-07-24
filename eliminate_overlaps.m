function [deformed, area_gap, area_overlap] =  eliminate_overlaps(curves, Trunk)
%% detect any curve segments that are along a gap region
%  

nCurves = length(curves);
deformed = cell(1, nCurves);
nh = 1;
handle_anchor = zeros(1, nh);
offsets = zeros(nh, 2);
iter = 1;
max_iter = 5;
n2static = 2;
min_dist_thr = 0.3;
dist_threshold = 0.01;
min_num_points_in_seg = 10;
[~, area_gap, area_overlap] = compute_gap_overlap_area(curves, Trunk);
ori_overlap = area_overlap;
% simplified polygon to check overlap
% polys = cell(1, nCurves);
% d_sample = 0.02;
% for i = 1 : nCurves
%     polys{i} = dpsimplify(curves{i}, d_sample);
% end
polys = curves;
while iter < max_iter % change offset length
    % 1. find curve intersection point to split the curve
    segIds = findCurveSegmentIds(curves, dist_threshold, min_num_points_in_seg);
    drawFeaturePoints(curves, segIds);
    for i = 1 : nCurves
        nseg = length(segIds{i}) - 1;
        deformed{i} = curves{i};
        % for each curve segment
        for s = 1 : nseg
            % take the previous and post points as static anchor
            st = segIds{i}(s);
            ed = segIds{i}(s + 1);
            activePids = st : ed;
            if length(activePids) < min_num_points_in_seg %n2static * 2
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
                min_dij = min_dist_thr;
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
                            inPoly = inpolygon(pi(1), pi(2), ...
                                    polys{j}(:, 1), polys{j}(:, 2));
                                if inPoly == 1
                                    inPoly = inpolygon(pj(1), pj(2), ...
                                    polys{i}(:, 1), polys{i}(:, 2));
                                end
                            if inPoly == 1
                                min_dij = d / 2;
                                min_offset_ij = d / 2 * dir;                          
                            end
                        end
                    end
                end % j curve
                if min_dij ~= min_dist_thr && min_dij > max_off_s
                    max_off_s = min_dij;
                    max_offset_s = min_offset_ij;
                    max_handleId = n2static + ip;
                end
            end % point ip at i curve s segment
            % deform     
            if max_off_s ~= 0
                handle_anchor(1) = max_handleId;   
                offsets(1, :) = max_offset_s;
                prevCurve = curves{i};
                curves{i}(activePids, :) = lap2D(sCurve, static_anchor, handle_anchor, offsets);
                [~, ~, cur_overlap] = compute_gap_overlap_area(curves, Trunk);
                if cur_overlap >= area_overlap
                    curves{i} = prevCurve;
                    translated = translateCurve(sCurve, max_offset_s);
                    s1 = activePids(2);
                    s2 = activePids(length(activePids) - 1);
                    curves{i}(s1 : s2, :) = translated(2 : length(translated) - 1, :);
                    [~, ~, cur_overlap] = compute_gap_overlap_area(curves, Trunk);
                end
                if cur_overlap >= area_overlap
                    curves{i} = prevCurve; %unchanged
                else
                    deformed{i}(s1 : s2, :) = curves{i}(s1 : s2, :);
                    area_overlap = cur_overlap;
                end
            end
        end % curve seg
    end
    if area_overlap < ori_overlap
        break;
    end
    iter = iter + 1;
    dist_threshold = dist_threshold / 2;
    min_num_points_in_seg = min_num_points_in_seg * 0.8;
%     d_sample = d_sample / 2;
%     for i = 1 : nCurves
%         polys{i} = dpsimplify(curves{i}, d_sample);
%     end
end
end