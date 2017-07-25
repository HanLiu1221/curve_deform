function [deformed, area_gap, area_overlap] =  eliminate_gaps(curves, Trunk, overlap_thr)
%% detect any curve segments that are along a gap region
%  

nCurves = length(curves);
deformed = cell(1, nCurves);
nh = 1;
handle_anchor = zeros(1, nh);
offsets = zeros(nh, 2);
iter = 1;
max_iter = 6;
n2static = 2;
dist_threshold = 0.01;
min_num_points_in_seg = 10;
[~, area_gap, area_overlap] = compute_gap_overlap_area(curves, Trunk);
ori_gap = area_gap;
ratio = 1.0;
while iter < max_iter % change offset length
    % 1. find curve intersection point to split the curve
    segIds = findCurveSegmentIds(curves, dist_threshold, min_num_points_in_seg);
    for i = 1 : nCurves
        nseg = length(segIds{i}) - 1;
        deformed{i} = curves{i};
        % for each curve segment
        for s = 1 : nseg
            % take the previous and post points as static anchor
            st = segIds{i}(s);
            ed = segIds{i}(s + 1);
            activePids = st : ed;
            if length(activePids) < n2static * 2
                continue;
            end
            % jump end points
            s1 = activePids(n2static + 1);
            s2 = activePids(length(activePids) - n2static);
            static_anchor = zeros(1, n2static * 2);
            for k = 1 : n2static
                static_anchor(k) = k;
                static_anchor(k + n2static) = length(activePids) - n2static + k;
            end
%             for k = 1 : n2static
%                 static_anchor(k) = activePids(k);
%                 static_anchor(k + n2static) = activePids(length(activePids) - n2static + k);
%             end
            %% only this curve segment
            sCurve = curves{i}(activePids, :);
            % points to compute
            ipoints = curves{i}(s1 : s2, :);
            max_off_s = 0;
            max_offset_s = [0, 0];
            max_handleId = length(activePids) / 2;
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
                            min_dij = d / 2 * ratio;
                            min_offset_ij = min_dij * dir;
                        end
                    end
                end % j curve
                if min_dij ~= intmax('int64') && min_dij > max_off_s
                    max_off_s = min_dij;
                    max_offset_s = min_offset_ij;
                    %max_handleId = s1 + ip - 1;
                    max_handleId = n2static + ip;
                end
            end % point ip at i curve s segment
            % deform     
            if max_off_s ~= 0
                handle_anchor(1) = max_handleId;   
                offsets(1, :) = max_offset_s;
                prevCurve = curves{i};
                %curves{i} = lap2D(curves{i}, static_anchor, handle_anchor, offsets);
                curves{i}(activePids, :) = lap2D(sCurve, static_anchor, handle_anchor, offsets);
                [~, cur_gap, cur_overlap] = compute_gap_overlap_area(curves, Trunk);
                if cur_overlap > overlap_thr || cur_gap >= area_gap
                    translated = translateCurve(sCurve, max_offset_s);
                    curves{i}(s1 : s2, :) = translated(n2static + 1 : length(activePids) - n2static,:);
                    [~, cur_gap, cur_overlap] = compute_gap_overlap_area(curves, Trunk);
                end
                if cur_overlap > overlap_thr || cur_gap >= area_gap
                    curves{i} = prevCurve; %unchanged
                else
                    deformed{i}(s1 : s2, :) = curves{i}(s1 : s2, :);
%                     area_overlap = cur_overlap;
%                     area_gap = cur_gap;
                end
            end
        end % curve seg
    end
    [~, area_gap, area_overlap] = compute_gap_overlap_area(curves, Trunk);
    if area_overlap <= overlap_thr && area_gap < ori_gap
        break;
    end
    iter = iter + 1;
    dist_threshold = dist_threshold / 2;
    %min_num_points_in_seg = min_num_points_in_seg * 0.8;
    ratio = ratio * 0.8;
end
end
