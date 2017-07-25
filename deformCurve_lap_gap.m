function [deformed, area_gap, area_overlap] = deformCurve_lap_gap(curves, Trunk, feaIds, overlap_thr)
%%

% n pieces of curves
nCurves = length(curves);

%% compute offset
n2static = 2;
% analyze each curve
nh = 1;
handle_anchor = zeros(1, nh);
offsets = zeros(nh, 2);
deformed = curves;
%% Iterate over each curve segment
[~, area_gap, area_overlap] = compute_gap_overlap_area(curves, Trunk);
for i = 1 : nCurves
    nseg = size(feaIds{i}, 1);
    deformed{i} = curves{i};
    % for each curve segment
    for s = 1 : nseg
        activePids = feaIds{i}(s, 1) : feaIds{i}(s, 3);
        if length(activePids) <= n2static * 2
            continue;
        end
        % jump end points
        s1 = activePids(1 + n2static);
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
        for ip = n2static : size(ipoints, 1) - n2static
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
                        min_dij = d / 2;
                        min_offset_ij = min_dij * dir;
                    end
                end
            end % all j curve
            if min_dij ~= intmax('int64') && min_dij > max_off_s
                max_off_s = min_dij;
                max_offset_s = min_offset_ij;
                max_handleId = s1 + ip - 1;
            end
        end % point ip at i curve s segment
        % deform     
        if max_off_s ~= 0
            handle_anchor(1) = max_handleId;   
            offsets(1, :) = max_offset_s;
            prevCurve = curves{i};
            curves{i} = lap2D(curves{i}, static_anchor, handle_anchor, offsets);
            [~, cur_gap, cur_overlap] = compute_gap_overlap_area(curves, Trunk);
%             if cur_overlap > overlap_thr || cur_gap >= area_gap
%                 curves{i} = prevCurve;
%                 curves{i}(s1 :s2, :) = translateCurve(curves{i}(s1 :s2, :), offsets(1, :));
%                 [~, cur_gap, cur_overlap] = compute_gap_overlap_area(curves, Trunk);
%             end
            if cur_overlap <= overlap_thr && cur_gap < area_gap
                deformed{i}(activePids, :) = curves{i}(activePids, :);
                %deformed{i}(s1 : s2, :) = curves(s1 : s2, :);
                area_gap = cur_gap;
                area_overlap = cur_overlap;
            else
                curves{i} = prevCurve;
            end
        end
    end % curve seg
end
end