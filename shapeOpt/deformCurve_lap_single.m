function deformed = deformCurve_lap_single(curves, Trunk, feaIds, mode)
% mode:
% 0 - eliminate overlaps
% 1 - eliminate gaps
% 2 - translate only

% n pieces of curves
nCurves = length(curves);
maxd = 0;
for i = 1 : length(Trunk)
    for j = 1 : length(Trunk)
        d = norm(Trunk(i, :) - Trunk(j, :));
        if d > maxd
            maxd = d;
        end
    end
end

%% compute offset
nthr = 6;
n2static = 3;
% analyze each curve
deformed = cell(1, nCurves);
nh = 1;
handle_anchor = zeros(1, nh);
offsets = zeros(nh, 2);
selectedCid = 0;
selectedSeg = 0;
max_off_s = 0;
max_offset_s = [0, 0];
max_handleId = 0;
% simplified polygon to check overlap
polys = cell(1, nCurves);
for i = 1 : nCurves
    polys{i} = dpsimplify(curves{i}, 0.02);
end
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
%         max_off_s = 0;
%         max_offset_s = [0, 0];
%         max_handleId = floor((s1 + s2)/2);
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
                        % check if it is overlap or gap
                        md = 1;
                        if mode == 0 || mode == 1
                            in = inpolygon(pi(1), pi(2), polys{j}(:, 1), polys{j}(:, 2));
                            md = in;
                            if mode == 1
                                md = 1 - in;
                            end
                        end
                        if md ~= 0
                            min_dij = d;
                            if mode == 2
                                min_offset_ij = d * dir;
                            else 
                                min_offset_ij = d / 2 * dir;
                            end
                        end
                    end
                end
            end % j curve
            if min_dij ~= intmax('int64') && min_dij > max_off_s
                max_off_s = min_dij;
                max_offset_s = min_offset_ij;
                max_handleId = feaIds{i}(s, 1) + ip - 1;
                selectedCid = i;
                selectedSeg = s;
            end
        end % point ip at i curve s segment
    end % curve seg
end
% deform a curve at a time
handle_anchor(1) = max_handleId;   
offsets(1, :) = max_offset_s;
static_anchor = zeros(1, n2static * 2);
segIds = feaIds{selectedCid}(selectedSeg, 1) : feaIds{selectedCid}(selectedSeg, 3);
s1 = segIds(1) + n2static;
s2 = segIds(length(segIds) - n2static);
for k = 1 : n2static
    static_anchor(k) = segIds(k);
    static_anchor(k + n2static) = segIds(length(segIds) - n2static + k);
end
iDeformed = lap2D(curves{selectedCid}, static_anchor, handle_anchor, offsets);
% update
deformed{selectedCid}(s1 : s2, :) = iDeformed(s1 : s2, :);
end