function deformed = deformCurve_lap_overlap_with_handle(curves, Trunk, feaIds, handleIds, mode)
%%
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
n2static = 2;
% analyze each curve
deformed = cell(1, nCurves);
nh = 1;
handle_anchor = zeros(1, nh);
offsets = zeros(nh, 2);
% simplified polygon to check overlap
polys = cell(1, nCurves);

if mode == 1
    d_sample = 0.02;
    for i = 1 : nCurves
    polys{i} = dpsimplify(curves{i}, 0.02);
    end
else % mode == 2
    polys = curves;
end

for i = 1 : nCurves
    nseg = size(feaIds{i}, 1);
    deformed{i} = curves{i};
    % for each curve segment
    for s = 1 : nseg
        % points & feature ids
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
        %% handle polylines or curves
        isLine = 0;
        if length(handleIds{i}) > 0
            isLine = 1;
        end
        % points to compute
        if isLine == 0
            ipoints = curves{i}(s1 : s2, :);
        else
            ipoints =  curves{i}(handleIds{i}, :);
        end
        [max_off_s, max_offset_s, handePointId] = computeOffset(curves, i, ipoints, polys, mode);
        if isLine == 1 && max_off_s == 0
            ipoints = curves{i}(s1 : s2, :);
            [max_off_s, max_offset_s, anyId] = computeOffset(curves, i, ipoints, polys, mode);
        end
        if isLine == 0
            max_handleId = feaIds{i}(s, 1) + handePointId - 1;
        else
            % poly lines, only for one case: the whole curve is a segment
            max_handleId = feaIds{i}(s, 1) + handleIds{i}(handePointId) - 1;
            % static anchor: previous and post handle ids
            if length(handleIds{i}) > 1
                static_anchor = zeros(1, 2);
                if handePointId > 1
                    static_anchor(1) = feaIds{i}(s, 1) + handleIds{i}(handePointId - 1) - 1;
                else
                    static_anchor(1) = feaIds{i}(s, 1);
                end
                if handePointId + 1 <= length(handleIds{i})
                    static_anchor(2) = feaIds{i}(s, 1) + handleIds{i}(handePointId + 1) - 1;
                else
                    static_anchor(2) = feaIds{i}(s, 3);
                end
            end
        end
        % laplacian deformation     
        if max_off_s ~= 0
            handle_anchor(1) = max_handleId;   
            offsets(1, :) = max_offset_s;
            if mode == 0 || mode == 1
                iDeformed = lap2D(curves{i}, static_anchor, handle_anchor, offsets);
            else
                iDeformed = translateCurve(curves{i}, max_offset_s);
            end
            deformed{i}(s1 : s2, :) = iDeformed(s1 : s2, :);
        end
    end % curve seg
end
end

function [max_off_s, max_offset_s, max_handleId] = computeOffset(curves, iCurv, ipoints, polys, mode)
max_off_s = 0;
max_offset_s = [0, 0];
max_handleId = length(ipoints) / 2;
for ip = 1 : size(ipoints, 1)
    pi = ipoints(ip, :);
    % find nearest neighboring point
    % find nearest neighboring point
    min_dij = intmax('int64');
    min_offset_ij = [0, 0];
    for j = 1 : length(curves)
        if iCurv == j
            continue;
        end
        jpoints = curves{j};
        for jp = 1 : size(jpoints, 1)
            pj = jpoints(jp, :);
            d = norm(pj - pi);
            dir = (pj - pi) / d;
            if d < min_dij
                % check if it is overlap or gap
                in = inpolygon(pi(1), pi(2), polys{j}(:, 1), polys{j}(:, 2));
                if in == 1
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
        max_handleId = ip;
    end
end % point ip at i curve s segment
end

function deformed = translateCurve(curve, offset)
deformed = zeros(size(curve));
for i = 1 : length(curve)
    deformed(i, :) = curve(i, :) + offset;
end
end