function deformed = deformCurve_lap_overlap(curves, Trunk, feaIds, mode)
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
n2static = 2;
min_dist_thr = 0.5;
% analyze each curve
deformed = cell(1, nCurves);
nh = 1;
handle_anchor = zeros(1, nh);
offsets = zeros(nh, 2);
% simplified polygon to check overlap
polys = cell(1, nCurves);

d_sample = 0.01;
for i = 1 : nCurves
    if mode == 1
        polys{i} = dpsimplify(curves{i}, d_sample);
    else
         polys{i} = curves{i};
    end
end

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
        max_off_val = 0;
        max_off_dir = [0, 0];
        max_handleId = floor((s1 + s2)/2);
        for ip = 1 : size(ipoints, 1)
            pi = ipoints(ip, :);
            % find nearest neighboring point
            min_dij = min_dist_thr;
            min_dir = [0, 0];
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
                        inPoly = inpolygon(pi(1), pi(2), ...
                                    polys{j}(:, 1), polys{j}(:, 2));
                        if inPoly == 1
                            min_dir = dir;
                            if mode == 2
                                min_dij = d / 2;
                            else 
                                min_dij = d / 2;
                            end
                        end
                    end
                end
            end % all j curve
            if min_dij ~= min_dist_thr && min_dij > max_off_val
                max_off_val = min_dij;
                max_off_dir = min_dir;
                max_handleId = s1 + ip - 1;
            end
        end % point ip at i curve s segment
        % laplacian deformation 
        
        if max_off_val ~= 0
            handle_anchor(1) = max_handleId;  
            offsets(1, :) = max_off_val * max_off_dir;
            if mode == 1
                iDeformed = lap2D(curves{i}, static_anchor, handle_anchor, offsets);
            else
                iDeformed = translateCurve(curves{i},  offsets(1, :));
            end
            deformed{i}(s1 : s2, :) = iDeformed(s1 : s2, :);
        end
%         itest = 0;
%         while max_off_val ~= 0 && itest < 5 
%             offsets(1, :) = max_off_dir * max_off_val;
%             if mode == 1
%                 iDeformed = lap2D(curves{i}, static_anchor, handle_anchor, offsets);
%             else
%                 iDeformed = translateCurve(curves{i}, max_off_dir);
%             end
%             if ifOutOfPoly(iDeformed, Trunk) == 0
%                 deformed{i}(s1 : s2, :) = iDeformed(s1 : s2, :);
%                 break;
%             end
%             max_off_val = max_off_val / 2;
%             itest = itest + 1;
%         end
    end % curve seg
end
end

function outOfPoly = ifOutOfPoly(seg, T)
nOut = 0;
for i = 1 : length(seg)
    inPoly = inpolygon(seg(i, 1), seg(i, 2), T(:, 1), T(:, 2));
    if inPoly == 0
        nOut = nOut + 1;
    end
end
if nOut < length(seg) * 0.1
    outOfPoly = 0;
else 
    outOfPoly = 1;
end
end