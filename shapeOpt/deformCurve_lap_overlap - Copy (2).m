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
min_dist_thr = 0.2;
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
    center = [0, 0];
    for j = 1 : length(curves{i})
        center = center + curves{i}(j, :);
    end
    center = center / length(curves{i});
    for s = 1 : nseg
        activePids = feaIds{i}(s, 1) : feaIds{i}(s, 3);
        if length(activePids) <= n2static * 2
            continue;
        end
        % jump end points
        t1 = 1 + n2static;
        t2 = length(activePids) - n2static;
        s1 = activePids(t1);
        s2 = activePids(t2);
        static_anchor = zeros(1, n2static * 2);
        for k = 1 : n2static
            static_anchor(k) = activePids(k);
            static_anchor(k + n2static) = activePids(length(activePids) - n2static + k);
        end
        % points to compute
        ipoints = curves{i}(activePids, :);%(s1 : s2, :);
        max_off_s = 0;
        max_offset_s = [0, 0];
        max_handleId = floor((s1 + s2)/2);
        for ip = t1 : 1 : t2
            pi = ipoints(ip, :);
            % find nearest neighboring point
            min_dij = min_dist_thr;
            min_offset_ij = [0, 0];
            % try normal as the direction
            ldir = ipoints(ip + 2, :) - ipoints(ip - 2, :);
            ndir = [-ldir(2), ldir(1)];
            ndir = ndir / norm(ndir);
            % overlap, the move dir should point into the region
            cdir = pi - center;
            cdir = cdir / norm(cdir);
            if dot(cdir, ndir) > 0
                ndir = -1.0 * ndir;
            end
            avg_off = [0, 0];
            n_off = 0;
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
                        %dir = ndir;
                        % check if it is overlap or gap
                        inPoly = inpolygon(pi(1), pi(2), ...
                                    polys{j}(:, 1), polys{j}(:, 2));
                        if inPoly == 1
                            if mode == 2
                                min_offset_ij = d * dir;
                                min_dij = d;
                            else 
                                min_offset_ij = d / 2 * dir;
                                min_dij = d / 2;
                            end
                            avg_off = avg_off + (pj - pi);
                            n_off = n_off + 1;
                        end
                    end
                end
            end % all j curve for ip point
            if min_dij ~= min_dist_thr && min_dij > max_off_s
                max_off_s = min_dij;
                max_offset_s = min_offset_ij;
                max_handleId = activePids(ip);
            end
        end % point ip at i curve s segment
        % laplacian deformation     
        if max_off_s ~= 0
            handle_anchor(1) = max_handleId;   
            offsets(1, :) = max_offset_s;
            if mode == 1
                iDeformed = lap2D(curves{i}, static_anchor, handle_anchor, offsets);
            else
                iDeformed = translateCurve(curves{i}, max_offset_s);
            end
            deformed{i}(s1 : s2, :) = iDeformed(s1 : s2, :);
        end
    end % curve seg
end
end