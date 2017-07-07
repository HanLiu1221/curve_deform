function deformed = deformCurve_lap(curve, Trunk, feaIds, mode)

% n pieces of curves
nCurves = length(curve);
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
R0 = 0.08;
thr = 0.001;
nthr = 6;
n2static = 3;
% analyze each curve
deformed = cell(1, nCurves);
for i = 1 : nCurves
    nseg = size(feaIds{i}, 1);
    deformed{i} = curve{i};
    currCurve = curve{i};
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
        ipoints = curve{i}(s1 : s2, :);
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
                jpoints = curve{j};
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
            if min_dij > max_off_s
                max_off_s = min_dij;
                max_offset_s = min_offset_ij;
                max_handleId = feaIds{i}(s, 1) + ip - 1;
            end
        end % point ip at i curve s segment
        % deform
        nh = 1;
        handle_anchor = zeros(1, nh);
        offsets = zeros(nh, 2);
        handle_anchor(1) = max_handleId;   
        offsets(1, :) = max_offset_s;
        iDeformed = lap2D(currCurve, static_anchor, handle_anchor, offsets);
        deformed{i}(s1 : s2, :) = iDeformed(s1 : s2, :);
    end % curve seg
end

end



function pnts = samplePntsOnLine(s, e, n_sample)
pnts = zeros(n_sample, 2);
d = norm(e - s);
dir = (e - s) / d;
unit = d / (n_sample - 1);
for i = 1 : n_sample
    pnts(i, :) = s + (i - 1) * unit * dir;
end
end