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
    % make a closed curve for lap neighboring
%     npnts = size(curve{i}, 1);
%     iline = samplePntsOnLine(curve{i}(npnts, :), curve{i}(1, :), npnts);
    %closedCurve = [curve{i}; iline];
    closedCurve = curve{i};
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
        offset = zeros(1, 2);
        nneig = 0;
        maxoff = 0;
        max_offset = [0, 0];
        max_handleId = floor((s1 + s2)/2);
        for j = 1 : nCurves
            if i == j
                continue;
            end
            jpoints = curve{j};
            for ip = 1 : size(ipoints, 1)
                pi = ipoints(ip, :);
                for jp = 1 : size(jpoints, 1)
                    pj = jpoints(jp, :);
                    d = norm(pj - pi);
                    dir = (pj - pi) / d;
                    if d < R0 && d > thr
                        d = d / maxd;
                        off = d * dir;
                        offset = offset + off;
                        nneig = nneig + 1;
                        if d > maxoff
                            maxoff = d;
                            max_handleId = feaIds{i}(s) + ip - 1;
                            max_offset = off;
                        end
                    end
                end % j curve
            end % i curve s segment
        end % compute offset from other curves 
        offset = offset / nneig;
        % deform
        nh = 1;
        handle_anchor = zeros(1, nh);
        offsets = zeros(nh, 2);
        if mode == 0
            handle_anchor(1) = feaIds{i}(s, 2);   
            offsets(1, :) = offset;
        else 
            handle_anchor(1) = max_handleId;
            offsets(1, :) = max_offset;
        end
        iDeformed = lap2D(closedCurve, static_anchor, handle_anchor, offsets);
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