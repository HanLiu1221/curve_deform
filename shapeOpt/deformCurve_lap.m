function deformed = deformCurve_lap(curve, Trunk, weights)

% n pieces of curves
nCurves = length(curve);

%% find feature points
feaIds = cell(1, nCurves);
figure;
for i = 1:nCurves
    simp = dpsimplify(curve{i}, 0.02);
    pids = findPointIndices(curve{i}, simp);
    feaIds{i} = pids;
    % draw
    plot(curve{i}(:,1), curve{i}(:,2), 'k-');
    hold on
    plot(simp(:,1), simp(:,2), 'r--o', 'LineWidth',2);
    hold on
end
legend('original polyline', 'simplified');

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
R0 = 0.1;
thr = 0.01;
nthr = 10;
n2static = 2;
% analyze each curve
deformed = cell(1, nCurves);
for i = 1 : nCurves
    nseg = length(feaIds{i}) - 1;
    deformed{i} = curve{i};
    % for each curve segment
    for s = 1 : nseg
        activePids = feaIds{i}(s) : feaIds{i}(s + 1);
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
        handleId = ceil((s1 + s2) / 2);
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
                            handleId = feaIds{i}(s) + ip - 1;
                            %offset = off;
                        end
                    end
                end % j curve
            end % i curve s segment
        end % compute offset from other curves 
        offset = offset / nneig;
        % deform
        nh = 1;
        handle_anchor = zeros(1, nh);
        handle_anchor(1) = handleId;
        offsets = zeros(nh, 2);
        offsets(1, :) = offset;
        iDeformed = lap2D(curve{i}, static_anchor, handle_anchor, offsets);
        deformed{i}(s1 : s2, :) = iDeformed(s1 : s2, :);
    end % curve seg
end

end

function ids = findPointIndices(curve, poly) 
ids = ones(size(poly, 1), 1);
for i = 1 : size(poly, 1)
    mind = 100000;
    for j = 1 : size(curve, 1)
        d = norm(curve(j, :) - poly(i, :));
        if d < mind
            mind = d;
            ids(i) = j;
        end
    end
end
end