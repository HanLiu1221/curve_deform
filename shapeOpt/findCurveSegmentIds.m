function segIds = findCurveSegmentIds(curves, dist_threshold, min_num_points_in_seg)
n = length(curves);
segIds = cell(1, n);
ids = zeros(1, n);
for i = 1 : n
    segIds{i} = zeros(1, length(curves{i}));
    ids(i) = 0;
end
for i = 1 : n - 1
    iCurve = curves{i};
    for j = i + 1 : n
        jCurve = curves{j};
        % find the contacts between any two curves
        % there may exist more than one contact between two curves
        for p = 1 : length(iCurve)
            for q = 1 : length(jCurve)
                d = norm(iCurve(p, :) - jCurve(q, :));
                if d < dist_threshold
                    ids(i) = ids(i) + 1;
                    segIds{i}(ids(i)) = p;
                    ids(j) = ids(j) + 1;
                    segIds{j}(ids(j)) = q;
                    %break;
                end
            end
        end
    end
end
for i = 1 : n
    segs = segIds{i}(1, 1 : ids(i));
    segs = sort(segs);
    tmp = zeros(1, length(segs) + 2);
    tmp(1) = segs(1);
    id = 1;
    for p = 2 : length(segs)
        if segs(p) - tmp(id) >= min_num_points_in_seg% * 2
            id = id + 1;
            tmp(id) = segs(p);
        end
    end
    id = id + 1;
    tmp(id) = length(curves{i});
    segIds{i} = tmp(1, 1 : id);
end
end