function optimizeBoundaryCurve(P, T_P, T_Q, scale, folder, reverse)

clc;
close all;
addpath(genpath('./'));

%% 1. sample points on the boundary curve
n = length(P);
d_sample = 0.04;
ctrlPnts = cell(n,1);
for i = 1:n
    [ctrlPnts{i},idx_CPs] = extract_control_points(P{i}, d_sample);   
end
% visualize
scale = 0;
oriControlPs = ctrlPnts;
tran_oriControlPs = transform_curves(ctrlPnts, T_Q, scale);
% show_curves(oriControlPs, tran_oriControlPs, 1);

%% 2. handle pieces that go beyond the trunk
controlPs = remove_outside_regions(ctrlPnts, T_P);
tran_controlPs = transform_curves(controlPs, T_Q, scale);
show_curves(controlPs, tran_controlPs, 0);

originalShapeFile = strcat(folder, 'iter_0.png');
saveas(gcf, originalShapeFile);

%% 3. iteratively eliminate overlaps
curves = controlPs;
for i = 1 : length(curves)
    curves{i} = removeReps(curves{i});
end
% 3.1 find feature points
nCurves = length(curves);
featIds_noSplit = cell(1, nCurves);
feaIds = cell(1, nCurves);
loaded = 0;
%featfile = 'data/featurePointIds.mat';
% if exist(featfile, 'file') == 2
%     load(featfile);
%     loaded = 1;
%     feaIds = featurePointIds;
% end
figure;
for i = 1:nCurves
    if loaded == 0
        feaIds{i} = getSplitCurvePointIds(curves{i}, reverse);
    end
    npnts = length(curves{i});
    featIds_noSplit{i} = [1, floor((1 + npnts) / 2), npnts];
    % draw
    plot(curves{i}(:,1), curves{i}(:,2), 'k-');
    hold on;
    % split points
    splitIds = [feaIds{i}(:, 1)', npnts];
    plot(curves{i}(splitIds,1), curves{i}(splitIds,2), 'r--o', 'LineWidth',2);
    hold on;
    for j = 1:length(splitIds)
        text(curves{i}(splitIds(j), 1), curves{i}(splitIds(j), 2), num2str(splitIds(j)));
        hold on;
    end
end
legend('original polyline', 'simplified');

cruveSegmentFile = strcat(folder, 'curve_seg.png');
saveas(gcf, cruveSegmentFile);

% 3.2 define a local region, using r, to compute the attraction force
[~, area_gap, area_overlap] = compute_gap_overlap_area(controlPs, T_P);
iter = 0;
overlap_thr = 1e-6;
max_iter = 12;
disp('===Start eleminating overlaps===');
total_iter = 0;
while area_overlap > overlap_thr
    if iter >= max_iter
        break;
    end
    %tic;
    prev = curves;
    if area_overlap > 0.005 
        curves = deformCurve_lap_simu(curves, T_P, featIds_noSplit, 0);
    elseif area_overlap > 0.0001
        curves = deformCurve_lap_simu(curves, T_P, feaIds, 0);
    else
        curves = deformCurve_lap_simu(curves, T_P, feaIds, 2);
    end
    %toc
    [~, area_gap, area_overlap] = compute_gap_overlap_area(curves, T_P);
    deformEnergy = 0; % getDeformationEnergy(curves, prev);
    iter = iter + 1;   
    total_iter = total_iter + 1;
    % visualize
    tran_curve = transform_curves(curves, T_Q, scale);
    show_curves(curves, tran_curve, 1);
    iterFile = strcat(folder, 'Iter_', num2str(total_iter), '_overlap', '.png');
    str = strcat('Iter ',  num2str(iter), ...
        ': overlap-area: ', num2str(area_overlap), ...
        ', gap-area: ', num2str(area_gap));
    text(-1, 1, str);
    saveas(gcf, iterFile);
    disp(strcat('===Iter', num2str(iter), '==='));
    disp(strcat('area of overlap: ', num2str(area_overlap)));
    disp(strcat('area of gap: ', num2str(area_gap)));
    disp(strcat('deformation energy: ', num2str(deformEnergy)));
end

%% 4. iteratively eliminate gaps
iter = 0;
gap_thr = 0.001;
disp('===Try to eliminate gaps by deforming each cuve===');
prev_gap = area_gap;
while area_gap > gap_thr && iter < max_iter
    [curves, area_gap, area_overlap] = ...
        deformCurve_lap_gap(curves, T_P, feaIds, overlap_thr);
    deformEnergy = 0; % getDeformationEnergy(curves, prev);
    iter = iter + 1;    
    if (area_overlap > overlap_thr)
        disp('Cause overlaps.');
        break;
    end
    if area_gap >= prev_gap
        disp('Gaps was not changed.');
        break;
    end
    % visualize
    tran_curve = transform_curves(curves, T_Q, scale);
    show_curves(curves, tran_curve, 1);
    total_iter = total_iter + 1;
    iterFile = strcat(folder, 'Iter_', num2str(total_iter), '_global_gap', '.png');
    str = strcat('Iter ',  num2str(iter), ...
        ': overlap-area: ', num2str(area_overlap), ...
        ', gap-area: ', num2str(area_gap));
    text(-1, 1, str);
    saveas(gcf, iterFile);
    disp(strcat('===Iter', num2str(iter), '==='));
    disp(strcat('area of overlap: ', num2str(area_overlap)));
    disp(strcat('area of gap: ', num2str(area_gap)));
    disp(strcat('deformation energy: ', num2str(deformEnergy)));
    prev_gap = area_gap;
end
disp('===Start eleminating gaps according to gap regions===');
iter = 0;
while area_gap > gap_thr && iter < max_iter
    [curves, area_gap, area_overlap] = eliminate_gaps(curves, T_P, overlap_thr);
    deformEnergy = 0; % getDeformationEnergy(curves, prev);
    iter = iter + 1;    
    if (area_overlap > overlap_thr)
        disp('Cause overlaps.');
        break;
    end
    if (area_gap >= prev_gap || prev_gap - area_gap < 1e-4)
        disp('Gaps was not changed.');
        break;
    end
    % visualize
    tran_curve = transform_curves(curves, T_Q, scale);
    show_curves(curves, tran_curve, 1);
    total_iter = total_iter + 1;
    iterFile = strcat(folder, 'Iter_', num2str(total_iter), '_gap', '.png');
    str = strcat('Iter ',  num2str(iter), ...
        ': overlap-area: ', num2str(area_overlap), ...
        ', gap-area: ', num2str(area_gap));
    text(-1, 1, str);
    saveas(gcf, iterFile);
    disp(strcat('===Iter', num2str(iter), '==='));
    disp(strcat('area of overlap: ', num2str(area_overlap)));
    disp(strcat('area of gap: ', num2str(area_gap)));
    disp(strcat('deformation energy: ', num2str(deformEnergy)));
    prev_gap = area_gap;
end

end

function cleanedPnts = removeReps(pnts)
tmp = pnts;
n = 1;
for i = 2 : length(pnts)
    d = norm(pnts(i, :) - tmp(n, :));
    if (d > 1e-3)
        tmp(n + 1, :) = pnts(i, :);
        n = n + 1;
    end
end
cleanedPnts = tmp(1 : n, :);
end