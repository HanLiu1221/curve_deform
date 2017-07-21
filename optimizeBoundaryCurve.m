function [P_o, P_g] = optimizeBoundaryCurve(P, T_P, T_Q, folder, reverse)

clc;
close all;
addpath(genpath('./'));

%% 2. handle pieces that go beyond the trunk
scale = 0;
cP = remove_outside_regions(P, T_P);
tran_cP = transform_curves(cP, T_Q, scale);
show_curves(cP, tran_cP, 0);

originalShapeFile = strcat(folder, 'iter_0.png');
saveas(gcf, originalShapeFile);

%% 3. iteratively eliminate overlaps
curves = cP;
for i = 1 : length(curves)
    curves{i} = removeReps(curves{i});
end
% 3.1 find feature points
nCurves = length(curves);
featIds_noSplit = cell(1, nCurves);
feaIds = cell(1, nCurves);
handleIds = cell(1, nCurves);
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
        [feaIds{i}, handleIds{i}] = getSplitCurvePointIds(curves{i}, reverse);
    end
    npnts = length(curves{i});
    featIds_noSplit{i} = [1, floor((1 + npnts) / 2), npnts];
    % draw
    plot(curves{i}(:,1), curves{i}(:,2), 'k-');
    hold on;
    % split points
    splitIds = [feaIds{i}(:, 1)', npnts];
    plot(curves{i}(splitIds, 1), curves{i}(splitIds, 2), 'r--o', 'LineWidth',2);
    hold on;
    if length(handleIds{i}) > 0
        plot(curves{i}(handleIds{i}, 1), curves{i}(handleIds{i}, 2), 'm*', 'LineWidth',2);
        hold on;
    end
    for j = 1:length(splitIds)
        text(curves{i}(splitIds(j), 1), curves{i}(splitIds(j), 2), num2str(splitIds(j)));
        hold on;
    end
end
legend('original polyline', 'simplified');

cruveSegmentFile = strcat(folder, 'curve_seg.png');
saveas(gcf, cruveSegmentFile);

% 3.2 define a local region, using r, to compute the attraction force
[~, area_gap, area_overlap] = compute_gap_overlap_area(cP, T_P);
iter = 0;
overlap_thr = 1e-6;
max_iter = 10;
disp('===Try to eliminate overlaps by deforming each cuve===');
total_iter = 0;
prev_overlap = area_overlap;
while area_overlap > overlap_thr
    if iter >= max_iter
        break;
    end
    %tic;
    if area_overlap > 0.005
        curves = deformCurve_lap_overlap(curves, T_P, featIds_noSplit, 1);
    elseif area_overlap > 0.0001
        curves = deformCurve_lap_overlap(curves, T_P, feaIds, 1);
    else
        curves = deformCurve_lap_overlap(curves, T_P, feaIds, 2);
    end
    %toc
    [~, area_gap, area_overlap] = compute_gap_overlap_area(curves, T_P);
%     if area_overlap >= prev_overlap
%         disp('overlap  not changed.');
%         break;
%     end
    deformEnergy = 0; % getDeformationEnergy(curves, prev);
    iter = iter + 1;   
    total_iter = total_iter + 1;
    % visualize
    tran_curve = transform_curves(curves, T_Q, scale);
    show_curves(curves, tran_curve, 1);
    iterFile = strcat(folder, 'Iter_', num2str(total_iter), '_global_overlap', '.png');
    str = strcat('Iter ',  num2str(iter), ...
        ': overlap-area: ', num2str(area_overlap), ...
        ', gap-area: ', num2str(area_gap));
    text(-1, 1, str);
    saveas(gcf, iterFile);
    disp(strcat('===Iter', num2str(iter), '==='));
    disp(strcat('area of overlap: ', num2str(area_overlap)));
    disp(strcat('area of gap: ', num2str(area_gap)));
    disp(strcat('deformation energy: ', num2str(deformEnergy)));
    prev_overlap = area_overlap;
end

disp('===Start eleminating overlaps===');
iter = 0;
while area_overlap > overlap_thr
    if iter >= max_iter
        break;
    end
    [curves, area_gap, area_overlap] = eliminate_overlaps(curves, T_P, overlap_thr);
    if area_overlap >= prev_overlap
        disp('overlap  not changed.');
        break;
    end
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
    prev_overlap = area_overlap;
end

P_o = curves;
P_g = curves;

if area_overlap > 1e-3
    return;
end
if area_overlap > overlap_thr
    overlap_thr = area_overlap;
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
    if area_gap == prev_gap
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
%     if (area_gap >= prev_gap)% || prev_gap - area_gap < 1e-4)
%         disp('Gaps was not changed.');
%         break;
%     end
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
P_g = curves;

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