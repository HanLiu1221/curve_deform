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
loaded = 0;
%featfile = 'data/featurePointIds.mat';
% if exist(featfile, 'file') == 2
%     load(featfile);
%     loaded = 1;
%     feaIds = featurePointIds;
% end

featurePointIds = cell(1, nCurves);
handlePointIds = cell(1, nCurves);
for i = 1:nCurves
    if loaded == 0
        [feaIds{i}, handlePointIds{i}] = getSplitCurvePointIds(curves{i}, reverse);
    end
    npnts = length(curves{i});
    featIds_noSplit{i} = [1, floor((1 + npnts) / 2), npnts];
    % split points
    splitIds = [feaIds{i}(:, 1)', npnts];
    featurePointIds{i} = splitIds;
end
drawFeaturePoints(curves, featurePointIds);

cruveSegmentFile = strcat(folder, 'curve_seg.png');
saveas(gcf, cruveSegmentFile);

% 3.2 define a local region, using r, to compute the attraction force
[~, area_gap, area_overlap] = compute_gap_overlap_area(cP, T_P);
iter = 0;
overlap_thr = 1e-6;
max_iter = 10;
disp('***Start Overlap Elimination***');
disp('===Overlap Stage 1: Deforming each cuve in low or mid resolution===');
total_iter = 0;
prev_overlap = area_overlap;
useSeg = 0;
resolution = '_global';
while area_overlap > overlap_thr
    if iter >= max_iter
        break;
    end
    %tic;
    prevCurves = curves;
    if area_overlap > 0.005 && useSeg == 0
        curves = deformCurve_lap_overlap(curves, T_P, featIds_noSplit, 1);
        resolution = '_global';
    elseif area_overlap > 0.0001
        curves = deformCurve_lap_overlap(curves, T_P, feaIds, 1);
        resolution = '_mid';
    else
        curves = deformCurve_lap_overlap(curves, T_P, feaIds, 2);
    end
    %toc
    [~, area_gap, area_overlap] = compute_gap_overlap_area(curves, T_P);
    if area_overlap >= prev_overlap
        disp('overlap not changed.');
        curves = prevCurves;
        if useSeg == 1
            break;
        end
        useSeg = 1;
    end
    deformEnergy = 0; % getDeformationEnergy(curves, prev);
    iter = iter + 1;   
    total_iter = total_iter + 1;
    % visualize
    tran_curve = transform_curves(curves, T_Q, scale);
    show_curves(curves, tran_curve, 1);    
    iterFile = strcat(folder, 'Iter_', num2str(total_iter), resolution, '_overlap.png');
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

disp('===Overlap Stage 2: Deforming each cuve in high resolution===');
iter = 0;
resolution = '_local';
while area_overlap > overlap_thr
    if iter >= max_iter
        break;
    end
    prevCurves = curves;
    [curves, area_gap, area_overlap] = eliminate_overlaps(curves, T_P);
    if area_overlap >= prev_overlap
        curves = prevCurves;
        disp('overlap  not changed.');
        break;
    end
    iter = iter + 1;   
    total_iter = total_iter + 1;
    % visualize
    tran_curve = transform_curves(curves, T_Q, scale);
    show_curves(curves, tran_curve, 1);
    iterFile = strcat(folder, 'Iter_', num2str(total_iter), resolution, '_overlap.png');
    str = strcat('Iter ',  num2str(iter), ...
        ': overlap-area: ', num2str(area_overlap), ...
        ', gap-area: ', num2str(area_gap));
    text(-1, 1, str);
    saveas(gcf, iterFile);
    disp(strcat('===Iter', num2str(iter), '==='));
    disp(strcat('area of overlap: ', num2str(area_overlap)));
    disp(strcat('area of gap: ', num2str(area_gap)));
    %disp(strcat('deformation energy: ', num2str(deformEnergy)));
    prev_overlap = area_overlap;
end

P_o = curves;
P_g = curves;

if area_overlap > 1e-3
    return;
end
overlap_thr = max(overlap_thr, area_overlap);
%% 4. iteratively eliminate gaps
iter = 0;
gap_thr = 0.0001;
disp('***Start Gap Elimination***');
disp('===Gap Stage 1: Deforming each cuve in mid resolution===');
prev_gap = area_gap;
% global is useless, since it's already in a full state from overlap stage
resolution = '_mid'; 
while area_gap > gap_thr && iter < max_iter
    preCurves = curves;
    [curves, area_gap, area_overlap] = ...
        deformCurve_lap_gap(curves, T_P, feaIds, overlap_thr);
    deformEnergy = 0; % getDeformationEnergy(curves, prev);
    iter = iter + 1;    
    if (area_overlap > overlap_thr)
        curves = preCurves;
        disp('Cause overlaps.');
        break;
    end
    if area_gap == prev_gap
        curves = preCurves;
        disp('Gaps was not changed.');
        break;
    end
    % visualize
    tran_curve = transform_curves(curves, T_Q, scale);
    show_curves(curves, tran_curve, 1);
    total_iter = total_iter + 1;
    iterFile = strcat(folder, 'Iter_', num2str(total_iter), resolution, '_gap.png');
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

disp('===Gap Stage 2: Deforming each cuve in high resolution===');
iter = 0;
max_iter = 15;
resolution = '_local';
while area_gap > gap_thr && iter < max_iter
    [curves, area_gap, area_overlap] = eliminate_gaps(curves, T_P, overlap_thr);
    deformEnergy = 0; % getDeformationEnergy(curves, prev);
    iter = iter + 1;    
    if (area_overlap > overlap_thr)
        disp('Cause overlaps.');
        break;
    end
    if (area_gap == prev_gap)% || prev_gap - area_gap < 1e-4)
        disp('Gaps was not changed.');
        break;
    end
    % visualize
    tran_curve = transform_curves(curves, T_Q, scale);
    show_curves(curves, tran_curve, 1);
    total_iter = total_iter + 1;
    iterFile = strcat(folder, 'Iter_', num2str(total_iter), resolution, '_gap.png');
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