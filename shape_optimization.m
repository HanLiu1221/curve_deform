%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundary deformation for a given shape, represented as a contour curve
% Input: shape P, with Trunk T_P; Q is the RIOT of P, with Trunk T_Q
% Output: 
% Han Liu, 06/20/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function shape_optimization() 

%% 0.0 add references in the subdirectories
clc;
close all;
clear;
addpath(genpath('./'));

%% 1 load curves
load('./input/curves4deform_noGap.mat');
id = 1;
P = curves4deform(id).curves1; % P
Q = curves4deform(id).curves2; % Q

% 1.1 get trunks
T_P = get_candidate_trunk(P); % T_P
T_Q = get_candidate_trunk(Q); % T_Q

% 1.2 get the transform
scale = 0;
tran_P = transform_curves(P, T_Q, scale);
tran_Q = transform_curves(Q, T_P, scale);

% show_curves(P, tran_P, 0, Q, tran_Q);
% show_curves(P, tran_P, 0);

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

%% 3. iteratively eliminate overlaps
%  3.1 define a local region, using r, to compute the attraction force
[~,~,area_overlap] = compute_gap_overlap_area(controlPs, T_P);
iter = 1;
thr = 1e-4;
curve = controlPs;
for i = 1 : length(curve)
    curve{i} = removeReps(curve{i});
end
% find feature points
%% find feature points
nCurves = length(curve);
featIds_noSplit = cell(1, nCurves);
featfile = 'data/featurePointIds.mat';
feaIds = cell(1, nCurves);
loaded = 0;
if exist(featfile, 'file') == 2
    load(featfile);
    loaded = 1;
    feaIds = featurePointIds;
end
figure;
for i = 1:nCurves
    if loaded == 0
        feaIds{i} = getSplitCurvePointIds(curve{i});
    end
    npnts = length(curve{i});
    featIds_noSplit{i} = [1, floor((1 + npnts) / 2), npnts];
    % draw
    plot(curve{i}(:,1), curve{i}(:,2), 'k-');
    hold on;
    % split points
    splitIds = [feaIds{i}(:, 1)', npnts];
    plot(curve{i}(splitIds,1), curve{i}(splitIds,2), 'r--o', 'LineWidth',2);
    hold on;
    for j = 1:length(splitIds)
        text(curve{i}(splitIds(j), 1), curve{i}(splitIds(j), 2), num2str(splitIds(j)));
        hold on;
    end
end
legend('original polyline', 'simplified');

while area_overlap > thr
    %tic;
    %curve = deformCurve_force(curve, T_P);
    area_overlap
    if area_overlap > 0.005 
        curve = deformCurve_lap(curve, T_P, featIds_noSplit, 1);
    elseif area_overlap > 0.001 
        curve = deformCurve_lap(curve, T_P, feaIds, 0);
    else
        curve = deformCurve_force(curve, T_P);
    end
    %toc
    [~,~,area_overlap] = compute_gap_overlap_area(curve, T_P);
    iter = iter + 1;
    
    % visualize
    tran_curve = transform_curves(curve, T_Q, scale);
    show_curves(curve, tran_curve, 1);
end

% C = controlPs;
% s = 0.99;
% while area_overlap > thr
%     for i = 1:n
%         newC = C{i};
%         newC = scale_poly(newC, s);
%         C{i}(2:end-1,:) = newC(2:end-1,:);
%         [~,~,area_overlap] = compute_gap_overlap_area(C, T_P);
%     end
%     % visualize
%     tC = transform_curves(C, T_Q, scale);
%     show_curves(C, tC);
% end

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