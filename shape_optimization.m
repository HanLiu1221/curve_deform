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
load('./input/curves4deform.mat');
P = curves4deform(90).curves1; % P
Q = curves4deform(90).curves2; % Q

% 1.1 get trunks
T_P = get_candidate_trunk(P); % T_P
T_Q = get_candidate_trunk(Q); % T_Q

% 1.2 get the transform
scale = 0;
tran_P = transform_curves(P, T_Q, scale);
tran_Q = transform_curves(Q, T_P, scale);

show_curves(P, tran_P, Q, tran_Q);
show_curves(P, tran_P);

%% 1. sample points on the boundary curve
n = length(P);
d_sample = 0.05;
ctrlPnts = cell(n,1);
for i = 1:n
    [ctrlPnts{i},idx_CPs] = extract_control_points(P{i}, d_sample);   
end
% visualize
scale = 0;
oriControlPs = ctrlPnts;
tran_oriControlPs = transform_curves(ctrlPnts, T_Q, scale);
show_curves(oriControlPs, tran_oriControlPs);

%% 2. handle pieces that go beyond the trunk
controlPs = remove_outside_regions(ctrlPnts, T_P);
tran_controlPs = transform_curves(controlPs, T_Q, scale);
show_curves(controlPs, tran_controlPs);

%% 3. iteratively eliminate overlaps
%  3.1 define a local region, using r, to compute the attraction force
[~,~,area_overlap] = compute_gap_overlap_area(controlPs, T_P);
iter = 1;
thr = 1e-4;
curve = controlPs;

while area_overlap > thr
    tic;
    curve = deformCurve(curve, T_P);
    toc
    [~,~,area_overlap] = compute_gap_overlap_area(curve, T_P);
    iter = iter + 1;
    
    % visualize
    tran_curve = transform_curves(curve, T_Q, scale);
    show_curves(curve, tran_curve);
end

end