
%% 0. add references in the subdirectories
clc;
close all;
clear;
addpath(genpath('./'));

%% 1. load curves
load('./input/curves4deform.mat');
P = curves4deform(90).curves1; % P
Q = curves4deform(90).curves2; % Q

%% 2. visulize curves
T_P = get_candidate_trunk(P); % T_P
T_Q = get_candidate_trunk(Q); % T_Q

% get the transform
scale = 0;
tran_P = transform_curves(P, T_Q, scale);
tran_Q = transform_curves(Q, T_P, scale);

show_curves(P, tran_P, Q, tran_Q);
show_curves(P, tran_P);

%% 3. detect feature points and insert sample points as control points 
curves = P; % deform P firstly
CT = T_P;
target_CT = T_Q;

n = length(curves);
d_sample = 0.05;
controlPs = cell(n,1);
for i = 1:n
    [controlPs{i},idx_CPs] = extract_control_points(curves{i},d_sample);   
end

oriControlPs = controlPs;
tran_oriControlPs = transform_curves(controlPs, target_CT, scale);
show_curves(oriControlPs, tran_oriControlPs);

%% 4. remove regions outside the trunk
controlPs = remove_outside_regions(controlPs,CT);

tran_controlPs = transform_curves(controlPs,target_CT,scale);
show_curves(controlPs, tran_controlPs);

%% 5. remove overlaps among pieces
global d_step iter;
d_step = 0.05;
iter = 1;
[~,~,area_overlap] = compute_gap_overlap_area(controlPs,CT);

while area_overlap > 1e-4
    newControlPs = cell(n,1);
    for i = 1:n
        newControlPs{i} = generate_new_control_points(controlPs{i},d_step);
    end  
    controlPs = reduce_overlaps(controlPs,CT,target_CT,area_overlap,newControlPs,oriControlPs); 
    [~,~,area_overlap] = compute_gap_overlap_area(controlPs,CT);
    iter = iter + 1;
end

%% 6. to do: reduce gaps among pieces

