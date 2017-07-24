%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundary deformation for a given shape, represented as a contour curve
% Input: shape P, with Trunk T_P; Q is the RIOT of P, with Trunk T_Q
% Output: 
% Han Liu, 06/20/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function shape_optimization_no_texture() 

%% 0.0 add references in the subdirectories
clc;
close all;
clear;
addpath(genpath('./'));

%% load curves
%load('./input/curves4deform_noGap.mat');
load('./input/curves4deform_6examples.mat');
for id = 1 : length(curves4deform)
    %% 1. get the original RIOT
    
    P = curves4deform(id).curves1; % P
    Q = curves4deform(id).curves2; % Q

    T_P = get_candidate_trunk(P); % T_P
    T_Q = get_candidate_trunk(Q); % T_Q

    scale = 0;
    tran_P = transform_curves(P, T_Q, scale);
    tran_Q = transform_curves(Q, T_P, scale);

    show_curves(P, tran_P, 0, Q, tran_Q);
    
    %% 2. file paths
    folder = strcat('shapeOpt/output/data_1/example_', num2str(id), '/');
    
    if ~exist(folder)
        mkdir(folder);
    else
        delete(strcat(folder, '*.png'));
    end
    riot_file = strcat(folder, 'riot_', num2str(id), '.png');
    saveas(gcf, riot_file);
    
    pfolder = strcat(folder, 'P/');
    qfolder = strcat(folder, 'Q/');
    if ~exist(pfolder)
        mkdir(pfolder);
    else
        delete(strcat(pfolder, '*.png'));
    end
    if ~exist(qfolder)
        mkdir(qfolder);
    else
        delete(strcat(qfolder, '*.png'));
    end
    
    %% 3. prepare data
    % 3.1. sample points on the boundary curve
    P = sampleControPoints(P);
    Q = sampleControPoints(Q);
    show_curves(P, tran_P, 0, Q, tran_Q);
    sample_points_file = strcat(folder, 'riot_', num2str(id), '_sp.png');
    saveas(gcf, sample_points_file);    
   
    
    %% 5. boundary curve optimization
    [P_o, P_g] = optimizeBoundaryCurve(P, T_P, T_Q, pfolder, 0);
    [Q_o, Q_g]  = optimizeBoundaryCurve(Q, T_Q, T_P, qfolder, 1);
    overlap_file = strcat(folder, 'riot_', num2str(id), '_no_overlap', '.png');
    drawARiot(P_o, Q_o, T_P, T_Q, overlap_file);
    gap_file = strcat(folder, 'riot_', num2str(id), '_min_gap', '.png');
    drawARiot(P_g, Q_g, T_P, T_Q, gap_file);
    
end
end

function drawARiot(P, Q, T_P, T_Q, filename)
figure;
tran_P = transform_curves(P, T_Q, 0);
tran_Q = transform_curves(Q, T_P, 0);
show_curves(P, tran_P, 0, Q, tran_Q);
saveas(gcf, filename);
end

function cP = sampleControPoints(P)
n = length(P);
d_sample = 0.04;
cP = cell(n,1);
for i = 1 : n
    [cP{i}, ~] = extract_control_points(P{i}, d_sample);   
end
end