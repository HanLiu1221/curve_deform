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

%% load curves
load('./input/curves4deform_noGap.mat');
for id = 1 : length(curves4deform)
    P = curves4deform(id).curves1; % P
    Q = curves4deform(id).curves2; % Q

    % 1. get trunks
    T_P = get_candidate_trunk(P); % T_P
    T_Q = get_candidate_trunk(Q); % T_Q

    % 2. get the transform
    scale = 0;
    tran_P = transform_curves(P, T_Q, scale);
    tran_Q = transform_curves(Q, T_P, scale);

    show_curves(P, tran_P, 0, Q, tran_Q);
    %show_curves(P, tran_P, 0);
    
    folder = strcat('shapeOpt/output/example_', num2str(id), '/');
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
    
    % 3. boundary curve optimization
    optimizeBoundaryCurve(P, T_P, T_Q, scale, pfolder, 0);
    optimizeBoundaryCurve(Q, T_Q, T_P, scale, qfolder, 1);
end
end