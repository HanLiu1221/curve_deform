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
%load('./input/curves4deform_noGap.mat');
for id = 1 : 6
    %% 1. get the original RIOT
    dataFolder = strcat('./input/riot_', num2str(id), '/');
    if ~exist(dataFolder)
        continue;
    end
    load(strcat(dataFolder, 'curves4deform_texture.mat'));
    P = curves4deform_texture.curves1;
    Q = curves4deform_texture.curves2;
    
%     P = curves4deform(id).curves1; % P
%     Q = curves4deform(id).curves2; % Q

    T_P = get_candidate_trunk(P); % T_P
    T_Q = get_candidate_trunk(Q); % T_Q

    scale = 0;
    tran_P = transform_curves(P, T_Q, scale);
    tran_Q = transform_curves(Q, T_P, scale);

    show_curves(P, tran_P, 0, Q, tran_Q);
    
    %% 2. file paths
    folder = strcat('shapeOpt/output/texture/example_', num2str(id), '/');
    
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
   
    %% 4. texture 
    p_path = strcat(dataFolder,  'p.tif');
    q_path = strcat(dataFolder,  'q.tif');
    p_image = imread(p_path);
    q_image = imread(q_path);
    S_p = [curves4deform_texture.scale1_x, curves4deform_texture.scale1_y];
    T_p = [curves4deform_texture.trans1_x, curves4deform_texture.trans1_y];
    S_q = [curves4deform_texture.scale2_x, curves4deform_texture.scale2_y];
    T_q = [curves4deform_texture.trans2_x, curves4deform_texture.trans2_y];
    %[Sp, Tp] = estimateTrans(pnts, [size(q_image, 1), size(q_image, 2)] * 1.0);
    tex_file_p = strcat(folder, 'ori_tex_p_', num2str(id), '.png');
    tex_file_q = strcat(folder, 'ori_tex_q_', num2str(id), '.png');
    mesh_file_p = strcat(folder, 'ori_mesh_p_', num2str(id), '.png');
    mesh_file_q = strcat(folder, 'ori_mesh_q_', num2str(id), '.png');
    [FP, VP, TVP, pnts_P] = computeTexture(tran_P, S_p, T_p, mesh_file_p);
    VF_P = findVertexRing(VP, FP);
    drawTexture(p_image, FP, VP, TVP, tex_file_p);
    [FQ, VQ, TVQ, pnts_Q] = computeTexture(tran_Q, S_q, T_q, mesh_file_q);
    VF_Q = findVertexRing(VQ, FQ);
    drawTexture(q_image, FQ, VQ, TVQ, tex_file_q);
    %set(h, 'CData', p_image, 'FaceColor','texturemap');
    
    %% 5. boundary curve optimization
    [P_o, P_g] = optimizeBoundaryCurve(P, T_P, T_Q, pfolder, 0);
    [Q_o, Q_g]  = optimizeBoundaryCurve(Q, T_Q, T_P, qfolder, 1);
    overlap_file = strcat(folder, 'riot_', num2str(id), '_no_overlap', '.png');
    drawARiot(P_o, Q_o, T_P, T_Q, overlap_file);
    gap_file = strcat(folder, 'riot_', num2str(id), '_min_gap', '.png');
    drawARiot(P_g, Q_g, T_P, T_Q, gap_file);
    
    % after overlap
    tran_P_o = transform_curves(P_o, T_Q, scale);
    tran_Q_o = transform_curves(Q_o, T_P, scale);
    show_curves(P_o, tran_P_o, 0, Q_o, tran_Q_o);
    tex_file_p_o = strcat(folder, 'overlap_tex_p_', num2str(id), '.png');
    tex_file_q_o = strcat(folder, 'overlap_tex_q_', num2str(id), '.png');
    mesh_file_p_o = strcat(folder, 'overlap_mesh_p_', num2str(id), '.png');
    mesh_file_q_o = strcat(folder, 'overlap_mesh_q_', num2str(id), '.png');
    [FPO, VPO, TVPO] = computeTexture(tran_P_o, S_p, T_p, mesh_file_p_o);
    drawTexture(p_image, FPO, VPO, TVPO, tex_file_p_o);
    [FQO, VQO, TVQO] = computeTexture(tran_Q_o, S_q, T_q, mesh_file_q_o);
    drawTexture(q_image, FQO, VQO, TVQO, tex_file_q_o);
    
    % after overlap
    tran_P_g = transform_curves(P_g, T_Q, scale);
    tran_Q_g = transform_curves(Q_g, T_P, scale);
    show_curves(P_g, tran_P_g, 0, Q_g, tran_Q_g);
    tex_file_p_g = strcat(folder, 'gap_tex_p_', num2str(id), '.png');
    tex_file_q_g = strcat(folder, 'gap_tex_q_', num2str(id), '.png');
    mesh_file_p_g = strcat(folder, 'gap_mesh_p_', num2str(id), '.png');
    mesh_file_q_g = strcat(folder, 'gap_mesh_q_', num2str(id), '.png');
    [FPG, VPG, TVPG] = computeTexture(tran_P_o, S_p, T_p, mesh_file_p_g);
    drawTexture(p_image, FPG, VPG, TVPG, tex_file_p_g);
    [FQG, VQG, TVQG] = computeTexture(tran_Q_o, S_q, T_q, mesh_file_q_g);
    drawTexture(q_image, FQG, VQG, TVQG, tex_file_q_g);
    
end
end

function VF = findVertexRing(V, F)
n = length(V);
VF = zeros(n, 6);
vids = ones(1, n);
for i = 1 : length(F)
    v1 = F(i, 1);
    v2 = F(i, 2);
    v3 = F(i, 3);
    % v1
    VF(v1, vids(v1)) = v2;
    vids(v1) = vids(v1) + 1;
    VF(v1, vids(v1)) = v3;
    vids(v1) = vids(v1) + 1;
    % v2
    VF(v2, vids(v2)) = v1;
    vids(v2) = vids(v2) + 1;
    VF(v2, vids(v2)) = v3;
    vids(v2) = vids(v2) + 1;
    % v3
    VF(v3, vids(v3)) = v1;
    vids(v3) = vids(v3) + 1;
    VF(v3, vids(v3)) = v2;
    vids(v3) = vids(v3) + 1;
end
end

function VM = findApproximateMap(tran_P, V)
% extract the boundary p 
npnt = 0;
    for i = 1 : length(tran_P)
        npnt = npnt + length(tran_P{i});
    end
    pnts = zeros(npnt, 2);
end

function [FF, VV, TV, pnts] = computeTexture(tran_P, S, T, mesh_file_p)
    %% 4. Texture test
    npnt = 0;
    for i = 1 : length(tran_P)
        npnt = npnt + length(tran_P{i});
    end
    pnts = zeros(npnt, 2);
    bbox = zeros(2, 2);
    idx = 1;
    for i = 1 : length(tran_P)
        for j = 1 : length(tran_P{i})
            pnts(idx, :) = tran_P{i}(j, :);
            bbox(1, 1) = min(bbox(1, 1), pnts(idx, 1));
            bbox(1, 2) = min(bbox(1, 2), pnts(idx, 2));
            bbox(2, 1) = max(bbox(2, 1), pnts(idx, 1));
            bbox(2, 2) = max(bbox(2, 2), pnts(idx, 2));
            idx = idx + 1;
        end
    end
%     idx = 1;
%     for i = 1 : length(tran_P)
%         pnts(idx : idx + length(tran_P{i}) - 1, :) = tran_P{i};
%         idx = idx + length(tran_P{i});
%     end

    % 4.1 triangulation
    % VV: n * 2 vertices
    % FF: nf * 3 faces
    figure;
    [VV, FF, h] = distmesh2d(@dpoly, @huniform, 0.1, bbox, pnts, pnts);
    axis equal;
    saveas(gcf, mesh_file_p);
%     tri = delaunayTriangulation(pc);
%     figure;
%     triplot(tri, pc(:, 1), pc(:, 2));

    TV = zeros(size(VV)); %texture coordinate
%     maxp = max([pnts(:, 1), pnts(:,2)]) * 1.0;
%     minp = min([pnts(:, 1), pnts(:,2)]) * 1.0;
%     imgSize = [size(q_image, 1), size(q_image, 2)] * 1.0;
    for i = 1 : length(TV)
        pos = VV(i, :);
        TV(i, :) = mapTextureCoord(pos, S, T);
        %TV(i, :) = [TV(i, 1), imgSize(2) - TV(i, 2)];
        TV(i, :) = [TV(i, 2), TV(i, 1)];
        %TV(i, :) = TV(i, :) ./ imgSize;
    end
end

function drawTexture(p_image, FF, VV, TV, tex_file_p)
%     figure;
%     imshow(p_image); hold on;
%     plot(TV(:, 1), TV(:, 2), 'b.');
%     axis equal;
    figure;
    nv = length(VV);
    V = [VV zeros(nv,1)];
    patcht(FF, V, FF, TV, p_image);
    axis equal; %axis off;
    saveas(gcf, tex_file_p, 'png');
end

function tp = mapTextureCoord(p, S, T)
tp = p.*S + T;
end

function [S, T] = estimateTrans(pnts, imSize)
maxp = max([pnts(:, 1), pnts(:,2)]) * 1.0;
minp = min([pnts(:, 1), pnts(:,2)]) * 1.0;
S = imSize / (maxp - minp);
T = imSize / 2 - (maxp + minp) / 2;
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