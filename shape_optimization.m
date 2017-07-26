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
%     P = sampleControPoints(P);
%     Q = sampleControPoints(Q);
    show_curves(P, tran_P, 0, Q, tran_Q);
    sample_points_file = strcat(folder, 'riot_', num2str(id), '_sp.png');
    saveas(gcf, sample_points_file);    
    
    %% 4. texture 
    p_path = strcat(dataFolder,  'p.tif');
    q_path = strcat(dataFolder,  'q.tif');
    p_image = imread(p_path);
    q_image = imread(q_path);
    S_p = [curves4deform_texture.scale2_x, curves4deform_texture.scale2_y];
    T_p = [curves4deform_texture.trans2_x, curves4deform_texture.trans2_y];
    S_q = [curves4deform_texture.scale1_x, curves4deform_texture.scale1_y];
    T_q = [curves4deform_texture.trans1_x, curves4deform_texture.trans1_y];
    %[Sp, Tp] = estimateTrans(pnts, [size(q_image, 1), size(q_image, 2)] * 1.0);
    tex_file_p = strcat(folder, 'ori_tex_p_', num2str(id), '.png');
    tex_file_q = strcat(folder, 'ori_tex_q_', num2str(id), '.png');
    mesh_file_p = strcat(folder, 'ori_mesh_p_', num2str(id), '.png');
    mesh_file_q = strcat(folder, 'ori_mesh_q_', num2str(id), '.png');
    
%     %% original points & texture
%     dC2 = curves4deform_texture.deformed_C2;%(length(curves4deform_texture.deformed_C2) : -1 : 1, :);
%     dC1 = curves4deform_texture.deformed_C1;%(length(curves4deform_texture.deformed_C1) : -1 : 1, :);
%     [dC2, id2] = reOrderCurve(dC2, tran_P{1}(1, :));
%     [dC1, id1] = reOrderCurve(dC1, tran_Q{1}(1, :));
%     C2 = curves4deform_texture.C2;%(length(curves4deform_texture.C2) : -1 : 1, :);
%     C1 = curves4deform_texture.C1;%(length(curves4deform_texture.C1) : -1 : 1, :);
%     oP = getCurves(C2, P);
%     oQ = getCurves(C1, Q);
%     
%     
%     [FP, VP, TVP, fixIds_P, handleIds_P] = ...
%         computeTexture(oP, S_p, T_p, mesh_file_p);
%     VF_P = findVertexRing(VP, FP);
%     drawTexture(p_image, FP, VP, TVP, tex_file_p);
%     [FQ, VQ, TVQ, fixIds_Q, handleIds_Q] = ...
%         computeTexture(oQ, S_q, T_q, mesh_file_q);
%     VF_Q = findVertexRing(VQ, FQ);
%     drawTexture(q_image, FQ, VQ, TVQ, tex_file_q);
%     
%     tex_file_p_ori = strcat(folder, 'original_tex_p_', num2str(id), '.png');
%     tex_file_q_ori = strcat(folder, 'original_tex_q_', num2str(id), '.png');
%     
% 
%     Off_P = computeOffset(VP, handleIds_P, tran_P);
%     VPD = lap2D_Tri(VP, VF_P, fixIds_P, handleIds_P, Off_P);
%     simpplot(VPD, FP);
%     drawTexture(p_image, FP, VPD, TVP, tex_file_p_ori);
%     Off_Q = computeOffset(VQ, handleIds_Q,tran_Q);
%     VQD = lap2D_Tri(VQ, VF_Q, fixIds_Q, handleIds_Q, Off_Q);
%     drawTexture(q_image, FQ, VQD, TVQ, tex_file_q_ori);
    
    [FP, VP, TVP, fixIds_P, handleIds_P] = ...
        computeTexture(tran_P, S_p, T_p, mesh_file_p);
    VF_P = findVertexRing(VP, FP);
    [FQ, VQ, TVQ, fixIds_Q, handleIds_Q] = ...
        computeTexture(tran_Q, S_q, T_q, mesh_file_q);
    VF_Q = findVertexRing(VQ, FQ);
    
    figure;
    subplot(1,2,1);
    drawTexture(q_image, FQ, VQ, TVQ, tex_file_q);
    subplot(1,2,2);
    drawTexture(p_image, FP, VP, TVP, tex_file_p);   
    tex_file = strcat(folder, 'original_texture_', num2str(id));
    saveas(gcf, tex_file, 'png');
    %set(h, 'CData', p_image, 'FaceColor','texturemap');
    
    mesh_file = strcat(folder, 'original_mesh_', num2str(id));
    figure;
    subplot(1, 2, 1);
    simpplot(VQ, FQ);
    subplot(1, 2, 2);
    simpplot(VP, FP);
    saveas(gcf, mesh_file, 'png');
    
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
    
%     Off_P_O = computeOffset(VP, handleIds_P, tran_P_o);
%     VPO = lap2D_Tri(VP, VF_P, fixIds_P, handleIds_P, Off_P_O);
%     Off_Q_O = computeOffset(VQ, handleIds_Q, tran_Q_o);
%     VQO = lap2D_Tri(VQ, VF_Q, fixIds_Q, handleIds_Q, Off_Q_O);

    %%test

    VPO_v = updatePoints(VP, handleIds_P, tran_P_o);
    VPO = lap2D_Tri(VPO_v, VF_P, handleIds_P, [], []);
    VQO_v = updatePoints(VQ, handleIds_Q, tran_Q_o);
    VQO = lap2D_Tri(VQO_v, VF_Q, handleIds_Q, [], []);
    
    VPO = updatePoints(VPO, handleIds_P, tran_P_o);
    VQO = updatePoints(VQO, handleIds_Q, tran_Q_o);
    figure;
    subplot(1,2,1);
    drawTexture(q_image, FQ, VQO, TVQ, tex_file_q_o);
    subplot(1,2,2);
    drawTexture(p_image, FP, VPO, TVP, tex_file_p_o);
    
    tex_file_o = strcat(folder, 'overlap_texture_', num2str(id));
    saveas(gcf, tex_file_o, 'png');
    
    mesh_file_g = strcat(folder, 'overlap_mesh_', num2str(id));
     figure;
     subplot(1, 2, 1);
     simpplot(VQO, FQ);
     subplot(1, 2, 2);
     simpplot(VPO, FP);
     saveas(gcf, mesh_file_g, 'png');

    % after overlap
    tran_P_g = transform_curves(P_g, T_Q, scale);
    tran_Q_g = transform_curves(Q_g, T_P, scale);
    show_curves(P_g, tran_P_g, 0, Q_g, tran_Q_g);
    tex_file_p_g = strcat(folder, 'gap_tex_p_', num2str(id), '.png');
    tex_file_q_g = strcat(folder, 'gap_tex_q_', num2str(id), '.png');    

%     Off_P_G = computeOffset(VP, handleIds_P, tran_P_g);
%     VPG = lap2D_Tri(VP, VF_P, fixIds_P, handleIds_P, Off_P_G);
%     Off_Q_G = computeOffset(VQ, handleIds_Q, tran_Q_g);
%     VQG = lap2D_Tri(VQ, VF_Q, fixIds_Q, handleIds_Q, Off_Q_G);
    VPG_v = updatePoints(VP, handleIds_P, tran_P_g);
    VPG = lap2D_Tri(VPG_v, VF_P, handleIds_P, [], []);
    VQG_v = updatePoints(VQ, handleIds_Q, tran_Q_g);
    VQG = lap2D_Tri(VQG_v, VF_Q, handleIds_Q, [], []);
    
    VPG = updatePoints(VPG, handleIds_P, tran_P_g);
    VQG = updatePoints(VQG, handleIds_Q, tran_Q_g);
    
    figure;
    subplot(1,2,1);
    drawTexture(q_image, FQ, VQG, TVQ, tex_file_q_g);
    subplot(1,2,2);
    drawTexture(p_image, FP, VPG, TVP, tex_file_p_g);
    
    tex_file_g = strcat(folder, 'gap_texture_', num2str(id));
    saveas(gcf, tex_file_g, 'png');
    
    mesh_file_g = strcat(folder, 'gap_mesh_', num2str(id));
    figure;
    subplot(1, 2, 1);
    simpplot(VQG, FQ);
    subplot(1, 2, 2);
    simpplot(VPG, FP);
    saveas(gcf, mesh_file_g, 'png');
    
end
end

function [rP, rid] = reOrderCurve(P, start)
idx = find(P(:,1) == start(1));
idy = find(P(:,2) == start(2));
id = intersect(idx, idy);
if isempty(id)
    if ~isempty(idx)
        id = idx;
    else
        id = idy;
    end
end
rid = id(1)
rP = [P(rid : length(P), :); P(1 : rid - 1, :)];
end

function curves = getCurves(P, C)
curves = cell(1, length(C));
id = 1;
figure;
for i = 1 : length(curves)
    curves{i} = zeros(size(C{i}));
    for j = 1 : length(C{i}) - 1
        curves{i}(j, :) = P(id, :);
        id = id + 1;
    end
    if i == length(curves)
        curves{i}(length(C{i}), :) = P(1, :);
    else
        curves{i}(length(C{i}), :) = P(id + 1, :);
    end
    plot(curves{i}(:, 1), curves{i}(:, 2),'k-'); hold on;
end
end

function VF = findVertexRing(V, F)
n = length(V);
VF = ones(n, 6);
VF = -1 * VF;
vids = ones(1, n);
for i = 1 : length(F)
    v1 = F(i, 1);
    v2 = F(i, 2);
    v3 = F(i, 3);
    % v1
    if isempty(find(VF(v1, :) == v2))
        VF(v1, vids(v1)) = v2;
        vids(v1) = vids(v1) + 1;
    end
    if isempty(find(VF(v1, :) == v3))
        VF(v1, vids(v1)) = v3;
        vids(v1) = vids(v1) + 1;
    end
    % v2
    if isempty(find(VF(v2, :) == v1))
        VF(v2, vids(v2)) = v1;
        vids(v2) = vids(v2) + 1;
    end
    if isempty(find(VF(v2, :) == v3))
        VF(v2, vids(v2)) = v3;
        vids(v2) = vids(v2) + 1;
    end
    % v3
    if isempty(find(VF(v3, :) == v1))
        VF(v3, vids(v3)) = v1;
        vids(v3) = vids(v3) + 1;
    end
    if isempty(find(VF(v3, :) == v2))
        VF(v3, vids(v3)) = v2;
        vids(v3) = vids(v3) + 1;
    end
end
VF(VF == 0) = -1;
end

function Off = computeOffset(V, handleId, curves)
Off = zeros(length(handleId), 2);
id = 1;
%figure;
for i = 1 : length(curves)
    for j = 2 : length(curves{i}) - 1
        Off(id, :) = curves{i}(j, :) - V(handleId(id), :);
%         plot(curves{i}(j, 1),curves{i}(j, 2),'b.'); hold on;
%         plot(V(handleId(id), 1),V(handleId(id), 2),'r.'); hold on;
        id = id + 1;
    end
end
end

function V0 = updatePoints(V, handleId, curves)
V0 = V;
id = 1;
%figure;
for i = 1 : length(curves)
    for j = 2 : length(curves{i}) - 1
        V0(handleId(id), :) = curves{i}(j, :);
        id = id + 1;
    end
end
end

function [F, V, TV, fixIds, handleIds] = ...
    computeTexture(tran_P, S, T, mesh_file_p)
    nboundaryPnts = 0;
    nFix = length(tran_P);
    fixIds = zeros(1, nFix);
    for i = 1 : length(tran_P)
        fixIds(i) = nboundaryPnts + 1;
        nboundaryPnts = nboundaryPnts + length(tran_P{i}) - 1;
    end
    pnts = zeros(nboundaryPnts, 2);
    handleIds = zeros(1, nboundaryPnts - length(tran_P));
    bbox = zeros(2, 2);
    idx = 1;
    hid = 1;
    for i = 1 : length(tran_P)
        handleIds(1, hid : hid + length(tran_P{i}) - 3) = ...
            idx + 1 : idx + length(tran_P{i}) - 2;
        hid = hid + length(tran_P{i}) - 2;
        for j = 1 : length(tran_P{i}) - 1
            pnts(idx, :) = tran_P{i}(j, :);
            bbox(1, 1) = min(bbox(1, 1), pnts(idx, 1));
            bbox(1, 2) = min(bbox(1, 2), pnts(idx, 2));
            bbox(2, 1) = max(bbox(2, 1), pnts(idx, 1));
            bbox(2, 2) = max(bbox(2, 2), pnts(idx, 2));
            idx = idx + 1;
        end
    end
    

    %% triangulation
    % VV: n * 2 vertices
    % FF: nf * 3 faces
    figure;
    % IDX: re-ordered points index
    [V, F, I1, I2] = distmesh2d(@dpoly, @huniform, 0.1, bbox, pnts, pnts);
    axis equal; axis on;
    saveas(gcf, mesh_file_p);

    TV = zeros(size(V)); %texture coordinate
    for i = 1 : length(TV)
        pos = V(i, :);
        TV(i, :) = mapTextureCoord(pos, S, T);
        TV(i, :) = [TV(i, 2), TV(i, 1)];
    end
    
    center = (bbox(1, :) + bbox(2, :)) / 2;
    fid = 1;
    fIds = zeros(1, length(V) - length(handleIds));
    for i = length(handleIds) : 1 : length(V)
        pos = V(i,:);
        if norm(pos - center) < 0.1
            fIds(fid) = i;
            fid = fid + 1;
        end
    end

    for i = 1 : length(fixIds)
        id = fixIds(i);
        fixIds(i) = I2(id); % in VV
    end
    
    % add fix points
    isInsidePoly = feval(@dpoly, V, V(fixIds, :)) < 0.001;
    insideIds = find(isInsidePoly == 1);
    fixIds = unique([fixIds insideIds']);
    
    for i = 1 : 1 : length(handleIds)
        id = handleIds(i);
        handleIds(i) = I2(id);
    end
    fixIds = setdiff(fixIds, handleIds);
%     intersect(fixIds,handleIds)
    figure; 
    plot(V(handleIds,1),V(handleIds,2),'k.'); hold on;
    plot(V(fixIds,1),V(fixIds,2),'r.');
end

function drawTexture(p_image, FF, VV, TV, tex_file_p)
%     figure;
%     imshow(p_image); hold on;
%     plot(TV(:, 1), TV(:, 2), 'b.');
%     axis equal;
    %figure;
    nv = length(VV);
    V = [VV zeros(nv,1)];
    patcht(FF, V, FF, TV, p_image);
    axis equal; axis off;
    %saveas(gcf, tex_file_p, 'png');
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