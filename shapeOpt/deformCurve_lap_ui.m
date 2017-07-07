function deformed = deformCurve_lap_ui()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundary deformation for a given shape, represented as a contour curve
% Input: shape P, with Trunk T_P; Q is the RIOT of P, with Trunk T_Q
% Output: 
% Han Liu, 06/20/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User Interface
%  1. click to set static anchors > 2
%  2. click to set handles [1, 2]
%  3. click the new position of the handle

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

%show_curves(P, tran_P, Q, tran_Q);
%show_curves(P, tran_P);

%% 1. handle pieces that go beyond the trunk
curves = remove_outside_regions(P, T_P);
Rcurves = transform_curves(P, T_Q, scale);
show_curves(curves, Rcurves, 1);
hold on;
%% user loop
userFinish = false;
while ~userFinish
    % start with anchor points
    set(gcf, 'CurrentCharacter', 's'); 
    prevkey = 's';
    userRun = false;
    func = 1;
    nMaxStatic = 6;
    % user set static anchors
    static_ids = zeros(1, nMaxStatic);
    handle_anchor = 0;
    newPos = zeros(1, 2);
    stId = 1;
    % deform one curve at a time
    curveId = 0;
    handlePnt = zeros(1, 2);
    % nHandle = 1 as default
    while ~userRun
        % accept user settings
        key = get(gcf, 'CurrentCharacter');
        switch key
            case 's'
                func = 1;
            case 'h'
                func = 2;
            case 'e'
                func = 3;
            case 'r'
                userRun = true;
            case 'q'
                userFinish = true;
            otherwise
                func = 1;
        end
        if prevkey ~= key
            prevkey = key;
            continue;
        end
        [x, y] = ginput(1);
        mousePos = [x, y];
        if func == 1 || func == 2
            [cid, pid] = findPointIds(curves, mousePos); 
            if cid == 0 || pid == 0
                continue;
            end
            if curveId == 0
                curveId = cid;
            end
            % edit a curve at a time
            if curveId ~= cid
                continue;
            end
        end
        subplot(1,2,1);
        if func == 1 && stId <= nMaxStatic
            if ~isempty(find(static_ids == pid))
                continue;
            end
            % static
            static_ids(stId) = pid;
            % if the user continues to set static points,
            % the previous points will be rewritten
            stId = stId + 1;
            plot(curves{cid}(pid, 1), curves{cid}(pid, 2), 'k*');
        elseif func == 2
            % handle
            handlePnt = [x, y];
            handle_anchor = pid;
            plot(curves{cid}(pid, 1), curves{cid}(pid, 2), 'r*');
        elseif func == 3
            newPos = mousePos;
            plot(x, y, 'ko');
        end
    end
    % run laplacian editing
    static_anchor = static_ids(1, 1 : stId - 1);
    static_anchor = sort(static_anchor);
    deformed = lap2D(curves{curveId}, static_anchor, handle_anchor, newPos - handlePnt);
    curves{curveId} = deformed;
    Rcurves = transform_curves(curves, T_Q, scale);
    show_curves(curves, Rcurves, 1);
end

end


function [cid, pid] = findPointIds(curves, pnt) 
% find the nearest point on a curve
cid = 0;
pid = 0;
mind = 1000;
for i = 1 : length(curves)
    for j = 1 : size(curves{i}, 1)
        d = norm(curves{i}(j, :) - pnt);
        if d < mind
            mind = d;
            cid = i;
            pid = j;
        end
    end
end
end