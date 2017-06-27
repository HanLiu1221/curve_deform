%% --copyRight @ Shuhua Li <sue142857@gmail.com> 05/09/2017--%%
%-function: extract control points to deform curve

function [controlPs,idx_CPs] = extract_control_points(curve,d_sample)


%% 1. detect feature points
t_curvature = 10;
n = size(curve,1);
Lines = zeros(n,2);
Lines(:,1) = (1:n)';
Lines(:,2) = [2:n 1]';
k=LineCurvature2D(curve,Lines);
% N=LineNormals2D(curve,Lines);
% figure,  hold on;
% plot([curve(:,1) curve(:,1)+k.*N(:,1)/100]',[curve(:,2) curve(:,2)+k.*N(:,2)/100]','g');
% plot([curve(Lines(:,1),1) curve(Lines(:,2),1)]',[curve(Lines(:,1),2) curve(Lines(:,2),2)]','b');
% plot(curve(:,1),curve(:,2),'r.');
% axis equal;

fea_ps = find(abs(k)>t_curvature);
%% 2. make sure two end points of the curve being included
if isempty(fea_ps) 
    fea_ps = 1;
elseif fea_ps(1) ~=1
    tmp = fea_ps;
    fea_ps(1) = 1;
    fea_ps(2:end+1) = tmp;
end
if fea_ps(end) ~=n
    fea_ps(end+1) = n;
end
%% 3. insert sample points 
s = sample_points_via_midpoint(curve,fea_ps,d_sample);
idx_CPs = s.index;
controlPs = curve(idx_CPs,:);
