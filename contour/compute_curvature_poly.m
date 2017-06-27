%-----copyRight(c) Shuhua Li<sue142857@gmail.com> 03.14.2017-----%
%function: compute angles at the vertices of the polygon
% -- input: 
%          -poly: arranged in counter-clockwise direction;
%           poly(1,:) == poly(end,:)
function [angles,types] = compute_curvature_poly(poly)

if ~ismember(poly(1,:),poly(end,:),'rows')
    poly(end+1,:) = poly(1,:);
end
n = size(poly,1)-1; % poly(1,:) == poly(end,:)
angles = zeros(n,1);
types = zeros(n,1);
for i = 1:n
    v2 = poly(i,:);
    if i == 1
        v1 = poly(end-1,:);
    else
        v1 = poly(i-1,:);
    end
    v3 = poly(i+1,:);
   [convex,angle] = check_convex_angle(v1,v2,v3);
   if ~convex
       angles(i) = 360 + angle;
       types(i) = 0; % concave
   else
       angles(i) = angle;
       types(i) = 1;
   end
end