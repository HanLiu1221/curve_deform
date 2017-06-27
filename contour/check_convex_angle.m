%-----copyRight(c) Shuhua Li<sue142857@gmail.com> 01.12.2017-----%

% v1, v2 and v3 are arranged in counter-clockwise direction
function [convex,angle] = check_convex_angle(v1,v2,v3)

convex = 0; 
e1 = v3 - v2; % edge<v2,v3>
e2 = v1- v2;  % edge<v2,v1>

% compute the angle in counter clockwise direction from V1 to V2
% the range of angle is from -180 to 180
angle = atan2d(e1(1)*e2(2)-e1(2)*e2(1),e1(1)*e2(1)+e1(2)*e2(2));
if angle > 0 && angle < 180
    convex = 1;
end