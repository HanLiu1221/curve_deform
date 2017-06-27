%-----copyRight(c) Shuhua Li<sue142857@gmail.com> 01.12.2017-----%

% v1, v2 and v3 are arranged in counter-clockwise direction
function [angle,angle_pi] = compute_angle(v1,v2,v3)
 
e1 = v3 - v2; % edge<v2,v3>
e2 = v1- v2;  % edge<v2,v1>

% compute the angle in counter clockwise direction from V1 to V2
% the range of angle is from -180 to 180
angle = atan2d(e1(1)*e2(2)-e1(2)*e2(1),e1(1)*e2(1)+e1(2)*e2(2));
angle_pi = pi*angle/180;