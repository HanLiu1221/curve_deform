%-----copyRight(c) Shuhua Li<sue142857@gmail.com> 03.15.2017-----%
%function: counter-clockwise rotation for point sets around the first point
% --input:
%         -theta: unit(pi)
function rotated_p = rotate_point(pc,p,theta)

debug = 0;

R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

n = size(p,1);
pc = pc';
center = repmat(pc,1,n);
p = p';

rotated_p = R*(p-center) + center;
rotated_p = rotated_p';

if debug
    p = p';
    plot(p(:,1), p(:,2), 'k-', rotated_p(:,1), rotated_p(:,2), 'r-', pc(1), pc(2), 'bo');
    axis equal;
end


