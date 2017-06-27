function labyrinth_test()

clc;
clear all;
close all;
T_P = [0, 0; 1, 0; 1, 1; 0, 1];
T_Q = [0, 0; 0, 1; 1, 1; 1, 0];
curve = cell(1, 4);

r = 0.5;
curve{1} = sampleCircle([0.5, 0], r, -pi / 2);
curve{2} = sampleCircle([1, 0.5], r, pi);
curve{3} = sampleCircle([0.5, 1], r, pi / 2);
curve{4} = sampleCircle([0, 0.5], r, 0);

figure;
for i = 1:4
    plot(curve{i}(:, 1), curve{i}(:, 2));
    axis equal;
    hold on;
end

[~,~,area_overlap] = compute_gap_overlap_area(curve, T_P);
iter = 1;
thr = 1e-4;
scale = 0;
while area_overlap > thr
    curve = deformCurve(curve, T_P);
    [~,~,area_overlap] = compute_gap_overlap_area(curve, T_P);
    iter = iter + 1;
    
    % visualize
    tran_curve = transform_curves(curve, T_Q, scale);
    show_curves(curve, tran_curve);
end

end

function points = sampleCircle(origin, r, s)

n = 20;
step =  pi / (n - 1);
points = zeros(n, 2);
for i = 0:n-1
    a = s + step * i;
    points(i+1, :) =  origin + [r * sin(a), r * cos(a)];
end

end