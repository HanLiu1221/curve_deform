%-----CopyRight(c) Shuhua Li<sue142857@gmail.com> 04.18.2017-----%
%--function: compute the area of regions outside the conjugate polygon, the
%area of overlaps and gaps inside the conjugate polygon.

function [area_out,area_gap,area_overlap] = compute_gap_overlap_area(curves,CT)

debug = 0;
scale = 100; % scale the area to 10000*area
n = length(curves);
%% 1. convert to clockwise direction and scale the area 
xMin = Inf; yMin = Inf;
xMax = -Inf; yMax = -Inf;

[x,y] = poly2cw(CT(:,1),CT(:,2)); 
x = scale*[x;x(1)]; y = scale*[y;y(1)];
CT = [x,y];
if min(x) < xMin
    xMin = min(x);
end
if min(y) < yMin
    yMin = min(y);
end
if max(x) > xMax
    xMax = max(x);
end
if max(y) > yMax
    yMax = max(y);
end

for i = 1:n
    [x,y] = poly2cw(curves{i}(:,1),curves{i}(:,2)); 
    x = scale*[x;x(1)]; y = scale*[y;y(1)];
    curves{i} = [x,y];
    if min(x) < xMin
        xMin = min(x);
    end
    if min(y) < yMin
        yMin = min(y);
    end
    if max(x) > xMax
        xMax = max(x);
    end
    if max(y) > yMax
        yMax = max(y);
    end
end

%% 2. compute area of the union of inside regions
xVec = xMin:xMax;
yVec = yMin:yMax;

BW_union = mpoly2mask(curves, xVec, yVec); % Union all objects
BW_CT = mpoly2mask({CT}, xVec, yVec); % Union all objects

BW_out = BW_union & (~BW_CT); 
area_out = bwarea(BW_out)/scale^2;

BW_gap = BW_CT & (~BW_union); 
area_gap = bwarea(BW_gap)/scale^2;

BW_inside_union = BW_union & BW_CT; 
area_inside_union = bwarea(BW_inside_union)/scale^2;

area_inPieces = zeros(n,1);
for i = 1:n
    BW_pi = mpoly2mask(curves(i), xVec, yVec); 
    BW_inside_pi = BW_pi & BW_CT; 
    area_inPieces(i) = bwarea(BW_inside_pi)/scale^2;
end
area_overlap = sum(area_inPieces) - area_inside_union; 
% debug
if debug
    bg = {'w','k'};
    colors = distinguishable_colors(n,bg);
    figure;
    subplot(2,2,1);
    plot([CT(:,1);CT(1,1)],[CT(:,2);CT(1,2)],'k.-','MarkerSize',20,'LineWidth',2);hold on;
    for i = 1:n
        plot(curves{i}(:,1),curves{i}(:,2),'--','color',colors(i,:),'LineWidth',2);hold on;
        patch([curves{i}(:,1);curves{i}(1,1)],[curves{i}(:,2);curves{i}(1,2)],colors(i,:),'FaceAlpha',0.3); hold on;
    end
    axis equal; axis off; title('outside-overlaps-gaps');
    
    subplot(2,2,2); imagesc(xVec,yVec,BW_union); 
    set(gca,'Ydir','Normal'); axis image; title('union');
    subplot(2,2,3); imagesc(xVec,yVec,BW_out); 
    set(gca,'Ydir','Normal');axis image; title('outside');
    subplot(2,2,4); imagesc(xVec,yVec,BW_gap); 
    set(gca,'Ydir','Normal');axis image; title('gap');
end
