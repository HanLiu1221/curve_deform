%% --copyRight @ Shuhua Li <sue142857@gmail.com> 05/24/2017--%%

function newP = scale_poly(P,scale)

debug = 0;
n = size(P,1);
cen = center_gravity(P(:,1),P(:,2));
newP = scale*(P - repmat(cen,n,1)) + repmat(cen,n,1);

if debug
    figure;
    plot(P(:,1),P(:,2),'b.-'); hold on;
    plot(newP(:,1),newP(:,2),'go-'); hold on;
    plot(cen(1),cen(2),'b.','MarkerSize',10);
    hold off;
    axis equal;
end