%% --copyRight @ Shuhua Li <sue142857@gmail.com> 05/10/2017--%%
% -function: for move each control point along m directions and obtain m new
% control points

function newControlPs = generate_new_control_points(controlPs,d_step)

n = size(controlPs,1);
m = 8; % the number of directions
newControlPs = cell(n,2*m);
%% 1. compute m directions
% normal direction
Lines = zeros(n,2);
Lines(:,1) = (1:n)';
Lines(:,2) = [2:n 1]';
N=LineNormals2D(controlPs,Lines);
% rotate normal direction to obtain m-1 directions
Ns = cell(n,m);
Ns(:,1) = num2cell(N,2);
for i = 1:n
    for j = 2:m
        theta = (j-1)*2*pi/m;
        Ns{i,j} = rotate_point([0 0],N(i,:),theta);
    end
end

% k=LineCurvature2D(controlPs,Lines);
% figure,  hold on;
% for i = 1:m
%     Ni = cell2mat(Ns(:,i));
%     plot([controlPs(:,1) controlPs(:,1)+k.*Ni(:,1)/100]',[controlPs(:,2) controlPs(:,2)+k.*Ni(:,2)/100]','g');
% end
% plot([controlPs(Lines(:,1),1) controlPs(Lines(:,2),1)]',[controlPs(Lines(:,1),2) controlPs(Lines(:,2),2)]','b');
% plot(controlPs(:,1),controlPs(:,2),'r.');
% axis equal;

%% 2. move control points along inward/outward directions to obtain new control points
for i = 2:n-1 % the first and last control points are fixed
    for j = 1:m
        newControlPs{i,j} = controlPs(i,:) + d_step*Ns{i,j};
    end
    for j = 1:m
        newControlPs{i,m+j} = controlPs(i,:) - d_step*Ns{i,j};
    end
end
% % the first and last control points are fixed
% newControlPs{1,1} = controlPs(1,:);
% newControlPs{n,1} = controlPs(n,:);