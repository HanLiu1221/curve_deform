%% --copyRight @ Shuhua Li <sue142857@gmail.com> 06/15/2017--%%
% - function: rotate/translate/scale(option) curves based on
% the target trunk, which provides end point correspondeces
% - input:
%         -curves: a set of curves
%         -CT: target trunk
%         -scale: option for scale 
% - output: 
%          -new_curves: curves after rigid transformations
%          -Con: contour closed by new curves

function new_curves = transform_curves(curves,CT,scale)
debug = 0;
n = length(curves);
new_curves = cell(n,1);
for i = 1:n
    p1 = curves{i};
    if i == n
        p2 = [CT(n,:);CT(1,:)];
    else
        p2 = CT(i:i+1,:);
    end
    
    if scale
        % scale
        s = norm(p2(1,:)-p2(2,:),2)/norm(p1(1,:)-p1(end,:),2);
        p1_scale = s*p1;
    else
        p1_scale = p1;
    end
    % translate
    num = size(p1,1);
    t = p2(1,:) - p1_scale(1,:);
    p1_trasl = p1_scale + repmat(t,num,1);
    % rotate
    v1 = p2(2,:);v2 = p2(1,:); v3 = p1_trasl(end,:);
    [~,rotatA] = compute_angle(v1,v2,v3);
    if rotatA < 0
        rotatA = 2*pi + rotatA;
    end
    p1_rotat = rotate_point(p1_trasl(1,:),p1_trasl,rotatA);
  
    new_curves{i} = p1_rotat;

    if debug
        figure; ax = cell(2,1);
        ax{1} = subplot(1,2,1);
        for j = 1:n
            plot(curves{j}(:,1),curves{j}(:,2),'b-','LineWidth',2);hold on;
            e = [curves{j}(1,:);curves{j}(end,:)];
            plot(e(:,1),e(:,2),'k.-','LineWidth',2,'MarkerSize',20);hold on;
        end
        axis equal;
        
        ax{2} = subplot(1,2,2);
        for j = 1:n
            if j < n
                plot(CT(j:j+1,1),CT(j:j+1,2),'k.-','LineWidth',2,'MarkerSize',20);hold on;
            else
                plot([CT(n,1) CT(1,1)],[CT(n,2) CT(1,2)],'k.-','LineWidth',2,'MarkerSize',20);hold on;
            end
            if ~isempty(new_curves{j})
                plot(new_curves{j}(:,1),new_curves{j}(:,2),'r-','LineWidth',2);hold on;
            end
        end
        axis equal;
        linkaxes([ax{1},ax{2}],'xy'); 
        hold off;
    end
    
end


