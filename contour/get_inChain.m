%-----copyRight(c) Shuhua Li<sue142857@gmail.com> 04.05.2017-----%
% function: rotate/translate/scale(option) each piece in outChain2 based on
% two end point correspondeces
 
function inChain1 = get_inChain(vs_T1,outChain2,scale)
debug = 0;
n = size(vs_T1,1);
inChain1 = cell(n,1);
for i = 1:n
    p1 = outChain2{i};
    if i == n
        p2 = [vs_T1(n,:);vs_T1(1,:)];
    else
        p2 = vs_T1(i:i+1,:);
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

    inChain1{i} = p1_rotat;

    if debug
        figure; ax = cell(2,1);
        ax{1} = subplot(1,2,1);
        for j = 1:n
            plot(outChain2{j}(:,1),outChain2{j}(:,2),'b-','LineWidth',2);hold on;
            e = [outChain2{j}(1,:);outChain2{j}(end,:)];
            plot(e(:,1),e(:,2),'k.-','LineWidth',2,'MarkerSize',20);hold on;
        end
        axis equal;
        
        ax{2} = subplot(1,2,2);
        for j = 1:n
            if j < n
                plot(vs_T1(j:j+1,1),vs_T1(j:j+1,2),'k.-','LineWidth',2,'MarkerSize',20);hold on;
            else
                plot([vs_T1(n,1) vs_T1(1,1)],[vs_T1(n,2) vs_T1(1,2)],'k.-','LineWidth',2,'MarkerSize',20);hold on;
            end
            if ~isempty(inChain1{j})
                plot(inChain1{j}(:,1),inChain1{j}(:,2),'r-','LineWidth',2);hold on;
            end
        end
        axis equal;
        linkaxes([ax{1},ax{2}],'xy'); 
        hold off;
    end
end


