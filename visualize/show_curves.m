%% -----copyRight(c) Shuhua Li<sue142857@gmail.com> 06.15.2017-----%
% -function: show correspondences btw inside curves and outside curves
function show_curves(curves,tran_curves,toHold,varargin)

curves1 = curves;
tran_curves1 = tran_curves;
if ~isempty(varargin)
    curves2 = varargin{1};
    tran_curves2 = varargin{2};
end
%% 
n = length(curves1);
bg = {'r','g','b','k'};
colors = distinguishable_colors(n,bg);
figure;%('Visible','off'); 
%%  show curves and tran_curves in the left and right subplots
ax=cell(2,1);
xMin = Inf; xMax = -Inf;
yMin = Inf; yMax = -Inf;
%% the first subplot
ax{1} = subplot(1,2,1); h = zeros(n,1); Leg = cell(n,1);
for i = 1:n
    %% curves1
    x1 = curves1{i}(:,1); y1 = curves1{i}(:,2);
    if min(x1)< xMin
        xMin = min(x1);
    end
    if min(y1)< yMin
        yMin = min(y1);
    end
    if max(x1)> xMax
        xMax = max(x1);
    end
    if max(y1)> yMax
        yMax = max(y1);
    end
    % fill the piece based on the edge color
    patch([x1;x1(1)],[y1;y1(1)],colors(i,:),'FaceAlpha',0.3); hold on;
    %text
    geom = center_gravity(x1,y1);
    xc = geom(1); yc = geom(2); 
    
    Str=['Q' num2str(i)];
    text(xc,yc,Str,...%'BackgroundColor',[1 1 1],...
        'Color','black',...
        'FontSize',10,...
        'VerticalAlignment','middle',...
        'HorizontalAlignment','center');
    % boundary
    plot(x1,y1,'r-','LineWidth',2);hold on;
    
    h(i) = plot([x1(1) x1(end)],[y1(1) y1(end)],'-','color',colors(i,:),'LineWidth',2);hold on;
    plot([x1(1) x1(end)],[y1(1) y1(end)],'k.','MarkerSize',20);hold on;
    %
    Leg{i} = ['Piece ' int2str(i)];
    
    %% tran_curves2
    if ~isempty(varargin)
        x1 = tran_curves2{i}(:,1); y1 = tran_curves2{i}(:,2);
        if min(x1)< xMin
            xMin = min(x1);
        end
        if min(y1)< yMin
            yMin = min(y1);
        end
        if max(x1)> xMax
            xMax = max(x1);
        end
        if max(y1)> yMax
            yMax = max(y1);
        end
        % fill the piece based on the edge color
        patch([x1;x1(1)],[y1;y1(1)],colors(i,:),'FaceAlpha',0.3); hold on; 
        %text
        geom = center_gravity(x1,y1);  
        xc = geom(1); yc = geom(2); 

        Str=['P' num2str(i)];
        text(xc,yc,Str,...%'BackgroundColor',[1 1 1],...
            'Color','black',...
            'FontSize',10,...
            'VerticalAlignment','middle',...
            'HorizontalAlignment','center');
        % boundary
        plot(x1,y1,'b-','LineWidth',2);hold on;
        h(i) = plot([x1(1) x1(end)],[y1(1) y1(end)],'-','color',colors(i,:),'LineWidth',2);hold on;
        plot([x1(1) x1(end)],[y1(1) y1(end)],'k.','MarkerSize',20);hold on;
    end
end
% legend1 = legend(h,Leg);clear Leg; 
% set(legend1,'Location','northeast');
axis equal;axis off;
if toHold == 1
    hold on;
else
    hold off;
end

%% the second subplot
ax{2} = subplot(1,2,2); h = zeros(n,1); Leg = cell(n,1);
for i = 1:n
    %% curves1
    x1 = tran_curves1{i}(:,1); y1 = tran_curves1{i}(:,2);
    if min(x1)< xMin
        xMin = min(x1);
    end
    if min(y1)< yMin
        yMin = min(y1);
    end
    if max(x1)> xMax
        xMax = max(x1);
    end
    if max(y1)> yMax
        yMax = max(y1);
    end
    % fill the piece based on the edge color
    patch([x1;x1(1)],[y1;y1(1)],colors(i,:),'FaceAlpha',0.3); hold on;
    %text
    geom = center_gravity(x1,y1);
    xc = geom(1); yc = geom(2); 
    
    Str=['Q' num2str(i)];
    text(xc,yc,Str,...%'BackgroundColor',[1 1 1],...
        'Color','black',...
        'FontSize',10,...
        'VerticalAlignment','middle',...
        'HorizontalAlignment','center');
    % boundary
    plot(x1,y1,'r-','LineWidth',2);hold on;
    
    h(i) = plot([x1(1) x1(end)],[y1(1) y1(end)],'-','color',colors(i,:),'LineWidth',2);hold on;
    plot([x1(1) x1(end)],[y1(1) y1(end)],'k.','MarkerSize',20);hold on;
    %
    Leg{i} = ['Piece ' int2str(i)];
    
    %% curves2
    if ~isempty(varargin)
        x1 = curves2{i}(:,1); y1 = curves2{i}(:,2);
        if min(x1)< xMin
            xMin = min(x1);
        end
        if min(y1)< yMin
            yMin = min(y1);
        end
        if max(x1)> xMax
            xMax = max(x1);
        end
        if max(y1)> yMax
            yMax = max(y1);
        end
        % fill the piece based on the edge color
        patch([x1;x1(1)],[y1;y1(1)],colors(i,:),'FaceAlpha',0.3); hold on; 
        %text
        geom = center_gravity(x1,y1);  
        xc = geom(1); yc = geom(2); 

        Str=['P' num2str(i)];
        text(xc,yc,Str,...%'BackgroundColor',[1 1 1],...
            'Color','black',...
            'FontSize',10,...
            'VerticalAlignment','middle',...
            'HorizontalAlignment','center');
        % boundary
        plot(x1,y1,'b-','LineWidth',2);hold on;
        h(i) = plot([x1(1) x1(end)],[y1(1) y1(end)],'-','color',colors(i,:),'LineWidth',2);hold on;
        plot([x1(1) x1(end)],[y1(1) y1(end)],'k.','MarkerSize',20);hold on;
    end
end
% legend1 = legend(h,Leg);clear Leg; 
% set(legend1,'Location','northeast');
axis equal;axis off;
if toHold == 1
    hold on;
else
    hold off;
end
%%
linkaxes([ax{1},ax{2}],'xy');
ax{2}.XLim = [xMin,xMax];
ax{2}.YLim = [yMin,yMax];

set(gcf,'color','w');
if toHold == 1
    hold on;
else
    hold off;
end

