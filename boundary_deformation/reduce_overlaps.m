%% --copyRight @ Shuhua Li <sue142857@gmail.com> 05/22/2017--%%
% - function: select the optimal new control point, which reduce overlaps
% and keep local features

function controlPs = reduce_overlaps(controlPs,CT,target_CT,area_overlap1,...
    newControlPs,oriControlPs)

otol = 1e-4;
debug = 1;
%% evaluate each new control point based the overlap area and local feature
n = length(controlPs);
opt_d_fea = Inf; % optimal feature distance 
for i = 1:n % for each curve
    oriCi = oriControlPs{i};
    fea = compute_localFea(oriCi);
    
    m = size(oriCi,1);
    numl = size(newControlPs{i},2);
    for j = 2:m-1   % for each sample point
        for k = 1:numl % for each direction
%             fprintf('%d %d %d;\n',i,j,k);
            cur_controlPs = controlPs;
            v = cur_controlPs{i}(j,:);
            newv = newControlPs{i}{j,k};
            cur_controlPs{i}(j,:) = newv;
            cur_Ci = cur_controlPs{i};
           
            % 1. check whether the point is moved outside the polygon
            check = check_inside_poly(newv,CT);
            if ~check
                continue;
            end
            % 2. check the distance to the trunk
            xv = [CT(:,1);CT(1,1)]; xv = xv';
            yv = [CT(:,2);CT(1,2)]; yv = yv';
            d_min = p_poly_dist(v(1), v(2), xv, yv);
            new_d_min = p_poly_dist(newv(1), newv(2), xv, yv);
            if ~(abs(d_min) < 0.01) && (abs(new_d_min) < 0.01)
                continue;
            end
            % 3. check self-intesection of new curve
            cur_Ci = cur_Ci';
            P = InterX(cur_Ci); 
            if ~isempty(P)
                continue;
            end
            % 4. compute new gap and overlap 
            [~,~,area_overlap2] = compute_gap_overlap_area(cur_controlPs,CT);
            % 5. compute new local feature  
            cur_Ci = cur_Ci';
            cur_fea = compute_localFea(cur_Ci);
            d_fea = compute_localFea_dist(fea,cur_fea);
            d_fea3 = d_fea(j-1:j+1);
            if area_overlap2 > area_overlap1 - otol % make sure the overlaps descrese 
                continue;
            end
            if sum(d_fea3) < opt_d_fea % keep the feature as more as possible
                opt_d_fea = sum(d_fea3);
                opt_d_fea3 = d_fea3;
                loc_opt = [i j k];
                area_overlap = area_overlap2;
            end  
        end
    end
end
if opt_d_fea == Inf
    return;
end
%% display
optP = newControlPs{loc_opt(1)}{loc_opt(2),loc_opt(3)};
global iter;
if debug %&& mod(iter,10)==0
   scale = 0;
   tran_controlPs = transform_curves(controlPs,target_CT,scale);
   new_controlPs = controlPs;
   new_controlPs{loc_opt(1)}(loc_opt(2),:) =  optP;
   tran_new_controlPs = transform_curves(new_controlPs,target_CT,scale);
   
   newSeg = new_controlPs{loc_opt(1)}(loc_opt(2)-1:loc_opt(2)+1,:);
   tran_newSeg = tran_new_controlPs{loc_opt(1)}(loc_opt(2)-1:loc_opt(2)+1,:);
   show_new_curves(controlPs,tran_controlPs,newSeg,tran_newSeg);
   str = ['iter' num2str(iter)];
   saveas(gcf,str,'png');
   close all;
end
controlPs{loc_opt(1)}(loc_opt(2),:) = optP;
%% improve_localFea
newC = improve_localFea(controlPs{loc_opt(1)},oriControlPs{loc_opt(1)},loc_opt(2),CT);
controlPs{loc_opt(1)} = newC;
