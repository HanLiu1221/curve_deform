%% --copyRight @ Shuhua Li <sue142857@gmail.com> 05/24/2017--%%

function controlPs = reduce_gaps(controlPs,T,area_gap1,newControlPs,oriControlPs)

otol = 1e-4;
debug = 1;
%% evaluate each new control point based the gap area and local feature
n = length(controlPs);
opt_d_fea = Inf;
for i = 1:n
    oriCi = oriControlPs{i};
    fea = compute_localFea(oriCi);
    
    m = size(oriCi,1);
    numl = size(newControlPs{i},2);
    for j = 2:m-1   
        for k = 1:numl
%             fprintf('%d %d %d;\n',i,j,k);
            cur_controlPs = controlPs;
            v = cur_controlPs{i}(j,:);
            newv = newControlPs{i}{j,k};
            cur_controlPs{i}(j,:) = newv;
            cur_Ci = cur_controlPs{i};
            % 0. check whether there are overlaps
            for kk = 1:n
                if i == kk
                    continue;
                end
                Ckk = controlPs{kk}; 
                check1 = check_inside_poly(cur_Ci(2:end-1,:),Ckk);
                check2 = check_inside_poly(Ckk(2:end-1,:),cur_Ci);
                if check1 || check2
                    break;
                end
            end
            if check1 || check2
                continue;
            end
            % 1. check whether the point is moved outside the polygon
            check = check_inside_poly(newv,T);
            if ~check
                continue;
            end
            % 2. check the distance to the trunk
            xv = [T(:,1);T(1,1)]; xv = xv';
            yv = [T(:,2);T(1,2)]; yv = yv';
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
            [~,area_gap2,~] = compute_gap_overlap_area(cur_controlPs,T);
            if area_gap2 > area_gap1 - otol
                continue;
            end
            % 5. compute new local feature
            cur_Ci = cur_Ci';
            cur_fea = compute_localFea(cur_Ci);
            d_fea = compute_localFea_dist(fea,cur_fea);
            d_fea3 = d_fea(j-1:j+1);
            if sum(d_fea3) < opt_d_fea
                opt_d_fea = sum(d_fea3);
                opt_d_fea3 = d_fea3;
                loc_opt = [i j k];
                area_gap = area_gap2;
%                 show_new_curve(controlPs,cur_Ci,loc_opt(1:2),0,area_gap-area_gap1,opt_d_fea3,'');
            end  
        end
    end
end
if opt_d_fea == Inf
    return;
end
%% display
global iter;
optP = newControlPs{loc_opt(1)}{loc_opt(2),loc_opt(3)};
if debug %&& mod(iter,10)==0
   cur_C = controlPs{loc_opt(1)};
   cur_C(loc_opt(2),:) =  optP;
   dirName = ['iter ' num2str(iter)];
   show_new_curve(controlPs,cur_C,loc_opt(1:2),area_gap-area_gap1,0,opt_d_fea3,dirName);
   close all;
end
controlPs{loc_opt(1)}(loc_opt(2),:) = optP;
%% improve_localFea
newC = improve_localFea(controlPs{loc_opt(1)},oriControlPs{loc_opt(1)},loc_opt(2),T,{controlPs,loc_opt(1)});
controlPs{loc_opt(1)} = newC;
