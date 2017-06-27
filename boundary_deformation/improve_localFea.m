%% --copyRight @ Shuhua Li <sue142857@gmail.com> 05/23/2017--%%

function newC = improve_localFea(C,oriC,p_fix,CT,varagin)

global d_step;
debug = 0;
if nargin > 4
    controlPs  = varagin{1};
    idx_C = varagin{2};
end
%% improve other point based on two measures
fea = compute_localFea(oriC);
cur_fea = compute_localFea(C);
d_fea = compute_localFea_dist(fea,cur_fea);
% show_localFea_dist(C,d_fea);
n = size(C,1);
while 1
    newvs = generate_new_control_points(C,d_step);
    m = size(newvs,2);
    loc = [];
    for i = 2:n-1
         if i == p_fix
             continue;
         end
         for j = 1:m
             % check 1
             check = check_inside_poly(newvs{i,j},CT);
             if ~check
                 continue;
             end
             curC = C;
             curC(i,:) = newvs{i,j};
             % check 2: check whether the point is moved inside other pieces
             if nargin > 4
                 for kk = 1:length(controlPs)
                    if idx_C == kk
                        continue;
                    end
                    Ckk = controlPs{kk}; 
                    check1 = check_inside_poly(curC(2:end-1,:),Ckk);
                    check2 = check_inside_poly(Ckk(2:end-1,:),curC);
                    if check1 || check2
                        break;
                    end
                end
                if check1 || check2
                    continue;
                end
             end
             
             cur_fea_ = compute_localFea(curC);
             d_fea_ = compute_localFea_dist(fea,cur_fea_);
             check1 = max(d_fea_) < max(d_fea);
             check2 = (max(d_fea_) == max(d_fea)) && (sum(d_fea_) < sum(d_fea));
             if check1 || check2
                 d_fea = d_fea_;
                 loc = [i,j];
             end
         end
    end
    if isempty(loc)
        break;
    end
    if debug
        curC = C;
        curC(loc(1),:) = newvs{loc(1),loc(2)};
%         show_localFea_dist(curC,d_fea);
    end
    C(loc(1),:) = newvs{loc(1),loc(2)};
end
newC = C;