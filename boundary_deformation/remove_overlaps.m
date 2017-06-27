%% --copyRight @ Shuhua Li <sue142857@gmail.com> 05/23/2017--%%


function controlPs = remove_overlaps(controlPs)

debug = 0;
n = length(controlPs);
scale = 0.99;
for i = 1:n
    for j = 1:n
        if i == j
            continue;
        end
        Ci = controlPs{i}; newCi = Ci;
        Cj = controlPs{j}; newCj = Cj;
        while 1
            check1 = check_inside_poly(Ci(2:end-1,:),Cj);
            check2 = check_inside_poly(Cj(2:end-1,:),Ci);
            if ~check1 && ~check2
                break;
            end
            newCi = scale_poly(newCi,scale);
            newCj = scale_poly(newCj,scale);
            Ci(2:end-1,:) = newCi(2:end-1,:);
            Cj(2:end-1,:) = newCj(2:end-1,:);
            if debug
                cur_controlPs = controlPs;
                cur_controlPs{i} = Ci;
                cur_controlPs{j} = Cj;
                show_curves(cur_controlPs);
            end
        end
        controlPs{i} = Ci;
        controlPs{j} = Cj;
        if debug
            show_curves(controlPs);
        end
    end
end