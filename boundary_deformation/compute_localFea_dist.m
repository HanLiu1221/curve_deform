

function d_fea = compute_localFea_dist(fea,cur_fea)

% weight = [1 1 0.5];
d = abs(fea - cur_fea);

% idx1 = (fea(:,3)>0 & fea(:,3) <180) & (~(cur_fea(:,3) < 180));
% idx2 = (fea(:,3)>180) & (~(cur_fea(:,3) > 180));
% d(idx1,3) = Inf;
% d(idx2,3) = Inf;
d(:,3) = 0.8*d(:,3)/360;

d_fea = sum(d,2);
