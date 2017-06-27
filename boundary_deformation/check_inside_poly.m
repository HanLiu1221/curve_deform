%% --copyRight @ Shuhua Li <sue142857@gmail.com> 05/27/2017--%%

function check = check_inside_poly(v,CT)

xv = [CT(:,1);CT(1,1)];
yv = [CT(:,2);CT(1,2)];
% check whether the point is moved outside the polygon
[in,on] = inpolygon(v(:,1),v(:,2),xv,yv);
check = max(in-on); % not include on the polygon

