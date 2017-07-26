%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adapted from: Matlab script for 2D Laplacian Editing
% 06/20/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function defCurv = lap2D_Tri(curve, VF, static_anchors, handle_anchors, offset)

n=length(curve);

% the Laplacian matrix (uniform weighting)
L = spdiags(ones(n,1),0,n,n) - spdiags(ones(n,1),1,n,n);
L = L+L';
L(1,n)= -1;
L(n,1) = -1;
L = L./2;


delta = L*curve;


% we want to construct the matrix of the system for v-primes
L_prime = [   L     zeros(n)   % the x-part
	       zeros(n)    L    ]; % the y-part
	      
for i=1:n
  % the neighbors of i are i-1 and i+1
  %ring = [(1 + mod(i-2,n)), i, (1 + mod(i,n))]; % i-1, i, i+1
  j = 1;
  while j <= length(VF(i, :)) && VF(i, j) ~= -1
      j = j + 1;
  end
  ring = VF(i, 1 : j - 1);
      
  V = curve(ring,:)';
  V = [V
       ones(1,length(ring))];


  % the coeff matrix for the system that solves for T
  %      s  a  t1
  % T = -a  s  t2
  %      0  0  1
  
  C = zeros(length(ring) * 2,4);
  % ... Fill C in
  for r=1:length(ring)
    C(r,:) =                [V(1,r)       V(2,r)  V(3,r)      0  ];
    C(length(ring)+r,:) =   [V(2,r)  (-1)*V(1,r)       0  V(3,r) ];

  end;
    
  Cinv = pinv(C);
  s =  Cinv(1,:);
  a =  Cinv(2,:);

  delta_i = delta(i,:)';
  delta_ix = delta_i(1);
  delta_iy = delta_i(2);
  
  % T*delta gives us an array of coefficients
  Tdelta = [delta_ix*s      + delta_iy*a 
	        delta_ix*(-1)*a + delta_iy*s];
	    
  
  % updating the weights in Lx_prime, Ly_prime, Lw_prime
  L_prime(i,[ring (ring + n)]) = L_prime(i,[ring (ring + n)]) +...
                                              (-1)*Tdelta(1,:);
  L_prime(i+n,[ring (ring + n)]) = L_prime(i+n,[ring (ring + n)]) +...
                                                (-1)*Tdelta(2,:);
end;


% additional constraints - we need at least 3

% weight for the constraints
w=1;

% building the least-squares system matrix
A_prime = L_prime;
rhs = zeros(2*n,1);
anch_pos = [];

anchors = [static_anchors handle_anchors];
for j=1:length(anchors)
  A_prime = [A_prime
	     w*((1:(2*n))==anchors(j))
	     w*((1:(2*n))==(anchors(j)+n))];
  rhs = [rhs
	 w*curve(anchors(j),1)
	 w*curve(anchors(j),2)];
  
  anch_pos = [anch_pos
	      curve(anchors(j),1:2)];
end;


% displaying the curve: static anchors are in black, the handle to be moved
% in red.
% figure;
% plot(curve(:,1),curve(:,2),'b-',...
%      curve(static_anchors,1),curve(static_anchors,2),'*k',...
%      curve(handle_anchors,1),curve(handle_anchors,2),'*r');
% axis equal;

% moving the handle
if ~isempty(handle_anchors)
handles = curve(handle_anchors, :);
new_handles = handles + offset;
%[x_input y_input] = ginput(1);
x_input = new_handles(:, 1); 
y_input = new_handles(:, 2);
lenr = length(rhs);
start = lenr - length(x_input) * 2 + 1;
rhs( start: 2 : lenr - 1) = x_input;
rhs( start + 1: 2 : lenr) = y_input;
anch_pos(length(anch_pos) - length(x_input) + 1 : length(anch_pos),:) = [x_input y_input];
end
% moving second handle vertex if exists
% if length(handle_anchors) > 1
%     [x_input y_input] = new_handles(2, :);
%     rhs((lenr-3):(lenr-2)) = [x_input y_input]';
%     anch_pos(length(anch_pos)-1,:) = [x_input y_input];
% end


% solving for v-primes
curve_col = A_prime\rhs;
defCurv = [curve_col(1:n) curve_col((n+1):(2*n))];

% figure;
% plot(curve(:,1),curve(:,2),'g-',...
%      defCurv(:,1),defCurv(:,2),'-b',...
%      defCurv(anchors,1),defCurv(anchors,2),'*k',...
%      anch_pos(:,1),anch_pos(:,2), 'mo',...
%      defCurv(anchors(length(anchors)),1),defCurv(anchors(length(anchors)),2),'*r');
% axis equal;
% title('Editing result');

end 



