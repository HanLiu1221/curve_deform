%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundary deformation for a given shape, represented as a contour curve
% Input: shape P, with Trunk T_P; Q is the RIOT of P, with Trunk T_Q
% Output: 
% Han Liu, 06/20/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function shape_optimization(P, T_P, Q, T_Q) 

%% 1. sample points on the boundary curve

%% 2. handle pieces that go beyond the trunk


%% 3. handle overlaps between any two pieces


% 3.1 detect overlaps

% 3.2 iteratively handle overlap one by one
% as default, each curve is a closed region by connecting two end points
% - the two end points serve as static anchors
% - the furtheres point on the curve in the overlapping region serve as
%   the handle anchor

end