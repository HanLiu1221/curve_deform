function Cen = center_gravity( x, y ) 
% This function computes the center gravity of a polygon whose coordinates
% are given with x and y

% This function is extracted from Polygeom function.
% Credit goes to 
% H.J. Sommer III - 02.05.14 - tested under MATLAB v5.2
%
% code available at:
%    http://www.me.psu.edu/sommer/me562/polygeom.m
% derivation of equations available at:
%    http://www.me.psu.edu/sommer/me562/polygeom.doc

% number of vertices
[ x ] = shiftdim( x );
[ y ] = shiftdim( y );
[ n ] = size( x,1 );

% temporarily shift data to mean of vertices for improved accuracy
xm = mean(x);
ym = mean(y);
x = x - xm*ones(n,1);
y = y - ym*ones(n,1);

% delta x and delta y
dx = x( [ 2:n 1 ] ) - x;
dy = y( [ 2:n 1 ] ) - y;

% summations for CW boundary integrals
A = sum( y.*dx - x.*dy )/2;
Axc = sum( 6*x.*y.*dx -3*x.*x.*dy +3*y.*dx.*dx +dx.*dx.*dy )/12;
Ayc = sum( 3*y.*y.*dx -6*x.*y.*dy -3*x.*dy.*dy -dx.*dy.*dy )/12;

% check for CCW versus CW boundary
if A < 0,
  A = -A;
  Axc = -Axc;
  Ayc = -Ayc;
end

% centroidal moments
xc = Axc / A;
yc = Ayc / A;

% replace mean of vertices
Cen = [xc + xm yc + ym];