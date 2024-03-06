% Function to return the mth derivative of nth order Legendre polynomial
% Inputs:
%   m - derivative order
%   n - Order
%   time_int [OPTIONAL] - these polynomials will be orthogonal over the
%   specified interval
%   val [OPTIONAL] - evaluate at this value
function P = int_legendre_poly(m,n,time_int,val)
    if(nargin < 3)
        time_int = [-1 1] ;
    end
    
    P = legendre_poly(n,time_int) ;
    for j=1:m
        P = int(P) ;
    end
    
    if(nargin == 4)
        syms x real
        P = subs(P, x, val) ;
    end
end