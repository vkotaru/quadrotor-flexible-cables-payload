% Function to return the nth order Legendre polynomial
% Inputs:
%   n - Order
%   time_int [OPTIONAL] - these polynomials will be orthogonal over the
%   specified interval
%   val [OPTIONAL] - evaluate at this value
function P = legendre_poly(n,time_int,val)
    syms x real ;
    if(n==0)
        % diff(..., 0) gives wrong value.
        if(nargin == 3)
            P = 1 ;
        else
            P = sym(1) ;
        end
        return ;
    else
        P = 1/(2^n * factorial(n)) * diff( (x^2-1)^n, n) ;
    end
    if(nargin > 1) % Time shift specified
        t0 = time_int(1) ; t1 = time_int(2) ;
        P = sqrt(2/(t1-t0)) * subs(P,x,(x-t0)/(t1-t0)+(x-t1)/(t1-t0)) ;
    end
    if(nargin == 3)
        P = subs(P,x,val) ;
    end
end