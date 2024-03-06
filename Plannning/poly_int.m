function[X] =  poly_int(x,t0,tf,var)
% Definite Integration of polynomials
% -----------------------------------
% Input:
%   x [1x(n+1)]- coefficients of the polynomials of order 'n' OR
%   x [1x1] - degree of the polynomial
%   t0 [1x1] - Initial time of Intergration
%   tf [1x1] - Final time
%   var [string] - conditions
%
% Output:
% X [1x1] - Final value of the integration. 
% ------------------------------------

    if strcmp(var,'degree')

        T0 = t0*ones(1,x+1);
        Tf = tf*ones(1,x+1);
        N = [1:x+1];
        X = ((Tf.^N)-(T0.^N)).*(N.^(-1));

    elseif strcmp(var,'coeff')

        T0 = t0*ones(1,size(x,2));
        Tf = tf*ones(1,size(x,2));
        N = [1:size(x,2)];
        X = ((Tf.^N)-(T0.^N)).*(N.^(-1)).*x;

    end

end