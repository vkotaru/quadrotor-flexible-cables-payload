function[dx] =  poly_diff(x,r)
% function to differentiate polynomials interms of coefficients

if r >= size(x,2)
    dx = 0;
else
    dx = x;
    for i = 1:r
        n = length(dx);
        dx = dx(2:n).*[1:(n-1)];
    end
end

end