for i = 1:length(t)
    temp_point = [sum(x(1:n1).*((t(i)*N).^powers)),...
                    sum(x((1*n1+1):2*n1).*((t(i)*N).^powers))...
                    sum(x((2*n1+1):3*n1).*((t(i)*N).^powers))];
    points = [points;temp_point];
end


 % velocity, acceleration and jerk/jolt constraints
        for p = 1%:(r-1)
            if j== 1
                dAeq(kk,n1*(k-1)+1:n1*k) = [zeros(1,p), (ti(j)*poly_diff(N,p)).^(powers(1:end-p))];
                k = k+1;
                kk = kk+1;
            elseif j ==m 
                dAeq(kk,n1*(k-2)+1:n1*k) = [zeros(1,p), (ti(j)*poly_diff(N,p)).^(powers(1:end-p)),...
                    zeros(1,p), -(ti(j)*poly_diff(N,p)).^(powers(1:end-p))];
                kk = kk+1;
                dAeq(kk,n1*(k-1)+1:n1*k) = [zeros(1,p), (ti(j)*poly_diff(N,p)).^(powers(1:end-p))];
                k = k+1;
                kk = kk+1;
            else 
                dAeq(kk,n1*(k-2)+1:n1*k) = [zeros(1,p), (ti(j)*poly_diff(N,p)).^(powers(1:end-p)),...
                    zeros(1,p), -(ti(j)*poly_diff(N,p)).^(powers(1:end-p))];
                k = k+1;
                kk = kk+1;
            end
        end

