function [x, fail] = R_lp_MP(anc, rangems, p, maxiter, gamma)
%Paper: A Message Passing Based Iterative Algorithm for Robust TOA Positioning in Impulsive Noise

%-Inputs
%anc - matrix including sensor positions 
%rangems - measured distance vector
%p - model parameter
%maxiter - max iteration number
%gamma - a small tolerance

%-Outputs
%x - location estimate
%fail - true if fail to converge, otherwise false

L = length(anc);
%x = -10 + 20*rand(2,1);
x = zeros(2,1);
w = zeros(L,1);
epsilon = zeros(2,L);
zeta = zeros(L,1);

fail = false;

%k: counter for iterations
k = 0;

while 1
    
    x_old = x;
    
    %update w, epsilon, and zeta
    for i = 1:L
        w(i) = abs((rangems(i)-norm(x_old - anc(:,i)))^2+0.000001)^((p-2)/2);
        epsilon(:,i) = (x_old - anc(:,i))/norm(x_old - anc(:,i));
        zeta(i) = norm(x_old - anc(:,i)) - epsilon(:,i)'*x_old;
    end
    
    %update x
    sum1 = zeros(2,2);
    sum2 = zeros(2,1);
    for i = 1:L
        sum1 = sum1 + (2*w(i))*(epsilon(:,i)*(epsilon(:,i)'));
        sum2 = sum2 + epsilon(:,i)*(2*w(i))*(rangems(i) - zeta(i));
    end
    
    if ((isnan(sum(sum(sum1)))) || (isnan(sum(sum(sum2)))))
        fail = true;
        break
    end
    
    x = (pinv(sum1))*sum2;
    
    if ((norm(x - x_old)/min(norm(x),norm(x_old))) < gamma) || ((k+1) == maxiter)
        if ((k+1) == maxiter)
            fprintf('have reached the maximum iteration number\n')
            fail = true;
        end
        break
    end
    
    k = k + 1;
    
end

fprintf('It takes %d iterations to meet a termination condition\n', k)

end

