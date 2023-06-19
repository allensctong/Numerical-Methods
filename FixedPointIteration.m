%% Fixed-Point Iteration

%% Input Information

G = @(x) [(cos(x(2) * x(3)) + 0.5) / 3;
          1 / 25 * sqrt(x(1)^2 + 0.3125) - 0.03;
          -1 / 20 * exp(-x(1) * x(2)) - (10 * pi - 3) / 60];

p0 = [1;
      1;
      1];

tol = 1e-5;                 % tolerance, 1e-4 = 10^{-4}

max_iter = 30;              % max number of iterations

%% Fixed-Point Iteration

i = 1;                      % iteration count

fprintf('i\tx1\t\t\tx2\t\t\tx3\t\tinf error\n');          % for display
fprintf('%d\t%f\t%f\t%f\n',0,p0(1),p0(2),p0(3));

while( i <= max_iter)

    % get p_i 
    p = G(p0);              % p_i 
    
    % calculate infty norm
    inf_error = max(abs(p - p0));

    % display information
    fprintf('%d\t%f\t%f\t%f\t%f\n',i,p(1),p(2),p(3),inf_error);    % displays iteration i, p_i
    
    % check stopping condition 
    if(inf_error < tol)
        break;
    end
    
    % increase iteration count
    i = i + 1;
    
    % prepare for next iteration
    p0 = p;
    
end
