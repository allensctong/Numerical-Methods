%% Newton's Method for Systems

%% Input Information

F = @(x) [3 * x(1)^2 - x(2)^2;
          3 * x(1) * x(2)^2 - x(1)^3 - 1];

J = @(x) [6 * x(1),                    -2 * x(2);
          3 * x(2)^2 - 3 * x(1)^2,     6 * x(1) * x(2)];

x = [1;
     1];

tol = 1e-6;                 % tolerance, 1e-6 = 10^{-6}

max_iter = 10;              % max number of iterations

%% Newton's Method for Systems

k = 1;                      % iteration count

fprintf('i\tx1\t\t\tx2\t\t\tinf error\n');          % for display
fprintf('%d\t%f\t%f\t\n',0,x(1),x(2));

while( k <= max_iter)

    % solve the system J(x)y = -F(x)
    y = J(x)\-F(x);

    x = x+y;
    
    % calculate infty norm
    inf_error = max(abs(y));

    % display information
    fprintf('%d\t%f\t%f\t%f\n',k,x(1),x(2),inf_error);    % displays iteration i, p_i
    
    % check stopping condition 
    if(inf_error < tol)
        break;
    end
    
    % increase iteration count
    k = k + 1;
    
end
