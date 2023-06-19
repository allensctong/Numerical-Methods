%% Broyden's Method

%% Input Information

F = @(x) [x(1) * (1-x(1)) + 4 * x(2) - 12;
          (x(1) - 2)^2 + (2 * x(2) - 3)^2 - 25];

J = @(x) [1 - 2 * x(1),         4;
          2 * x(1) - 4,         4 * x(2) - 6];

x = [2.5;
     4];

tol = 1e-6;                 % tolerance, 1e-6 = 10^{-6}

max_iter = 10;              % max number of iterations

%% Broyden's Method

fprintf('i\tx1\t\tx2\t\t\tinf error\n');          % for display
fprintf('%d\t%f\t%f\t\n',0,x(1),x(2));

A0 = J(x);
v = F(x);

A = inv(A0);    % inv in general is really expensive

s = -A*v;
x = x+s;
k = 2;

fprintf('%d\t%f\t%f\t\n',1,x(1),x(2));

while( k <= max_iter)

    w = v;
    v = F(x);
    y = v-w;

    z = -A*y;
    
    p = -transpose(s)*z;

    u = transpose(s)*A;

    A = A +(1/p)*(s+z)*u;

    s = -A*v;

    x = x+s;
    
    % calculate infty norm
    inf_error = max(abs(s));

    % display information
    fprintf('%d\t%f\t%f\t%f\t\n',k,x(1),x(2),inf_error);    % displays iteration i, p_i
    
    % check stopping condition 
    if(inf_error < tol)
        break;
    end
    
    % increase iteration count
    k = k + 1;
    
end
