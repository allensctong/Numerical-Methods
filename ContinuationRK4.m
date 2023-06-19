%% Set up the problem parameters
f1 = @(x) x(1)^2 - x(2)^2 + 2 * x(2);
f2 = @(x) 2 * x(1) + x(2)^2 - 6;

F = @(x) [f1(x); f2(x)];
           

J = @(x) [ 2 * x(1), -2 * x(2) + 2;
           2,        2 * x(2)     ];
       
       
 
N = 1;   

x0 = [0; 0];
x = continuation_rk4(x0, F, J, N);

% x0 = [1; 1];
% x = continuation_rk4(x0, F, J, N)

% x0 = [3; -2];
% x = continuation_rk4(x0, F, J, N)

%% Function implementing the continuation using RK4

function x = continuation_rk4(x0, F, J, N)

    h = 1/N;
    b = - h * F(x0);
    
    fprintf('i\t\tx1\t\tx2\n')
    fprintf('%d\t%.9f\t%.9f\t\n', 0, x0(1), x0(2))
    for i = 1:N
       A = J(x0);
       k1 = linsolve(A,b);

       A = J(x0+1/2*k1);
       k2 = linsolve(A,b);

       A = J(x0+1/2*k2);
       k3 = linsolve(A,b);

       A = J(x0+k3);
       k4 = linsolve(A,b);
       
       x0 = x0 + (k1+2*k2+2*k3+k4)/6;

       fprintf('%d\t%.9f\t%.9f\t\n', i,x0(1),x0(2))

       
    end
    x = x0;
end
