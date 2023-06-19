%% Set up the problem parameters

% Defining functions for BVP
f = @(x,y,y_prime) -exp(-2*y);
dfy = @(x,y,y_prime) 2*exp(-2*y);
dfy_prime = @(x,y,y_prime) 0;


x_start = 1;
x_end = 2;
alpha = 0;
beta = log(2);

N = 9;
TOL = 10^(-4);
MaxIter = 10;
[w,x] = nonlinfindif_easy(x_start,x_end, alpha, beta, f, dfy, dfy_prime, N, MaxIter, TOL);

% Actual solution
y = log(x);



%% Function implementing the nonlinear finite difference algorithm with Newton's method (easier)
function [w,x] = nonlinfindif_easy(x_start,x_end, alpha, beta, f, dfy, dfy_prime, N, MaxIter, TOL)
% This returns the values w_1, ..., w_N approximating the solution y to the
% BVP y' = f(x,y,y') on [a,b] with y(a) = alpha, y(b) = beta.
% The values of w_0 = alpha and w_N+1 = beta are not returned.

    h=(x_end - x_start)/(N+1);
    x = linspace(x_start + h, x_end - h, N); 
    w = zeros(N,1);
    b = zeros(N,1);
    J = zeros(N,N);
    
    tick = 1;
    
    
    for i = 1:N
        w(i) = alpha + i * (beta-alpha)*h/(x_start-x_end);
    end
    
    
    
    while tick <= MaxIter
        
        % Set the first rows of J and row of b
        J(1,2) = -1 + (h/2) * dfy_prime(x(1), w(1), (w(2) - alpha)/(2*h));
        J(1,1) = 2 + h^2 * dfy(x(1), w(1), (w(2) - alpha)/(2*h));
        b(1) = -(2*w(1) - w(2) - alpha + h^2 * f(x(1), w(1), (w(2) - alpha)/(2*h)));
        
        % Set the last rows of J and b
        J(N,N) = 2 + h^2 * dfy(x(N), w(N), (beta - w(N-1))/(2*h));
        J(N,N-1) = -1 - (h/2) * dfy_prime(x(N), w(N), (beta - w(N-1))/(2*h));
        b(N) = -1*(-w(N-1) + 2* w(N) - beta + h^2 * f(x(N), w(N), (beta - w(N-1))/(2*h)));
        
        % Set the remaining entries of J and b
        for i = 2:(N-1)
            J(i,i-1) = -1 - (h/2) * dfy_prime(x(i), w(i), (w(i+1) - w(i-1))/(2*h));
            J(i,i) = 2 + h^2 * dfy(x(i), w(i), (w(i+1) - w(i-1))/(2*h));
            J(i,i+1) = -1 + (h/2) * dfy_prime(x(i), w(i), (w(i+1) - w(i-1))/(2*h));
            
            b(i) = -1*(-w(i-1) + 2* w(i) - w(i+1) + h^2 * f(x(i), w(i), (w(i+1) - w(i-1))/(2*h)));
        end
        
        % Solve the linear system Jv = b and set w = w + v
        v = linsolve(J,b);
        w = w + v;
        
        if norm(v) <= TOL
            break
        end
        
        tick = tick+1;
        disp(tick);
        
    end

        

end