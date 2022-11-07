%% Higher Order Taylor Method

%% Inputs

a = 1;          % left endpoint
b = 2;          % right endpoint
h = 0.05;       % stepsize
N = (b-a)/h;    % the number of steps
alpha = -1;     % initial y value

f = @(t,y) 1/t^2 - y/t - y^2;        % as in dy/dt = f(t,y);
df = @(t) -2/t^3;
d2f = @(t) 6/t^4;
d3f = @(t) -24/t^5;

%% Order 2

t = zeros(1,N+1);       % stores all the t values
w = zeros(1,N+1);       % stores all the approximation values for order 2

t(1) = a;               % initial value
w(1) = alpha;

%% Approximations (w)
for i=1:N
    w(i+1) = w(i) + h*f(t(i),w(i)) + (h^2/2)*df(t(i)); % using current value to approximate next value
    t(i+1) = a + i*h; % increment t to get to the next value
end

%% Order 4

u = zeros(1,N+1);       % stores all the approximation values for order 4

u(1) = alpha;

%% Approximations (u)
for i=1:N
    u(i+1) = u(i) + h*f(t(i),u(i)) + (h^2/2)*df(t(i)) + (h^3/factorial(3))*d2f(t(i)) + (h^4/factorial(4))*d3f(t(i));
    t(i+1) = a + i*h;
end

%% Plot the approximation

figure()
plot(t,w,'*-')
hold on;            % so we can plot multiple things on the same graph
plot(t,u,'o-')

%% Plot the exact solution

y = @(t) -1/t;

num_plot = 100;     % need a lot of plotting points to get a smooth graph!

t_plot = zeros(1,num_plot+1);
y_plot = zeros(1,num_plot+1);
h_plot = (b-a)/num_plot;

for i=1:num_plot+1
    t_plot(i) = a + (i-1) * h_plot;
    y_plot(i) = y(t_plot(i));
end

plot(t_plot,y_plot,'-')
title("Order 2 & Order 4 Methods for solving y' = 1/t^2 - y/t - y^2, 1 \leq t \leq 2")
legend("Order 2","Order 4","Exact Solution")

%% Compute the actual errors, error bound, and print information

error = zeros(1,N+1);
error_order_4 = zeros(1,N+1);
%error_bound = zeros(1,N+1);
%M = 1;
%L = 1;
fprintf('i\tt_i\t\tw_i\t\tu_i\t\ty(t_i)\t\t|y(t_i) - w_i)|\t|y(t_i) - u(i)|\n')

for i=1:N+1
    error(i) = abs( y(t(i)) - w(i) );                 % | y(t_i) - w_i |
    error_order_4(i) = abs( y(t(i)) - u(i) );                 % | y(t_i) - u_i |
    %error_bound(i) = (h*M/2*L) * (exp(L*(t(i)-a)) - 1 );
    fprintf('%d\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\n',i,t(i),w(i),u(i),y(t(i)),error(i),error_order_4(i))
end

%% Plot the error

figure()
semilogy(t,error,'*-')
hold on
semilogy(t,error_order_4,'o-')
title("Error using Euler Method to solve y' = 1/t^2 - y/t - y^2, 1 \leq t \leq 2")
legend("Order 2", "Order 4")















