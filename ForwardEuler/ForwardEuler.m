%% Forward Euler

%% Inputs

a = 0;          % left endpoint
b = 1;          % right endpoint
h = 0.1;        % stepsize
N = (b-a)/h;    % the number of steps
alpha = 0;      % initial y value

f = @(t,y) t*exp(3*t) - 2*y;        % as in dy/dt = f(t,y);

%% Forward Euler

t = zeros(1,N+1);       % stores all the t values
w = zeros(1,N+1);       % stores all the approximation values

t(1) = a;
w(1) = alpha;

for i=1:N
    w(i+1) = w(i) + h*f(t(i),w(i));
    t(i+1) = a + i*h;
end

%% Plot the approximation

figure()
plot(t,w,'*-')
hold on;            % so we can plot multiple functions on the same graph

%% Plot the exact solution

y = @(t) (1/5)*t*exp(3*t) - (1/25)*exp(3*t) + (1/25)*exp(-2*t);

num_plot = 100;     % need a lot of plotting points to get a smooth graph!

t_plot = zeros(1,num_plot+1);
y_plot = zeros(1,num_plot+1);
h_plot = (b-a)/num_plot;

for i=1:num_plot+1
    t_plot(i) = a + (i-1) * h_plot;
    y_plot(i) = y(t_plot(i));
end

plot(t_plot,y_plot)
title("Forward Euler Method for solving y' = te^{3t} - 2y, 0 \leq t \leq 1")
legend("Approximation","Exact Solution")
