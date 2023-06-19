%% Runge-Kutta Method for Systems of Differential Equations

%% Input information
a = 0;          % left endpoint
b = 1;          % right endpoint
m = 2;          % number of equations
h = 0.2;        % stepsize
N = (b-a)/h;    % number of subintervals
alpha1 = 1;     % initial conditions
alpha2 = 1;

f1 = @(t,u1,u2) 3*u1 + 2*u2 - (2*t^2 + 1)*exp(2*t);
f2 = @(t,u1,u2) 4*u1 + u2 + (t^2+2*t-4)*exp(2*t);

% exact solutions
u1 = @(t) 1/3*exp(5*t)-1/3*exp(-t)+exp(2*t);
u2 = @(t) 1/3*exp(5*t)+2/3*exp(-t)+t^2*exp(2*t);

%% Do the method

t = a;

w1 = alpha1;
w2 = alpha2;

% output starting information
fprintf('t \t\t\t w1 \t\t u1 \t\t w2 \t\t u2\n')                    % header
fprintf('%f \t %f \t %f \t %f \t %f \n',t,w1,u1(t),w2,u2(t))        % initial information

for i=1:N

    k(1,1) = h * f1(t, w1, w2); % w1, w2 are in place of u1, u2
    k(1,2) = h * f2(t, w1, w2);

    k(2,1) = h * f1(t + h/2, w1 + k(1,1)/2, w2 + k(1,2)/2);
    k(2,2) = h * f2(t + h/2, w1 + k(1,1)/2, w2 + k(1,2)/2);

    k(3,1) = h * f1(t + h/2, w1 + k(2,1)/2, w2 + k(2,2)/2);
    k(3,2) = h * f2(t + h/2, w1 + k(2,1)/2, w2 + k(2,2)/2);

    k(4,1) = h * f1(t + h, w1 + k(3,1), w2 + k(3,2));
    k(4,2) = h * f2(t + h, w1 + k(3,1), w2 + k(3,2));

    w1 = w1 + (k(1,1) + 2*k(2,1) + 2*k(3,1) + k(4,1))/6;
    w2 = w2 + (k(1,2) + 2*k(2,2) + 2*k(3,2) + k(4,2))/6;

    t = a + i*h;

    fprintf('%f \t %f \t %f \t %f \t %f \n',t,w1,u1(t),w2,u2(t))
end





