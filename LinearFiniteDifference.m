%% Nonlinear Shooting With Newton's Method

%% Input Information
A = 0;          % left endpoint
B = 1;          % right endpoint
alpha = 2;      % boundary condition at left endpoint
beta = 1;       % boundary condition at right endpoint
N = 9;         % number of subintervals

% f = @(x,y,y_prime) y_prime*cos(x)-y*log(y);
% partialf_partialy = @(x,y,y_prime) -log(y)-1;
% partialf_partialy_prime = @(x,y,y_prime) cos(x);

p = @(x) -3;
q = @(x) 2;
r = @(x) 2*x+3;


%% Do the method

% step 1
h = (B-A)/(N+1);
x = A + h;
a = zeros(1,N+1);
b = zeros(1,N+1);
c = zeros(1,N+1);
d = zeros(1,N+1);
l = zeros(1,N+1);
u = zeros(1,N+1);
z = zeros(1,N+1);
a(1) = 2 + h^2 * q(x);
b(1) = -1 + h/2 * p(x);
d(1) = -h^2 * r(x) + (1 + h/2 * p(x)) * alpha;

fprintf('x \t\t\t w \n')

% step 2
for i = 2:N-1
    x = A + i*h;

    a(i) = 2 + h^2*q(x);
    b(i) = -1 + h/2 * p(x);
    c(i) = -1 - h/2 * p(x);
    d(i) = -h^2*r(x);

end

% step 3
x = B-h;
a(N) = 2 + h^2*q(x);
c(N) = -1 - h/2*p(x);
d(N) = -h^2*r(x)+(1-(h/2)*p(x))*beta;

% step 4
l(1) = a(1);
u(1) = b(1)/a(1);
z(1) = d(1)/l(1);

% step 5
for i=2:N-1
    l(i) = a(i) - c(i)*u(i-1);
    u(i) = b(i)/l(i);
    z(i) = (d(i)-c(i)*z(i-1))/l(i);
end

% step 6
l(N) = a(N)-c(N)*u(N-1);
z(N) = (d(N)-c(N)*z(N-1))/l(N);

% step 7
% modified indices
% for solving tridiagonal
w(1) = alpha; % will be overwritten
w(N+1) = beta;
w(N) = z(N);

% step 8
for i = N-1:-1:1 % for reverse indexing (put -1 in between so reverse works properly)
    w(i)=z(i)-u(i)*w(i+1);
end

% step 9
fprintf('%f \t %f \n', A, alpha)
for i = 1:N+1
    x = A+i*h;
    fprintf('%f \t %f \n',x,w(i))
end

% step 10


