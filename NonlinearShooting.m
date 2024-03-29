%% Nonlinear Shooting With Newton's Method

%% Input Information
a = 1;          % left endpoint
b = 2;          % right endpoint
alpha = 0;      % boundary condition at left endpoint
beta = log(2);  % boundary condition at right endpoint
N = 10;         % number of subintervals
tol = 1e-4;     % tolerance
M = 100;        % maximum number of iterations (not mentioned here)

f = @(x,y,y_prime) -exp(-2*y);
partialf_partialy = @(x,y,y_prime) 2*exp(-2*y);
partialf_partialy_prime = @(x,y,y_prime) 0;

% exact solutions
y = @(x) log(x);

%% Do the method

h = (b-a)/N;
j = 1;
TK = (beta - alpha)/(b-a);

fprintf('x \t\t\t w1 \t\t w2 \t\t y(x) \n')

while(j <= M)
    w(1,1) = alpha;
    w(2,1) = TK;
    u1 = 0;
    u2 = 1;

    for i=2:N+1
        x = a + (i-2)*h;

        k(1,1) = h * w(2,i-1);
        k(1,2) = h * f( x, w(1,i-1), w(2,i-1) );

        k(2,1) = h * ( w(2,i-1) + k(1,2)/2 );
        k(2,2) = h * f( x + h/2, w(1,i-1) + k(1,1)/2, w(2,i-1) + k(1,2)/2 );

        k(3,1) = h * ( w(2,i-1) + k(2,2)/2 );
        k(3,2) = h * f( x + h/2, w(1,i-1) + k(2,1)/2, w(2,i-1) + k(2,2)/2 );

        k(4,1) = h * ( w(2,i-1) + k(3,2) );
        k(4,2) = h * f( x + h, w(1,i-1) + k(3,1), w(2,i-1) + k(3,2) );

        w(1,i) = w(1,i-1) + ( k(1,1) + 2*k(2,1) + 2*k(3,1) + k(4,1) )/6;
        w(2,i) = w(2,i-1) + ( k(1,2) + 2*k(2,2) + 2*k(3,2) + k(4,2) )/6;

        k_prime(1,1) = h * u2;
        k_prime(1,2) = h * ( partialf_partialy( x, w(1,i-1), w(2,i-1) ) * u1 ...
            + partialf_partialy_prime( x, w(1,i-1), w(2,i-1) )*u2 );

        k_prime(2,1) = h * ( u2 + k_prime(1,2)/2 );
        k_prime(2,2) = h * ( partialf_partialy( x + h/2, w(1,i-1), w(2,i-1) ) * ( u1 + k_prime(1,1)/2 ) ...
            + partialf_partialy_prime( x + h/2, w(1,i-1), w(2,i-1) ) * ( u2 + k_prime(1,2)/2 ) );

        k_prime(3,1) = h * ( u2 + k_prime(2,2)/2 );
        k_prime(3,2) = h * ( partialf_partialy( x + h/2, w(1,i-1), w(2,i-1) ) * ( u1 + k_prime(2,1)/2 ) ...
            + partialf_partialy_prime( x + h/2, w(1,i-1), w(2,i-1) ) * ( u2 + k_prime(2,2)/2 ) );

        k_prime(4,1) = h * ( u2 + k_prime(3,2) );
        k_prime(4,2) = h * ( partialf_partialy( x + h, w(1,i-1), w(2,i-1) ) * ( u1 + k_prime(3,1) ) ...
            + partialf_partialy_prime( x + h, w(1,i-1), w(2,i-1) ) * ( u2 + k_prime(3,2) ) );

        u1 = u1 + ( k_prime(1,1) + 2*k_prime(2,1) + 2*k_prime(3,1) + k_prime(4,1) )/6;
        u2 = u2 + ( k_prime(1,2) + 2*k_prime(2,2) + 2*k_prime(3,2) + k_prime(4,2) )/6;
    end

    if(abs(w(1,N+1) - beta) <= tol)
        for i = 1:N+1
            x = a + (i-1) * h;
            fprintf('%f \t %f \t %f \t %f \n',x,w(1,i),w(2,i),y(x))
        end

        fprintf('The procedure is complete for j = %d \n',j)
        break;
    end

    TK = TK - ( w(1,N+1) - beta )/u1;

    j = j+1;
end






