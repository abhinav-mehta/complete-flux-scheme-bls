clear all;
clc;
sigma = 0.05; r = 0.2; K = 1;
X = 1; T = 100;
N = 50; M = 500;
x(N) = X; x(1) = 0;
t(M) = T; t(1) = 0;
dx = X/(N-1); dt = T/(M-1);
for i = 2:N
    x(i) = x(i-1) + dx;
end
for i = 2:M 
    t(i) = t(i-1) + dt;
end
V(N,M) = 0; 
for j = 1:M
    for i = 1:N
        V(i,j) = 0;
    end
end
for i = 2:N 
    V(i,M) = max(-x(i) + K,0);
end
for i = 1:M
    V(1,i) = K*exp(-r*(T-t(i)));
end
l(N) = 0;
l(1) = 2*(sigma^2 - r)/(sigma^2*x(2));
for i = 2:N 
    l(i) = 2*(sigma^2 - r)/(sigma^2*x(i)) ;
end
e(N)=0;
for i = 1:N
    e(i) = -(sigma^2*x(i)^2)/2 ;
end
p(N-1) = 0; alpha(N-1) = 0;
beta(N-1) =0; gamma(N-1) =0; 
for i = 1:N-1 
    p(i) = dx*(l(i)/2 + l(i+1)/2);
    alpha(i) = ((N-1)/X)*b((-1)*p(i))*(2/(l(i)+l(i+1)))*(w((-1)*p(i))*l(i) + w(p(i))*l(i+1))*(w((-1)*p(i))*e(i) + w(p(i))*e(i+1));
    beta(i) = ((N-1)/X)*b(p(i))*(2/(l(i)+l(i+1)))*(w((-1)*p(i))*l(i) + w(p(i))*l(i+1))*(w((-1)*p(i))*e(i) + w(p(i))*e(i+1));
    gamma(i) = max((0.5 - w(p(i))) , 0);
end
A(N-2,N-2) = alpha(N-1);
for i = 1:N-3 
    A(i,i) = alpha(i+1) + beta(i+1); 
    A(i,i+1) = (-1) * beta(i+1); 
    A(i+1,i) = (-1) * alpha(i+1); 
end
B(N-2,N-2) = 1 - gamma(N-1);
for i = 1:N-3 
    B(i,i) = 1 - gamma(i+1);
    B(i+1,i) = gamma(i+1) ; 
end
B = B*dx;
b_vector(N-2,M) = 0;
for i = 1:M 
    b_vector(1,i) = (-dx*gamma(1)*(K*r*exp(-r*(T-t(i)))) + alpha(1)*V(1,i) + dx*gamma(1)*(-(sigma^2 - r)*V(1,i)));
    %%b_vector((N-2),i) = beta(N-1)*V(N,i);
end
theta = 0.5; 
A_matrix = B - dt*(1-theta)*A + (r-sigma^2)*dt*(1-theta)*B;
A_matrix_inv = inv(A_matrix);
B_matrix(N-2,M-2) = 0;
source(N-2,1) = 0;
for i = (M-1) : -1 : 1
    source(:,1) = -(sigma^2 - r)*V(2:N-1,i+1);
    B_matrix(:,i) = -dt*((1-theta)*b_vector(:,i) + theta*b_vector(:,i+1)) -dt*theta*B*source(:,1) + B*V(2:(N-1),i+1) + dt*theta*A*V(2:(N-1),i+1);
    V(2:(N-1), i) = A_matrix\B_matrix(:,i); 
end

V2(N,M) = 0;
tm(M)= 0;
for i=1:N
    for j = 1: M
        tm(j) = abs(T-t(j));
        [call,put] = blsprice(x(i), K, r, tm(j), sigma, 0);
        V2(i,j) = put;   
    end 
end