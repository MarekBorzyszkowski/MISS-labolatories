mu=0.01; tau=100; b=0.5; beta=b/2/mu; L=1; M=100; N=2000; x=linspace(0,L,M)';
v=sqrt(tau/mu); dx=x(2)-x(1); dt=dx/v; p=(v*dt/dx)^2;
q=1+beta*dt;
u=1-beta*dt;
s=sin(3*pi*x+pi/4);
g=x-x.^2; 
h=zeros(1,N)- 3*pi*sqrt(2)/2;
r=zeros(1,N) + sqrt(2)/2; 

f(:,1)=s;
f(1,2)=r(2);
f(2:M-1,2)=p/2*(f(3:M,1)-2*f(2:M-1,1)+ f(1:M-2,1))+f(2:M-1,1)+u*dt*g(2:M-1);
f(M,2)=p*(f(M-1,1)- f(M,1)+dx*h(1)) + f(M,1) +u*dt*g(1);

for n=2:N-1
    f(1,n+1)=r(n+1);
    f(2:M-1,n+1)=p/q*(f(3:M,n)-2*f(2:M-1,n)+f(1:M-2,n))+2/q*f(2:M-1,n)-u/q*f(2:M-1,n-1);
    f(M,n+1)=2*p/q*(f(M-1,n)-f(M,n)+dx*h(n))+2/q*f(M,n)-u/q*f(M,n-1);
end

surf(1:N, 1:M, f)
shading interp;
xlabel('N'); ylabel('M'); zlabel('y(x, t)');
title('Drgania struny z odwróconymi warunkami brzegowymi');

y_eq = -3*pi*sqrt(2)/2*x + sqrt(2)/2;

figure;
plot(x, f(:, end), 'b-', 'DisplayName', 'Ostatnia iteracja symulacji');
hold on;
plot(x, y_eq, 'r-', 'DisplayName', 'Rozwiązanie równowagowe');
xlabel('x'); ylabel('y');
legend show;
title('Ostatnia iteracja vs rozwiązanie równowagowe');
