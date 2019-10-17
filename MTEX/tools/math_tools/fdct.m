function C = fdct(x)
    
N = length(x);
y = zeros(1,2*N);
y(1:N) = fliplr(x);
y(N+1:2*N) = x;
Y = fft(y);
k=0:N-1;
C = real(exp(-1j.* pi.*k./(2*N)).*Y(1:2:end));

return


%%


x = zeros(1024,1);
x(4) = 1;
y = idctt(x);

plot(sqrt(1024/2)*y)


a = dctt(y);


a(1:10)


%%

fun = @(x) cos(x);

N = 1024;
x = ((0:(N-1))+0.5)/N*pi;

A = dctt(fun(x)) / sqrt(N/2);

A(1:10)



%% 

k = 10;
N = 2^k;
j = 2;
%fun = @(x) ones(size(x));
fun = @(x) x.^1;
%fun = @(x) [zeros(1,j) 1] * legendre0(j,x);

x = cos(((0:(N-1))+0.5)/N*pi);
B = dctt(fun(x)) / sqrt(N/2);

B(1:10)

%plan = fptmex('init',k+1,0);
%fptmex('precompute',plan,S2Kernel.alpha,S2Kernel.beta,S2Kernel.gamma,0);

%A = fptmex('trafo',plan,B.',N-1,0);

%fptmex('finalize',plan);

n = 0:(N-1);
A = cheb2leg(B./sqrt(2));
      
A(1:10)



S2K = S2Kernel(A)
plot(S2K)
hold on
xx = linspace(-1,1);
plot(xx,fun(xx))
hold off

%%
fun = @(x) x.^5;
t =chebfun(fun);
B = chebcoeffs(t);
A = cheb2leg(B)

S2K = S2Kernel(A)
plot(S2K)
hold on
xx = linspace(-1,1);
plot(xx,fun(xx))
hold off


%%


N=5
bma = 2;
c=zeros(N,2);

c(1:2:N,1)=(2./[1 1-(2:2:(N-1)).^2 ])';
c(2,2)=1;
f=real(ifft([c(1:N,:);c((N-1):-1:2,:)]));
w=bma * ([f(1,1); 2*f(2:(N-1),1); f(N,1)])/2
x=0.5*(N*bma*f(1:N,2))






