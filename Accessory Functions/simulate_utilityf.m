ufunc = @(x, theta)  -6*(x-theta).^2;
ufunc2 = @(x, theta)  -6.*(x-theta).^2-150;

ufunc3 = @(x, theta)  -6.*(x-theta).^2+0.5*(300-(x-80).^2);

x=[0:100];
y=ufunc2(x,30);
y2=ufunc(x,50);
y3=ufunc3(x,80);

plot(x,[y;y2;y3])