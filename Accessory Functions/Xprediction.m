F = makedist('Beta',1,1);
n=20;
theta = (random(F,n,1));


fmin = @(x) (1-cdf(F,x)).^(n-1).*pdf(F,x).*n;
exmin = integral(@(x) x.*fmin(x),0,1)
fmax = @(x)(cdf(F,x)).^(n-1).*pdf(F,x).*n;
exmax = integral(@(x) x.*fmax(x),0,1)

fk= @(x,k) (factorial(n)./(factorial(k-1).*factorial(n-k))).*(cdf(F,x).^(k-1)).*((1-cdf(F,x)).^(n-k)).*pdf(F,x);
exmax2 = integral(@(x) x.*fk(x,n),0,1)


fkjoint = @(x,y,j,k)  (factorial(n)./(factorial(j-1).*factorial(k-j-1).*factorial(n-k))).* ...
    (cdf(F,x).^(j-1)).*((cdf(F,y)-cdf(F,x)).^(k-j-1)).* ...
    ((1-cdf(F,y)).^(n-k)).*pdf(F,x).*pdf(F,y);

ymax = @(x) x;
test=integral2(@(x,y) fkjoint(x,y,19,20),0,1,0,ymax)

x=linspace(0,1)';
xy= x*x'

asdf=fkjoint(x,x,9,10)