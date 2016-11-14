function w = w6j(j1,j2,j3,j4,j5,j6)

% W = W6J(j1,j2,j3,j4,j5,j6)
%
% Calculates Wigner 6-j Symbol using Racah W-coefficients
%       { j1 j2 j3 }
%       { j4 j5 j6 }
%
% (j1,j2,j3),(j1,j5,j6),(j4,j2,j6) and (j4,j3,j5) must satisfy triangle
% condition |j2 - j3|<=j1<=j2 + j3
%
% All square roots use logarithms to cope with large factorials
%
% Equations (7.25-7) in 'Angular Momentum:Techniques in Quantum Mechanics'
% by V. Devanathan, Kluwer Academic Publishers (1999).
%
%J. Pritchard Durham University 2009

if(triangle(j1,j2,j3)==1)
    error(sprintf('j1,j2,j3 must satisfy triangle relation:\n\n\t\t|j2 - j3|<=j1<=j2 + j3'));
elseif(triangle(j1,j5,j6)==1)
    error(sprintf('j1,j5,j6 must satisfy triangle relation:\n\n\t\t|j5 - j6|<=j1<=j5 + j6'));
elseif(triangle(j4,j2,j6)==1)
    error(sprintf('j4,j2,j6 must satisfy triangle relation:\n\n\t\t|j2 - j6|<=j4<=j2 + j6'))
elseif(triangle(j4,j3,j5)==1)
    error(sprintf('j4,j3,j4 must satisfy triangle relation:\n\n\t\t|j3 - j5|<=j4<=j3 + j5'));
else
    w = (-1)^(j1+j2+j4+j5)*racah(j1,j2,j5,j4,j3,j6);
end

%Reduce times need to type factorial!
function y=f(x)
y=factorial(x);

%Racah Coefficient
function W = racah(a,b,c,d,e,f)
W = tri(a,b,e)*tri(c,d,e)*tri(a,c,f)*tri(b,d,f);

xstart = max([-1,(a+b+e),(c+d+e),(a+c+f),(b+d+f)]); 
xstop = min([(a+b+c+d),(a+d+e+f),(b+c+e+f)]);
xsum=0;
for(x=xstart:xstop)
%     xsum = xsum + ((-1)^(x+a+b+c+d)*f(x+1)/...
%         (f(x-a-b-e)*f(x-c-d-e)*f(x-a-c-f)*f(x-b-d-f)...
%         *f(a+b+c+d-x)*f(a+d+e+f-x)*f(b+c+e+f-x)));
    xsum = xsum + (-1)^(x+a+b+c+d)...
        *exp(lgf(x+1)-lgf(x-a-b-e)-lgf(x-c-d-e)...
        -lgf(x-a-c-f)-lgf(x-b-d-f)-lgf(a+b+c+d-x)...
        -lgf(a+d+e+f-x)-lgf(b+c+e+f-x));
end
W = W.*xsum; 

%Triangle Coefficient
function t=tri(a,b,c)
% t = sqrt(f(a+b-c)*f(a-b+c)*f(-a+b+c)/f(a+b+c+1));
t = exp(0.5*(lgf(a+b-c)+lgf(a-b+c)...
    +lgf(-a+b+c)-lgf(a+b+c+1)));

%Triangle relations test |a-b|<=c<=a+b
function test=triangle(a,b,c)
test=(a<abs(b-c))||(a>(b+c));

%Stirlings approximation ln(n!) = nln(n)-n+0.5ln(2pin)
function y=lgf(x)
if(x<170)
    y=log(factorial(x));
else
    y=x*log(x)-x+0.5*log(2*pi*x)+1/(12*x)-1/(360*x^3)+1/(1260*x^5)...
        -1/(1680*x^7)+1/(1188*x^9);
end