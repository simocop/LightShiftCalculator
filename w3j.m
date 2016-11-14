function w=w3j(j1,m1,j2,m2,j3,m3)

% W = W3J(j1,m1,j2,m2,j3,m3)
%
% j1,j2,j3 must satisfy triangle condition |j2 - j3|<=j1<=j2 + j3
%
% Wigner 3-j symbol is evaluated using equation found in 'Angular 
% Momentum: An Illustrated guide to Rotational Symmetries for
% Physical Systems', W. J. Thompson
%
%J. Pritchard Durham University 2009

%Check Triangular relation
if(((j3<abs(j1-j2))||(j3>(j1+j2))))
    error(sprintf('Addition of angular momentum requires triangle relation\n\t|j1-j2|<=j3<=j1+j2'));
%Evaluate w3j
else
    if(m1+m2+m3~=0)
        w=0;
    elseif(j2==0)
        w=1;
    else
%         w=(-1)^(j1-j2-m3)*...
%             sqrt((f(j3+j1-j2)*f(j3-j1+j2)*f(j1+j2-j3)*f(j3-m3)*f(j3+m3))/...
%             (f(j1+j2+j3+1)*f(j1-m1)*f(j1+m1)*f(j2-m2)*f(j2+m2)))...
%             *ksum(j1,m1,j2,m2,j3,m3);
        w=(-1)^(j1-j2-m3)*...
            exp(0.5*(lgf(j3+j1-j2)+lgf(j3-j1+j2)+lgf(j1+j2-j3)...
            +lgf(j3-m3)+lgf(j3+m3) - lgf(j1+j2+j3+1)...
            -lgf(j1-m1)-lgf(j1+m1)-lgf(j2-m2)-lgf(j2+m2)))...
            *ksum(j1,m1,j2,m2,j3,m3);

    end
end

%Summation performed for all values of k which give non-negative fs
function s=ksum(j1,m1,j2,m2,j3,m3)
s=0;
kmin=max([0,-j1+j2-m3]);
kmax=min([j3-j1+j2,j3-m3]);
for(k=kmin:kmax)
%     s =s+(-1)^(k+j2+m2)*f(j2+j3+m1-k)*f(j1-m1+k)/...
%         (f(k)*f(j3-j1+j2-k)*f(j3-m3-k)*f(k+j1-j2+m3));
    s = s +(-1)^(k+j2+m2)*exp(lgf(j2+j3+m1-k)+lgf(j1-m1+k)...
        -lgf(k)-lgf(j3-j1+j2-k)-lgf(j3-m3-k)-lgf(k+j1-j2+m3));
end

%Stirlings approximation ln(n!) = nln(n)-n+0.5ln(2pin)
function y=lgf(x)
if(x<170)
    y=log(factorial(x));
else
    y=x*log(x)-x+0.5*log(2*pi*x)+1/(12*x)-1/(360*x^3)+1/(1260*x^5)...
        -1/(1680*x^7)+1/(1188*x^9);
end