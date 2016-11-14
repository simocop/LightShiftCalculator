function dz = MakeDipoleMatrix(S,element,delta)
% Calculates dipole transition matrix for light parallel with the
% quantisation axis. This can be rotated using the rotation matrices to
% find the other Cartesian components.
% Simon Coop
% Last edited 26/09/2016

if nargin == 1
    element = 1;
    delta = 0;
end

load Rb87_data.mat

DipMat(element) = (1 + delta)*DipMat(element);

ns = length(S(:,1));
dz = zeros(ns,ns);

for i1 = 1:ns
    for i2 = 1:ns
        l = S(i1,:);
        lp = S(i2,:);
        
        n = l(2); 
        L = l(3);
        J = l(4);
        F = l(5);
        M = l(6);
        
        np = lp(2); % p for prime
        Lp = lp(3);
        Jp = lp(4);
        Fp = lp(5);
        Mp = lp(6);
        
        if DipMat(n,L+1,J+1/2,np,Lp+1,Jp+1/2) ...
                && abs(Mp - M) == 0 ...
                && abs(Fp - F) <= 1 ...
                && abs(Lp - L) <= 1 ...
                && abs(Jp - J) <= 1 ...
                && ~(M == 0 && (Fp - F) == 0)
                
            g = GeoFactor(I,J,F,M,Jp,Fp,Mp);
            dz(i1,i2) = g*DipMat(n,L+1,J+1/2,np,Lp+1,Jp+1/2);
        end
    end
end
assert(nnz(dz.*conj(dz)') == 0);
dz = (dz + conj(dz)');

end

function p = GeoFactor(I,J,F,M,Jp,Fp,Mp)
t1 = (-1)^(M+J+I);
t2 = sqrt((2*F+1)*(2*Fp+1));
t3 = w3j(Fp,Mp,1,0,F,-M); % q = 0
t4 = w6j(J,Jp,1,Fp,F,I);
p = t1*t2*t3*t4;
end