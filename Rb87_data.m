%% Generate variables for Stark shift calculation
% Simon Coop, ICFO
% Last edited 31/01/2017

clear

if exist('Rb87_data.mat')
    delete Rb87_data.mat
end

%% Constants
h = 6.626070040e-34; % Js
ec = 1.6021766208e-19; % C
a0 = 5.2917721067e-11; % m, Bohr radius
c = 299792458; % m/s
e0 = 8.8541878176e-12;
I = 3/2; % Nuclear spin of Rb87

%%

nmax = 6;
Lmax = 3;
Jmax = 7/2;

%% Energies from NIST ASD

f = NaN(nmax,Lmax+1,Jmax+1/2);

n = 5,L = 0, J = 1/2;
f(n,L+1,J+1/2) = 0;

n = 5,L = 1, J = 1/2;
f(n,L+1,J+1/2) = 12578.95;

n = 5,L = 1, J = 3/2;
f(n,L+1,J+1/2) = 12816.545;

n = 4,L = 2, J = 3/2;
f(n,L+1,J+1/2) = 19355.649;

n = 4,L = 2, J = 5/2;
f(n,L+1,J+1/2) = 19355.203;

n = 6,L = 0, J = 1/2;
f(n,L+1,J+1/2) = 20132.51;

n = 4,L = 3, J = 5/2;
f(n,L+1,J+1/2) = 26792.118;

n = 4,L = 3, J = 7/2;
f(n,L+1,J+1/2) = 26792.092;

f = f*100*c;

%% Manual entry of dipole matrix elements

DipMat = zeros(nmax,Lmax+1,Jmax+1/2,nmax,Lmax+1,Jmax+1/2);

%%% From PRA 83 052508 2011 %%%
n = 5,L = 0, J = 1/2, np = 5, Lp = 1, Jp = 1/2;
DipMat(n,L+1,J+1/2,np,Lp+1,Jp+1/2) = 4.253;

n = 5,L = 0, J = 1/2, np = 5, Lp = 1, Jp = 3/2;
DipMat(n,L+1,J+1/2,np,Lp+1,Jp+1/2) = 6.003;

n = 5,L = 1, J = 1/2, np = 6, Lp = 0, Jp = 1/2;
DipMat(n,L+1,J+1/2,np,Lp+1,Jp+1/2) = 4.145;

n = 5,L = 1, J = 3/2, np = 6, Lp = 0, Jp = 1/2;
DipMat(n,L+1,J+1/2,np,Lp+1,Jp+1/2) = 6.047;

n = 5,L = 1, J = 1/2, np = 4, Lp = 2, Jp = 3/2;
DipMat(n,L+1,J+1/2,np,Lp+1,Jp+1/2) = 8.037;

n = 5,L = 1, J = 3/2, np = 4, Lp = 2, Jp = 3/2;
DipMat(n,L+1,J+1/2,np,Lp+1,Jp+1/2) = 3.628;

n = 5,L = 1, J = 3/2, np = 4, Lp = 2, Jp = 5/2;
DipMat(n,L+1,J+1/2,np,Lp+1,Jp+1/2) = 10.889;

%%% Directly from M. Safronova of U. Delaware. %%%
n = 4,L = 2, J = 3/2, np = 4, Lp = 3, Jp = 5/2;
DipMat(n,L+1,J+1/2,np,Lp+1,Jp+1/2) =  10.4185;

n = 4,L = 2, J = 5/2, np = 4, Lp = 3, Jp = 5/2;
DipMat(n,L+1,J+1/2,np,Lp+1,Jp+1/2) =  2.7861;

n = 4,L = 2, J = 5/2, np = 4, Lp = 3, Jp = 7/2;
DipMat(n,L+1,J+1/2,np,Lp+1,Jp+1/2) =  12.4601;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DipMat = DipMat*a0*ec;

%% Hyperfine A & B contants
% Experimental & calc'd values from PRA 83 052508 2011
A = ones(nmax,Lmax+1,Jmax+1/2)/1e6; % non-zero to prevent degeneracies

A(5,1,1) = 3417.341;
A(5,2,1) = 406.2;
A(5,2,2) = 84.845;
A(6,1,1) = 807.66;
A(4,3,2) = 25.1;
A(4,3,3) = -16.9;
A(5,3,2) = 16.57;
A(5,3,3) = -7.44;

A = A*1e6;

B = ones(nmax,Lmax+1,Jmax+1/2)/1e6;

B(5,2,2) = 12.52;
B(4,3,2) = 2.23; % calc'd
B(4,3,3) = 4.149; % from OPTICS LETTERS / Vol. 32, Iss. 19, p. 2810 / October 1, 2007

B = B*1e6;

%% Linewidths
% linewidth(n,L,J)
% from PRA 83 052508 2011
lifetime = zeros(nmax,Lmax+1,Jmax+1/2);

lifetime(5,1,1) = 0;
lifetime(5,1,2) = 0;
lifetime(6,1,1) = 45.4;
lifetime(7,1,1) = 88.3;
lifetime(8,1,1) = 161.9;
lifetime(9,1,1) = 271.7;
lifetime(10,1,1) = 426;

lifetime(5,2,1) = 27.4;
lifetime(5,2,2) = 26.0;
lifetime(6,2,1) = 122.5;
lifetime(6,2,2) = 112.4;
lifetime(7,2,1) = 277.8;
lifetime(7,2,2) = 255.2;
lifetime(8,2,1) = 501.0;
lifetime(8,2,2) = 464.2;
lifetime(9,2,1) = 809.3;
lifetime(9,2,2) = 753.8;

lifetime(4,3,2) = 83.0;
lifetime(4,3,3) = 89.4;
lifetime(5,3,2) = 240.3;
lifetime(5,3,3) = 231.6;
lifetime(6,3,2) = 258.0;
lifetime(6,3,3) = 247.4;
lifetime(7,3,2) = 339.5;
lifetime(7,3,3) = 327.0;
lifetime(8,3,2) = 468.8;
lifetime(8,3,3) = 452.7;

lifetime(4,4,3) = 60.7;
lifetime(4,4,4) = 60.8;
lifetime(5,4,3) = 109.3;
lifetime(5,4,4) = 109.4;

linewidth = 1./(lifetime*1e-9);
linewidth(isinf(linewidth)) = 0;

save Rb87_data