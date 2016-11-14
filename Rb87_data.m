%% Generate variables for Stark shift calculation
% Simon Coop, ICFO
% Last edited 26/09/2016

clear

if exist('Rb87_data.mat')
    delete Rb87_data.mat
end

%% Constants
h = 6.626070040e-34; % Js
e = 1.6021766208e-19; % C
a0 = 5.2917721067e-11; % m, Bohr radius
c = 299792458; % m/s
e0 = 8.8541878176e-12;
I = 3/2; % Nuclear spin of Rb87

%% Maximum value of n and L to include in the list of states.
% The data file has data up to n = 11 and L = 3. 
max_desired_n = 6;
max_desired_l = 3;

%% Read in dipole matrix elements into a big matrix
alldata = importdata('Rbdata.txt',' ',19);
data = alldata.data;
nmax = 0;
Lmax = 0;
Jmax = 0;
for ii = 1:length(data(:,1))
    
    if ~isnan(data(ii,1)) && data(ii,7) < 0
        n = data(ii,1);
        np = data(ii,3);
        
        switch data(ii,2)
            case -1
                L = 0;
                J = 1/2;
            case 1
                L = 1;
                J = 1/2;
            case -2
                L = 1;
                J = 3/2;
            case 2
                L = 2;
                J = 3/2;
            case -3
                L = 2;
                J = 5/2;
            case 3
                L = 3;
                J = 5/2;
            case -4
                L = 3;
                J = 7/2;
        end
        
        switch data(ii,4)
            case -1
                Lp = 0;
                Jp = 1/2;
            case 1
                Lp = 1;
                Jp = 1/2;
            case -2
                Lp = 1;
                Jp = 3/2;
            case 2
                Lp = 2;
                Jp = 3/2;
            case -3
                Lp = 2;
                Jp = 5/2;
            case 3
                Lp = 3;
                Jp = 5/2;
            case -4
                Lp = 3;
                Jp = 7/2;
        end
        
        if (L <= max_desired_l && Lp <= max_desired_l) ...
                && (n <= max_desired_n && np <= max_desired_n)
            
            f1 = data(ii,6)*100*c; % convert wavenumbers to Hz
            f2 = data(ii,7)*100*c;
            
            f(n,L+1,J+1/2) = f1;
            f(np,Lp+1,Jp+1/2) = f2;
            
            DipMat(n,L+1,J+1/2,np,Lp+1,Jp+1/2) = data(ii,9)*e*a0;
            
            if n > nmax
                nmax = n;
            end
            if L > Lmax
                Lmax = L;
            end
            if J > Jmax
                Jmax = J;
            end
            
            if np > nmax
                nmax = np;
            end
            if Lp > Lmax
                Lmax = Lp;
            end
            if Jp > Jmax
                Jmax = Jp;
            end
        end
    end
end

if nmax == 11
    DipMat(nmax,:,:,:,:,:) = 0;
end

f(f == 0) = NaN;

%% Hyperfine A & B contants
% Experimental & calc'd values from PRA 83 052508 2011
A = ones(nmax,Lmax+1,Jmax+1/2)/1e6; % non-zero to prevent degeneracies

A(5,1,1) = 3417.341;
A(5,2,1) = 406.2;
A(5,2,2) = 84.845;
A(6,1,1) = 807.66;
A(6,2,1) = 132.56;
A(6,2,2) = 27.7;
A(7,1,1) = 319.759;
A(7,2,1) = 59.32;
A(7,2,2) = 12.57;
A(8,1,1) = 159.2;
A(8,2,1) = 32.12;
A(8,2,2) = 6.57;
A(9,1,1) = 90.9;
A(9,2,1) = 19.59; % calc'd
A(9,2,2) = 4.02; % calc'd
A(4,3,2) = 25.1;
A(4,3,3) = -16.9;
A(5,3,2) = 16.57;
A(5,3,3) = -7.44;
A(6,3,2) = 8.76;
A(6,3,3) = -3.3;
A(7,3,2) = 5.04;
A(7,3,3) = -1.78;
A(8,3,2) = 3.13;
A(8,3,3) = -1.06;
A(9,3,2) = 2.07;
A(9,3,3) = -0.69;

A = A*1e6;

B = ones(nmax,Lmax+1,Jmax+1/2)/1e6;

B(5,2,2) = 12.52;
B(6,2,2) = 3.953;
B(7,2,2) = 1.762;
B(8,2,2) = 0.935;
B(9,2,2) = 0.55;
B(4,3,2) = 2.23; % calc'd
B(4,3,3) = 4.149; % from OPTICS LETTERS / Vol. 32, No. 19 / October 1, 2007
B(5,3,2) = 0.913; % calc'd
B(5,3,3) = 1.29; % calc'd
B(6,3,2) = 0.53;
B(6,3,3) = 0.623; % calc'd
B(7,3,2) = 0.26;
B(7,3,3) = 0.343; % calc'd
B(8,3,2) = 0.17;
B(8,3,3) = 0.207;% calc'd
B(9,3,2) = 0.11;
B(9,3,3) = 0.133; % calc'd

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