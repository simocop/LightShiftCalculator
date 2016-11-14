% MATLAB script to calculate light shifts in Rb87 using a result from
% Floquet's theorem. 
% Simon Coop ICFO
% 26/09/2016

clear

tic

% Load Rb87 data
Rb87_data

% Make a vector of all states with energies and quantum numbers
% S(1,:) = [energy1 n1 L1 J1 F1 M1];
% S(2,:) = [energy2 n2 L2 J2 F2 M2];
% etc

S = MakeStateVector(0);

% Include 4D, 4F, 5S, 5P, and 6S states only.
S = S([1:80 145:200],:);
   
dz = MakeDipoleMatrix(S);

[Fx,Fy] = MakeRotationMatrix(S);
ns = length(S);

figure(1)
clf
PlotAllLevels(S)

%% Quantisation axis is in the z-direction
% Vx, Vy, Vz are the Cartesian components of the interaction part of the 
% Hamiltonian
dx = expm(-1i*Fy*pi/2)*dz*expm(1i*Fy*pi/2);
dy = expm(1i*Fx*pi/2)*dz*expm(-1i*Fx*pi/2);

V(:,:,1) = dx;
V(:,:,2) = dy;
V(:,:,3) = dx;
V(:,:,4) = dy;

p1529 = linspace(0,1.5e7,50);
p1560 = 3.0e9;

for I = 1:length(p1529)
    I
    
    r1 = 51;
    r2 = 50;
    w1529 = 1529.27e-9;
    w1560 = r1/r2*w1529;
    lambda = [w1560 w1560 w1529 w1529];
    power = [p1560 p1560 p1529(I) p1529(I)]/2; % divide by 2 as have two "waves" at each wavelength.
    steps_per_cycle = 200;
    period = lambda(3)/c;
    dt = period/steps_per_cycle;
    t0 = r1*period;
    nsteps = r1*steps_per_cycle;
    phase = [0 0.1 0 pi];
    
    [floquet_shifts,floquet_energies_t] = FloquetShiftCalc(S,V,lambda,power,nsteps,dt,phase);
    floquet_energies(:,I) = floquet_energies_t;
    
end


%%

e = S(:,1);
n = S(:,2);
L = S(:,3);
J = S(:,4);
F = S(:,5);
M = S(:,6);

refE = e(26);

lvlinds = find(n == 5 & L == 1 & J == 3/2);

shifts_f = (floquet_energies(lvlinds,:) - refE)/1e6

figure(2)
plot(p1529,sort(shifts_f),'b','LineWidth',2)
set(gcf,'Color','w')
xlabel('1529 intensity (Wm^{-2})')
ylabel('Detuning from 2 ^\rightarrow{3} transition (MHz)')
ylim([-700 0])
