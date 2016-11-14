function [FloquetShifts,energies,evecs_srt] = FloquetShiftCalc(S,V,wavelength,intensity,nsteps,dt,phase)

% Finds energies of dressed states of an Rb87 atom subject to an 
% oscillating electric field. 
% S is list of free-atom states made with MakeStateVector()
% V is a m x n x o matrix of dipole matrix elements, where the third 
% dimension o corresponds to the number of polarisations and/or frequencies
% used.
% WAVELENGTH is a vector containing the wavelength of each polarisation
% and/or frequency used. 
% INTENSITY is a vector containing the intensity of each polarisation/wavelength.
% NSTEPS is the number of steps in the calculation, more steps means a more
% accurate calculation. For light around 1560 nm it seems pretty well
% converged above about 200 steps. 
% DT is the period of the total electric field. 
% PHASE is a vector with the relative phase of each of the 
% Note that dimensions of wavelength, intens, (third dimension of) V, and 
% phase should match. For example if you want to do a multi-frequency 
% simulation, you could have wavelength = [w1 w2], and phase = [0 0]
% However if you want circularly polarised light, have wavelength = [w1 w1]
% and phase = [0 pi/2].

load Rb87_data.mat

[~,~,np] = size(V);

assert(np == length(wavelength) && ...
    np == length(intensity) && ...
    np == length(phase))

t0 = dt*nsteps;

ns = length(S); % Number of states used in calc
U = eye(ns);
H0 = diag(S(:,1));
E0 = sqrt(2*intensity/(c*e0));

flaser = c./wavelength;

prev_pc = 0;

Vh = V/h;

coeff = -1i*2*pi*dt;

for ii = 1:nsteps
    H = sparse(zeros(size(U)));
    for jj = 1:np
        E = E0(jj)*sin(2*pi*flaser(jj)*(ii-1)*dt + phase(jj));
        H = H + E*Vh(:,:,jj);
    end
    H = coeff*(H0 + H);
    U = expm(H)*U;
    
    % Print out progress on command line
    pc = round(100*ii/nsteps);
    if pc ~= prev_pc
        disp([num2str(round(ii/nsteps*100)) ' %'])
        prev_pc = pc;
    end
end

ee = S(:,1); % eigenenergies of unperturbed eigenstates (in Hz)

[evecs,evals] = eig(U);

evals = diag(evals);

[~,princ_comp] = max(abs(evecs));

[~,perm] = sort(princ_comp);
evals_srt = evals(perm);
evecs_srt = evecs(:,perm);
ee_srt = ee(perm);

y = unique(princ_comp);
keep_inds = find(y == princ_comp(perm(1:length(y))));
keep_inds(end) = [];
if nargout > 1
    keep_inds = 1:length(S);
end

shift_arg = evals_srt(keep_inds).*exp(1i.*2*pi.*ee(keep_inds).*t0); % note ee in Hz

FloquetShifts = -angle(shift_arg)/(2*pi*t0);

% If need all energies
energies = -angle(evals(perm).*exp(1i.*2*pi.*ee.*t0))/(2*pi*t0) + ee;