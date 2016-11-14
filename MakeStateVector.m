function S = MakeStateVector(degeneracy_shift)

% Make vector of energy levels in 87Rb

if nargin == 0
    degeneracy_shift = 0;
end

S = [];

load Rb87_data.mat

for n = 4:nmax
    for L = 0:Lmax
        for J = abs(L - 1/2) : L + 1/2
            for F = abs(J - I) : J + I
                e = f(n,L+1,J+1/2);
                if ~isnan(e)
                    de = GetHyperfineShiftHz(n,I,L,J,F,A,B);
                    et = e + de; % energy of level in Hz
                    for M = F:-1:-F
                        DGB = degeneracy_shift*M;
                        S = [S; et+DGB n L J F M e+de];
                        % add the extra e+de for later plotting of energy
                        % levels (don't want to plot literally every one,
                        % just down to F-levels)
                    end
                end
            end
        end
    end
end

S = sortrows(S);

end

function V_hfs = GetHyperfineShiftHz(n,I,L,J,F,A,B)
g =  F*(F+1) - I*(I+1) - J*(J+1);
num = 1.5*g*(g+1) - 2*I*(I+1)*J*(J+1);
if num == 0
    V_hfs = 0.5*A(n,L+1,J+1/2)*g;
else
    den = 2*I*(2*I-1)*2*J*(2*J-1);
    V_hfs = 0.5*A(n,L+1,J+1/2)*g + B(n,L+1,J+1/2)*num/den;
end

end
