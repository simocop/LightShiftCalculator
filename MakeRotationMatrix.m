function [Fx,Fy] = MakeRotationMatrix(S)
% Creates the rotation matrices for a given list of F and M numbers
% Simon Coop
% Last edited 26/09/2016

F = S(:,5); % assuming S has been made with MakeStateVector();
M = S(:,6);
ns = length(F);

Fmat = repmat(F,1,ns);
Mmat = repmat(M,1,ns);

[Ij,Ii] = meshgrid(1:ns);

Fp = sqrt((Fmat - Mmat + 1).*(Fmat + Mmat)).*(Ii == Ij + 1);
Fm = conj(Fp');

Fx = 1/2*(Fp + Fm);
Fy = -1i/2*(Fp - Fm);

end
