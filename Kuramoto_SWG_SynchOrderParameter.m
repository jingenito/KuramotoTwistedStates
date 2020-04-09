%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates h_inf order parameter for the synchronized
% state of the Small World Graph model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = Kuramoto_SWG_SynchOrderParameter(u,beta,kappa,p,r)
q = getMinEigenvalueIndex(p,r);
val = sqrt(kappa / beta) * exp(1i * 2 * pi * q .* u);
