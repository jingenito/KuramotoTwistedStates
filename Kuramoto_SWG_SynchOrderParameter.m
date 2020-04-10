%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates h_inf order parameter for the synchronized
% state of the Small World Graph model.
%
% Outputs a matrix where the first column consists of 
% the positive eigenvalues and the second column 
% consists of the negative eigenvalues.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = Kuramoto_SWG_SynchOrderParameter(u,p,r)
q = getMinEigenvalueIndex(p,r);
val(:,1) = exp(1i * 2 * pi * q .* u)';
val(:,2) = exp(-1i * 2 * pi * q .* u)';
