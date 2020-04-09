%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input a vector and an adjaceny matrix, and output the 
%% modified complex order parameter.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = Kuramoto_SWG_OrderParameter(x,W)

val = (exp(1i.*x) * W) / length(x);