%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input the small world graph parameters and output
% the index of the minimum eigenvalue.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function min = getMinEigenvalueIndex(p,r)
min = 2*r + p - 4*r*p;

for k = 2:6
    tmp = (pi * (k-1))^(-1) * (1-2*p) * sin(2*pi*(k-1)*r);
    if tmp < min 
        min = tmp;
    end
end