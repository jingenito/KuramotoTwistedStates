%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2nd order Kuramoto Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = kuramoto_2(u,a,w,k,n,G)

pos = u(1:n);
vel = u(n+1:end);

u_mat = repmat(pos,1,n);
F = (k/n)*sum(G.*sin(u_mat' - u_mat),2);
output = [vel; -a*vel + w + F];
