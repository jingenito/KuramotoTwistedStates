%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create a small world graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W_out = sw_graph(n,p,r)
% p = .2; %probability of forming a far edge, 0 < p < 1/2
% r = .3; %nearest-neighbor range 0 < r < 1/2

W_out = rand(n,n);
for i = 1:n %y-index
    x_value = (i-1)/n; %Translate matrix indices to values on unit square
    
    for j = 1:n %x-index
        y_value = (j-1)/n;
        
        if min(abs(2*pi*x_value-2*pi*y_value),2*pi-abs(2*pi*x_value-2*pi*y_value))<2*pi*r %Check if neighbors; note that matrix indices
            %run the opposite direction for the y-variable
            %Prob(W_out(i,j) < 1-p) = 1-p
            if W_out(i,j) < 1-p
                W_out(i,j) = 1;
            else
                W_out(i,j) = 0;
            end
        else
            if W_out(i,j) < p
                W_out(i,j) = 1;
            else
                W_out(i,j) = 0;
            end
        end
    end
end
W_upper = triu(W_out,1);

W_out = diag(zeros(n,1))+W_upper+W_upper';