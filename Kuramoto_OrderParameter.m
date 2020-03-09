%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input a vector and output the modulous of the complex order 
%% parameter for the Kuramoto Model.
%%
%% Using j as the compex number, sqrt(-1), to not be confused with the loop 
%% counter i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mod = Kuramoto_OrderParameter(x)
z = 0; %will hold sum of all the complex nubmers
for i=1:length(x)
    z = z + exp(1j*x(i)); 
end

z = z / length(x); %take the average of the sum
mod = abs(z); %output the modulous