%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input a vector and an adjaceny matrix, and output the 
%% modified complex order parameter.
%%
%% Using j as the compex number, sqrt(-1), to not be confused with the loop 
%% counter i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = Kuramoto_SWG_OrderParameter(x,W)
z = 0; %will hold sum of all the complex nubmers
for i=1:length(x)
    for j=1:length(x)
        z = z + W(i,j)*exp(1j*x(i)); 
    end
end

z = z / length(x)^2; %take the average of the sum
val = abs(z); %output the modulous