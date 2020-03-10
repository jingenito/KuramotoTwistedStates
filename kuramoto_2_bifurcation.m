%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create a bifurcation diagram of the 2nd Order Kuramoto Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear         % clear any variables
clf           % clears any figures already up

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 50; %number of oscillators
w = randn(n,1); %Random internal frequencies chosen from normal distribution
u_int = rand(n,1)*2*pi; %Random initial conditions
u_prime_int = randn(n,1); %random initial velocity conditions

%going to use the same connections for each (K,a) pair
G = sw_graph(n,.2,.4);   %Adjacency matrix of network connections

KVec = linspace(0,12,100);
a = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop over all (K,a) pairs and track the long term behavior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = []; %vector to track long term behavior
for i=1:length(KVec)
    [t,u]=ode45(@(t,y) kuramoto_2(y,a,w,KVec(i),n,G),[0,50],[u_int; u_prime_int]);
    
    r = u(length(t), 1:n); %get the theta vector 
    Z(i) = Kuramoto_OrderParameter(r); %caclulate the complex order parameter
    
    %set initial conditions to the previous solution
    u_int = r;
    u_prime_int = u(length(t), n+1:end);
end

Z1 = [];
for i=length(KVec):-1:1
    [t,u]=ode45(@(t,y) kuramoto_2(y,a,w,KVec(i),n,G),[0,50],[u_int; u_prime_int]);
    
    r = u(length(t), 1:n); %get the theta vector 
    Z1(i) = Kuramoto_OrderParameter(r); %caclulate the complex order parameter

end

plot(KVec,Z,'o', KVec,Z1,'o')