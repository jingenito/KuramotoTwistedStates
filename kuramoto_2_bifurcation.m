%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create a bifurcation diagram of the 2nd Order Kuramoto Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear         % clear any variables
clf           % clears any figures already up

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 1000; %number of oscillators
w = randn(n,1); %Random internal frequencies chosen from normal distribution
u_int = rand(n,1)*2*pi; %Random initial conditions
u_prime_int = randn(n,1); %random initial velocity conditions

%initialize parameters for synchronization vector calculation
p = 0.2;
r = 0.4;

%going to use the same connections for each (K,a) pair
G = sw_graph(n,p,r);   %Adjacency matrix of network connections

KVec = linspace(-40,0,100);
a = 10; %inertia term

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop over all (K,a) pairs and track the long term behavior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%storing result to calculate the l2 norm on each iteration
twistedStateVec = exp(1i .* linspace(0,1,n));
h_sync_vecs = Kuramoto_SWG_SynchOrderParameter(twistedStateVec,p,r); %transposing once for l-2 norm
h_sync1 = h_sync_vecs(:,1); %positive eigenvalues
h_sync2 = h_sync_vecs(:,2); %negative eigenvalues
h_sync = (1/n) * sqrt(h_sync1.^2 + h_sync2.^2); %calculte the l2 norm
%projected in both twisted state directions

Z = zeros(1,length(KVec)); %preallocating memory for optimization
for i=1:length(KVec)
    [t,u]=ode45(@(t,y) kuramoto_2(y,a,w,KVec(i),n,G),[0,50],[u_int; u_prime_int]);
    
    theta = u(length(t), 1:n); %get the theta vector 
    h = Kuramoto_SWG_OrderParameter(theta,G); %caclulate the complex order parameter
    Z(i) = abs(h * h_sync) / length(h)^2;
    
    %set initial conditions to the previous solution
    u_int = theta;
    u_prime_int = u(length(t), n+1:end);
end

% Z1 = zeros(1,length(KVec)); %preallocating memory for optimization
% for i=length(KVec):-1:1
%     [t,u]=ode45(@(t,y) kuramoto_2(y,a,w,KVec(i),n,G),[0,50],[u_int; u_prime_int]);
%     
%     theta = u(length(t), 1:n); %get the theta vector 
%     h = Kuramoto_SWG_OrderParameter(theta,G); %caclulate the complex order parameter
%     Z1(i) = abs(h * h_sync) / length(h)^2;
% end

plot(KVec,Z,'.','Color','b')
% hold on
% plot(KVec,Z1,'.','Color','r')
% hold off