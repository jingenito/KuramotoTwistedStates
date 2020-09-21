%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2nd Order Kuramoto Model
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

k = -100; %Coupling strength
a = 0.3; %alpha term on the first derivative

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Graph connectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G = sw_graph(n,.2,.4);   %Adjacency matrix of network connections
osc_list = 1:n; %for visualization of connectivity graph
sw_nodes = [cos(2*pi*osc_list/n)',sin(2*pi*osc_list/n)']; %graph node locations

figure(1)
gplot(G,sw_nodes,'.-'); %plot graph

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve Kuramoto Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Unit circle (to later plot)
figure(2)
x = linspace(0,2*pi,n);
x1 = cos(x);
x2 = sin(x);


%Plot initial conditions
subplot(1,2,1) %plot on plane         
plot(u_int,'.')
axis([1 n 0 2*pi])
    
subplot(1,2,2) %plot on unit circle
plot(cos(u_int),sin(u_int),'.',x1,x2)

pause()
 
sol_end_point = 200;
%Integrate model
disp('Solving ODE')
% ode45
[t,u]=ode45(@(t,y) kuramoto_2(y,a,w,k,n,G),[0,sol_end_point],[u_int; u_prime_int]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make movie of solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Interpolate time points for even movie 
t0 = 0:5:sol_end_point;
u0 = interp1(t,u,t0);

t= t0;
u = mod(u0, 2 * pi);
frame = 1;

figure(2)
for i=1:length(t)

    subplot(1,3,1)       
    plot(u(i,(1:n)),'.')
    axis([1 n 0 2*pi])
    
    subplot(1,3,2)
    plot(cos(u(i,(1:n))),sin(u(i,(1:n))),'.',x1,x2)
    
    subplot(1,3,3)
    plot(u(i,(n+1:end)),'.')
    axis([1 n -3/abs(a) 3/abs(a)])
    
    F(frame) = getframe(gcf);
    frame = frame+ 1;
    
    pause(0)
    clf
   
end

v=VideoWriter('kuramoto.avi');
open(v);
writeVideo(v,F);
close(v);
