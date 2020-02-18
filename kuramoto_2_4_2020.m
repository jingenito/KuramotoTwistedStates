%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Simulation of Kuramoto system %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                           % 
%                                                           %
% Simulate u'_i = w_i + K/N \sum_{j=1}^N G_{ij}sin(u_j-u_i) %
% Coupling matrix G - refers to W_graph_struc.m code        % 
%                                                           %
% Method: uses ode45 to solve the ODE system                % 
%                                                           % 
% Output: Creates movie of evolution of u over time,        %
%         represented both on plane and on unit circle      %
%                                                           % 
%                                                           % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear         % clear any variables
clf           % clears any figures already up

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 50; %number of oscillators
w = randn(n,1); %Random internal frequencies chosen from normal distribution
u_int = rand(n,1)*2*pi; %Random initial conditions

k = 20; %Coupling strength


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
x = linspace(0,2*pi,100);
x1 = cos(x);
x2=sin(x);


%Plot initial conditions
    subplot(1,2,1) %plot on plane         
        plot(u_int,'.')
        axis([1 n 0 2*pi])
    
        subplot(1,2,2) %plot on unit circle
            plot(cos(u_int),sin(u_int),'.',x1,x2)
            
              %  j= j+1;
                
            pause()
 
%Integrate model
disp('Solving ODE')
% ode45
[t,u]=ode45(@(t,y) kuramoto(y,w,k,n,G),[0,500],u_int);


disp('Plotting solution')

%Interpolate time points for even movie 

t0 = 0:5:500;
u0 = interp1(t,u,t0);

t= t0;
u=u0;

u=mod(u,2*pi); %Return to unit circle


%r = abs((1/n)*(sum(exp(1i*u),2)));   % order parameter



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make movie of solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

frame = 1;

figure(2)
for i=1:length(t)
   % if mod(i,5)==0
        %subplot(1,2,1)         % plot order parameter
        %    plot(t(1:i),r(1:i))

        subplot(1,2,1)        % plot circle
             plot(u(i,:),'.')
             axis([1 n 0 2*pi])
      
             subplot(1,2,2)
            plot(cos(u(i,:)),sin(u(i,:)),'.',x1,x2)
            
             F(frame) = getframe(gcf);
                frame = frame+ 1;
              %  j= j+1;
                
            pause(0)
        clf
    %end
end

                
v=VideoWriter('kuramoto.avi');
open(v);
writeVideo(v,F);
close(v);  

function output = kuramoto(u,w,k,n,G)

u_mat = repmat(u,1,n);
output = w + (k/n)*sum(G.*sin(u_mat' - u_mat),2);

end

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
end