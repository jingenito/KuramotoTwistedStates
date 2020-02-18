clear         % clear any variables
clf           % clears any figures already up

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 50; %number of oscillators
w = randn(n,1); %Random internal frequencies chosen from normal distribution
u_int = rand(n,1)*2*pi; %Random initial conditions

k = 20; %Coupling strength
a = 10; %alpha term on the first derivative

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

pause()
 
%Integrate model
disp('Solving ODE')
% ode45
[t,u]=ode45(@(t,y) kuramoto_2([y(2),y(1)],a,w,k,n,G),[0,500],u_int)


function output = kuramoto_2(u,a,w,k,n,G)

u_mat = repmat(u(1),1,n);
output = [u(2); -a*u(1) + w + (k/n)*sum(G.*sin(u_mat' - u_mat),2)];

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