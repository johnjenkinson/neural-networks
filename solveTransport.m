% call: solveTransport.m
% John Jenkinson, UTSA ECE 2014/10/13
%
% Numerical solution to the transport equation
% u_x + u_t = 0 on the region 
%{x\in[0,1],t\in[0,10]} with boundary conditions 
% u(x,0)=sin(2*pi*x) and u(0,t)=sin(10*pi*t), 
% where u is the solution. The time derivative
% is approximated by the forward difference
% quotient, [u_i(n+1)-u_i(n)]/n*dt and the
% space derivative is approximated by the
% backward difference quotient, 
% [u_i(n)-u_{i-1}(n)]/i*dx.
%
% Example:
% gridX=[0,1];  Solution interval in space
% gridT=[0,10]; Solution interval in time
% x=gridX(1):dx:gridX(2); 
% t=gridT(1):dt:gridT(2);
% dt=0.01;     time step
% dx=dt;       space step
% u0x=sin(2*pi*x);  boundary condition in space
% u0t=sin(10*pi*t); boundary condition in time
%
function[u]=solveTransport(gridX,gridT,...
    dx,dt,u0x,u0t,x,t)

% Create grid
X=length(x);
T=length(t);

% Initialize solution
u=zeros(T,X);
u(1,:)=u0x;
u(:,1)=u0t;

% Numerical method
k=dt/dx;
a=(1-k);
for n=1:(T-1)
    for i=2:X
        u(n+1,i)=k*u(n,i-1)+a*u(n,i);
    end
end
end
