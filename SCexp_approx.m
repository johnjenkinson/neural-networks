% call: SCexp_approx.m
% John Jenkinson UTSA ECE Nov. 2014.
%
% Forward Euler solution to Eccles' space-clamped eq.,
% u_t + u = delta(t-T1), where the delta function
% is approximated by a sequence of Gaussian functions
% of decreasing variance.
%
% Inputs are enumerated below:
% dt - time step; 
% x  - vector of grid points;
% ic - initial condition;
% b,c,d - parameters for the Gaussian function,
% a*exp( -(x-b)^2 / 2*c^2 ) + d, where
% b = expected value;
% c = standard deviation
% d - value that the function asymptotically approaches
% Note: in practice d=0.
%
% Output of the function is a matrix ulimit where
% each row is the solution for a value of c.

function[ulimit]=SCexp_approx(dt,x,ic,b,c,d)

    u=zeros(1,numel(x));
    ulimit=zeros(numel(c),numel(x));
    u(1)=ic;
    uc=(1-dt);
    xend=end(x)-dt;
           
    for v=1:numel(c);
        a=1/( c(v)*sqrt(2*pi) );
        for k=1:dt:xend
            if( k>=(b+5*c) )
                e=0;
            else
            e=a*exp( -((k-b)^2 / (2*c(v)^2)) )+d;
            u(k+1)=uc*u(k)+dt*e;
        end
        ulimit(v,:)=u;
    end
    
end
