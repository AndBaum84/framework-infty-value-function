function [ g ] = driftfunc( x , ATilde , muL , muH , lambda )
% Reputatinal Drift Function (in the absence of a signal)
% is an autonomous ordinary differential equation (ODE): dx/dt=g(x)

g = lambda*(ATilde-x)-(muH-muL)*x.*(1-x);


end

