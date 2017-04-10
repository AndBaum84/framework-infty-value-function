function [ Q , G , g] = chebfun( weights , a , b , c , lambda )
%creating the chby polynom in concrete form
%   Inputs:
%   weights:    the cheby weights calculated by the fun chebweights
%   a & b:      start- end end-points of the true grid [a,b]=[0,1]
%   Output:
%   G:          cheby function ready to evaluate on the true grid [0,1]
%   g:          derivative of G
%   C:          cheby function - c/lambda

% p=zeros(length(weights));
for i=1:length(weights)
    p(i)= weights(i)*(2*((sym('x')-a)/(b-a))-1)^(i-1);
end
p=sum(p);

syms f(x);
f(x)=p;
g = diff(f(x));
G = matlabFunction(f);
g = matlabFunction(g);

weightsnew=weights;
weightsnew(1,1)=weightsnew(1,1)-(c/lambda);
% q=zeros(length(weights));
for t=1:length(weightsnew)
    q(t)= weightsnew(t)*(2*((sym('x')-a)/(b-a))-1)^(t-1);
end
l=sum(q);

syms C(x);
C(x)=l;

Q = matlabFunction(C);
end

