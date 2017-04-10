function [c , PolyWeights] = chebweights( spline, n , X , d)
%Calculating Chebyweights with Cubic spline function to find Cheby nodes
%   Inputs:
%   X:          vector of single value on the true grid [0,1] for the 
%               approximation with cubic spline.
%   n:          Degree of Polynomial Approximation
%   d:          Vector of Dividends at each X(i)
%   Outputs:
%   c:          Coefficients of the cheby
%   Polyweights:weights of the cheby polynom

chebyX = -cos(((2*(0:n)+1)*pi)/(2*n+2));
xx = (chebyX + 1)*0.5;       % Transforming the Cheby Nodes to the intervall [0 1]
yy = spline(X,d,xx);         % Evaluating the polynomials at the cubic splines

% Defining the cheby coefficients c
c=zeros(n+1,1);
c(1)=sum(yy)/(n+1);
for k=1:n
c(k+1)=2*sum(yy.*cos(k*acos(chebyX)))/(n+1);
end

% Defining the Cheby Polynomial Matrix T
T=zeros(n+1,n+1);
T(1,1)=1; 
T(2,2)=1;
for i=3:(n+1),
    T(i,2:i)=T(i,2:i)+2*T(i-1,1:(i-1));
    T(i,1:(i-2))=T(i,1:(i-2))-T(i-2,1:(i-2));
end

PolyWeights=sum(repmat(c,1,n+1).*T,1);