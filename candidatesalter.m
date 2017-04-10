function [ candidates , IOC ] = candidatesalter( spline , X , d , thresh)
% Find candidates with intersection analysis with cubic splines
%   Output:
%           candidates, can be a vector all candidates near intersection
%           IOC, are there any candidates

dnew=zeros(length(d),1);
for t=1:length(d)
    dnew(t,1)=d(t,1)-thresh;
end

t1=1;
candidates=[];
for i=1:length(X)-1
    a1=X(i,1);
    a2=X(i+1,1);
    y1=spline(X,dnew,a1);
    y2=spline(X,dnew,a2);
    c=y1*y2;
    if c<0
        candidates(t1,1)=i;
        candidates(t1,2)=X(i,1);
        t1=t1+1;
    end
end

s=isempty(candidates);
if s==1
    IOC=0;
else
    IOC=1;
    if length(candidates(:,1))>1
        warning('potentially more than one solution')
    end
end

end

