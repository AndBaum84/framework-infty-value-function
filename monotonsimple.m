function [ monotosimple , monotol , monotoh] = monotonsimple( X , d , l , h)
% Check monotonicity of value-curves H-type, L-type, and value of quality
%   self explaining

maxalter=length(X)-2;
honi=zeros(maxalter,1);
honil=zeros(maxalter,1);
honih=zeros(maxalter,1);
for i1=1:length(X)-2
    a3=d(i1,1)-d(i1+1,1);
    a4=d(i1+1,1)-d(i1+2,1);
    al1=l(i1,1)-l(i1+1,1);
    al2=l(i1+1,1)-l(i1+2,1);
    ah1=h(i1,1)-h(i1+1,1);
    ah2=h(i1+1,1)-h(i1+2,1);
    c2=a3*a4;
    cl=al1*al2;
    ch=ah1*ah2;
    if c2>=0
        honi(i1,1)=1;
    else
        honi(i1,1)=0;
    end
    if cl>=0
        honil(i1,1)=1;
    else
        honil(i1,1)=0;
    end
    if ch>=0
        honih(i1,1)=1;
    else
        honih(i1,1)=0;
    end
end

checkalter=sum(honi(:,1));
checkalterl=sum(honil(:,1));
checkalterh=sum(honih(:,1));

if checkalter==maxalter
    monotosimple=1;
else
    monotosimple=0;
end
if checkalterl==maxalter
    monotol=1;
else
    monotol=0;
end
if checkalterh==maxalter
    monotoh=1;
else
    monotoh=0;
end

end

