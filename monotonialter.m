function [ monoto , change , monotosimple , monotol , monotoh , strictmonoto , strictmonotol , strictmonotoh] = monotonialter( f , X , d , l , h)
% Check (strict) monotonicity during Bisection method of L-type, H-type, and value
% of quality


max=length(X)-1;
moni=zeros(max,1);
t1=1;
change = [];
for i=1:length(X)-1
    a1=X(i,1);
    a2=X(i+1,1);
    b1=f(a1);
    b2=f(a2);
    c1=b1*b2;
    if c1>=0
        moni(i,1)=1;
    else
        moni(i,1)=0;
        change(t1,1)=i;
        change(t1,2)=X(i,1);
        t1=t1+1;
    end
end

check=sum(moni);

if check==max
    monoto=1;
else
    monoto=0;
%     warning('Dividend function is not monoton on [0,1] and may have more than one solution to lambda*D(xstar)=c')
end

maxalter=length(X)-2;
honi=zeros(maxalter,1);
honil=zeros(maxalter,1);
honih=zeros(maxalter,1);
stricti=zeros(maxalter,1);
strictil=zeros(maxalter,1);
strictih=zeros(maxalter,1);
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
    if c2>0
        stricti(i1,1)=1;
    else
        stricti(i1,1)=0;
    end
    if cl>0
        strictil(i1,1)=1;
    else
        strictil(i1,1)=0;
    end
    if ch>0
        strictih(i1,1)=1;
    else
        strictih(i1,1)=0;
    end
end

checkalter=sum(honi(:,1));
checkalterl=sum(honil(:,1));
checkalterh=sum(honih(:,1));
checkstricti=sum(stricti(:,1));
checkstrictil=sum(strictil(:,1));
checkstrictih=sum(strictih(:,1));

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
if checkstricti==maxalter
    strictmonoto=1;
else
    strictmonoto=0;
end
if checkstrictil==maxalter
    strictmonotol=1;
else
    strictmonotol=0;
end
if checkstrictih==maxalter
    strictmonotoh=1;
else
    strictmonotoh=0;
end

end

