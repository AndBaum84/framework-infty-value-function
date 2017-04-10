function [ AD , D , d , condD ] = DPGNshort( xstar , r , muH , lambda , l , h , Gws , Gfs , Gfw , X , IOsgrid , IOallback , IOcond , precision,ATilde,c)
%Evaluating the linear system of equations for L-type
%   Output:     AD, coefficient matrix
%               d,  Value-function results
%               condD condition number of the coefficinet matrix

if xstar==0
    G=Gfs(:,1);
elseif xstar==1
    G=Gfw(:,1);
else
    G=Gws(:,1);
end

if IOallback==1
    
elseif IOallback==0
    
else
    error('IOallback must be either 0 or 1');
end

if length(G)~=length(X)
    error('vector dimension of X and G must be the same!')
elseif length(l)~=length(X)
    error('vector dimension of X and l must be the same!')
elseif length(h)~=length(X)
    error('vector dimension of X and h must be the same!')
end

griddigits1=numel(num2str(xstar));
griddigits2=numel(num2str(precision));
sgrid2digits=max(griddigits1,griddigits2);


AD=zeros(length(X),length(X));
D=zeros(length(X),1);
t1=1;
if IOsgrid==1
    error('does not work with IOsgrid=1');
else
    while t1<=length(X)
        ATi=ATilde(t1,1);
       if t1==1
           splus=round(X(t1+1,1)-X(t1,1),sgrid2digits);
           D(t1,1)=(muH*(h(length(X),1)-h(t1,1))-ATi*c)*splus;
           AD(t1,t1)=((r+lambda-(ATi*lambda))*splus)+G(t1,1);
           AD(t1,t1+1)=(-1)*G(t1,1);
       elseif t1==length(X)
           sminus=round(X(t1,1)-X(t1-1,1),sgrid2digits);
           D(t1,1)=(muH*(h(length(X),1)-h(t1,1))-ATi*c)*sminus;
           AD(t1,t1-1)=G(t1,1);
           AD(t1,t1)=((r+lambda-(ATi*lambda))*sminus)-G(t1,1);
       else
           a1=round(X(t1,1)-X(t1-1,1),sgrid2digits);
           a2=round(X(t1+1,1)-X(t1,1),sgrid2digits);
           if a1==a2
               symmet=round(X(t1+1,1)-X(t1-1,1),sgrid2digits);
               D(t1,1)=(muH*(h(length(X),1)-h(t1,1))-ATi*c)*symmet;
               AD(t1,t1-1)=G(t1,1);
               AD(t1,t1)=(r+lambda-(ATi*lambda))*symmet;
               AD(t1,t1+1)=(-1)*G(t1,1);
           else
               if xstar==0
                   sminus=round(X(t1,1)-X(t1-1,1),sgrid2digits);
                   D(t1,1)=(muH*(h(length(X),1)-h(t1,1))-ATi*c)*sminus;
                   AD(t1,t1-1)=G(t1,1);
                   AD(t1,t1)=((r+lambda-(ATi*lambda))*sminus)-G(t1,1);
               elseif xstar==1
                   splus=round(X(t1+1,1)-X(t1,1),sgrid2digits);
                   D(t1,1)=(muH*(h(length(X),1)-h(t1,1))-ATi*c)*splus;
                   AD(t1,t1)=((r+lambda-(ATi*lambda))*splus)+G(t1,1);
                   AD(t1,t1+1)=(-1)*G(t1,1);
               else
                   check=Gws(t1,1);
                   if check<0
                       sminus=round(X(t1,1)-X(t1-1,1),sgrid2digits);
                       D(t1,1)=(muH*(h(length(X),1)-h(t1,1))-ATi*c)*sminus;
                       AD(t1,t1-1)=G(t1,1);
                       AD(t1,t1)=((r+lambda-(ATi*lambda))*sminus)-G(t1,1);
                   else
                       splus=round(X(t1+1,1)-X(t1,1),sgrid2digits);
                       D(t1,1)=(muH*(h(length(X),1)-h(t1,1))-ATi*c)*splus;
                       AD(t1,t1)=((r+lambda-(ATi*lambda))*splus)+G(t1,1);
                       AD(t1,t1+1)=(-1)*G(t1,1);
                   end
               end
           end
       end
       t1=t1+1;
    end
end

if IOcond==1
    condD=cond(D);
else
    condD=0;
end

d=AD\D;


end


