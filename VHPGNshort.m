function [ AH , C , h , condH ] = VHPGNshort(  xstar, roh , muH , lambda , c , r , G , ATilde , X , l , IOsgrid , IOallback , IOcond , precision)
%Evaluating the linear system of equations for H-type
%   Output:     AH, coefficient matrix
%               h,  Value-function results
%               condH condition number of the coefficinet matrix

if length(G)~=length(X)
    error('vector dimension of X and G must be the same!')
elseif length(l)~=length(X)
    error('vector dimension of X and l must be the same!')
elseif length(ATilde)~=length(X)
    error('vector dimension of X and h must be the same!')
end
if IOallback==1
    
elseif IOallback==0
    
else
    error('IOallback must be either 0 or 1');
end
griddigits1=numel(num2str(xstar));
griddigits2=numel(num2str(precision));
sgrid2digits=max(griddigits1,griddigits2);

   AH=zeros(length(X),length(X));
   C=zeros(length(X),1);
   t1=1;
   if IOsgrid==1
        error('does not work with IOsgrid=1');
   else
       while t1<=length(X)
           Ati=ATilde(t1,1);
           if t1==1
               splus=round(X(t1+1,1)-X(t1,1),sgrid2digits);
               C(t1,1)=(X(t1,1)+(lambda*l(t1,1)*(1-Ati))-(Ati*c))*splus;
               AH(t1,t1)=((roh-(Ati*lambda))*splus)+G(t1,1);
               AH(t1,t1+1)=(-1)*G(t1,1);
               AH(t1,length(X))=(-1)*muH*splus;
           elseif t1==length(X)-1
               symmet=round(X(t1+1,1)-X(t1-1,1),sgrid2digits);
               C(t1,1)=(X(t1,1)+(lambda*l(t1,1)*(1-Ati))-(Ati*c))*symmet;
               AH(t1,t1-1)=G(t1,1);
               AH(t1,t1)=(roh-(Ati*lambda))*symmet;
               AH(t1,t1+1)=(-1)*((muH*symmet)+G(t1,1));
           elseif t1==length(X)
               sminus=round(X(t1,1)-X(t1-1,1),sgrid2digits);
               C(t1,1)=(X(t1,1)+(lambda*l(t1,1)*(1-Ati))-(Ati*c))*sminus;
               AH(t1,t1-1)=G(t1,1);
               AH(t1,t1)=((r+((1-Ati)*lambda))*sminus)-G(t1,1);
           else
               a1=round(X(t1,1)-X(t1-1,1),sgrid2digits);
               a2=round(X(t1+1,1)-X(t1,1),sgrid2digits);
               if a1==a2
                   symmet=round(X(t1+1,1)-X(t1-1,1),sgrid2digits);
                   C(t1,1)=(X(t1,1)+(lambda*l(t1,1)*(1-Ati))-(Ati*c))*symmet;
                   AH(t1,t1-1)=G(t1,1);
                   AH(t1,t1)=(roh-(Ati*lambda))*symmet;
                   AH(t1,t1+1)=(-1)*G(t1,1);
                   AH(t1,length(X))=(-1)*(muH*symmet);
               else
                   if xstar==0
                       sminus=round(X(t1,1)-X(t1-1,1),sgrid2digits);
                       C(t1,1)=(X(t1,1)+(lambda*l(t1,1)*(1-Ati))-(Ati*c))*sminus;
                       AH(t1,t1-1)=G(t1,1);
                       AH(t1,t1)=((roh-(Ati*lambda))*sminus)-G(t1,1);
                       AH(t1,length(X))=(-1)*(muH*sminus);
                   elseif xstar==1
                       splus=round(X(t1+1,1)-X(t1,1),sgrid2digits);
                       C(t1,1)=(X(t1,1)+(lambda*l(t1,1)*(1-Ati))-(Ati*c))*splus;
                       AH(t1,t1)=((roh-(Ati*lambda))*splus)+G(t1,1);
                       AH(t1,t1+1)=(-1)*G(t1,1);
                       AH(t1,length(X))=(-1)*muH*splus;
                   else
                       check=G(t1,1);
                       if check<0
                           sminus=round(X(t1,1)-X(t1-1,1),sgrid2digits);
                           C(t1,1)=(X(t1,1)+(lambda*l(t1,1)*(1-Ati))-(Ati*c))*sminus;
                           AH(t1,t1-1)=G(t1,1);
                           AH(t1,t1)=((roh-(Ati*lambda))*sminus)-G(t1,1);
                           AH(t1,length(X))=(-1)*(muH*sminus);
                       else
                           splus=round(X(t1+1,1)-X(t1,1),sgrid2digits);
                           C(t1,1)=(X(t1,1)+(lambda*l(t1,1)*(1-Ati))-(Ati*c))*splus;
                           AH(t1,t1)=((roh-(Ati*lambda))*splus)+G(t1,1);
                           AH(t1,t1+1)=(-1)*G(t1,1);
                           AH(t1,length(X))=(-1)*muH*splus;
                       end
                   end
               end
           end
           t1=t1+1;
       end 
   end
   
if IOcond==1
    condH=cond(AH);
else
    condH=0;
end

h=AH\C;
   
end

