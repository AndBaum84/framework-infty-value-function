function [ AL , B , l , condL ] = VLPGNshort(xstar,r , G , X , IOsgrid , IOallback , IOcond , precision)
%Evaluating the linear system of equations for L-type
%   Output:     AL, coefficient matrix
%               l,  Value-function results
%               condL condition number of the coefficinet matrix

if length(G)~=length(X)
    error('vector dimension of X and G must be the same!')
end
if IOallback==1
    
elseif IOallback==0
    
else
    error('IOallback must be either 0 or 1');
end

griddigits1=numel(num2str(xstar));
griddigits2=numel(num2str(precision));
sgrid2digits=max(griddigits1,griddigits2);

AL=zeros(length(X),length(X));
B=zeros(length(X),1);
t1=1;
       while t1<=length(X)
           if IOsgrid==1
                error('does not work with IOsgrid=1');
           else
              if t1==1
                  splus=round(X(t1+1,1)-X(t1,1),sgrid2digits);
                  B(t1,1)=X(t1,1)*splus;
                  AL(t1,t1)=G(t1,1)+(r*splus);
                  AL(t1,t1+1)=(-1)*G(t1,1);
              elseif t1==length(X)
                  sminus=round(X(t1,1)-X(t1-1,1),sgrid2digits);
                  B(t1,1)=X(t1,1)*sminus;
                  AL(t1,t1-1)=G(t1,1);
                  AL(t1,t1)=(r*sminus)-G(t1,1);
              else
                  a1=round(X(t1,1)-X(t1-1,1),sgrid2digits);
                  a2=round(X(t1+1,1)-X(t1,1),sgrid2digits);
                  if a1==a2
                      symmet=round(X(t1+1,1)-X(t1-1,1),sgrid2digits);
                      B(t1,1)=X(t1,1)*symmet;
                      AL(t1,t1-1)=G(t1,1);
                      AL(t1,t1)=r*symmet;
                      AL(t1,t1+1)=(-1)*G(t1,1);
                  else
                      if xstar==0
                          sminus=round(X(t1,1)-X(t1-1,1),sgrid2digits);
                          B(t1,1)=X(t1,1)*sminus;
                          AL(t1,t1-1)=G(t1,1);
                          AL(t1,t1)=(r*sminus)-G(t1,1);
                      elseif xstar==1
                          splus=round(X(t1+1,1)-X(t1,1),sgrid2digits);
                          B(t1,1)=X(t1,1)*splus;
                          AL(t1,t1)=G(t1,1)+(r*splus);
                          AL(t1,t1+1)=(-1)*G(t1,1);
                      else
                          check=G(t1,1);
                          if check<0
                              sminus=round(X(t1,1)-X(t1-1,1),sgrid2digits);
                              B(t1,1)=X(t1,1)*sminus;
                              AL(t1,t1-1)=G(t1,1);
                              AL(t1,t1)=(r*sminus)-G(t1,1);
                          else
                              splus=round(X(t1+1,1)-X(t1,1),sgrid2digits);
                              B(t1,1)=X(t1,1)*splus;
                              AL(t1,t1)=G(t1,1)+(r*splus);
                              AL(t1,t1+1)=(-1)*G(t1,1);
                          end
                      end
                  end
              end 
           end
           t1=t1+1;
       end

    if length(B)~=length(X)
        error('vector dimensions of B and X must be the same!')
    end
    
if IOcond==1    
    condL=cond(AL);
else
    condL=0;
end

l=AL\B;
       
end

