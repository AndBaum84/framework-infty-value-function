function [ Gws , ATildews , Gfs , ATildefs , Gfw , ATildefw ] = GPGNshort( driftfunc , xstar , xg , X , muL , muH , lambda , wherexstar , precision)
%
%   Gws => Output of driftfunction in a Work-Shirk setup given xstar
%   Gfs => Output of driftfunction in a Full-Shirk setup given xstar=0
%   Gfw => Output of driftfunction in a Full-Work setup given xstar=1
%   ATilde => Vector of believed/actual investments given xstar

griddigits1=numel(num2str(xstar));
griddigits2=numel(num2str(precision));
sgrid2digits=max(griddigits1,griddigits2);

Gws=zeros(length(X),1);
ATildews=zeros(length(X),1);
Gfs=zeros(length(X),1);
ATildefs=zeros(length(X),1);
Gfw=zeros(length(X),1);
ATildefw=zeros(length(X),1);
       
    if xstar==0
       t1=1;
       while t1<=length(X)
           aws01=round(X(t1,1),sgrid2digits);
           if t1==1
               Gws(t1,1)=driftfunc(aws01 , 0 , muL , muH , lambda);
               ATildews(t1,1)=0;
           else
               Gws(t1,1)=driftfunc(aws01, 1 , muL , muH , lambda);
               ATildews(t1,1)=1;
           end
           t1=t1+1;
       end
       t2=1;
       while t2<=length(X)
           afs1=round(X(t2,1),sgrid2digits);
           Gfs(t2,1)=driftfunc( afs1 , 0 , muL , muH , lambda );
           ATildefs(t2,1)=0;
           t2=t2+1; 
       end
    elseif xstar==1
       t3=1;
       while t3<=length(X)
            aws11=round(X(t3,1),sgrid2digits);
            if t3==length(X)
                Gws(t3,1)=driftfunc(aws11 , 1 , muL , muH , lambda);
                ATildews(t3,1)=1;
            else
                Gws(t3,1)=driftfunc(aws11 , 0 , muL , muH , lambda);
                ATildews(t3,1)=0;
            end
            t3=t3+1;
       end 
       t4=1;
       while t4<=length(X)
            afw1=round(X(t4,1),sgrid2digits);
            Gfw(t4,1)=driftfunc( afw1 , 1 , muL , muH , lambda);
            ATildefw(t4,1)=1;
            t4=t4+1;
       end
    else
       if xstar<xg
           t5=1;
           while t5<=length(X)
                aws21=round(X(t5,1),sgrid2digits);
                if t5<wherexstar
                    Gws(t5,1)=driftfunc(aws21 , 1 , muL , muH , lambda);
                    ATildews(t5,1)=1;
                elseif t5==wherexstar
                    aws22=aws21*(1+((muH/lambda)*(1-aws21)));
%                     Gws(t5,1)=driftfunc(aws21 , aws22 , muL , muH , lambda);
                    Gws(t5,1)=0;
                    ATildews(t5,1)=aws22;
                else
                    Gws(t5,1)=driftfunc(aws21 , 0 , muL , muH , lambda);
                    ATildews(t5,1)=0;
                end
                t5=t5+1;
           end
       else
           t6=1;
           while t6<=length(X)
                aws31=round(X(t6,1),sgrid2digits);
                if t6<=wherexstar
                    Gws(t6,1)=driftfunc( aws31 , 1 , muL , muH , lambda);
                    ATildews(t6,1)=1;
                else
                    Gws(t6,1)=driftfunc( aws31 , 0 , muL , muH , lambda);
                    ATildews(t6,1)=0;
                end
                t6=t6+1;
           end
       end
    end

end

