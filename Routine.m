clc;
clear;
close all;
% System Parameters as vectors
LAMBDA=[0.4]';                       
MUH=[0.7]';                        
COST=[0.2]';             
DISCOUNT=[0.17]';       
muL=0;

% Save Parameter vectors
% csvwrite('COST.dat',COST)
% csvwrite('DISCOUNT.dat',DISCOUNT)
% csvwrite('LAMBDA.dat',LAMBDA)
% csvwrite('MUH.dat',MUH)

% Grid for approximation
IOxstar=1;                  % should xstar be on the grid?
IOsgrid=0;                  % do you want to have small grid points around xstar? => Y=1
IOallback=0;
IOcond=1;                   % should the condition numbers be evaluated?
lgrid = 0.01;               % length of step for xstar

digitsgrid=numel(num2str(lgrid));  % precision indicator for xstar-grid

% xstar and startvalue for finding xstar which holds for lambda*D(xstar)=c
xstar=0;                    % this is also the startvalue for the approximation routine

% Chebyshev approximation parameters
n=30;                       % degree of Chebyshev polynomials
a=0;                        % beginning of main interval for Cheby nodes
b=1;                        % ending of main interval for Cheby nodes

% find indifference point x=xstar with dynamic xstar
gridpoints=1000;            % number of grid-point in z
precision=1/gridpoints;     % precision indicator for evalutaion of the value functions
epsilon=0.0001;             % error tolerance for Bisection Method
descend=1/2;                % interval updating of the Bisection Method

XSTAR=[];
SOLUTION=[];
t1=1;                       % Step counter for saving candidates over the whole routine
time=1;                     % Step counter for finding candidates -> only for inspection
time3=1;                    % Step counter for saving solutions over the whole routine

tic
for j4=1:length(MUH)
    muH=MUH(j4,1);
    for j3=1:length(LAMBDA)
        lambda=LAMBDA(j3,1);
        for j2=1:length(DISCOUNT)
            r=DISCOUNT(j2,1);
            for j1=1:length(COST)
                c=COST(j1,1);
                xstar=0;
                roh=r+muH+lambda;
                thresh=c/lambda;
                xg=lambda/muH;
                XSEMI=[];
                tsemi=1;
                while xstar<=1
                %% Evaluating HJB-Equations for VH, VL, and D
                    [X,wherexstar]=gridfunshort( xstar , gridpoints , IOsgrid , IOxstar);
                    [Gws,ATildews,Gfs,ATildefs,Gfw,ATildefw]=GPGNshort(@driftfunc,xstar,xg,X,muL,muH,lambda,wherexstar,precision);
                    if xstar==0
                        [AL,B,l,condL]=VLPGNshort(xstar,r,Gfs,X,IOsgrid,IOallback,IOcond,precision);
                        [AH,C,h,condH]=VHPGNshort(xstar,roh,muH,lambda,c,r,Gfs,ATildefs,X,l,IOsgrid,IOallback,IOcond,precision);
                        [AD,D,d,condD]=DPGNshort(xstar,r,muH,lambda,l,h,Gws,Gfs,Gfw,X,IOsgrid,IOallback,IOcond,precision,ATildefs,c);
                    elseif xstar==1
                        [AL,B,l,condL]=VLPGNshort(xstar,r,Gfw,X,IOsgrid,IOallback,IOcond,precision);
                        [AH,C,h,condH]=VHPGNshort(xstar,roh,muH,lambda,c,r,Gfw,ATildefw,X,l,IOsgrid,IOallback,IOcond,precision);
                        [AD,D,d,condD]=DPGNshort(xstar,r,muH,lambda,l,h,Gws,Gfs,Gfw,X,IOsgrid,IOallback,IOcond,precision,ATildefw,c);
                    else
                        [AL,B,l,condL]=VLPGNshort(xstar,r,Gws,X,IOsgrid,IOallback,IOcond,precision);
                        [AH,C,h,condH]=VHPGNshort(xstar,roh,muH,lambda,c,r,Gws,ATildews,X,l,IOsgrid,IOallback,IOcond,precision);
                        [AD,D,d,condD]=DPGNshort(xstar,r,muH,lambda,l,h,Gws,Gfs,Gfw,X,IOsgrid,IOallback,IOcond,precision,ATildews,c);
                    end
                
                    % find candidates
                    [Candidates,IOC]=candidatesalter(@spline,X,d,thresh);

                    % check monotonicity
                    [monoto,monotol,monotoh]=monotonsimple(X,d,l,h);

                    % find indifferent points for specific xstar
                    if IOC~=0
                        totinnersteps=length(Candidates(:,2));
                        for step=1:totinnersteps
                            sol=Candidates(step,2);
                            XSTAR(t1,1)=xstar;
                            XSEMI(tsemi,1)=xstar;
                            XSTAR(t1,2)=Candidates(step,2);
                            XSEMI(tsemi,2)=Candidates(step,2);
                            XSTAR(t1,3)=spline(X,d,Candidates(step,2));
                            XSTAR(t1,4)=monoto;
                            XSTAR(t1,5)=monotol;
                            XSTAR(t1,6)=monotoh;
                            XSTAR(t1,7)=lambda;
                            XSTAR(t1,8)=muH;
                            XSTAR(t1,9)=c;
                            XSTAR(t1,10)=r;
                            XSTAR(t1,16)=1;
                            XSTAR(t1,18)=condD;
                            XSTAR(t1,19)=condL;
                            XSTAR(t1,20)=condH;
                            t1=t1+1;
                            tsemi=tsemi+1;
                        end
                    end
                    xstar=round(xstar+lgrid,digitsgrid);
                    time=time+1;
                    clear candidates sol totinnersteps monoto change Q F f weights xcandi approxcheb approxspline approxerr
                    clear l h d faked condD Gws ATildews Gfs ATildefs Gfw ATildefw X wherexstar int1 int2
                end

                %% Find optimum
                if IOC==1
                XSEMI=sortrows(XSEMI);
                y45=XSEMI(:,1);
                ycheb=XSEMI(:,2);
                ydiff=ycheb-y45;
                candiopt=[];
                ta2=1;
                for ta1=1:length(ydiff)-1
                    a1=ydiff(ta1,1);
                    a2=ydiff(ta1+1);
                    ba1=a1*a2;
                    if ba1<0
                        candiopt(ta2,1)=ta1;
                        candiopt(ta2,2)=XSEMI(ta1,1);
                        candiopt(ta2,3)=XSEMI(ta1,2);
                        candiopt(ta2,4)=1;
                        ta2=ta2+1;
                    end
                end
                clear ta1 ta2
                emptycandiopt=isempty(candiopt);
                if emptycandiopt==0
                    lengthcandiopt=length(candiopt(:,1));
                    if lengthcandiopt>1
                        warning('More than one optimal point!')
                    end
                else
                    for ta3=1:length(ydiff)
                        a3=ydiff(ta3,1);
                        if a3>=0
                            altercandi=1;
                        else
                            altercandi=0;
                        end
                    end
                    if altercandi==1
                        ydiffabs=abs(ydiff);
                        [H,I]=min(ydiffabs);
                        candiopt(1,1)=I;
                        candiopt(1,2)=XSEMI(I,1);
                        candiopt(1,3)=XSEMI(I,2);
                        candiopt(1,4)=1;
                    end
                    if altercandi==0
                        ydiffabs=abs(ydiff);
                        [H,I]=min(ydiffabs);
                        if H<=lgrid && XSEMI(I,1)~=1
                            candiopt(1,1)=I;
                            candiopt(1,2)=XSEMI(I,1);
                            candiopt(1,3)=XSEMI(I,2);
                            candiopt(1,4)=0;
                        end
                    end
                end
                clear y45 ycheb ydiff ta1 ta2 a1 a2 xstar
                emptycandiopt=isempty(candiopt);
                if emptycandiopt==0
                    %% new setup
                    intervalsize=lgrid;
                    IOsgrid=0;
                    SOL=zeros(1,2);
                    time2=1;                % Step counter for saving Solution for one specific system setting
                    for i=1:length(candiopt(:,1))
                        xstar=candiopt(i,2);
                        origin=candiopt(i,2)-candiopt(i,3);
                        diffsoluabs=1;
                        while diffsoluabs>epsilon && time2<=40 && xstar<=1 && xstar>=0
                            [X,wherexstar]=gridfunshort(xstar,gridpoints,IOsgrid,IOxstar);
                            [Gws,ATildews,Gfs,ATildefs,Gfw,ATildefw]=GPGNshort(@driftfunc,xstar,xg,X,muL,muH,lambda,wherexstar,precision);
                            if xstar==0
                                [AL,B,l,condL]=VLPGNshort(xstar,r,Gfs,X,IOsgrid,IOallback,IOcond,precision);
                                [AH,C,h,condH]=VHPGNshort(xstar,roh,muH,lambda,c,r,Gfs,ATildefs,X,l,IOsgrid,IOallback,IOcond,precision);
                                [AD,D,d,condD]=DPGNshort(xstar,r,muH,lambda,l,h,Gws,Gfs,Gfw,X,IOsgrid,IOallback,IOcond,precision,ATildefs,c);
                            elseif xstar==1
                                [AL,B,l,condL]=VLPGNshort(xstar,r,Gfw,X,IOsgrid,IOallback,IOcond,precision);
                                [AH,C,h,condH]=VHPGNshort(xstar,roh,muH,lambda,c,r,Gfw,ATildefw,X,l,IOsgrid,IOallback,IOcond,precision);
                                [AD,D,d,condD]=DPGNshort(xstar,r,muH,lambda,l,h,Gws,Gfs,Gfw,X,IOsgrid,IOallback,IOcond,precision,ATildefw,c);
                            else
                                [AL,B,l,condL]=VLPGNshort(xstar,r,Gws,X,IOsgrid,IOallback,IOcond,precision);
                                [AH,C,h,condH]=VHPGNshort(xstar,roh,muH,lambda,c,r,Gws,ATildews,X,l,IOsgrid,IOallback,IOcond,precision);
                                [AD,D,d,condD]=DPGNshort(xstar,r,muH,lambda,l,h,Gws,Gfs,Gfw,X,IOsgrid,IOallback,IOcond,precision,ATildews,c);
                            end
                            [~,weights]=chebweights(@spline,n,X,d);
                            [Q,F,f]=chebfun(weights,a,b,c,lambda);

                            checkxstar=xstar;
                            if checkxstar<=0-lgrid
                                intv=[X(wherexstar,1),X(wherexstar+1,1),1]';
                            elseif checkxstar>=1-lgrid
                                intv=[0,X(wherexstar-1,1),X(wherexstar,1)]';
                            else
                                intv=[0,X(wherexstar-1,1),X(wherexstar,1),X(wherexstar+1,1),1]';
                            end
                            approxcheb=F(intv);
                            approxspline=spline(X,d,intv);
                            approxerr=abs(approxspline)-abs(approxcheb);
                            % check monotonicity
                            [monoto,change,monotosimple,monotol,monotoh,strictmonoto,strictmonotol,strictmonotoh]=monotonialter(f,X,d,l,h);
                            % Newton's method
                            [sol,~,ierr,~]=nsold(xstar,Q,[10^-12,12^-12]);
                            % Save candidates
                            XSTAR(t1,1)=xstar;
                            XSTAR(t1,2)=sol;
                            XSTAR(t1,3)=F(sol);
                            XSTAR(t1,4)=monotosimple;
                            XSTAR(t1,5)=monotol;
                            XSTAR(t1,6)=monotoh;
                            XSTAR(t1,7)=lambda;
                            XSTAR(t1,8)=muH;
                            XSTAR(t1,9)=c;
                            XSTAR(t1,10)=r;
                            XSTAR(t1,16)=2;
                            XSTAR(t1,17)=monoto;
                            XSTAR(t1,18)=condD;
                            XSTAR(t1,19)=condL;
                            XSTAR(t1,20)=condH;
                            XSTAR(t1,21)=ierr;
                            XSTAR(t1,22)=strictmonoto;
                            XSTAR(t1,23)=strictmonotol;
                            XSTAR(t1,24)=strictmonotoh;
                            if checkxstar>0+lgrid && checkxstar<1-lgrid
                                XSTAR(t1,11)=approxerr(1,1);
                                XSTAR(t1,12)=approxerr(2,1);
                                XSTAR(t1,13)=approxerr(3,1);
                                XSTAR(t1,14)=approxerr(4,1);
                                XSTAR(t1,15)=approxerr(5,1);
                            else
                                if checkxstar>=0+lgrid
                                    XSTAR(t1,13)=approxerr(1,1);
                                    XSTAR(t1,14)=approxerr(2,1);
                                    XSTAR(t1,15)=approxerr(3,1);
                                elseif checkxstar<=1-lgrid
                                    XSTAR(t1,11)=approxerr(1,1);
                                    XSTAR(t1,12)=approxerr(2,1);
                                    XSTAR(t1,13)=approxerr(3,1);
                                end
                            end
                            t1=t1+1;
                            % Bisection method
                            if candiopt(i,4)==1
                                if origin<0
                                    diffsolu=xstar-sol;
                                    if diffsolu<0
                                        a1=xstar;
                                        a2=xstar+intervalsize;
                                        xstar=abs(a1+a2)*descend;
                                        intervalsize=abs(a1-xstar);
                                    else
                                        a1=xstar-intervalsize;
                                        a2=xstar;
                                        xstar=abs(a1+a2)*descend;
                                        intervalsize=abs(a1-xstar);
                                    end
                                else
                                    diffsolu=xstar-sol;
                                    if diffsolu<0
                                        a1=xstar-intervalsize;
                                        a2=xstar;
                                        xstar=abs(a1+a2)*descend;
                                        intervalsize=abs(a1-xstar);
                                    else
                                        a1=xstar;
                                        a2=xstar+intervalsize;
                                        xstar=abs(a1+a2)*descend;
                                        intervalsize=abs(a1-xstar);
                                    end
                                end
                            else
                                if origin<0
                                    diffsolu=xstar-sol;
                                    if diffsolu<0
                                        a1=xstar;
                                        a2=xstar-intervalsize;
                                        xstar=abs(a1+a2)*descend;
                                        intervalsize=abs(a1-xstar);
                                    else
                                        a1=xstar+intervalsize;
                                        a2=xstar;
                                        xstar=abs(a1+a2)*descend;
                                        intervalsize=abs(a1-xstar);
                                    end
                                else
                                    diffsolu=xstar-sol;
                                    if diffsolu<0
                                        a1=xstar+intervalsize;
                                        a2=xstar;
                                        xstar=abs(a1+a2)*descend;
                                        intervalsize=abs(a1-xstar);
                                    else
                                        a1=xstar;
                                        a2=xstar-intervalsize;
                                        xstar=abs(a1+a2)*descend;
                                        intervalsize=abs(a1-xstar);
                                    end
                                end
                            end
                            diffsoluabs=abs(diffsolu);
                            % save solution
                            SOL(time2,1)=xstar;
                            SOL(time2,2)=diffsolu;
                            SOL(time2,3)=diffsoluabs;
                            SOL(time2,4)=lambda;
                            SOL(time2,5)=muH;
                            SOL(time2,6)=c;
                            SOL(time2,7)=r;
                            SOL(time2,8)=monoto;
                            SOL(time2,9)=condD;
                            SOL(time2,10)=condL;
                            SOL(time2,11)=condH;
                            SOL(time2,17)=ierr;
                            SOL(time2,18)=monotosimple;
                            SOL(time2,19)=monotol;
                            SOL(time2,20)=monotoh;
                            if checkxstar>0+lgrid && checkxstar<1-lgrid
                                SOL(time2,12)=approxerr(1,1);
                                SOL(time2,13)=approxerr(2,1);
                                SOL(time2,14)=approxerr(3,1);
                                SOL(time2,15)=approxerr(4,1);
                                SOL(time2,16)=approxerr(5,1);
                            else
                                if checkxstar>=0+lgrid
                                    SOL(time2,14)=approxerr(1,1);
                                    SOL(time2,15)=approxerr(2,1);
                                    SOL(time2,16)=approxerr(3,1);
                                elseif checkxstar<=1-lgrid
                                    SOL(time2,12)=approxerr(1,1);
                                    SOL(time2,13)=approxerr(2,1);
                                    SOL(time2,14)=approxerr(3,1);
                                end
                            end
                            time2=time2+1;
                            clear candidates sol totinnersteps monoto change weights xcandi approxcheb approxspline
                            clear l h d faked condD Gws ATildews Gfs ATildefs Gfw ATildefw X wherexstar int1 int2 intv
                            clear approxerr checkxstar
                        end
                    end
                    [~,Idiff]=min(SOL(:,3));
                    % Save Solution for each system-setting
                    SOLUTION(time3,1)=SOL(Idiff,1);
                    SOLUTION(time3,2)=F(SOL(Idiff,1));
                    SOLUTION(time3,3)=SOL(Idiff,4);
                    SOLUTION(time3,4)=SOL(Idiff,5);
                    SOLUTION(time3,5)=SOL(Idiff,6);
                    SOLUTION(time3,6)=SOL(Idiff,7);
                    SOLUTION(time3,7)=SOL(Idiff,8);
                    SOLUTION(time3,8)=SOL(Idiff,9);
                    SOLUTION(time3,9)=SOL(Idiff,10);
                    SOLUTION(time3,10)=SOL(Idiff,11);
                    SOLUTION(time3,16)=SOL(Idiff,17);
                    SOLUTION(time3,17)=SOL(Idiff,18);
                    SOLUTION(time3,18)=SOL(Idiff,19);
                    SOLUTION(time3,19)=SOL(Idiff,20);
                    SOLUTION(time3,20)=SOL(Idiff,2);
                    if SOL(Idiff,1)>0+lgrid && SOL(Idiff,1)<1-lgrid
                        SOLUTION(time3,11)=SOL(Idiff,12);
                        SOLUTION(time3,12)=SOL(Idiff,13);
                        SOLUTION(time3,13)=SOL(Idiff,14);
                        SOLUTION(time3,14)=SOL(Idiff,15);
                        SOLUTION(time3,15)=SOL(Idiff,16);
                    else
                        if SOL(Idiff,1)>=0+lgrid
                            SOLUTION(time3,13)=SOL(Idiff,14);
                            SOLUTION(time3,14)=SOL(Idiff,15);
                            SOLUTION(time3,15)=SOL(Idiff,16);
                        elseif SOL(Idiff,1)<=1-lgrid
                            SOLUTION(time3,11)=SOL(Idiff,12);
                            SOLUTION(time3,12)=SOL(Idiff,13);
                            SOLUTION(time3,13)=SOL(Idiff,14);
                        end
                    end

                    time3=time3+1;
                    XSTAR=sortrows(XSTAR);
                    clear IOC intv intervalsize Idiff a1 a2 condH condL AD AH AL B ba1
                    clear C Candidates candipot cu2 D diffsolu diffsoluabs emptycandiopt ia2 ic2 lengthcandiopt
                    clear monotoh monotol monotosimple origin SOL time2 Q F f
                end
                end
                clear XSEMI tsemi
            end
        end
    end
end
ending=toc;
% save output 
% csvwrite('XSTAR.csv',XSTAR)
% csvwrite('SOLUTION.csv',SOLUTION)