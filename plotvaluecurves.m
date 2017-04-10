clc;
clear;
close all;
% System Parameters
lambda = 0.4;
muH = 0.7;                                              
muL = 0;
c = 0.2;
r = 0.17;                                                
%% Perfect Good News (PGN) Value function approximation by HJB equations Routine
% In case of PGN, only three possible equilibrium dynamics are possible: 
% 1) work-shirk equilibrium
%   i)  xstar<xg => ATilde(xstar)=xstar*(1+((muH/lambda)*(1-xstar)))
%   ii) xstar>xg => ATilde(xstar)=1
% 2) full shirk equilibrium
%   ATilde(xstar)=0 (all ATilde(x)=0
%
% In equilibrium i) the cutoff is convergent whereas in 2) the cutoff is
% permeable.
% 
% Furthermore, if 1=>muH/lambda => the equilibrium is unique
% hence, for xstar element of (0,1) -> (1/xstar)>1=>muH/lambda
% and therefore, the equilibrium is permeable and unique.
% But for xstar element of (0,1) -> (1/xstar)>muH/lambda>1
% the equilibrium is just permeable but maybe not unique.

roh=r+muH+lambda;
thresh=c/lambda;

% Grid for approximation
IOxstar=1;                  % should xstar be on the grid
IOsgrid=0;                  % do you want to have small grid points around xstar => YN=1
IOallback=0;                % should always be 0
IOcond=1;
lgrid=0.01;                 % length of main grid
sgridn=100;                 % how many small intervalls you want to have around xstar (must be even)
digitsgrid=numel(num2str(lgrid));

gridpoints=1000;
gridintervals=gridpoints;
precision=1/gridpoints;

% xstar and startvalue for finding xstar which holds for lambda*D(xstar)=c
xstar=0.3;                  % this is also the startvalue for the approximation
% Xstar=[0:lgrid:1]';
xg=lambda/muH;

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
    [AD,D,d,condD]=DPGNshort(xstar,r,muH,lambda,l,h,Gws,Gfs,Gfw,X,IOsgrid,IOallback,IOcond,precision,ATildefw);
else
    [AL,B,l,condL]=VLPGNshort(xstar,r,Gws,X,IOsgrid,IOallback,IOcond,precision);
    [AH,C,h,condH]=VHPGNshort(xstar,roh,muH,lambda,c,r,Gws,ATildews,X,l,IOsgrid,IOallback,IOcond,precision);
    [AD,D,d,condD]=DPGNshort(xstar,r,muH,lambda,l,h,Gws,Gfs,Gfw,X,IOsgrid,IOallback,IOcond,precision,ATildews,c);
end

% csvwrite('d1000.dat',d)
% csvwrite('X1000.dat',X)
% csvwrite('cond1000.dat',[condL,condH,condD])

%% Print the Value- and Dividend-Function 
% Calculate VH-VL
faked=h-l;
% Define Thresholdvalue for investment insentives
Thresh=zeros(length(X),1);
for t1=1:length(X)
    Thresh(t1,1)=thresh;
end
% print value functions
%%%%%%%%%%%%
% writing styles: sansserif
%%%%%%%%%%%%
GraphicFontSize=20;
%%%%%%%%%%%%

figure(1)
plot(X,h,'--k')
xlabel('$x$','FontSize',GraphicFontSize,'FontName','sansserif','Interpreter','latex');
ylabel('$V_{H}^{x^{*}}(x)$ , $V_{L}^{x^{*}}(x)$ , $D^{x^{*}}(x)$','FontSize',GraphicFontSize,'FontName','sansserif','Interpreter','latex')
ax=gca;
ax.FontSize = GraphicFontSize-5;
ax.FontName = 'SansSerif';
ax.TickLabelInterpreter='latex';
hold on
plot(X,l,'-.k')
hold on
plot(X,d,'-k')
hold on
plot(X,Thresh,':k')
hold off
legend({'$V_{H}^{x^{*}}(x)$','$V_{L}^{x^{*}}(x)$','$D^{x^{*}}(x)$','$c/ \lambda$'},'Location','eastoutside','FontSize',GraphicFontSize-5,'FontName','sansserif','Interpreter','latex','EdgeColor',[1 1 1])
box off

%% Approximation of Dividend Function
n=30;       %30 is for many cases a good choice
a=0;
b=1;

[coef,weights] = chebweights(@spline,n,X,d);
[y,xpoly]=chebapprox(weights,X,n,a,b);

% create concrete function of the chebyshev-polynomial where you can use
% the normal grid on [0,1]

[Q,F,f]=chebfun(weights,a,b,c,lambda);

% approximation error at x0 , xstar-5 ,xstar , xstar+5 , x1
if xstar~=0 && xstar~=1
    intv=[0,X(wherexstar-5,1),X(wherexstar,1),X(wherexstar+5,1),1]';
else
    intv=[0,1]';

end
approxcheb=F(intv);
approxspline=spline(X,d,intv);
approxerr=abs(approxspline)-abs(approxcheb);

% Check monotonicity

% [monoto , change]=monotoni(f,X);
[monoto,change,monotosimple,monotol,monotoh]=monotonialter(f,X,d,l,h);

% %% Plot Cheby approximation

figure(2)
plot(X,y,'-.k')
xlabel('$x$','FontSize',GraphicFontSize,'FontName','sansserif','Interpreter','latex');
ylabel('$D^{x^{*}}(x)$','FontSize',GraphicFontSize,'FontName','sansserif','Interpreter','latex')
ax=gca;
ax.FontSize = GraphicFontSize-5;
ax.FontName = 'SansSerif';
ax.TickLabelInterpreter='latex';
xlim([0.4 0.8]);
% grid on
% grid minor
hold on
plot(X,d,'-k')
% hold on
% plot(X,F(X),'--r')
hold off
legend({'$D^{x^{*}}_{cheby}(x)$','$D^{x^{*}}_{cubic}(x)$'},'Location','eastoutside','FontSize',GraphicFontSize-5,'FontName','sansserif','Interpreter','latex','EdgeColor',[1 1 1])

%% Find indifference points with nsold.m

% find candidates
[candidates,IOC]=candidatesalter(@spline,X,d,thresh);
if IOC==1
    totinnersteps=length(candidates(:,2));
    sol=zeros(totinnersteps,1);
    for step=1:totinnersteps
        xcandi=candidates(step,2);
%         [sol,it_hist,ierr,x_hist]=nsold(xcandi,Q,[10^-8,10^-8]);
    end
else
    warning('no x found under this configuration which fulfills indifference condition')
end



 