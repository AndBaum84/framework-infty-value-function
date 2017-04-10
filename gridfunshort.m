function [ X , wherexstar ] = gridfunshort( xstar , gridintervals , IOsgrid , IOxstar)
%Generating the full grid vector.
%   The function will be approximated on the main intervall of grid1.
%   There are critical points like xstar in [0,1].
%   Around those critical point the grid will be smaller. 
%   Work-shirk: xstar element of (0,1).
%   Full work:  xstar is at the position of x=1=xstar.
%   Full shirk: xstar is at the position of x=0=xstar.
%   Output:
%   X:          grid for approximation
%   wherexstar: index of xstar in the grid (X) vector


% Check inputs
if IOxstar==1
        
    elseif IOxstar==0
        if IOsgrid==1
            warning('xstar does not exist on the grid!')
        elseif IOsgrid==0
            warning('Even though you don`t want to, xstar does exist on the grid!')
        end
    else
        error('IOxstar must be either 0 or 1')
end

griddigits=numel(num2str(1/gridintervals));
xstardigits=numel(num2str(xstar))-2;
xstarinterval=10^(-xstardigits);
griddigitstrue=numel(num2str(1/gridintervals))-2;
xinterval=10^(-griddigitstrue);

if xstar~=0 && xstar~=1 && xstardigits>griddigitstrue
    if xstar-xstarinterval<=xinterval
        Aroundxstar=[xstar,xstar+xstarinterval];
        candi=round(xstar+xstarinterval,griddigitstrue);
        if xstar-candi<=0
            afteraroundxstar=candi;
        else
            afteraroundxstar=candi+xinterval;
        end
        Xgrid=[0,Aroundxstar,afteraroundxstar:xinterval:1]';
    elseif xstar+xstarinterval>=1-xinterval
        Aroundxstar=[xstar-xstarinterval,xstar];
        candi=round(xstar-xstarinterval,griddigitstrue);
        if xstar-candi<=0
            beforearoundxstar=candi-xinterval;
        else
            beforearoundxstar=candi;
        end
        Xgrid=[0:xinterval:beforearoundxstar,Aroundxstar,1]';
    else
        Aroundxstar=[xstar-xstarinterval,xstar,xstar+xstarinterval];
        candibefore=round(xstar-xstarinterval,griddigitstrue);
        if xstar-candibefore<=0
            beforearoundxstar=candibefore-xinterval;
        else
            beforearoundxstar=candibefore;
        end
        candiafter=round(xstar+xstarinterval,griddigitstrue);
        if xstar-candiafter<=0
            afteraroundxstar=candiafter;
        else
            afteraroundxstar=candiafter+xinterval;
        end
        Xgrid=[0:xinterval:beforearoundxstar,Aroundxstar,afteraroundxstar:xinterval:1]';
    end   
    %Xgrid=unique(Xgrid,'stable');
else
    lgrid=1/gridintervals;
    Xgrid=[0:lgrid:1]';
end
    

if IOsgrid==1
    error('does not work with IOsgrid=1');
elseif IOsgrid==0
    X=Xgrid;
else
    error('YN must be either 0 or 1')
end


if IOxstar==1
    take=max([xstardigits,griddigitstrue]);
    wherexstar=find(round(X,take)==round(xstar,take));
    if length(wherexstar)>1
        clear wherexstar
        wherexstar=find(round(X,take+2)==round(xstar,take+2));
    end
    if length(wherexstar)>1
        clear wherexstar
        wherexstar=find(round(X,take+2)==round(xstar,take+3));
    end
elseif IOxstar==0
    b2=(b-a)/2;
    wherexstar=a+b2;
end

end