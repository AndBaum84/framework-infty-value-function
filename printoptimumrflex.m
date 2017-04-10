%% Print candidated with optimum for one parameter 'c' changing

clear all;
clc;

%%%%%%%%%%%%%
printall=1;         % print all points or only candidates
%%%%%%%%%%%%%

XSTAR=load('bigXSTAR_04_07_flex_flex.dat');
% XSTAR2=load('bigXSTAR_04_07_02_flexsmall.dat');
% 
% XSTAR=[XSTAR(:,1:20);XSTAR2(:,1:20)];

DISCOUNT=unique(XSTAR(:,10));
howmanyr=length(DISCOUNT);

if printall==0
    selectcandonly=find(XSTAR(:,16)==1);
    XSTAR=XSTAR(selectcandonly,:);
end
    
selectpos=find(XSTAR(:,2)>=0);
XSTAR=XSTAR(selectpos,:);

%% at which cost
c=0.2;

selectc=find(XSTAR(:,9)==c);
XSTAR=XSTAR(selectc,:);

select1=find(XSTAR(:,10)==DISCOUNT(2));
XSTAR1=XSTAR(select1(1:end-35),1:2);


select2=find(XSTAR(:,10)==DISCOUNT(5));
XSTAR2=XSTAR(select2,1:2);

select3=find(XSTAR(:,10)==DISCOUNT(6));
XSTAR3=XSTAR(select3,1:2);

select4=find(XSTAR(:,10)==DISCOUNT(7));
XSTAR4=XSTAR(select4,1:2);

select5=find(XSTAR(:,10)==DISCOUNT(10));
XSTAR5=XSTAR(select5,1:2);

% select6=find(XSTAR(:,10)==DISCOUNT(10));
% XSTAR6=XSTAR(select6,1:2);

line45=[0:0.01:1]';

%%%%%%%%%%%%%%
GraphicFontSize=20;
%%%%%%%%%%%%%%

figure(1)
scatter(XSTAR1(:,1),XSTAR1(:,2),'o','MarkerEdgeColor',[0 0 0])
xlabel('$x^{*}$','FontSize',GraphicFontSize,'FontName','sansserif','Interpreter','latex');
ylabel('$\hat{x}^{x^{*}}$','FontSize',GraphicFontSize,'FontName','sansserif','Interpreter','latex')
ax=gca;
ax.FontSize = 15;
ax.FontName = 'SansSerif';
ax.TickLabelInterpreter='latex';
hold on
scatter(XSTAR2(:,1),XSTAR2(:,2),'o','MarkerEdgeColor',[0.2 0.2 0.2])
hold on
scatter(XSTAR3(:,1),XSTAR3(:,2),'o','MarkerEdgeColor',[0.4 0.4 0.4])
hold on
scatter(XSTAR4(:,1),XSTAR4(:,2),'o','MarkerEdgeColor',[0.5 0.5 0.5])
hold on
scatter(XSTAR5(:,1),XSTAR5(:,2),'o','MarkerEdgeColor',[0.6 0.6 0.6])
hold on
% scatter(XSTAR6(:,1),XSTAR6(:,2),'o','MarkerEdgeColor',[0.6 0.6 0.6])
% hold on
plot(line45,line45,':','Color',[0 0 0])
hold off
legend({'$r=0.05$','$r=0.17$','$r=0.2$','$r=0.22$','$r=0.35$'},'Location','eastoutside','FontSize',GraphicFontSize-5,'FontName','sansserif','Interpreter','latex','EdgeColor',[1 1 1])
box off