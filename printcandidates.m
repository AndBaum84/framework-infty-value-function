%% Print candidates only for one parameter setting
clear all;
clc;

XSTAR=load('bigXSTAR_04_07_02_017.dat');
line45=[0:0.01:1]';

% select c and r
c=0.2;
r=0.17;

selectc=find(XSTAR(:,9)==c);
XSTAR=XSTAR(selectc,:);
selectr=find(XSTAR(:,10)==r);
XSTAR=XSTAR(selectr,:);
selectcandonly=find(XSTAR(:,16)==1);
XSTARcandonly=XSTAR(selectcandonly,1:2);
selctoptonly=find(XSTAR(:,16)==2);
XSTARoptonly=XSTAR(selctoptonly,1:2);



%%%%%%%%%%%%%%
GraphicFontSize=20;
%%%%%%%%%%%%%%

figure(1)
scatter(XSTARcandonly(:,1),XSTARcandonly(:,2),'o','MarkerEdgeColor',[0.501960813999176 0.501960813999176 0.501960813999176])
xlabel('$x^{*}$','FontSize',GraphicFontSize,'FontName','sansserif','Interpreter','latex');
ylabel('$\hat{x}^{x^{*}}$','FontSize',GraphicFontSize,'FontName','sansserif','Interpreter','latex')
ax=gca;
ax.FontSize = 15;
ax.FontName = 'SansSerif';
ax.TickLabelInterpreter='latex';
hold on
scatter(XSTARoptonly(:,1),XSTARoptonly(:,2),'o','MarkerEdgeColor',[0 0 0])
hold on
plot(line45,line45,':','Color',[0 0 0])
hold off
legend({'Newton-method','Newton- and bisection-method','$x^{*}=\hat{x}^{x^{*}}$'},'Location','eastoutside','FontSize',GraphicFontSize-5,'FontName','sansserif','Interpreter','latex','EdgeColor',[1 1 1])
box off
            