%% Print value of quality under several gridpoints
lambda=0.4;
c=0.2;
thresh=c/lambda;

%%%%%%%%%%%%
GraphicFontSize=20;
%%%%%%%%%%%%

d10=load('d10.dat');
X10=load('X10.dat');
d100=load('d100.dat');
X100=load('X100.dat');
d1000=load('d1000.dat');
X1000=load('X1000.dat');
d10000=load('d10000.dat');
X10000=load('X10000.dat');

Thresh=zeros(length(X1000),1);
for t1=1:length(X1000)
    Thresh(t1,1)=thresh;
end

figure(1)
plot(X10,d10,'--k','LineWidth',0.3)
xlabel('$x$','FontSize',GraphicFontSize,'FontName','sansserif','Interpreter','latex');
ylabel('$D^{x^{*}}(x)$','FontSize',GraphicFontSize,'FontName','sansserif','Interpreter','latex')
ax=gca;
ax.FontSize = GraphicFontSize-5;
ax.FontName = 'SansSerif';
ax.TickLabelInterpreter='latex';
hold on
plot(X100,d100,'-.k','LineWidth',0.3)
hold on
plot(X1000,d1000,'-k','LineWidth',0.3)
hold on
plot(X10000,d10000,':k','LineWidth',0.3)
% hold on
% plot(X1000,Thresh,'Color',[0.5 0.5 0.5])
hold off
legend1=legend({'$10$','$100$','$1000$','$10000$'},'Location','eastoutside','FontSize',GraphicFontSize-5,'FontName','sansserif','Interpreter','latex','EdgeColor',[1 1 1]);
title(legend1,'Number of grid-points:','FontSize',GraphicFontSize-5,'Interpreter','latex');
box off
xlim([0.25,0.4])
ylim([0.49,0.51])