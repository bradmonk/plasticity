%
%
%
randn('state',0);
tmp = randn(size(x));
x = 1:0.1:100;
y = slidingavg(tmp,1000);
z1 = 0.05* exp(-x/50);
z2 = 0.08* exp(-x/50);

confplot(x,y,z1,z2,'Color',[1 0 0],'LineWidth',2);
grid on; box off;
xlabel('time [s]','FontName','Helvetica','FontSize',30);
ylabel('[a.u.]','FontName','Helvetica','FontSize',30);
legend('Estimated value','95% confidence boundaries');
set(gca,'FontName','Helvetica','FontSize',20,'YLim',[-0.2 0.15]);

title(title('CONFPLOT - {\copyright} 2002 Michele Giugliano, PhD');
%print('testpic.jpg','-djpeg90','-r300')