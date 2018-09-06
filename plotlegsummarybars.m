function plotlegsummarybars(data,figind,figname,starind,yminmax)

colorblindD  =[204 101 0]/256;
colorblindApT=[5 165 255]/256;
colorblindMpN=[0 59 179]/256;
time=[-.6 -.35 -.1 0 .1 .35 .6];

try close(figind); end
figure(figind);
[ph] = plot(time, data);
xlim([-.8 .8])
ylim(yminmax)
hold on
set(ph(1),'Color',colorblindApT,'Marker','^','MarkerSize',15,'LineWidth',3,'LineStyle',':')
set(ph(2),'Color',colorblindMpN,'Marker','v','MarkerSize',15,'LineWidth',3,'LineStyle',':')
set(ph(3),'Color',colorblindD,'LineWidth',3,'LineStyle',':')
ab=plot(time(starind), data(starind,3),'Color',colorblindD,'Marker','*','MarkerSize',15,'LineWidth',3,'LineStyle','none');
% set(ph(3),'Color',colorblindD,'MarkerIndices',starind,'Marker','*','MarkerSize',15,'LineWidth',3,'LineStyle',':')
aa=plot(time(2:end-1), data(2:end-1,:));
nonstar=setdiff(1:7,starind);
for nn=1:length(nonstar)
  ac(nn)=plot(time(nonstar(nn)), data(nonstar(nn),3),'Color',colorblindD,'Marker','o','MarkerSize',8,'LineWidth',3,'LineStyle','none');
end
set(aa(1),'Color',colorblindApT,'LineWidth',3,'LineStyle','-')
set(aa(2),'Color',colorblindMpN,'LineWidth',3,'LineStyle','-')
set(aa(3),'Color',colorblindD,'LineWidth',3,'LineStyle','-')
ab=plot(xlim, [0 0],'Color',colorblindD);
set(gca,'FontSize',15)
set(gca,'LineWidth',3)
set(gca,'Xtick',time)
set(gca,'XTickLabel',{ '-500' '-70' '-20' '0' '20' '70' '500'})
xlabel('Multisensory Asynchrony (ms)')
legend('A+T','MS+N','Difference')
% ylabel('ERP (uV)')
print(figind,figname,'-painters','-depsc')

