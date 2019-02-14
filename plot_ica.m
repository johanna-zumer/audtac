function plot_ica(cf,runicaout,statplot,xlimlim,timeadd,timepts,llind,name,fdir)
% function plot_ica(cf,runicaout,statplot,xlimlim,timeadd,timepts,llind,name,fdir)

close all
timepts=timepts+timeadd;
figure(21);lh=plot(timepts,cf(:,1),'k');
ylim([-max(abs(cf(:,1)))-.02 max(abs(cf(:,1)))+.02]);
xlim(xlimlim);
hold on; plot(timepts,zeros(1,length(cf(:,1))),'k');
set(lh(1),'linewidth',3);
set(lh(1).Parent,'FontSize',30);
set(gca,'XTick',[-.6:.1:1.1])
set(gca,'XTickLabel',{ '' '-0.5' '' ' ' '' ' ' '0' ' ' '' ' ' '' '0.5 ' '' ' ' ''  ' ' '1.0' ' '})
print(21,[fdir name '_time_ICcomp' num2str(llind) '.eps'],'-painters','-depsc')

figure(22);bar(runicaout.topo(:,llind));ylim([-3 3])
xticklabels({'AT500' 'AT70' 'AT20' 'AT0' 'TA20' 'TA70' 'TA500'})
set(gca,'XTickLabelRotation',45)
set(gca,'FontSize',30);
aa=get(gca,'Children');
set(aa,'FaceColor','black')
print(22,[fdir name '_bar_ICcomp' num2str(llind) '.eps'],'-painters','-depsc')

figure(23);cfg=[];cfg.latency=1;cfg.zlim='maxabs';cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
print(23,[fdir name '_topo_ICcomp' num2str(llind) '.eps'],'-painters','-depsc')

end

