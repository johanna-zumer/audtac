function [aa1,aaa1,vvv1]=pca_masked(mask_use,fulldata_reshape,labels,pcuse,timepts,xlimlim,fdir,name)

[s1,s2]=size(mask_use);
timemask=find(nanmean(mask_use,1));
chanmask=find(nanmean(mask_use,2));
blobmask=zeros(s1,s2);
blobmask(chanmask,timemask)=1;

maskind=find(reshape(mask_use,[1 s1*s2]));
blobind=find(reshape(blobmask,[1 s1*s2]));
maskeddata=fulldata_reshape(:,maskind);
[aa1,bb1,vv1]=svd(maskeddata,'econ');
for ct=1:length(maskind),bb7(ct,:)=regress(maskeddata(:,ct),aa1(:,1:7));end
[mx,mind]=max(abs(bb7)');
fakedata=nan(1,s1*s2);
fakedata(maskind)=mind;
figure;imagesc(reshape(fakedata,[s1 s2]));caxis([0 7]);
figure;imagesc(reshape(fakedata,[s1 s2]));caxis([0 4]);


projdata=inv(aa1)*fulldata_reshape;
% or use ppca (probablistically fills in missing data)
projblobdata=inv(aa1)*fulldata_reshape(:,blobind);


if isempty(pcuse)
  pcuse=input('which PC to use? ');
end
[aaa1,bbb1,vvv1]=svd(reshape(projdata(pcuse,:),[s1 s2]));
[aaa1b,bbb1b,vvv1b]=svd(reshape(projblobdata(pcuse,:),[length(chanmask) length(timemask)]));

figure(21);lh=plot(timepts,vvv1(:,1),'k');
ylim([-max(abs(vvv1(:,1)))-.02 max(abs(vvv1(:,1)))+.02]);
xlim(xlimlim);
hold on; plot(timepts,zeros(1,length(vvv1(:,1))),'k');
set(lh(1),'linewidth',3);
set(lh(1).Parent,'FontSize',30);
set(gca,'XTick',[-.6:.1:1.1])
set(gca,'XTickLabel',{ '' '-0.5' '' ' ' '' ' ' '0' ' ' '' ' ' '' '0.5 ' '' ' ' ''  ' ' '1.0' ' '})
print(21,[fdir 'PC_time_' name '.eps'],'-painters','-depsc')

statplot.avg=aaa1(:,1);
statplot.time=1;
statplot.dimord='chan_time';
statplot.label=labels;
figure(23);cfg=[];cfg.latency=1;cfg.layout='eeg1010';cfg.zlim='maxabs';ft_topoplotER(cfg,statplot);
print(23,[fdir 'PC_topo_' name '.eps'],'-painters','-depsc')

figure(22);bar(aa1(:,pcuse));ylim([-max(abs(aa1(:,pcuse)))-.02 max(abs(aa1(:,pcuse)))+.02]);
xticklabels({'AT500' 'AT70' 'AT20' 'AT0' 'TA20' 'TA70' 'TA500'})
set(gca,'XTickLabelRotation',45)
set(gca,'FontSize',30);
aa=get(gca,'Children');
set(aa,'FaceColor','black')
print(22,[fdir 'PC_bar_' name '.eps'],'-painters','-depsc')
