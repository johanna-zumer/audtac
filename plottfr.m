function tfrsave = plottfr(tfr_diff,tfr_tpa,tfr_mspn,xlimlim,purpose,stat,band,param,ll,tacbasemax,plotallmask,chanplot,stattimwin,soades,fdir,statuse)
% function tfrsave = plottfr(tfr_diff,tfr_tpa,tfr_mspn,xlimlim,purpose,stat,band,param,ll,tacbasemax,plotallmask,chanplot,stattimwin,soades,fdir,statuse)

close all

colorblindD  =[204 101 0]/256;
colorblindApT=[5 165 255]/256;
colorblindMpN=[0 59 179]/256;
colorblindT  =[255 51 166]/256;
colorblindA  =[0 158 115]/256;
colorblindM  =[0 0 0]/256;
colorblindN  =[128 128 128]/256;

cfg=[];
tmpA=ft_selectdata(cfg,tfr_diff);
tmpA.mask=logical(zeros(size(tmpA.powspctrm)));
tmpA.mask(:,:,dsearchn(tmpA.time',stat.time(1)):dsearchn(tmpA.time',stat.time(end)))=stat.(statuse);
tmpsA=ft_selectdata(cfg,tfr_mspn);
tmpsA.mask=logical(zeros(size(tmpsA.powspctrm)));
tmpsA.mask(:,:,dsearchn(tmpsA.time',stat.time(1)):dsearchn(tmpsA.time',stat.time(end)))=stat.(statuse);
tmpuA=ft_selectdata(cfg,tfr_tpa);
tmpuA.mask=logical(zeros(size(tmpuA.powspctrm)));
tmpuA.mask(:,:,dsearchn(tmpuA.time',stat.time(1)):dsearchn(tmpuA.time',stat.time(end)))=stat.(statuse);

pow{1}=tmpuA;
pow{2}=tmpsA;
pow{3}=tmpA;

% fix hard setting of 49:53?????
baset=pow{2};
[ll pow{1}.time(49:53)]
baset.(param)=repmat(nanmean(cat(3,nanmean(pow{1}.(param)(:,:,49:53),3),nanmean(pow{2}.(param)(:,:,49:53),3)),3),[1 1 length(pow{1}.time)]);
baset.mask=repmat(nanmean(pow{2}.mask(:,:,49:53),3),[1 1 length(pow{1}.time)]);  % use 'nul' as baseline for all three conditions for plotting
cfg=[];
cfg.parameter={param 'mask'};
cfg.operation='subtract';
powb{1}=ft_math(cfg,pow{1},baset);
powb{2}=ft_math(cfg,pow{2},baset);
powb{3}=pow{3};

for pp=1:3
  cfg=[];
  if ll==1
    cfg.latency=[tacbasemax(ll) 1.2];
  elseif ll==3
    cfg.latency=[tacbasemax(ll) 1.3];
  elseif ll==4
    cfg.latency=[tacbasemax(ll) 1.3];
  elseif ll==5
    cfg.latency=[tacbasemax(ll) 1.3];
  elseif ll==6
    cfg.latency=[tacbasemax(ll) 1.3+.02];
  elseif ll==7
    cfg.latency=[tacbasemax(ll) 1.3+.07];
  elseif ll==9
    cfg.latency=[tacbasemax(ll) 1.3+.2];
  end
  powb{pp}=ft_selectdata(cfg,powb{pp});
  
  powb{pp}.avg=squeeze(powb{pp}.(param));
  powb{pp}.mask=squeeze(powb{pp}.mask);
  powb{pp}.dimord='chan_time';
  powb{pp}=rmfield(powb{pp},'freq');
  powb{pp}=rmfield(powb{pp},param);
  if plotallmask
    powb{pp}.mask=repmat(ceil(mean(powb{pp}.mask,1)),[size(powb{pp}.mask,1) 1]);
  end
end % pp

for cg=1:length(chanplot)
% for cg=1:2
  cfg=[];
  cfg.parameter='avg';
  cfg.layout='elec1010.lay';
  cfg.ylim='maxabs';
%   switch param
%     case 'powspctrm'
%       switch band
%         case 'theta'
%           cfg.ylim=[-2.5 2.5];
%         case 'alpha'
%           cfg.ylim=[-4.5 4.5];
%         case 'beta'
%           cfg.ylim=[-.8 .8];
%       end
%     case 'plvabs'
%       cfg.ylim=[-.1 .4];
%   end
  cfg.linewidth=3;
  cfg.xlim=xlimlim;
  cfg.channel=chanplot{cg};
  cfg.graphcolor=[colorblindApT; colorblindMpN; colorblindD];
  cfg.interactive='no';
  cfg.maskparameter='mask';
  cfg.maskstyle='box'; % default
  figure(ll+10*(cg+1))
  cfg = ft_singleplotER(cfg,powb{1},powb{2},powb{3});
  hold on;plot(powb{1}.time,0,'k');
  cfg.ylim=get(gca,'YLim');
  set(gca,'XTick',[-.6:.1:1.8])
  set(gca,'XTickLabel',{ '' '-0.5' '' ' ' '' ' ' '0' ' ' '' ' ' '' '0.5 ' '' ' ' ''  ' ' '1.0'  '' ' ' ''  ' ' '1.5' '' ' ' '' })
  set(gca,'FontSize',30)
  title([])
  plot([stattimwin(1) stattimwin(1)],cfg.ylim,'k--','LineWidth',6)
  plot([stattimwin(2) stattimwin(2)],cfg.ylim,'k--','LineWidth',6)
  plot([0 0],cfg.ylim,'Color',colorblindT,'LineWidth',6)
  plot([soades(ll) soades(ll)],cfg.ylim,'Color',colorblindA,'LineWidth',6)
  plot(cfg.xlim,[0 0],'Color','k')
  axis([cfg.xlim(1) cfg.xlim(2) cfg.ylim(1) cfg.ylim(2)])
  if cg==3
    legend('A+T','MS+N','Difference')
  end
end % cg


print(ll+20,[fdir purpose '_tfr_' param '_tacPaud_MSpN_diff_chan1_' num2str(ll) band '_' statuse '.png'],'-dpng')
try 
  print(ll+30,[fdir purpose '_tfr_' param '_tacPaud_MSpN_diff_chan2_' num2str(ll) band '_' statuse '.png'],'-dpng')
end
try 
  print(ll+40,[fdir purpose '_tfr_' param '_tacPaud_MSpN_diff_chan3_' num2str(ll) band '_' statuse '.png'],'-dpng')
end
print(ll+20,[fdir purpose '_tfr_' param '_tacPaud_MSpN_diff_chan1_' num2str(ll) band '_' statuse '.eps'],'-depsc')
try
  print(ll+30,[fdir purpose '_tfr_' param '_tacPaud_MSpN_diff_chan2_' num2str(ll) band '_' statuse '.eps'],'-depsc')
end
try
  print(ll+40,[fdir purpose '_tfr_' param '_tacPaud_MSpN_diff_chan3_' num2str(ll) band '_' statuse '.eps'],'-depsc')
end




selchan1 = match_str(powb{3}.label, chanplot{1});
tfrsave.time=powb{3}.time(nearest(powb{3}.time, xlimlim(1)):nearest(powb{3}.time, xlimlim(2)));
tfrsave.courseFC=mean(powb{3}.avg(selchan1,nearest(powb{3}.time,xlimlim(1)):nearest(powb{3}.time, xlimlim(2))),1);
try
  selchan2 = match_str(powb{3}.label, chanplot{2});
  tfrsave.courseOP=mean(powb{3}.avg(selchan2,nearest(powb{3}.time,xlimlim(1)):nearest(powb{3}.time, xlimlim(2))),1);
end

%% plot topoplot
cfg=[];
cfg.parameter=param;
cfg.layout='elec1010.lay';
cfg.maskalpha=0.5;
cfg.highlight='on';
cfg.comment='auto';
cfg.colorbar='yes';
masktime_tmp=find(squeeze(any(mean(stat.(statuse)(:,1,:),2),1)));
difftimes=diff(find(squeeze(any(mean(stat.(statuse)(:,1,:),2),1))));
clear masktime
if any(difftimes>1)
  finddifftimes=find(difftimes>1);
  for dd=1:length(finddifftimes)+1
    if dd==1
      masktime{dd}=masktime_tmp(1):masktime_tmp(find(difftimes>1));
    elseif dd==length(finddifftimes)+1
      masktime{dd}=masktime_tmp(finddifftimes(dd-1)+1):masktime_tmp(end);
    else
      masktime{dd}=masktime_tmp(finddifftimes(dd-1)+1):masktime_tmp(finddifftimes(dd));
    end
  end
else
  masktime{1}=masktime_tmp;
end
if ll==1
  switch band
    case 'theta'
      switch param
        case 'powspctrm'
          switch purpose
            case {'ica' 'within'}
              masktime_tmp=masktime{1};
              clear masktime
              masktime{1}=masktime_tmp(1:45);
              masktime{2}=masktime_tmp(46:end);
          end
      end
  end
end
switch band
  case 'theta'
    switch param
      case 'powspctrm'
        switch purpose
          case {'ica' 'within'}
            masktime_tmp=masktime{1};
            clear masktime
            masktime{1}=masktime_tmp(1:45);
            masktime{2}=masktime_tmp(46:end);
          case 'graphabs'
            if ll==7
              masktime{1}=11:48;
            elseif ll==3
              masktime{1}; % keep as is
            end
        end
      case 'plvabs'
        switch purpose
          case 'graphabs'
            if ll==7
              masktime{1}=masktime{1}(1:27); % this reduces to .07-.33 (or effectively 0-.26) to match AT70
            elseif ll==3
              masktime{1}; % keep as is
            end
        end
    end
end
if ll==9
  baseline2=[-.5 -.5+.08];
else
  baseline2=[tacbasemax(1) tacbasemax(1)+.08];
end
for dd=1:length(masktime)
  if ~isempty(masktime{dd})
    cfg.xlim=[stat.time(masktime{dd}(1)) stat.time(masktime{dd}(end))];
    cfg.highlightchannel=stat.label(find(ceil(mean(mean(stat.(statuse)(:,1,dsearchn(stat.time',cfg.xlim(1)):dsearchn(stat.time',cfg.xlim(2)) ),2),3))));
    cfg.baseline=baseline2;
    cfg.highlightsize=25;
    cfg.highlightsymbol='.';
    switch param
      case 'powspctrm'
        switch band
          case 'theta'
            cfg.zlim=[-1 1];
          case 'alpha'
            cfg.zlim=[-4 4];
          case 'beta'
            cfg.zlim=[-1.1 1.1];
        end
      case 'plvabs'
        switch band
          case 'theta'
            cfg.zlim=[-.15 .15];
          case 'alpha'
            cfg.zlim=[-.1 .1];
          case 'beta'
            cfg.zlim=[-.1 .1];
        end
    end
    figure(10*ll+4);
    ft_topoplotTFR(cfg,tmpA);
    print(10*ll+4,[fdir 'topoDiff_' param '_' band '_' num2str(ll) '_time' num2str(dd) '_' statuse '.eps'],'-depsc2')
    print(10*ll+4,[fdir 'topoDiff_' param '_' band '_' num2str(ll) '_time' num2str(dd) '_' statuse '.png'],'-dpng')
    tfrsave.topo{dd}=mean(tmpA.(param)(:,1,nearest(tmpA.time,cfg.xlim(1)):nearest(tmpA.time,cfg.xlim(2))),3);
    switch purpose
      case 'within'
        switch param
          case 'powspctrm'
            switch band
              case 'theta'
                cfg.zlim=[-1.5 1.5];
              case 'alpha'
                cfg.zlim=[-4 4];
              case 'beta'
                cfg.zlim=[-1.1 1.1];
            end
          case 'plvabs'
            switch band
              case 'theta'
                cfg.zlim=[0 .4];
              case 'alpha'
                cfg.zlim=[-.1 .1];
              case 'beta'
                cfg.zlim=[-.1 .1];
            end
        end
        figure(10*ll+4+100);
        ft_topoplotTFR(cfg,tmpuA);
        print(10*ll+4+100,[fdir 'topoU_' param '_' band '_' num2str(ll) '_time' num2str(dd) '_' statuse '.eps'],'-depsc2')
        print(10*ll+4+100,[fdir 'topoU_' param '_' band '_' num2str(ll) '_time' num2str(dd) '_' statuse '.png'],'-dpng')
        tfrsave.topoU{dd}=mean(tmpuA.(param)(:,1,nearest(tmpuA.time,cfg.xlim(1)):nearest(tmpuA.time,cfg.xlim(2))),3);
        figure(10*ll+4+200);
        ft_topoplotTFR(cfg,tmpsA);
        print(10*ll+4+200,[fdir 'topoM_' param '_' band '_' num2str(ll) '_time' num2str(dd) '_' statuse '.eps'],'-depsc2')
        print(10*ll+4+200,[fdir 'topoM_' param '_' band '_' num2str(ll) '_time' num2str(dd) '_' statuse '.png'],'-dpng')
        tfrsave.topoM{dd}=mean(tmpsA.(param)(:,1,nearest(tmpsA.time,cfg.xlim(1)):nearest(tmpsA.time,cfg.xlim(2))),3);
    end
  end
end

end

