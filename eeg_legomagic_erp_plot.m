eeg_legomagic_preamble

%% Plotting above ERP sensor results


% First, plotting the final stats/ difference of conditions
% Below, plotting conditions on their own.

clearvars -except sub *dir
% chanplot{1}={'Fz' 'Cz' 'F1' 'F2' 'FC1' 'FC2' 'C1' 'C2'}; % frontocentral
% chanplot{2}={'CP5' 'P5' 'PO3' 'POz' 'PO4' 'P6' 'CP6' 'Pz' 'P3' 'P4'}; % occipital
% chanplot{3}={'FC6' 'FC4' 'F4' 'F6' 'F8' 'FT8' 'T8' 'C6' 'C4'}; % right frontotemporal
% chanlabel{1}='Frontocentral electrodes';
% chanlabel{2}='Occipital-parietal electrodes';
% chanlabel{3}='Right frontotemporal electrodes';
chanplot{1}={'Fz' 'Cz' 'F1' 'F2' 'FC1' 'FC2' 'C1' 'C2'}; % frontocentral
chanplot{2}={'CP5' 'POz' 'Pz' 'P3' 'P4' 'C4' 'O1' 'O2' 'P7' 'PO7'}; % occipital
chanplot{3}={'C1' 'Cz' 'C2' 'CP1' 'CPz' 'CP2' 'P1' 'Pz' 'P2'}; % Centred on CPz
chanplot{4}={'P3' 'P1' 'Pz' 'P2' 'P4' 'PO3' 'POz' 'PO4' 'O1' 'O2' 'Oz'};
chanlabel{1}='Frontocentral electrodes';
chanlabel{2}='Occipital-parietal electrodes';
iterflag=1;
sleep=0;
tt=3;
statwinorig=0;  % =1 means the latency range of .1 to .45 (first way);  =0 means 0 to .5 (final /better way)
ftver=0;
mcseed=13;
usetr=0;
plotallmask=1; % =0 means each subplot will use significane mask relevant for those sensors & time plotted only
% =1 means each subplot will use mask relevant for all sensors even those not plotted.

if sleep
  iter=11;
  ss=12;
  trialkc=0;
else
  iter=27;  % previous result; final way (easier!)
  %   iter=31;  % smart sampling
  ss=10;
  trialkc=-1;
end

if iterflag
  if sleep
    try
      load([edir 'tlock_statmc_sleep' num2str(sleep) '_ss' num2str(ss) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '_statwinorig' num2str(statwinorig) '.mat']);
      load([edir 'tlock_grind_sleep' num2str(sleep) '_ss' num2str(ss) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '_statwinorig' num2str(statwinorig) '.mat']);
    catch
      load([edir 'tlock_statmc_sleep' num2str(sleep) '_ss' num2str(ss) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat']);
      load([edir 'tlock_grind_sleep' num2str(sleep) '_ss' num2str(ss) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat']);
    end
  else
    try
%       load([edir 'tlock_statmc_sleep' num2str(sleep) '_iter' num2str(iter) '_usetr' num2str(usetr) '_statwinorig' num2str(statwinorig) '_ftver' num2str(ftver) 'mcseed' num2str(mcseed) '.mat']);
%       load([edir 'tlock_grind_sleep' num2str(sleep) '_iter' num2str(iter) '_usetr' num2str(usetr) '_statwinorig' num2str(statwinorig) '_ftver' num2str(ftver) '.mat']);
      load([edir 'tlock_statmc_sleep' num2str(sleep) '_iter' num2str(iter) '_statwinorig' num2str(statwinorig) '_ftver' num2str(ftver) 'mcseed' num2str(mcseed) '.mat']);
      load([edir 'tlock_grind_sleep' num2str(sleep) '_iter' num2str(iter) '_statwinorig' num2str(statwinorig) '_ftver' num2str(ftver) '.mat']);
    catch
      try
        load([edir 'tlock_statmc_sleep' num2str(sleep) '_iter' num2str(iter) '_statwinorig' num2str(statwinorig) '.mat']);
        load([edir 'tlock_grind_sleep' num2str(sleep) '_iter' num2str(iter) '_statwinorig' num2str(statwinorig) '.mat']);
      catch
        error('must use one of above')
        load([edir 'tlock_statmc_sleep' num2str(sleep) '_iter' num2str(iter) '.mat']);
        load([edir 'tlock_grind_sleep' num2str(sleep) '_iter' num2str(iter) '.mat']);
      end
    end
  end
else
  load([edir 'tlock_statmc_sleep' num2str(sleep) '.mat']);
  load([edir 'tlock_grind_sleep' num2str(sleep) '.mat']);
end

soalist=[1 3 4 5 6 7 9];
timwinstatflag=0; % =1 for using only stat.time, =0 for using fixed time e.g. [-0.5 1.0];
soades=[-.5 nan -.07 -.02 0 .02 .07 nan .5];
tacbasemax=min(-.15,soades-.15);

scalediff=1;

coloruse=varycolor(10);
% optimsied to be maximally apart and distinct
% 1  TacPAud
% 2  Nul
% 3  MS_synch
% 4  Tac
% 5  diff_TPAMSPN
% 6  diff_synchAsynch
% 7  MS
% 8  MS_asynch
% 9  Aud
% 10 MSpN

% Inspired by http://jfly.iam.u-tokyo.ac.jp/color/index.html (for 8 colours all okay for colour-blind)
colorblindD  =[204 101 0]/256;
colorblindApT=[5 165 255]/256;
colorblindMpN=[0 59 179]/256;
colorblindT  =[255 51 166]/256;
colorblindA  =[0 158 115]/256;
colorblindM  =[0 0 0]/256;
colorblindN  =[128 128 128]/256;

% colorblind(1,:)=[0 0 0];
% colorblind(2,:)=[.9 .6 0];
% colorblind(3,:)=[.35 .7 .9];
% colorblind(4,:)=[0 .6 .5];
% colorblind(5,:)=[.95 .9 .25];
% colorblind(6,:)=[0 .45 .7];
% colorblind(7,:)=[.8 .4 0];
% colorblind(8,:)=[.8 .6 .7];

incondflag=0; % =1 for within-condition, =0 for across-condition (only showing <=70)
if incondflag==1 % for within-condition
  timwin=[-0.5 1];
  topozlim=[-5 5];
  topozlimdiff=[-5 5];
else % for across-condition (only showing <=70)
  timwin=[-0.05 .57];  % this will cause an error for ll=9, but not a problem as we don't need that one for ICA
  topozlim=[-5 5];
  topozlimdiff=[-3 3];
end

timeadd=max(0,soades);

close all
clear tmp*
clear grind_TPA_MSPN
for ll=soalist
% for ll=[6 7 9]
  cfg=[];
  cfg.operation='subtract';
  cfg.parameter='individual';
  grind_TPA_MSPN{ll,tt,ss}=ft_math(cfg,grind_tacPaud_save{ll,tt,ss},grind_tacMSpN_save{ll,tt,ss});
  
  cfg=[];
  if timwinstatflag==1
    cfg.latency=[statt_mc{ll,tt,ss}.time(1) statt_mc{ll,tt,ss}.time(end)];
  elseif timwinstatflag==0
    cfg.latency=timwin+timeadd(ll);
    stattimwin=[statt_mc{ll,tt,ss}.time(1) statt_mc{ll,tt,ss}.time(end)];
  end
  cfg.channel=statt_mc{ll,tt,ss}.label;
  
  tmpn2=ft_selectdata(cfg,grind_nultlock_save{ll,tt,ss});
  tmpn2.dimord='chan_time';
  tmpn2.avg=squeeze(mean(tmpn2.individual,1));
  tmpn2=rmfield(tmpn2,'individual');
  
  tmpt4=ft_selectdata(cfg,grind_tactlock_save{ll,tt,ss});
  tmpt4.dimord='chan_time';
  tmpt4.avg=squeeze(mean(tmpt4.individual,1));
  tmpt4=rmfield(tmpt4,'individual');
  
  tmpa9=ft_selectdata(cfg,grind_audtlock_save{ll,tt,ss});
  tmpa9.dimord='chan_time';
  tmpa9.avg=squeeze(mean(tmpa9.individual,1));
  tmpa9=rmfield(tmpa9,'individual');
  
  tmpb7=ft_selectdata(cfg,grind_MStlock_save{ll,tt,ss});
  tmpb7.dimord='chan_time';
  tmpb7.avg=squeeze(mean(tmpb7.individual,1));
  tmpb7=rmfield(tmpb7,'individual');
  
  tmpu1=ft_selectdata(cfg,grind_tacPaud_save{ll,tt,ss});
  tmpu1.dimord='chan_time';
  tmpu1.avg=squeeze(mean(tmpu1.individual,1));
  tmpu1=rmfield(tmpu1,'individual');
  tmpu1.mask=statt_mc{ll,tt,ss}.mask;
  
  tmpm10=ft_selectdata(cfg,grind_tacMSpN_save{ll,tt,ss});
  tmpm10.dimord='chan_time';
  tmpm10.avg=squeeze(mean(tmpm10.individual,1));
  tmpm10=rmfield(tmpm10,'individual');
  tmpm10.mask=statt_mc{ll,tt,ss}.mask;
  
  tmpd5=ft_selectdata(cfg,grind_TPA_MSPN{ll,tt,ss});
  tmpd5.dimord='chan_time';
  tmpd5.avg=scalediff*squeeze(mean(tmpd5.individual,1));
  tmpd5=rmfield(tmpd5,'individual');
  tmpd5.mask=statt_mc{ll,tt,ss}.mask;
  
  if timwinstatflag==0
    %       tmpmask=mean(statt_mc{ll,tt,ss}.mask(match_str(tmpd5.label,cfg.channel),:),1);
    %       tmpd5.mask=zeros(1,length(tmpd5.time));
    %       tmpd5.mask(dsearchn(tmpd5.time',stattimwin(1)):dsearchn(tmpd5.time',stattimwin(end)))=tmpmask;
    tmpmask=statt_mc{ll,tt,ss}.mask;
    tmpu1.mask=zeros(size(tmpu1.avg,1),length(tmpu1.time));
    tmpu1.mask(:,dsearchn(tmpu1.time',stattimwin(1)):dsearchn(tmpu1.time',stattimwin(end)))=tmpmask;
    tmpm10.mask=zeros(size(tmpm10.avg,1),length(tmpm10.time));
    tmpm10.mask(:,dsearchn(tmpm10.time',stattimwin(1)):dsearchn(tmpm10.time',stattimwin(end)))=tmpmask;
    tmpd5.mask=zeros(size(tmpd5.avg,1),length(tmpd5.time));
    tmpd5.mask(:,dsearchn(tmpd5.time',stattimwin(1)):dsearchn(tmpd5.time',stattimwin(end)))=tmpmask;
  end
  
  if plotallmask
    tmpu1.mask=repmat(ceil(mean(tmpu1.mask,1)),[size(tmpu1.mask,1) 1]);
    tmpm10.mask=repmat(ceil(mean(tmpm10.mask,1)),[size(tmpm10.mask,1) 1]);
    tmpd5.mask=repmat(ceil(mean(tmpd5.mask,1)),[size(tmpd5.mask,1) 1]);
  end
  
  
  cfg=[];
  if ll==1
    cfg.latency=[tacbasemax(ll) .65];
  elseif ll==3
    cfg.latency=[tacbasemax(ll) .65];
  elseif ll==4
    cfg.latency=[tacbasemax(ll) .65];
  elseif ll==5
    cfg.latency=[tacbasemax(ll) .65];
  elseif ll==6
    cfg.latency=[tacbasemax(ll) .65+.02];
  elseif ll==7
    cfg.latency=[tacbasemax(ll) .65+.07];
  elseif ll==9
      cfg.latency=[tacbasemax(ll) .65+.50];
  end
  tmpu1=ft_selectdata(cfg,tmpu1);
  tmpm10=ft_selectdata(cfg,tmpm10);
  tmpd5=ft_selectdata(cfg,tmpd5);
        

  
  %   cfg=[];
  %   if timwinstatflag==1
  %     cfg.latency=[statt_mc{ll,tt,ss}.time(1) statt_mc{ll,tt,ss}.time(end)];
  %   elseif timwinstatflag==0
  %     cfg.latency=[-0.5 1];
  %     stattimwin=[statt_mc{ll,tt,ss}.time(1) statt_mc{ll,tt,ss}.time(end)];
  %   end
  %   tmp=ft_selectdata(cfg,grave_TPA_MSPN{ll,tt,ss});
  %   tmp.mask=statt_mc{ll,tt,ss}.mask;
  %   tmpm=ft_selectdata(cfg,grind_tacMSpN_save{ll,tt,ss});
  %   tmpu=ft_selectdata(cfg,grind_tacPaud_save{ll,tt,ss});
  %   tmpm.mask=statt_mc{ll,tt,ss}.mask;
  %   tmpu.mask=statt_mc{ll,tt,ss}.mask;
  %
  %   tmpm.avg=squeeze(mean(tmpm.individual,1));
  %   tmpu.avg=squeeze(mean(tmpu.individual,1));
  %   tmpm.var=squeeze(var(tmpm.individual,1));
  %   tmpu.var=squeeze(var(tmpu.individual,1));
  %   tmpm.dimord='chan_time';
  %   tmpu.dimord='chan_time';
  %   tmpm=rmfield(tmpm,'individual');
  %   tmpu=rmfield(tmpu,'individual');
  
  %   figure(10+ll);
  %   for cg=1:length(chanplot)
  %     cfg=[];
  %     cfg.parameter='avg';
  %     cfg.layout='elec1010.lay';
  %     cfg.ylim=[-3 7];
  %     cfg.linewidth=3;
  %     cfg.xlim=timwin;
  %     cfg.channel=chanplot{cg};
  %     cfg.graphcolor=coloruse([4 9 1],:);
  %     subplot(3,1,cg)
  %     ft_singleplotER(cfg,tmp4,tmp9,tmp1);
  %   end
  %   legend('Tactile','Auditory','Sum Unisensory')
  %   print(10+ll,[fdir 'erp_TacAudTacPAud_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  %
  %   figure(20+ll);
  %   for cg=1:length(chanplot)
  %     cfg=[];
  %     cfg.parameter='avg';
  %     cfg.layout='elec1010.lay';
  %     cfg.ylim=[-3 8];
  %     cfg.linewidth=3;
  %     cfg.xlim=timwin;
  %     cfg.channel=chanplot{cg};
  %     cfg.graphcolor=coloruse([2 7 10],:);
  %     subplot(3,1,cg)
  %     ft_singleplotER(cfg,tmp2,tmp7,tmp10);
  %   end
  %   legend('Null','Multisensory','MultSens + Null')
  %   print(20+ll,[fdir 'erp_NulMSandMSPN_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  
  %   figure(ll);
  %   for cg=1:length(chanplot)+1
%   for cg=1:length(chanplot)
  for cg=1:3
    cfg=[];
    cfg.parameter='avg';
    cfg.layout='elec1010.lay';
    cfg.ylim=[-3 7];
    cfg.linewidth=3;
    cfg.xlim=timwin+timeadd(ll);
    if cg>length(chanplot)
      cfg.channel=tmpd5.label(any(tmpd5.mask,2));
    else
      cfg.channel=chanplot{cg};
    end
%     cfg.graphcolor=coloruse([2 4 9 7],:);
    cfg.graphcolor=[colorblindN; colorblindT; colorblindA; colorblindM];
    cfg.interactive='no';
    %     subplot(4,1,cg)
    figure(ll+10*(cg-1))
    ft_singleplotER(cfg,tmpn2,tmpt4,tmpa9,tmpb7);
    hold on;plot(tmpn2.time,0,'k');
    set(gca,'XTick',[-.5:.1:1])
    set(gca,'XTickLabel',{'-0.5' '' ' ' '' ' ' '0' ' ' '' ' ' '' '0.5 ' '' ' ' ''  ' ' '1.0'})
    set(gca,'FontSize',30)
    title([])
%     plot([0 0],cfg.ylim,'Color',coloruse(4,:),'LineWidth',6)
%     plot([soades(ll) soades(ll)],cfg.ylim,'Color',coloruse(9,:),'LineWidth',6)
    plot([0 0],cfg.ylim,'Color',colorblindT,'LineWidth',6)
    plot([soades(ll) soades(ll)],cfg.ylim,'Color',colorblindA,'LineWidth',6)
    plot(cfg.xlim,[0 0],'Color','k')
    axis([timwin(1)+timeadd(ll)-0.05 timwin(2)+timeadd(ll) cfg.ylim(1) cfg.ylim(2)])
    pbaspect([2 1 1]);
    if cg==3
      legend('Null','Tactile','Auditory','Multisensory')
    end
  end
  %   print(ll,[fdir 'erp_presum_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  print(ll,[fdir 'erp_presum_FC_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  print(ll+10,[fdir 'erp_presum_OP_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  print(ll,[fdir 'erp_presum_FC_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-painters','-depsc')
  print(ll+10,[fdir 'erp_presum_OP_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-painters','-depsc')
  print(ll+20,[fdir 'erp_presum_SigLabels_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-painters','-depsc')
  close(ll+20)
  
  %   figure(30+ll);
  %     for cg=length(chanplot)+1
%   for cg=1:length(chanplot)
  for cg=1:3
    cfg=[];
    cfg.parameter='avg';
    cfg.layout='elec1010.lay';
    cfg.ylim=[-5 8];
    cfg.linewidth=3;
    cfg.xlim=[timwin(1)+timeadd(ll)-.05 timwin(2)+timeadd(ll)+.05]; 
    if cg>length(chanplot)
      cfg.channel=tmpd5.label(any(tmpd5.mask,2));
    else
      cfg.channel=chanplot{cg};
    end
%     cfg.graphcolor=coloruse([1 10 5],:);
    cfg.graphcolor=[colorblindApT; colorblindMpN; colorblindD];
    cfg.interactive='no';
    cfg.maskparameter='mask';
    cfg.maskstyle='box'; % default
    %     subplot(4,1,cg)
    %     subplot(1,2,cg)
    figure(ll+10*(cg+1))
    ft_singleplotER(cfg,tmpu1,tmpm10,tmpd5);
    hold on;plot(tmpn2.time,0,'k');
    set(gca,'XTick',[-.6:.1:1.1])
    set(gca,'XTickLabel',{ '' '-0.5' '' ' ' '' ' ' '0' ' ' '' ' ' '' '0.5 ' '' ' ' ''  ' ' '1.0' ' '})
    set(gca,'FontSize',30)
    title([])
    plot([stattimwin(1) stattimwin(1)],cfg.ylim,'k--','LineWidth',6)
    plot([stattimwin(2) stattimwin(2)],cfg.ylim,'k--','LineWidth',6)
%     plot([0 0],cfg.ylim,'Color',coloruse(4,:),'LineWidth',6)
%     plot([soades(ll) soades(ll)],cfg.ylim,'Color',coloruse(9,:),'LineWidth',6)
    plot([0 0],cfg.ylim,'Color',colorblindT,'LineWidth',6)
    plot([soades(ll) soades(ll)],cfg.ylim,'Color',colorblindA,'LineWidth',6)
    plot(cfg.xlim,[0 0],'Color','k')
    axis([cfg.xlim(1) cfg.xlim(2) cfg.ylim(1) cfg.ylim(2)])
    if cg==3
      legend('mask','A+T','MS+N','Difference')
    end
  end
  try
    print(ll+20,[fdir 'erp_tacPaud_MSpN_diff_FC_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    print(ll+30,[fdir 'erp_tacPaud_MSpN_diff_OP_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    print(ll+20,[fdir 'erp_tacPaud_MSpN_diff_FC_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-painters','-depsc')
    print(ll+30,[fdir 'erp_tacPaud_MSpN_diff_OP_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-painters','-depsc')
  end
  try
    if ~isempty(tmpd5.label(any(tmpd5.mask,2)))
      print(ll+40,[fdir 'erp_tacPaud_MSpN_diff_SigLabels_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-painters','-depsc')
    end
  end
  
  xlimlim=[timwin(1)+timeadd(ll)-.05 timwin(2)+timeadd(ll)+.05];
  selchan1 = match_str(tmpd5.label, chanplot{1});
  selchan2 = match_str(tmpd5.label, chanplot{2});
  erpsave{ll}.time=tmpd5.time(nearest(tmpd5.time, xlimlim(1)):nearest(tmpd5.time, xlimlim(2)));
  erpsave{ll}.courseFC=mean(tmpd5.avg(selchan1,nearest(tmpd5.time,xlimlim(1)):nearest(tmpd5.time, xlimlim(2))),1);
  erpsave{ll}.courseOP=mean(tmpd5.avg(selchan2,nearest(tmpd5.time,xlimlim(1)):nearest(tmpd5.time, xlimlim(2))),1);

  
  %   cfg=[];
  %   cfg.avgoverchan='yes';
  %   cfg.channel=chanplot{1};
  %   tmp1=ft_selectdata(cfg,tmp);
  %   if timwinstatflag==0
  %     tmp1mask=mean(tmp.mask(match_str(tmp.label,cfg.channel),:),1);
  %     tmp1.mask=zeros(1,length(tmp1.time));
  %     tmp1.mask(dsearchn(tmp1.time',stattimwin(1)):dsearchn(tmp1.time',stattimwin(end)))=tmp1mask;
  %   end
  %   tmp1.mask=logical(ceil(tmp1.mask));
  %   tmpm1=ft_selectdata(cfg,tmpm);
  %   tmpm1.mask=logical(ceil(tmpm1.mask));
  %   tmpu1=ft_selectdata(cfg,tmpu);
  %   tmpu1.mask=logical(ceil(tmpu1.mask));
  %
  %   figure(10*ll);
  %   cfg=[];
  %   cfg.parameter='avg';
  %   cfg.layout='elec1010.lay';
  %   cfg.maskparameter='mask';
  %   cfg.maskstyle='box'; % default
  %   cfg.ylim=[-3 7];
  %   cfg.linewidth=3;
  %   ft_singleplotER(cfg,tmpu1,tmpm1,tmp1);
  %   print(10*ll,[fdir 'erp_final_FC_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  %
  %   cfg=[];
  %   cfg.avgoverchan='yes';
  %   cfg.channel=chanplot{2};
  %   tmp1=ft_selectdata(cfg,tmp);
  %   if timwinstatflag==0
  %     tmp1mask=mean(tmp.mask(match_str(tmp.label,cfg.channel),:),1);
  %     tmp1.mask=zeros(1,length(tmp1.time));
  %     tmp1.mask(dsearchn(tmp1.time',stattimwin(1)):dsearchn(tmp1.time',stattimwin(end)))=tmp1mask;
  %   end
  %   tmp1.mask=logical(ceil(tmp1.mask));
  %   tmpm1=ft_selectdata(cfg,tmpm);
  %   tmpm1.mask=logical(ceil(tmpm1.mask));
  %   tmpu1=ft_selectdata(cfg,tmpu);
  %   tmpu1.mask=logical(ceil(tmpu1.mask));
  %
  %   figure(10*ll+1);
  %   cfg=[];
  %   cfg.parameter='avg';
  %   cfg.layout='elec1010.lay';
  %   cfg.maskparameter='mask';
  %   cfg.maskstyle='box'; % default
  %   cfg.ylim=[-3 3];
  %   cfg.linewidth=3;
  %   ft_singleplotER(cfg,tmpu1,tmpm1,tmp1);
  %   print(10*ll+1,[fdir 'erp_final_OP_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  
  %   if [tt==2 && any(ll==[3 4 5 6])] || [tt==3 && any(ll==[4 6 7])]
  if any(statt_mc{ll,tt,ss}.mask(:))
    masktime=find(any(tmpd5.mask,1));
    cfg=[];
    cfg.parameter='avg';
    cfg.layout='elec1010.lay';
    cfg.maskalpha=0.5;
    cfg.zlim=topozlim;
    cfg.highlight='on';
    cfg.highlightsize=12;
    cfg.xlim=[tmpd5.time(masktime(1)) tmpd5.time(masktime(end))];
    cfg.comment='no';
    %     if ll==3
    %       cfg.xlim=[.1 .35];
    %     elseif ll==4
    %       cfg.xlim=[.06 .42];
    %     elseif ll==5
    %       cfg.xlim=[-.04 .36];
    %     elseif ll==6
    %       cfg.xlim=[.1 .44];
    %     end
%     sigchannels=tmpd5.label(find(ceil(mean(tmpd5.mask(:,dsearchn(tmpd5.time',cfg.xlim(1)):dsearchn(tmpd5.time',cfg.xlim(2))),2))));
    sigchannels=tmpd5.label(find(ceil(mean(statt_mc{ll,tt,ss}.mask,2))));
    cfg.highlightchannel=sigchannels;
    figure(100+ll);
    ft_topoplotER(cfg,tmpu1);
    print(100+ll,[fdir 'erp_topoU_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    print(100+ll,[fdir 'erp_topoU_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-painters','-depsc')
    figure(110+ll);
    ft_topoplotER(cfg,tmpm10);
    print(110+ll,[fdir 'erp_topoM_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    print(110+ll,[fdir 'erp_topoM_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-painters','-depsc')
    figure(120+ll);
    cfg.zlim=topozlimdiff;
    ft_topoplotER(cfg,tmpd5);
    print(120+ll,[fdir 'erp_topoDiff_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    print(120+ll,[fdir 'erp_topoDiff_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-painters','-depsc')
    erpsave{ll}.topo{1}=mean(tmpd5.avg(:,nearest(tmpd5.time,cfg.xlim(1)):nearest(tmpd5.time,cfg.xlim(2))),2);
    
    %     figure(50+ll);
    %     %   for cg=1:length(chanplot)
    %     cfg=[];
    %     cfg.parameter='avg';
    %     cfg.layout='elec1010.lay';
    %     cfg.ylim=[-3 8];
    %     cfg.linewidth=3;
    %     cfg.xlim=timwin;
    %     cfg.channel=sigchannels;
    %     cfg.graphcolor=coloruse([1 10 5],:);
    %     cfg.interactive='no';
    %     cfg.maskparameter='mask';
    %     cfg.maskstyle='box'; % default
    %     subplot(3,1,3)
    %     if timwinstatflag==0
    %       %       tmpmask=mean(statt_mc{ll,tt,ss}.mask(match_str(tmpd5.label,cfg.channel),:),1);
    %       %       tmpd5.mask=zeros(1,length(tmpd5.time));
    %       %       tmpd5.mask(dsearchn(tmpd5.time',stattimwin(1)):dsearchn(tmpd5.time',stattimwin(end)))=tmpmask;
    %       tmpmask=statt_mc{ll,tt,ss}.mask;
    %       tmpd5.mask=zeros(size(tmpd5.avg,1),length(tmpd5.time));
    %       tmpd5.mask(:,dsearchn(tmpd5.time',stattimwin(1)):dsearchn(tmpd5.time',stattimwin(end)))=tmpmask;
    %     end
    %     ft_singleplotER(cfg,tmpu1,tmpm10,tmpd5);
    %     hold on;plot(tmpn2.time,0,'k');
    %     set(gca,'XTick',[cfg.xlim(1):.1:cfg.xlim(2)])
    %     plot([stattimwin(1) stattimwin(1)],cfg.ylim,'k--')
    %     plot([stattimwin(2) stattimwin(2)],cfg.ylim,'k--')
    %     legend('Sum Unisensory','MultSens + Null','SumUnisens - MultsensNull')
    %     print(50+ll,[fdir 'erp_tacPaud_MSpN_diff_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  end
end
for ll=soalist,erpsave{ll}.label=statt_allmc{ll,tt,ss}.label;end
save erpsave.mat erpsave*

load erpsave.mat
load icasave.mat

timeadd=max(0,soades);
topocorrERP=nan(7,9);
courseFCcorrERP=nan(7,9);
courseOPcorrERP=nan(7,9);
for ll=soalist
  for llind=1:7
    if isfield(erpsave{ll},'topo')
      chanuse=match_str(icasaveERP{llind}.label,erpsave{ll}.label);
      topocorrERP(llind,ll)=corr(icasaveERP{llind}.topo(chanuse),erpsave{ll}.topo{1});
    end
    time1=nearest(erpsave{ll}.time,icasaveERP{llind}.time(1)+timeadd(ll));
    time2=nearest(erpsave{ll}.time,icasaveERP{llind}.time(end)+timeadd(ll));
    courseFCcorrERP(llind,ll)=corr(icasaveERP{llind}.course,erpsave{ll}.courseFC(time1:time2)');
    courseOPcorrERP(llind,ll)=corr(icasaveERP{llind}.course,erpsave{ll}.courseOP(time1:time2)');
  end
end
figure;imagescc(topocorrERP)
figure;imagescc(courseFCcorrERP)
figure;imagescc(courseOPcorrERP)

figure;imagesc(abs(topocorrERP));caxis([0 1]);colormap('gray')
figure;imagesc(abs(courseFCcorrERP));caxis([0 1]);colormap('gray')
figure;imagesc(abs(courseOPcorrERP));caxis([0 1]);colormap('gray')

abs(topocorrERP)>.45 & (abs(courseFCcorrERP)>.7 | abs(courseOPcorrERP)>.7)

figure;imagesc((abs(topocorrERP)+max(abs(courseFCcorrERP),abs(courseOPcorrERP)))/2);caxis([0 1]);colormap('gray')
figure;imagesc([(abs(topocorrERP)+max(abs(courseFCcorrERP),abs(courseOPcorrERP)))/2]>.55)

% plotting topo for peak times of interest, across all conditions
for ll=soalist
  cfg=[];
  cfg.operation='subtract';
  cfg.parameter='individual';
  grind_TPA_MSPN{ll,tt,ss}=ft_math(cfg,grind_tacPaud_save{ll,tt,ss},grind_tacMSpN_save{ll,tt,ss});

  cfg=[];
  if timwinstatflag==1
    cfg.latency=[statt_mc{ll,tt,ss}.time(1) statt_mc{ll,tt,ss}.time(end)];
  elseif timwinstatflag==0
    cfg.latency=timwin;
    stattimwin=[statt_mc{ll,tt,ss}.time(1) statt_mc{ll,tt,ss}.time(end)];
  end
  stattimwin=[statt_mc{ll,tt,ss}.time(1) statt_mc{ll,tt,ss}.time(end)];
  cfg.channel=statt_mc{ll,tt,ss}.label;  
  tmpd5=ft_selectdata(cfg,grind_TPA_MSPN{ll,tt,ss});
  tmpd5.dimord='chan_time';
  tmpd5.avg=scalediff*squeeze(mean(tmpd5.individual,1));
  tmpd5=rmfield(tmpd5,'individual');

  cfg=[];
  cfg.parameter='avg';
  cfg.layout='elec1010.lay';
  cfg.zlim=[-3 3];
  cfg.xlim=[stattimwin(1)+.18 stattimwin(1)+.22];
  cfg.highlight          = 'on';
  cfg.highlightchannel   =  chanplot{1};
  cfg.highlightsymbol    = '*';
  cfg.highlightsize      = 10;
  cfg.comment='no';
  figure(120+ll);
  ft_topoplotER(cfg,tmpd5);
  print(120+ll,[fdir 'erp_topoDiff_200ms_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  print(120+ll,[fdir 'erp_topoDiff_200ms_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-painters','-depsc')
  
  cfg.xlim=[stattimwin(1)+.105 stattimwin(1)+.145];
  cfg.zlim=[-3 3];
  cfg.highlightchannel   =  chanplot{3};
  figure(130+ll);
  ft_topoplotER(cfg,tmpd5);
  print(130+ll,[fdir 'erp_topoDiff_125ms_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  print(130+ll,[fdir 'erp_topoDiff_125ms_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-painters','-depsc')

%   cfg.xlim=[stattimwin(1)+.08 stattimwin(1)+.12];
%   cfg.zlim=[-3 3];
%   cfg.highlightchannel   =  chanplot{1};
%   figure(140+ll);
%   ft_topoplotER(cfg,tmpd5);
%   print(140+ll,[fdir 'erp_topoDiff_100ms_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
%   print(140+ll,[fdir 'erp_topoDiff_100ms_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-painters','-depsc')

  cfg.xlim=[stattimwin(1)+.38 stattimwin(1)+.42];
  cfg.zlim=[-3 3];
  cfg.highlightchannel   =  chanplot{4};
  figure(150+ll);
  ft_topoplotER(cfg,tmpd5);
  print(150+ll,[fdir 'erp_topoDiff_400ms_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  print(150+ll,[fdir 'erp_topoDiff_400ms_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-painters','-depsc')
  
end

for ll=soalist
  cfg=[];
  cfg.operation='subtract';
  cfg.parameter='individual';
  grind_TPA_MSPN{ll,tt,ss}=ft_math(cfg,grind_tacPaud_save{ll,tt,ss},grind_tacMSpN_save{ll,tt,ss});
end
chpl{1}=match_str(grind_TPA_MSPN{1,3,10}.label,chanplot{1});
chpl{2}=match_str(grind_TPA_MSPN{1,3,10}.label,chanplot{2});
chpl{3}=match_str(grind_TPA_MSPN{1,3,10}.label,chanplot{3});
chpl{4}=match_str(grind_TPA_MSPN{1,3,10}.label,chanplot{4});
for ll=soalist
  stattimwin=[statt_mc{ll,3,10}.time(1) statt_mc{ll,3,10}.time(end)];
%     erpU125(ll,:)=nanmean(nanmean(grind_tacPaud_save{ll,3,10}.individual(:,chpl{3},dsearchn(grind_tacPaud_save{ll,3,10}.time',stattimwin(1)+.1):dsearchn(grind_tacPaud_save{ll,3,10}.time',stattimwin(1)+.15)),2),3);
%   erpU100(ll,:)=nanmean(nanmean(grind_tacPaud_save{ll,3,10}.individual(:,chpl{1},dsearchn(grind_tacPaud_save{ll,3,10}.time',stattimwin(1)+.08):dsearchn(grind_tacPaud_save{ll,3,10}.time',stattimwin(1)+.12)),2),3);
  erpU125(ll,:)=nanmean(nanmean(grind_tacPaud_save{ll,3,10}.individual(:,chpl{3},dsearchn(grind_tacPaud_save{ll,3,10}.time',stattimwin(1)+.105):dsearchn(grind_tacPaud_save{ll,3,10}.time',stattimwin(1)+.145)),2),3);
  erpU200(ll,:)=nanmean(nanmean(grind_tacPaud_save{ll,3,10}.individual(:,chpl{1},dsearchn(grind_tacPaud_save{ll,3,10}.time',stattimwin(1)+.18):dsearchn(grind_tacPaud_save{ll,3,10}.time',stattimwin(1)+.22)),2),3);
  erpU400(ll,:)=nanmean(nanmean(grind_tacPaud_save{ll,3,10}.individual(:,chpl{4},dsearchn(grind_tacPaud_save{ll,3,10}.time',stattimwin(1)+.38):dsearchn(grind_tacPaud_save{ll,3,10}.time',stattimwin(1)+.42)),2),3);
  
%     erpM125(ll,:)=nanmean(nanmean(grind_tacMSpN_save{ll,3,10}.individual(:,chpl{3},dsearchn(grind_tacMSpN_save{ll,3,10}.time',stattimwin(1)+.1):dsearchn(grind_tacMSpN_save{ll,3,10}.time',stattimwin(1)+.15)),2),3);
%   erpM100(ll,:)=nanmean(nanmean(grind_tacMSpN_save{ll,3,10}.individual(:,chpl{1},dsearchn(grind_tacMSpN_save{ll,3,10}.time',stattimwin(1)+.08):dsearchn(grind_tacMSpN_save{ll,3,10}.time',stattimwin(1)+.12)),2),3);
  erpM125(ll,:)=nanmean(nanmean(grind_tacMSpN_save{ll,3,10}.individual(:,chpl{3},dsearchn(grind_tacMSpN_save{ll,3,10}.time',stattimwin(1)+.105):dsearchn(grind_tacMSpN_save{ll,3,10}.time',stattimwin(1)+.145)),2),3);
  erpM200(ll,:)=nanmean(nanmean(grind_tacMSpN_save{ll,3,10}.individual(:,chpl{1},dsearchn(grind_tacMSpN_save{ll,3,10}.time',stattimwin(1)+.18):dsearchn(grind_tacMSpN_save{ll,3,10}.time',stattimwin(1)+.22)),2),3);
  erpM400(ll,:)=nanmean(nanmean(grind_tacMSpN_save{ll,3,10}.individual(:,chpl{4},dsearchn(grind_tacMSpN_save{ll,3,10}.time',stattimwin(1)+.38):dsearchn(grind_tacMSpN_save{ll,3,10}.time',stattimwin(1)+.42)),2),3);
  
%     erpD125(ll,:)=nanmean(nanmean(grind_TPA_MSPN{ll,3,10}.individual(:,chpl{3},dsearchn(grind_TPA_MSPN{ll,3,10}.time',stattimwin(1)+.1):dsearchn(grind_TPA_MSPN{ll,3,10}.time',stattimwin(1)+.15)),2),3);
%   erpD100(ll,:)=nanmean(nanmean(grind_TPA_MSPN{ll,3,10}.individual(:,chpl{1},dsearchn(grind_TPA_MSPN{ll,3,10}.time',stattimwin(1)+.08):dsearchn(grind_TPA_MSPN{ll,3,10}.time',stattimwin(1)+.12)),2),3);
  erpD125(ll,:)=nanmean(nanmean(grind_TPA_MSPN{ll,3,10}.individual(:,chpl{3},dsearchn(grind_TPA_MSPN{ll,3,10}.time',stattimwin(1)+.105):dsearchn(grind_TPA_MSPN{ll,3,10}.time',stattimwin(1)+.145)),2),3);
  erpD200(ll,:)=nanmean(nanmean(grind_TPA_MSPN{ll,3,10}.individual(:,chpl{1},dsearchn(grind_TPA_MSPN{ll,3,10}.time',stattimwin(1)+.18):dsearchn(grind_TPA_MSPN{ll,3,10}.time',stattimwin(1)+.22)),2),3);
  erpD400(ll,:)=nanmean(nanmean(grind_TPA_MSPN{ll,3,10}.individual(:,chpl{4},dsearchn(grind_TPA_MSPN{ll,3,10}.time',stattimwin(1)+.38):dsearchn(grind_TPA_MSPN{ll,3,10}.time',stattimwin(1)+.42)),2),3);
end


% % Line plots for summary figure
% try close(55); end
% figure(55);
% x=[-.6 -.35 -.1 0 .1 .35 .6];
% [ph] = plot(x, [mean(erpU100(soalist,:),2), mean(erpM100(soalist,:),2), mean(erpD100(soalist,:),2) ]);
% xlim([-.8 .8])
% ylim([-2.6 3.6])
% set(ph(1),'Color',colorblindApT,'Marker','*','MarkerSize',15,'LineWidth',3,'LineStyle',':')
% set(ph(2),'Color',colorblindMpN,'Marker','*','MarkerSize',15,'LineWidth',3,'LineStyle',':')
% set(ph(3),'Color',colorblindD,'Marker','*','MarkerSize',15,'LineWidth',3,'LineStyle',':')
% set(gca,'FontSize',15)
% set(gca,'LineWidth',3)
% set(gca,'Xtick',x)
% set(gca,'XTickLabel',{ '-500' '-70' '-20' '0' '20' '70' '500'})
% xlabel('Multisensory Asynchrony (ms)')
% legend('A+T','MS+N','Difference')
% ylabel('ERP (uV)')
% print(55,[fdir 'erp_condDiff_100ms_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-painters','-depsc')

data=[mean(erpU125(soalist,:),2), mean(erpM125(soalist,:),2), mean(erpD125(soalist,:),2) ];
starind=4;
yminmax=[-.6 3.6];
figind=56;
figname=[fdir 'erp_condDiff_125ms.eps'];
plotlegsummarybars(data,figind,figname,starind,yminmax)

data=[mean(erpU200(soalist,:),2), mean(erpM200(soalist,:),2), mean(erpD200(soalist,:),2) ];
starind=[2 6];
figind=57;
yminmax=[-.6 7];
figname=[fdir 'erp_condDiff_200ms.eps'];
plotlegsummarybars(data,figind,figname,starind,yminmax)

data=[mean(erpU400(soalist,:),2), mean(erpM400(soalist,:),2), mean(erpD400(soalist,:),2) ];
starind=[3 5];
figind=58;
yminmax=[-1.6 0.7];
figname=[fdir 'erp_condDiff_400ms.eps'];
plotlegsummarybars(data,figind,figname,starind,yminmax)

% try close(56); end;
% figure(56);
% x=[-.6 -.35 -.1 0 .1 .35 .6];
% [ph] = plot(x, [mean(erpU125(soalist,:),2), mean(erpM125(soalist,:),2), mean(erpD125(soalist,:),2) ]);
% xlim([-.8 .8])
% ylim([-.6 3.6])
% hold on
% set(ph(1),'Color',colorblindApT,'Marker','^','MarkerSize',15,'LineWidth',3,'LineStyle',':')
% set(ph(2),'Color',colorblindMpN,'Marker','v','MarkerSize',15,'LineWidth',3,'LineStyle',':')
% set(ph(3),'Color',colorblindD,'Marker','o','MarkerSize',15,'LineWidth',3,'LineStyle',':')
% set(ph(3),'Color',colorblindD,'MarkerIndices',4,'Marker','*','MarkerSize',15,'LineWidth',3,'LineStyle',':')
% aa=plot(x(2:end-1), [mean(erpU125(3:7,:),2), mean(erpM125(3:7,:),2), mean(erpD125(3:7,:),2) ]);
% set(aa(1),'Color',colorblindApT,'LineWidth',3,'LineStyle','-')
% set(aa(2),'Color',colorblindMpN,'LineWidth',3,'LineStyle','-')
% set(aa(3),'Color',colorblindD,'LineWidth',3,'LineStyle','-')
% ab=plot(xlim, [0 0],'Color',colorblindD);
% set(gca,'FontSize',15)
% set(gca,'LineWidth',3)
% set(gca,'Xtick',x)
% set(gca,'XTickLabel',{ '-500' '-70' '-20' '0' '20' '70' '500'})
% xlabel('Multisensory Asynchrony (ms)')
% legend('A+T','MS+N','Difference')
% ylabel('ERP (uV)')
% print(56,[fdir 'erp_condDiff_125ms_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-painters','-depsc')
% 
% 
% try close(57); end;
% figure(57);
% x=[-.6 -.35 -.1 0 .1 .35 .6];
% [ph] = plot(x, [mean(erpU200(soalist,:),2), mean(erpM200(soalist,:),2), mean(erpD200(soalist,:),2) ]);
% xlim([-.8 .8])
% ylim([-.6 7])
% set(ph(1),'Color',colorblindApT,'Marker','*','MarkerSize',15,'LineWidth',3,'LineStyle',':')
% set(ph(2),'Color',colorblindMpN,'Marker','*','MarkerSize',15,'LineWidth',3,'LineStyle',':')
% set(ph(3),'Color',colorblindD,'Marker','*','MarkerSize',15,'LineWidth',3,'LineStyle',':')
% set(gca,'FontSize',15)
% set(gca,'LineWidth',3)
% set(gca,'Xtick',x)
% set(gca,'XTickLabel',{ '-500' '-70' '-20' '0' '20' '70' '500'})
% xlabel('Multisensory Asynchrony (ms)')
% legend('A+T','MS+N','Difference')
% ylabel('ERP (uV)')
% print(57,[fdir 'erp_condDiff_200ms_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-painters','-depsc')
% 
% 
% try close(58); end;
% figure(58);
% x=[-.6 -.35 -.1 0 .1 .35 .6];
% [ph] = plot(x, [mean(erpU400(soalist,:),2), mean(erpM400(soalist,:),2), mean(erpD400(soalist,:),2) ]);
% xlim([-.8 .8])
% ylim([-1.6 0.7])
% set(ph(1),'Color',colorblindApT,'Marker','*','MarkerSize',15,'LineWidth',3,'LineStyle',':')
% set(ph(2),'Color',colorblindMpN,'Marker','*','MarkerSize',15,'LineWidth',3,'LineStyle',':')
% set(ph(3),'Color',colorblindD,'Marker','*','MarkerSize',15,'LineWidth',3,'LineStyle',':')
% set(gca,'FontSize',15)
% set(gca,'LineWidth',3)
% set(gca,'Xtick',x)
% set(gca,'XTickLabel',{ '-500' '-70' '-20' '0' '20' '70' '500'})
% xlabel('Multisensory Asynchrony (ms)')
% legend('A+T','MS+N','Difference')
% ylabel('ERP (uV)')
% print(58,[fdir 'erp_condDiff_400ms_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-painters','-depsc')
% 


%% 
% OLD color scheme & plotyy
% Line plots for summary figure
try close(55); end;
figure(55);
x=[-.6 -.35 -.1 0 .1 .35 .6];
[ph, h1, h2] = plotyy(x, mean(erpU100(soalist,:),2), x, mean(erpD100(soalist,:),2));
line(x, mean(erpM100(soalist,:),2), 'Parent', ph(1));
line([-.8 .8], [0 0], 'Parent', ph(2));
xlim(ph(2), [-.8 .8])
set(ph(2).Children(1),'Color',colorblind(7,:))
set(h2,'Color',coloruse([5 ],:),'Marker','*','MarkerSize',15,'LineWidth',3,'LineStyle',':');
set(ph(1).Children(1),'Color',colorblind(6,:),'Marker','^','MarkerSize',15,'LineWidth',3,'LineStyle',':')
set(ph(1).Children(2),'Color',colorblind(3,:),'Marker','v','MarkerSize',15,'LineWidth',3,'LineStyle',':')
set(ph(1),'YColor',colorblind(6,:));
set(ph(2),'YColor',colorblind(7,:));
set(ph,'FontSize',15)
set(ph,'LineWidth',3)
xlim(ph(2), [-.8 .8])
xlim(ph(1), [-.8 .8])
% ylim(ph(1), [0 3.6])
% ylim(ph(2), [-.1 .8])
ylim(ph(1), [-2.6 3.6])
ylim(ph(2), [-2.6 3.6])
set(gca,'Xtick',x)
set(gca,'XTickLabel',{ '-500' '-70' '-20' '0' '20' '70' '500'})
xlabel('Multisensory Asynchrony (ms)')
legend('A+T','MS+N','Difference')
ylabel('ERP (uV)')
print(55,[fdir 'erp_condDiff_100ms_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-painters','-depsc')

try close(56); end;
figure(56);
x=[-.6 -.35 -.1 0 .1 .35 .6];
[ph, h1, h2] = plotyy(x, mean(erpU125(soalist,:),2), x, mean(erpD125(soalist,:),2));
line(x, mean(erpM125(soalist,:),2), 'Parent', ph(1));
line([-.8 .8], [0 0], 'Parent', ph(2));
xlim(ph(2), [-.8 .8])
set(ph(2).Children(1),'Color',colorblind(7,:))
set(h2,'Color',coloruse([5 ],:),'Marker','*','MarkerSize',15,'LineWidth',3,'LineStyle',':');
set(ph(1).Children(1),'Color',colorblind(6,:),'Marker','^','MarkerSize',15,'LineWidth',3,'LineStyle',':')
set(ph(1).Children(2),'Color',colorblind(3,:),'Marker','v','MarkerSize',15,'LineWidth',3,'LineStyle',':')
set(ph(1),'YColor',colorblind(6,:));
set(ph(2),'YColor',colorblind(7,:));
set(ph,'FontSize',15)
set(ph,'LineWidth',3)
xlim(ph(2), [-.8 .8])
xlim(ph(1), [-.8 .8])
% ylim(ph(1), [0 3.6])
% ylim(ph(2), [-.1 .8])
ylim(ph(1), [-.6 3.6])
ylim(ph(2), [-.6 3.6])
set(gca,'Xtick',x)
set(gca,'XTickLabel',{ '-500' '-70' '-20' '0' '20' '70' '500'})
xlabel('Multisensory Asynchrony (ms)')
legend('A+T','MS+N','Difference')
ylabel('ERP (uV)')
print(56,[fdir 'erp_condDiff_125ms_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-painters','-depsc')


try close(57); end;
figure(57);
x=[-.6 -.35 -.1 0 .1 .35 .6];
[ph, h1, h2] = plotyy(x, mean(erpU200(soalist,:),2), x, mean(erpD200(soalist,:),2));
line(x, mean(erpM200(soalist,:),2), 'Parent', ph(1));
line([-.8 .8], [0 0], 'Parent', ph(2));
set(h2,'Color',coloruse([5 ],:),'Marker','*','MarkerSize',15,'LineWidth',3,'LineStyle',':');
set(ph(2).Children(1),'Color',colorblind(7,:))
set(ph(1).Children(1),'Color',colorblind(6,:),'Marker','^','MarkerSize',15,'LineWidth',3,'LineStyle',':')
set(ph(1).Children(2),'Color',colorblind(3,:),'Marker','v','MarkerSize',15,'LineWidth',3,'LineStyle',':')
set(ph(1),'YColor',colorblind(6,:));
set(ph(2),'YColor',colorblind(7,:));
set(ph,'FontSize',15)
set(ph,'LineWidth',3)
xlim(ph(2), [-.8 .8])
xlim(ph(1), [-.8 .8])
% ylim(ph(1), [1.5 7])
% ylim(ph(2), [-.6 1.5])
ylim(ph(1), [-.6 7])
ylim(ph(2), [-.6 7])
set(gca,'Xtick',x)
set(gca,'XTickLabel',{ '-500' '-70' '-20' '0' '20' '70' '500'})
xlabel('Multisensory Asynchrony (ms)')
legend('A+T','MS+N','Difference')
ylabel('ERP (uV)')
print(57,[fdir 'erp_condDiff_200ms_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-painters','-depsc')


try close(58); end;
figure(58);
x=[-.6 -.35 -.1 0 .1 .35 .6];
[ph, h1, h2] = plotyy(x, mean(erpU400(soalist,:),2), x, mean(erpD400(soalist,:),2));
line(x, mean(erpM400(soalist,:),2), 'Parent', ph(1));
line([-.8 .8], [0 0], 'Parent', ph(2));
set(h2,'Color',coloruse([5 ],:),'Marker','*','MarkerSize',15,'LineWidth',3,'LineStyle',':');
set(ph(2).Children(1),'Color',colorblind(7,:))
set(ph(1).Children(1),'Color',colorblind(6,:),'Marker','^','MarkerSize',15,'LineWidth',3,'LineStyle',':')
set(ph(1).Children(2),'Color',colorblind(3,:),'Marker','v','MarkerSize',15,'LineWidth',3,'LineStyle',':')
set(ph(1),'YColor',colorblind(6,:));
set(ph(2),'YColor',colorblind(7,:));
set(ph,'FontSize',15)
set(ph,'LineWidth',3)
xlim(ph(2), [-.8 .8])
xlim(ph(1), [-.8 .8])
% ylim(ph(1), [-1.6 0.7])
% ylim(ph(2), [-1.6 0.1])
ylim(ph(1), [-1.6 0.7])
ylim(ph(2), [-1.6 0.7])
set(gca,'Xtick',x)
set(gca,'XTickLabel',{ '-500' '-70' '-20' '0' '20' '70' '500'})
xlabel('Multisensory Asynchrony (ms)')
legend('A+T','MS+N','Difference')
ylabel('ERP (uV)')
print(58,[fdir 'erp_condDiff_400ms_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-painters','-depsc')













% plotting stat pvalues in front-to-back channel order
cfg=[];cfg.layout='elec1010';layout=ft_prepare_layout(cfg);
for lb=1:length(statt_mc{5,3,10}.label),layoutind(lb)=match_str(layout.label,statt_mc{5,3,10}.label(lb));end
[sortval,sortind]=sort(layoutind);
for ll=soalist
  figure(ll);imagesc(statt_mc{ll,3,10}.time,1:61,statt_mc{ll,3,10}.prob(sortind,:));caxis([0 1]);colorbar
end
for ll=soalist
  figure(ll);
  h=imagesc(statt_mc{ll,3,10}.time,1:61,statt_mc{ll,3,10}.stat(sortind,:));caxis([-4 4]);colorbar
  set(h,'AlphaData',0.5)
  hold on
  f=imagesc(statt_mc{ll,3,10}.time,1:61,statt_mc{ll,3,10}.stat(sortind,:));caxis([-4 4]);colorbar
  set(f,'AlphaData',statt_mc{ll,3,10}.mask(sortind,:))
  hold off
  fighand=get(ll,'Children');
  %   set(fighand(1),'FontSize',22)
  set(gca,'FontSize',25)
  set(gca,'YTick',[10 30 50])
  set(gca,'YTickLabel',{'Front' 'Cent.' 'Post.'})
  set(gca,'YTickLabelRotation',90)
  print(ll,[fdir 'erp_maskedTstat_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-painters','-depsc')
  print(ll,[fdir 'erp_maskedTstat_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
end


% correlation following Cecere et al. 2017 of topo RSA (tRSA)
% first average over all time for one overall map
llcnt=0;lllcnt=0;trsa=zeros(7,7);tcov=zeros(7,7);
for ll=soalist
  llcnt=llcnt+1;
  lllcnt=0;
  for lll=soalist
    lllcnt=lllcnt+1;
    trsa(llcnt,lllcnt)=corr(statt_mc{ll,3,10}.stat(:),statt_mc{lll,3,10}.stat(:));
    tmpcov=cov(statt_mc{ll,3,10}.stat(:),statt_mc{lll,3,10}.stat(:));
    tcov(llcnt,lllcnt)=tmpcov(1,2);
  end
end
figure;imagesc(trsa);caxis([-1 1]);colorbar;
figure;imagesc(tcov);caxis([-max(abs(tcov(:))) max(abs(tcov(:)))]);colorbar;
% % then get spatial correlation for sliding time window averages
% llcnt=0;lllcnt=0;trsat=zeros(7,7,50);
% for ll=soalist
%   llcnt=llcnt+1;
%   lllcnt=0;
%   for lll=soalist
%     lllcnt=lllcnt+1;
%     for tt=1:50 % every 10 ms
%       stat1=mean(statt_mc{ll,3,10}.stat(:,tt*10-9:tt*10+1,:),2);
%       stat2=mean(statt_mc{lll,3,10}.stat(:,tt*10-9:tt*10+1,:),2);
%       trsat(llcnt,lllcnt,tt)=corr(stat1,stat2);
%     end
%   end
% end
% figure;for aa=1:7,for bb=1:7,subplot(7,7,(aa-1)*7+bb);imagesc(squeeze(trsat(aa,bb,:))');caxis([-1 1]);end;end
% every 20 ms
llcnt=0;lllcnt=0;trsat=zeros(7,7,25);covt=zeros(7,7,25);
for ll=soalist
  llcnt=llcnt+1;
  lllcnt=0;
  for lll=soalist
    lllcnt=lllcnt+1;
    for tt=1:25 % every 20 ms
      stat1=mean(statt_mc{ll,3,10}.stat(:,tt*20-19:tt*20+1,:),2);
      stat2=mean(statt_mc{lll,3,10}.stat(:,tt*20-19:tt*20+1,:),2);
      trsat(llcnt,lllcnt,tt)=corr(stat1,stat2,'type','Spearman');
      tmpcov=cov(stat1,stat2);
      covt(llcnt,lllcnt,tt)=tmpcov(1,2);
      ed(llcnt,lllcnt,tt)=norm(stat1-stat2); % Euclidean distance
    end
  end
end
figure;for aa=1:7,for bb=1:7,subplot(7,7,(aa-1)*7+bb);imagesc(squeeze(trsat(aa,bb,:))');caxis([-1 1]);end;end
figure;for aa=1:7,for bb=1:7,subplot(7,7,(aa-1)*7+bb);imagesc(squeeze(covt(aa,bb,:))');caxis([-max(abs(covt(:)))/2 max(abs(covt(:)))/2]);end;end
figure;for aa=1:7,for bb=1:7,subplot(7,7,(aa-1)*7+bb);plot(10:20:490,squeeze(trsat(aa,bb,:))');ylim([-1 1]);end;end
figure;for aa=1:7,for bb=1:7,subplot(7,7,(aa-1)*7+bb);plot(10:20:490,squeeze(covt(aa,bb,:))');ylim([-max(abs(covt(:)))/2 max(abs(covt(:)))/2]);end;end

figure;for tt=1:25,subplot(5,5,tt);imagesc(squeeze(trsat(:,:,tt)));caxis([-1 1]);title(num2str(tt*20-10));end;
figure;for tt=1:25,subplot(5,5,tt);imagesc(squeeze(ed(:,:,tt)));caxis([0 15]);title(num2str(tt*20-10));end;

% % every 50 ms
% llcnt=0;lllcnt=0;trsat=zeros(7,7,10);
% for ll=soalist
%   llcnt=llcnt+1;
%   lllcnt=0;
%   for lll=soalist
%     lllcnt=lllcnt+1;
%     for tt=1:10 % every 50 ms
%       stat1=mean(statt_mc{ll,3,10}.stat(:,tt*50-49:tt*50+1,:),2);
%       stat2=mean(statt_mc{lll,3,10}.stat(:,tt*50-49:tt*50+1,:),2);
%       trsat(llcnt,lllcnt,tt)=corr(stat1,stat2);
%     end
%   end
% end
% figure;for aa=1:7,for bb=1:7,subplot(7,7,(aa-1)*7+bb);imagesc(squeeze(trsat(aa,bb,:))');caxis([-1 1]);end;end

%% obsolete
close all;
% funny 'temporal' (in)congruency contrast
for ll=soalist
  if ll<5
    cfg=[];
    cfg.operation='subtract';
    cfg.parameter='individual';
    grind_diffMSsynch{ll,tt,ss}=ft_math(cfg,grind_tMSsynch_save{ll,tt,ss},grind_tMSasynch_save{ll,tt,ss});
    cfg.channel=statt_synch{ll,tt,ss}.label;
    
    tmpd6=ft_selectdata(cfg,grind_diffMSsynch{ll,tt,ss});
    tmpd6.dimord='chan_time';
    tmpd6.avg=scalediff*squeeze(mean(tmpd6.individual,1));
    tmpd6=rmfield(tmpd6,'individual');
    tmpd6.mask=statt_synch{ll,tt,ss}.mask;
    
    tmps3=ft_selectdata(cfg,grind_tMSsynch_save{ll,tt,ss});
    tmps3.dimord='chan_time';
    tmps3.avg=squeeze(mean(tmps3.individual,1));
    tmps3=rmfield(tmps3,'individual');
    tmps3.mask=statt_synch{ll,tt,ss}.mask;
    
    tmpa8=ft_selectdata(cfg,grind_tMSasynch_save{ll,tt,ss});
    tmpa8.dimord='chan_time';
    tmpa8.avg=squeeze(mean(tmpa8.individual,1));
    tmpa8=rmfield(tmpa8,'individual');
    tmpa8.mask=statt_synch{ll,tt,ss}.mask;
    
    stattimwin=[statt_synch{ll,tt,ss}.time(1) statt_synch{ll,tt,ss}.time(end)];
    stattime=statt_synch{ll,tt,ss}.time;
    
    
    if timwinstatflag==0
      %       tmpmask=mean(statt_mc{ll,tt,ss}.mask(match_str(tmpd5.label,cfg.channel),:),1);
      %       tmpd5.mask=zeros(1,length(tmpd5.time));
      %       tmpd5.mask(dsearchn(tmpd5.time',stattimwin(1)):dsearchn(tmpd5.time',stattimwin(end)))=tmpmask;
      if stattimwin(end)>tmps3.time(end)
        stattimwin(end)=tmps3.time(end);
      end
      tmpmask=statt_synch{ll,tt,ss}.mask;
      tmps3.mask=zeros(size(tmps3.avg,1),length(tmps3.time));
      tmps3.mask(:,dsearchn(tmps3.time',stattimwin(1)):dsearchn(tmps3.time',stattimwin(end)))=tmpmask(:,1:dsearchn(stattime',stattimwin(end)));
      tmpa8.mask=zeros(size(tmpa8.avg,1),length(tmpa8.time));
      tmpa8.mask(:,dsearchn(tmpa8.time',stattimwin(1)):dsearchn(tmpa8.time',stattimwin(end)))=tmpmask(:,1:dsearchn(stattime',stattimwin(end)));
      tmpd6.mask=zeros(size(tmpd6.avg,1),length(tmpd6.time));
      tmpd6.mask(:,dsearchn(tmpd6.time',stattimwin(1)):dsearchn(tmpd6.time',stattimwin(end)))=tmpmask(:,1:dsearchn(stattime',stattimwin(end)));
    end
    
    %     cfg=[];
    %     if timwinstatflag==1
    %       cfg.latency=stattimwin;
    %     elseif timwinstatflag==0
    %       cfg.latency=[-0.5 1];
    %     end
    %     tmp=ft_selectdata(cfg,grave_TMSs_TMSa{ll,tt,ss});
    %     tmp.mask=statt_synch{ll,tt,ss}.mask;
    %     tmps=ft_selectdata(cfg,grind_tMSsynch_save{ll,tt,ss});
    %     tmpa=ft_selectdata(cfg,grind_tMSasynch_save{ll,tt,ss});
    %     tmps.mask=statt_synch{ll,tt,ss}.mask;
    %     tmpa.mask=statt_synch{ll,tt,ss}.mask;
    %
    %     tmps.avg=squeeze(mean(tmps.individual,1));
    %     tmpa.avg=squeeze(mean(tmpa.individual,1));
    %     tmps.var=squeeze(var(tmps.individual,1));
    %     tmpa.var=squeeze(var(tmpa.individual,1));
    %     tmps.dimord='chan_time';
    %     tmpa.dimord='chan_time';
    %     tmps=rmfield(tmps,'individual');
    %     tmpa=rmfield(tmpa,'individual');
    %
    %     cfg=[];
    %     cfg.avgoverchan='yes';
    %     cfg.channel=chanplot{1};
    %     tmp1=ft_selectdata(cfg,tmp);
    %     tmps1=ft_selectdata(cfg,tmps);
    %     tmpa1=ft_selectdata(cfg,tmpa);
    %     if timwinstatflag==0
    %       tmp1mask=mean(tmp.mask(match_str(tmp.label,cfg.channel),:),1);
    %       tmp1.mask=zeros(1,length(tmp1.time));
    %       tmp1.mask(dsearchn(tmp1.time',stattimwin(1)):dsearchn(tmp1.time',stattimwin(end)))=tmp1mask;
    %       tmps1mask=mean(tmps.mask(match_str(tmps.label,cfg.channel),:),1);
    %       tmps1.mask=zeros(1,length(tmps1.time));
    %       tmps1.mask(dsearchn(tmps1.time',stattimwin(1)):dsearchn(tmps1.time',stattimwin(end)))=tmps1mask;
    %       tmpa1mask=mean(tmpa.mask(match_str(tmpa.label,cfg.channel),:),1);
    %       tmpa1.mask=zeros(1,length(tmpa1.time));
    %       tmpa1.mask(dsearchn(tmpa1.time',stattimwin(1)):dsearchn(tmpa1.time',stattimwin(end)))=tmpa1mask;
    %     end
    %     tmp1.mask=logical(ceil(tmp1.mask));
    %     tmps1.mask=logical(ceil(tmps1.mask));
    %     tmpa1.mask=logical(ceil(tmpa1.mask));
    
    %     figure(40+ll);
    %     for cg=1:length(chanplot)+1
    for cg=1:length(chanplot)
      cfg=[];
      cfg.parameter='avg';
      cfg.layout='elec1010.lay';
      cfg.ylim=[-5 14];
      cfg.linewidth=3;
      cfg.xlim=timwin;
      if cg>length(chanplot)
        cfg.channel=tmpd6.label(any(tmpd6.mask,2));
      else
        cfg.channel=chanplot{cg};
      end
      cfg.graphcolor=coloruse([3 8 6],:);
      cfg.interactive='no';
      cfg.maskparameter='mask';
      cfg.maskstyle='box'; % default
      %       subplot(4,1,cg)
      %       subplot(1,2,cg)
      if cg==1
        figure(ll+40);
      elseif cg==2
        figure(ll+45);
      end
      ft_singleplotER(cfg,tmps3,tmpa8,tmpd6);
      hold on;plot(tmps3.time,0,'k');
      set(gca,'XTick',[cfg.xlim(1):.1:cfg.xlim(2)])
      plot([stattimwin(1) stattimwin(1)],cfg.ylim,'k--')
      plot([stattimwin(2) stattimwin(2)],cfg.ylim,'k--')
      plot([0 0],cfg.ylim,'Color',coloruse(4,:))
      plot([soades(10-ll) soades(10-ll)],cfg.ylim,'Color',coloruse(9,:))
      axis([-0.55 1 cfg.ylim(1) cfg.ylim(2)])
      if cg==1
        legend('MultSens Synchronous','MultSens Aynchronous','MSsynch - MSasynch')
      end
    end
    %     print(40+ll,[fdir 'erp_MSsynch_asynch_diff_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    print(ll+40,[fdir 'erp_MSsynch_asynch_diff_FC_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    print(ll+45,[fdir 'erp_MSsynch_asynch_diff_OP_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    
    %   figure(40+ll);
    %     for cg=1:length(chanplot)
    %       cfg=[];
    %       cfg.parameter='avg';
    %       cfg.layout='elec1010.lay';
    %       cfg.ylim=[-5 14];
    %       cfg.linewidth=3;
    %       cfg.xlim=timwin;
    %       cfg.channel=chanplot{cg};
    %       cfg.graphcolor=coloruse([3 8 6],:);
    %       subplot(3,1,cg)
    %       ft_singleplotER(cfg,tmp3,tmp8,tmp6);
    %     end
    %     legend('MultSens Synchronous','MultSens Aynchronous','MSsynch - MSasynch')
    %     print(40+ll,[fdir 'erp_MSsynch_asynch_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    
    
    %     figure(10*ll);
    %     cfg=[];
    %     cfg.parameter='avg';
    %     cfg.layout='elec1010.lay';
    %     cfg.maskparameter='mask';
    %     cfg.maskstyle='box'; % default
    %     cfg.ylim=[-5 14];
    %     cfg.linewidth=3;
    %     ft_singleplotER(cfg,tmps1,tmpa1,tmp1);
    %     print(10*ll,[fdir 'erp_synch_final_FC_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    
    %     cfg=[];
    %     cfg.avgoverchan='yes';
    %     cfg.channel=chanplot{2};
    %     tmp1=ft_selectdata(cfg,tmp);
    %     tmps1=ft_selectdata(cfg,tmps);
    %     tmpa1=ft_selectdata(cfg,tmpa);
    %     if timwinstatflag==0
    %       tmp1mask=mean(tmp.mask(match_str(tmp.label,cfg.channel),:),1);
    %       tmp1.mask=zeros(1,length(tmp1.time));
    %       tmp1.mask(dsearchn(tmp1.time',stattimwin(1)):dsearchn(tmp1.time',stattimwin(end)))=tmp1mask;
    %       tmps1mask=mean(tmps.mask(match_str(tmps.label,cfg.channel),:),1);
    %       tmps1.mask=zeros(1,length(tmps1.time));
    %       tmps1.mask(dsearchn(tmps1.time',stattimwin(1)):dsearchn(tmps1.time',stattimwin(end)))=tmps1mask;
    %       tmpa1mask=mean(tmpa.mask(match_str(tmpa.label,cfg.channel),:),1);
    %       tmpa1.mask=zeros(1,length(tmpa1.time));
    %       tmpa1.mask(dsearchn(tmpa1.time',stattimwin(1)):dsearchn(tmpa1.time',stattimwin(end)))=tmpa1mask;
    %     end
    %     tmp1.mask=logical(ceil(tmp1.mask));
    %     tmps1.mask=logical(ceil(tmps1.mask));
    %     tmpa1.mask=logical(ceil(tmpa1.mask));
    
    %     figure(10*ll+1);
    %     cfg=[];
    %     cfg.parameter='avg';
    %     cfg.layout='elec1010.lay';
    %     cfg.maskparameter='mask';
    %     cfg.maskstyle='box'; % default
    %     cfg.ylim=[-5 3];
    %     cfg.linewidth=3;
    %     ft_singleplotER(cfg,tmps1,tmpa1,tmp1);
    %     print(10*ll+1,[fdir 'erp_synch_final_OP_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    
    %     if [tt==2 && any(ll==[1 4])] || [tt==3 && any(ll==[1 4])]
    if [tt==2 && any(ll==[1 4])] || [tt==3 && any(ll==[1])]
      %       masktime=find(any(tmp.mask,1));
      masktime=find(any(tmpd6.mask,1));
      if find(diff(masktime)>2)
        error('stats gave different answer');
      end
      cfg=[];
      cfg.parameter='avg';
      cfg.layout='elec1010.lay';
      cfg.maskalpha=0.5;
      cfg.zlim=[-6 6];
      cfg.highlight='on';
      cfg.highlightsize=12;
      cfg.xlim=[tmpd6.time(masktime(1)) tmpd6.time(masktime(end))];
      %       cfg.highlightchannel=tmp.label(find(ceil(mean(tmp.mask(:,dsearchn(stattime',cfg.xlim(1)):dsearchn(stattime',cfg.xlim(2))),2))));
      %       cfg.highlightchannel=tmpd6.label(find(ceil(mean(tmpd6.mask(:,dsearchn(stattime',cfg.xlim(1)):dsearchn(stattime',cfg.xlim(2))),2))));
      cfg.highlightchannel=tmpd6.label(find(ceil(mean(tmpd6.mask(:,dsearchn(tmpd6.time',cfg.xlim(1)):dsearchn(tmpd6.time',cfg.xlim(2))),2))));
      figure(50+ll);
      ft_topoplotER(cfg,tmps3);
      print(50+ll,[fdir 'erp_synch_topoS_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
      figure(60+ll);
      ft_topoplotER(cfg,tmpa8);
      print(60+ll,[fdir 'erp_synch_topoA_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
      figure(70+ll);
      ft_topoplotER(cfg,tmpd6);
      print(70+ll,[fdir 'erp_synch_topoDiff_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    end
    if [tt==2 && any(ll==[3])] || [tt==3 && any(ll==[3]) && ss~=11]
      %       masktime=find(any(tmp.mask,1));
      masktime=find(any(tmpd6.mask,1));
      % two different effects going on in these 2 different time windows
      if length(find(diff(masktime)>2))~=1
        error('stats gave something different');
      end
      masktime1=masktime(1:find(diff(masktime)>2));
      masktime2=masktime(find(diff(masktime)>2)+1:end);
      cfg=[];
      cfg.parameter='avg';
      cfg.layout='elec1010.lay';
      cfg.maskalpha=0.5;
      cfg.zlim=[-6 6];
      cfg.highlight='on';
      cfg.highlightsize=12;
      %       cfg.xlim=[tmp.time(masktime(1)) tmp.time(masktime(end))];
      cfg.xlim=[tmpd6.time(masktime1(1)) tmpd6.time(masktime1(end))];
      %       cfg.highlightchannel=tmp.label(find(ceil(mean(tmp.mask(:,dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2))),2))));
      %       cfg.highlightchannel=tmp.label(find(ceil(mean(tmp.mask(:,dsearchn(stattime',cfg.xlim(1)):dsearchn(stattime',cfg.xlim(2))),2))));
      cfg.highlightchannel=tmpd6.label(find(ceil(mean(tmpd6.mask(:,dsearchn(tmpd6.time',cfg.xlim(1)):dsearchn(tmpd6.time',cfg.xlim(2))),2))));
      figure(50+ll);
      ft_topoplotER(cfg,tmps3);
      print(50+ll,[fdir 'erp_synchm1_topoS_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
      figure(60+ll);
      ft_topoplotER(cfg,tmpa8);
      print(60+ll,[fdir 'erp_synchm1_topoA_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
      figure(70+ll);
      ft_topoplotER(cfg,tmpd6);
      print(70+ll,[fdir 'erp_synchm1_topoDiff_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
      cfg.xlim=[tmpd6.time(masktime2(1)) tmpd6.time(masktime2(end))];
      %       cfg.highlightchannel=tmp.label(find(ceil(mean(tmp.mask(:,dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2))),2))));
      cfg.highlightchannel=tmpd6.label(find(ceil(mean(tmpd6.mask(:,dsearchn(tmpd6.time',cfg.xlim(1)):dsearchn(tmpd6.time',cfg.xlim(2))),2))));
      figure(80+ll);
      ft_topoplotER(cfg,tmps3);
      print(80+ll,[fdir 'erp_synchm2_topoS_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
      figure(90+ll);
      ft_topoplotER(cfg,tmpa8);
      print(90+ll,[fdir 'erp_synchm2_topoA_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
      figure(100+ll);
      ft_topoplotER(cfg,tmpd6);
      print(100+ll,[fdir 'erp_synchm2_topoDiff_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    end
  end % ll
  
end

% % Not useful, as it plots one topo for every data sample (thus ever
% 1000hz)
% for ll=soalist
%   cfg=[];
%   cfg.layout='elec1010.lay';
%   ft_clusterplot(cfg,statt_mc{ll,tt,ss});
% end



