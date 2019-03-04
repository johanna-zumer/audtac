eeg_legomagic_preamble

%%  Plotting TF image results


chanplot{1}={'Fz' 'Cz' 'F1' 'F2' 'FC1' 'FC2' 'C1' 'C2'}; % frontocentral
% chanplot{2}={'CP5' 'P5' 'PO3' 'POz' 'PO4' 'P6' 'CP6' 'Pz' 'P3' 'P4'}; % occipital
chanplot{2}={'CP5' 'POz' 'Pz' 'P3' 'P4' 'C4' 'O1' 'O2' 'P7' 'PO7'}; % occipital; % same as for ERP
% chanplot{3}={'FC6' 'FC4' 'F4' 'F6' 'F8' 'FT8' 'T8' 'C6' 'C4'}; % right frontotemporal
chanlabel{1}='Frontocentral electrodes';
chanlabel{2}='Occipital-parietal electrodes';
% chanlabel{3}='Right frontotemporal electrodes';
soalist=[1 3 4 5 6 7 9];
tt=3;
fftaddflag=0;
synchasynch=0;
tacbasemax=min(-.15,soades-.15);
baseline2=[tacbasemax(1) tacbasemax(1)+.08];

% load([edir 'statsgrave_TFR_cond' num2str(ll) num2str(tt) num2str(ss) '.mat']);
% stdir='/home/zumerj/audtac/eeg_data/statsgrave_oldstattimwin/';

itcdepsampflag=1;
%
sleep=1;
if sleep
  ss=12;
  trialkc=-1;
  iter=11;
  usetr=1;
  mcseed=13;
else
  ss=10;
  %   iter=27;
  %   usetr=1;
  iter=31;
  usetr=3;
  mcseed=13;
  trialkc=-1;
end


pre30plot=0; % awkward.... leave for 0 for new iter27.
commentson=0; % sometimes we want comments on, to see xlim & ylim used; for final figures turn off.
plotplvflag=1;
plottfrflag=1;
combval=1;
adda=3;

% for all ll, print at least a TFR even if nothing significant.
for ll=soalist
% for ll=[1]
  close all
  clear grave* stat*
  if pre30plot
    load([edir 'statsgrave_TFR_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '.mat']);
    %   load([stdir 'statsgrave_TFR_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '.mat']);
  else
    if itcdepsampflag
      load([edir 'statsgrave_TFR_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '_usetr' num2str(usetr) '_mcseed' num2str(mcseed) '_itc' '.mat']);
    else
      load([edir 'statsgrave_TFR_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '_usetr' num2str(usetr) '_mcseed' num2str(mcseed) '.mat']);
    end
  end
  
%   for combval=1
    close all
    clear adda
    % power
    if plottfrflag
%       if fftaddflag
%         cfg=[];
%         cfg.latency=[stattl_mc_fftadd.time(1) stattl_mc_fftadd.time(end)];
%         tmp=ft_selectdata(cfg,gravelo_TPA_MSPN_fftadd);
%         tmp.mask=stattl_mc_fftadd.mask;
%         tmps=ft_selectdata(cfg,gravelo_tacMSpN_fftadd);
%         tmps.mask=stattl_mc_fftadd.mask;
%         tmpu=ft_selectdata(cfg,gravelo_tacPaud_fftadd);
%         tmpu.mask=stattl_mc_fftadd.mask;
%       else
%         if combval==1
          cfg=[];
          cfg.latency=[stattl_mc_comb1.time(1) stattl_mc_comb1.time(end)];
          tmp=ft_selectdata(cfg,gravelo_TPA_MSPN_comb1);
          tmp.mask=stattl_mc_comb1.mask;
          tmps=ft_selectdata(cfg,gravelo_tacMSpN_comb1);
          tmps.mask=stattl_mc_comb1.mask;
          tmpu=ft_selectdata(cfg,gravelo_tacPaud_comb1);
          tmpu.mask=stattl_mc_comb1.mask;
          
          cfg=[];
          tmpA=ft_selectdata(cfg,gravelo_TPA_MSPN_comb1);
          tmpA.mask=logical(ones(size(tmpA.powspctrm)));
          tmpA.mask(:,:,dsearchn(tmpA.time',stattl_mc_comb1.time(1)):dsearchn(tmpA.time',stattl_mc_comb1.time(end)))=stattl_mc_comb1.mask;
          tmpsA=ft_selectdata(cfg,gravelo_tacMSpN_comb1);
          tmpsA.mask=logical(ones(size(tmpsA.powspctrm)));
          tmpsA.mask(:,:,dsearchn(tmpsA.time',stattl_mc_comb1.time(1)):dsearchn(tmpsA.time',stattl_mc_comb1.time(end)))=stattl_mc_comb1.mask;
          tmpuA=ft_selectdata(cfg,gravelo_tacPaud_comb1);
          tmpuA.mask=logical(ones(size(tmpuA.powspctrm)));
          tmpuA.mask(:,:,dsearchn(tmpuA.time',stattl_mc_comb1.time(1)):dsearchn(tmpuA.time',stattl_mc_comb1.time(end)))=stattl_mc_comb1.mask;
%         elseif combval==2
%           cfg=[];
%           cfg.latency=[stattl_mc_comb2.time(1) stattl_mc_comb2.time(end)];
%           tmp=ft_selectdata(cfg,gravelo_TPA_MSPN_comb2);
%           tmp.mask=stattl_mc_comb2.mask;
%           tmps=ft_selectdata(cfg,gravelo_tacMSpN_comb2);
%           tmps.mask=stattl_mc_comb2.mask;
%           tmpu=ft_selectdata(cfg,gravelo_tacPaud_comb2);
%           tmpu.mask=stattl_mc_comb2.mask;
%         end
        %         cfg=[];
        %         tmpn=ft_selectdata(cfg,gravelo_tNulAlone_comb );
        %         tmpt=ft_selectdata(cfg,gravelo_tTacAlone_comb );
        %         tmpa=ft_selectdata(cfg,gravelo_tAudAlone_comb );
        %         tmpm=ft_selectdata(cfg,gravelo_tMSAlone_comb );
        
        %       gravelo_tNulAlone_comb.mask=logical(ones(size(gravelo_tNulAlone_comb.powspctrm)));
        %       gravelo_tTacAlone_comb.mask=logical(ones(size(gravelo_tTacAlone_comb.powspctrm)));
        %       gravelo_tAudAlone_comb.mask=logical(ones(size(gravelo_tAudAlone_comb.powspctrm)));
        %       gravelo_tMSAlone_comb.mask= logical(ones(size(gravelo_tMSAlone_comb.powspctrm)));
%       end
      
      %     if 0
      %       cfg=[];
      %       cfg.avgoverchan='yes';
      %       if ~isempty(tmp.label(find(ceil(mean(mean(tmp.mask(:,:,:),2),3)))))
      %         cfg.channel=tmp.label(find(ceil(mean(mean(tmp.mask(:,:,:),2),3))));
      %       else
      %         cfg.channel='all';
      %       end
      %       tmp1=ft_selectdata(cfg,tmp);
      %       tmp1.mask=logical(ceil(tmp1.mask));
      %
      %       figure(10*ll);
      %       cfg=[];
      %       cfg.parameter='powspctrm';
      %       cfg.layout='elec1010.lay';
      %       cfg.maskparameter='mask';
      %       cfg.maskalpha=0.5;
      %       cfg.zlim='maxabs';
      %       ft_singleplotTFR(cfg,tmp1);
      %       print(10*ll,[fdir 'tfrlo_final_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
      %     end
      
      
      % plot TFR of powspctrm, for each suplot: u and m and difference.
      zlim=[-3 3];
      baseline1=[tacbasemax(ll) tacbasemax(ll)+.08];
      
      if 0
        pow{1}=gravelo_tMSAlone_comb;
        pow{2}=tmpu;
        pow{3}=tmps;
        pow{4}=tmp;
        
        chansel=chanplot{1};
        figinds=10*ll;
        figstrings=[];
        figstrings{1}=[fdir 'tfrlo_crop_FC_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        tfr_subchannel_3cond_plot_pow(pow,chansel,zlim,figinds,figstrings,baseline1)
        
        chansel=chanplot{2};
        figinds=10*ll+1;
        figstrings=[];
        figstrings{1}=[fdir 'tfrlo_crop_OP_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        tfr_subchannel_3cond_plot_pow(pow,chansel,zlim,figinds,figstrings,baseline1)
        
        chansel=[chanplot{1} chanplot{2}];
        figinds=10*ll+2;
        figstrings=[];
        figstrings{1}=[fdir 'tfrlo_crop_FCOP_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        tfr_subchannel_3cond_plot_pow(pow,chansel,zlim,figinds,figstrings,baseline1)
      end
      
      %     pow{1}=gravelo_tMSAlone_comb;
      %     pow{2}=tmpuA;
      %     pow{3}=tmpsA;
      %     pow{4}=tmpA;
      pow{1}=tmpuA;
      pow{2}=tmpsA;
      pow{3}=tmpA;
      
      % baseline correct first, then remove baseline area for plotting
      for pp=1:3
        
        cfg=[];
        cfg.baseline=baseline2;
        cfg.baselinetype='absolute';
        cfg.parameter={'powspctrm' 'mask'};
        pow{pp}=ft_freqbaseline(cfg, pow{pp})
        
        cfg=[];
        if ll==1 | ll==3 | ll==4 | ll==5
          cfg.latency=[-0.5 1.3];
        elseif ll==6
          cfg.latency=[-0.5+.02 1.3+.02];
        elseif ll==7
          cfg.latency=[-0.5+.07 1.3+.07];
        elseif ll==9
          cfg.latency=[-0.5+.50 1.3+.50];
        end
        pow{pp}=ft_selectdata(cfg,pow{pp});
      end
      
      %%
      
      chansel=chanplot{1};
      figinds=10*ll;
      figstrings=[];
      figstrings{1}=[fdir 'tfrlo_all_FC_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
      tfr_subchannel_3cond_plot_pow(pow,chansel,zlim,figinds,figstrings,baseline2)
      %     tfr_subchannel_3cond_plotRows_pow(pow,chansel,zlim,figinds,figstrings,baseline2)
      
      chansel=chanplot{2};
      figinds=10*ll+1;
      figstrings=[];
      figstrings{1}=[fdir 'tfrlo_all_OP_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
      tfr_subchannel_3cond_plot_pow(pow,chansel,zlim,figinds,figstrings,baseline2)
      
      chansel=[chanplot{1} chanplot{2}];
      figinds=10*ll+2;
      figstrings=[];
      figstrings{1}=[fdir 'tfrlo_all_FCOP_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
      tfr_subchannel_3cond_plot_pow(pow,chansel,zlim,figinds,figstrings,baseline2)
      
      masktime=find(squeeze(any(mean(tmp.mask(:,1:end,:),2),1)));
      if ~isempty(masktime)
        sigchan=tmp.label(find(ceil(mean(mean(tmp.mask(:,:,dsearchn(tmp.time',tmp.time(masktime(1))):dsearchn(tmp.time',tmp.time(masktime(end)))),2),3))));
        chansel=sigchan;
        figinds=10*ll+9;
        figstrings=[];
        figstrings{1}=[fdir 'tfrlo_all_allsig_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        tfr_subchannel_3cond_plot_pow(pow,chansel,zlim,figinds,figstrings,baseline2)
        
        clear bk fr
        for sc=1:length(sigchan),bk(sc)=~isempty(cell2mat(strfind(sigchan(sc),'P'))) | ~isempty(cell2mat(strfind(sigchan(sc),'O')));end
        for sc=1:length(sigchan),fr(sc)=~isempty(cell2mat(strfind(sigchan(sc),'F'))) | ~isempty(cell2mat(strfind(sigchan(sc),'A')));end
        chansel=sigchan(find(fr));
        figinds=10*ll+9+50;
        figstrings=[];
        figstrings{1}=[fdir 'tfrlo_all_sigFC_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        tfr_subchannel_3cond_plot_pow(pow,chansel,zlim,figinds,figstrings,baseline2)
        
        chansel=sigchan(find(bk));
        figinds=10*ll+9+60;
        figstrings=[];
        figstrings{1}=[fdir 'tfrlo_all_sigOP_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        tfr_subchannel_3cond_plot_pow(pow,chansel,zlim,figinds,figstrings,baseline2)
      end
      
      
      % all channels
      if 0
        chansel=tmp.label;
        figinds=10*ll+8;
        figstrings=[];
        figstrings{1}=[fdir 'tfrlo_all_allchan_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        tfr_subchannel_3cond_plot_pow(pow,chansel,zlim,figinds,figstrings,baseline2)
      end
      
      %     % tried this out but doesn't look right; baselining seems to make sense for display prior to contrast
      %     zlim=[0 10];
      %     baseline='no';
      
      
      % old way
      %     chansel=chanplot{1};
      %     figinds=10*ll;
      %     figstrings=[];
      %     figstrings{1}=[fdir 'tfrlo_final_FC_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.png'];
      %     baseline=[tmp.time(1) tmp.time(9)];
      %     tfr_subchannel_3cond_plot_pow(tmp,tmpu,tmps,chansel,zlim,figinds,figstrings,baseline)
      %
      %     chansel=chanplot{2};
      %     figinds=10*ll+1;
      %     figstrings=[];
      %     figstrings{1}=[fdir 'tfrlo_final_OP_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.png'];
      %     baseline=[tmp.time(1) tmp.time(9)];
      %     tfr_subchannel_3cond_plot_pow(tmp,tmpu,tmps,chansel,zlim,figinds,figstrings,baseline)
      
      
      % theta topo
      skipplot=0;
      if ~isempty(find(squeeze(any(mean(tmp.mask(:,1:2,:),2),1))))
        cfg=[];
        cfg.parameter='powspctrm';
        cfg.layout='elec1010.lay';
        cfg.maskalpha=0.5;
        if sleep==0 && ll==1 && combval==1 && iter==31 && usetr==3
          cfg.ylim=[4 8.5]; % theta here includes 8 hz (nothing in alpha above 8)
        elseif sleep==0 && ll==3 && combval==1 && iter==31 && usetr==3
          cfg.ylim=[4 8.5]; % theta here includes 8 hz (nothing in alpha above 8)
        elseif sleep==0 && ll==1 && combval==1 && iter==31 && usetr==2
          cfg.ylim=[4 10.5]; % theta here includes 8 hz (nothing in alpha above 8)
        elseif sleep==1 && ll==9 && trialkc==-1
          skipplot=1; % combine with alpha
        else
          cfg.ylim=[4 6.5];  % ll==7 & ll==4
        end
        cfg.zlim=[-1.4 1.4];
        cfg.highlight='on';
        if commentson
          cfg.comment='auto';
        else
          cfg.comment='no';
        end
        %       masktime=find(squeeze(any(mean(tmp.mask(:,1:2,:),2),1)));
        %       cfg.xlim=[tmp.time(masktime(1)) tmp.time(masktime(end))];
        masktime_tmp=find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),:),2),1)));
        difftimes=diff(find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),:),2),1))));
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
        if sleep==0 && ll==1 && combval==1 && iter==31 && usetr==3
          masktime_tmp=masktime{1};
          clear masktime
          masktime{1}=masktime_tmp(1:55);
          masktime{2}=masktime_tmp(56:end);
        end
        for dd=1:length(masktime)
          %     if ll==1
          %       cfg.xlim=[.1 .35];
          %     elseif ll==3
          %       cfg.xlim=[.06 .42];
          %     elseif ll==5
          %       cfg.xlim=[-.04 .36];
          %     elseif ll==7
          %       cfg.xlim=[.1 .44];
          %     end
          cfg.xlim=[tmp.time(masktime{dd}(1)) tmp.time(masktime{dd}(end))];
          cfg.highlightchannel=tmp.label(find(ceil(mean(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2)) ),2),3))));
          %           cfg.highlightchannel=tmp.label(find(ceil(mean(mean(tmp.mask(:,1:2,dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2))),2),3))));
          %       cfg.baseline=[tmp1.time(1) tmp1.time(9)];
          cfg.baseline=baseline2;
          cfg.highlightsize=25;
          cfg.highlightsymbol='.';
          figure(10*ll+8);
          ft_topoplotTFR(cfg,tmpuA);
          print(10*ll+8,[fdir 'tfrlo_topoU_theta_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '_time' num2str(dd) '.eps'],'-depsc2')
          figure(10*ll+3);
          ft_topoplotTFR(cfg,tmpsA);
          print(10*ll+3,[fdir 'tfrlo_topoM_theta_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '_time' num2str(dd) '.eps'],'-depsc2')
          cfg.zlim=[-1.4 1.4];
          %     cfg.zlim='maxabs';
          %     cfg.baseline='no';
          figure(10*ll+4);
          ft_topoplotTFR(cfg,tmpA);
          print(10*ll+4,[fdir 'tfrlo_topoDiff_theta_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '_time' num2str(dd) '.eps'],'-depsc2')
        end
        
        
        if 0; %~isempty(setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2})))
          %         chansel=setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2}));
          chansel=cfg.highlightchannel;
          zlim=[-3 3];
          figinds=10*ll+1000;
          figstrings=[];
          figstrings{1}=[fdir 'tfrlo_final_Xtheta_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
          baseline=[tmp.time(1) tmp.time(9)];
          tfr_subchannel_3cond_plot_pow(tmp,tmpu,tmps,chansel,zlim,figinds,figstrings,baseline)
        end
        
      end
      
      % alpha topo
      skipplot=0;
      if ~isempty(find(squeeze(any(mean(tmp.mask(:,3:5,:),2),1))))
        cfg=[];
        if sleep==0 && ll==3 && combval==1  && pre30plot==1 % ??
          cfg.ylim=[10 12];
          %       cfg.xlim=[0 .1];
          %       elseif ll==7
          %         cfg.ylim=[8 10];
          %       cfg.xlim=[.16 .3];
        elseif sleep==0 && ll==1 && combval==1 && pre30plot==1
          cfg.ylim=[8 14];
        elseif sleep==0 && ll==1 && combval==1 && iter==31 && usetr==3
          skipplot=1; % theta goes up to 8Hz
        elseif sleep==0 && ll==3 && combval==1 && iter==31 && usetr==3
          %         cfg.ylim=[8 13];
          skipplot=1; % theta goes up to 8Hz;  beta goes down to 12 Hz
        elseif sleep==0 && ll==4 && combval==1 && iter==31 && usetr==3
          cfg.ylim=[10 18]; % final; merge with beta
        elseif sleep==0 && ll==5 && combval==1 && iter==31 && usetr==3
          cfg.ylim=[10 24]; % combine with beta
        elseif sleep==0 && ll==7 && combval==1 && iter==31 && usetr==3
          cfg.ylim=[8 12];  % and then ignore bits of higher beta
        elseif sleep==0 && ll==4 && combval==1 && iter==32 && usetr==3
          cfg.ylim=[8 12];
        elseif sleep==0 && ll==1 && combval==1 && iter==31 && usetr==2
          skipplot=1; % theta goes up to 10Hz
        elseif sleep==0 && ll==4 && combval==1 && iter==31 && usetr==2
          cfg.ylim=[8 12];
        elseif sleep==0 && ll==7 && combval==1 && iter==31 && usetr==2
          cfg.ylim=[8 12];
        elseif sleep==0 && ll==1 && iter==27 && usetr==1
          %         cfg.ylim=[8 14];
          cfg.ylim=[8 28]; % combine with beta (skip beta below)
        elseif sleep==0 && ll==3 && iter==27 && usetr==1
          %         cfg.ylim=[8 12];
          cfg.ylim=[8 30]; % combine with beta (skip beta below)
        elseif sleep==0 && ll==4 && iter==27 && usetr==1
          %         cfg.ylim=[10 12];
          cfg.ylim=[10 24]; % combine with beta (skip beta below)
        elseif sleep==0 && ll==5 && iter==27 && usetr==1
          %         cfg.ylim=[8 12];
          cfg.ylim=[8 30]; % combine with beta (skip beta below)
        elseif sleep==0 && ll==7 && iter==27 && usetr==1
          skipplot=1; % blur up from theta
        elseif sleep==1 && ll==9 && trialkc==-1
          cfg.ylim=[6 14]; % blur from theta to alpha to low beta, just 1 effect
        elseif sleep==1 && ll==5 && trialkc==0
          skipplot=1; % blur down from beta
        elseif sleep==1 && ll==9 && trialkc==0
          skipplot=1; % blur down from beta
        else
          disp('get ylim right per ll alpha')
          tmp.freq(find(mean(mean(tmp.mask,1),3)))
          keyboard
        end
        if skipplot==0
          if sleep==0 && ll==4 && iter==31
            masktime_tmp=find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',14),:),2),1))); % 14 Hz
            difftimes=diff(find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',14),:),2),1))));
          elseif sleep==0 && ll==5 && iter==31
            masktime_tmp=find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',18),:),2),1)));
            difftimes=diff(find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',18),:),2),1))));
          elseif sleep==0 && ll==7 && iter==31
            masktime_tmp=find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',10),:),2),1)));
            difftimes=diff(find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',10),:),2),1))));
          else
            masktime_tmp=find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),:),2),1)));
            difftimes=diff(find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),:),2),1))));
          end
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
          cfg.parameter='powspctrm';
          cfg.layout='elec1010.lay';
          cfg.maskalpha=0.5;
          %       cfg.zlim=[-1.4 1.4];
          cfg.highlight='on';
          cfg.baseline=baseline2;
          if commentson
            cfg.comment='auto';
          else
            cfg.comment='no';
          end
          for dd=1:length(masktime)
            cfg.xlim=[tmp.time(masktime{dd}(1)) tmp.time(masktime{dd}(end))];
            cfg.highlightchannel=tmp.label(find(ceil(mean(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2)) ),2),3))));
            if sleep==0 && ll==4 && iter==31 % blend with beta
              cfg.zlim=[-4 4];
            elseif sleep==0 && ll==5 && iter==31 % blend with beta
              cfg.zlim=[-4 4];
            else
              cfg.zlim=[-9 9];
            end
            cfg.highlightsize=25;
            cfg.highlightsymbol='.';
            figure(10*ll+5);
            ft_topoplotTFR(cfg,tmpuA);
            print(10*ll+5,[fdir 'tfrlo_topoU_alpha_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '_time' num2str(dd) '.eps'],'-depsc2')
            figure(10*ll+6);
            ft_topoplotTFR(cfg,tmpsA);
            print(10*ll+6,[fdir 'tfrlo_topoM_alpha_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '_time' num2str(dd) '.eps'],'-depsc2')
            if sleep==0 && ll==4 && iter==31 % blend with beta
              cfg.zlim=[-4 4];
            elseif sleep==0 && ll==5 && iter==31 % blend with beta
              cfg.zlim=[-4 4];
            elseif sleep==1
              cfg.zlim=[-3 3];
            else
              cfg.zlim=[-9 9];
            end
            figure(10*ll+7);
            ft_topoplotTFR(cfg,tmpA);
            print(10*ll+7,[fdir 'tfrlo_topoDiff_alpha_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '_time' num2str(dd) '.eps'],'-depsc2')
          end
          
          
          
          if 0 %~isempty(setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2})))
            chansel=setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2}));
            zlim=[-3 3];
            figinds=10*ll+1000+1;
            figstrings=[];
            figstrings{1}=[fdir 'tfrlo_final_Xalpha_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
            baseline=[tmp.time(1) tmp.time(9)];
            tfr_subchannel_3cond_plot_pow(tmp,tmpu,tmps,chansel,zlim,figinds,figstrings,baseline)
          end
        end % skipplot
      end
      
      % beta topo
      skipplot=0;
      if ~isempty(find(squeeze(any(mean(tmp.mask(:,6:10,:),2),1))))
        cfg=[];
        if sleep==0 && ll==5 && pre30plot==1
          cfg.ylim=[14 30];
          %       cfg.xlim=[.2 .28];
        elseif sleep==0 && ll==1 && combval==1 && pre30plot==1
          cfg.ylim=[13 15];
        elseif sleep==0 && ll==3 && combval==1 && pre30plot==1
          cfg.ylim=[14 28];
        elseif sleep==0 && ll==6 && combval==1 && pre30plot==1
          cfg.ylim=[16 28];
        elseif sleep==0 && ll==3 && combval==1 && iter==31 && usetr==3
          cfg.ylim=[14 30]; % final
        elseif sleep==0 && ll==4 && combval==1 && iter==31 && usetr==3
          %         cfg.ylim=[12 30];
          skipplot=1;  % final
        elseif sleep==0 && ll==5 && combval==1 && iter==31 && usetr==3
          skipplot=1; % merge with alpha
        elseif sleep==0 && ll==7 && combval==1 && iter==31 && usetr==3
          %         cfg.ylim=[14 22];
          skipplot=1;  % final
        elseif sleep==0 && ll==4 && combval==1 && iter==32 && usetr==3
          cfg.ylim=[14 22];
        elseif sleep==0 && ll==4 && combval==1 && iter==31 && usetr==2
          cfg.ylim=[14 30];
        elseif sleep==0 && ll==7 && combval==1 && iter==31 && usetr==2
          cfg.ylim=[14 18];
        elseif sleep==0 && ll==1 && iter==27 && usetr==1
          %         cfg.ylim=[16 28];
          skipplot=1; % combine with alpha above
        elseif sleep==0 && ll==3 && iter==27 && usetr==1
          %         cfg.ylim=[14 30];
          skipplot=1; % combine with alpha above
        elseif sleep==0 && ll==4 && iter==27 && usetr==1
          %         cfg.ylim=[14 24];
          skipplot=1; % combine with alpha above
        elseif sleep==0 && ll==5 && iter==27 && usetr==1
          %         cfg.ylim=[14 30];
          skipplot=1; % combine with alpha above
        elseif sleep==0 && ll==6 && iter==27 && usetr==1
          cfg.ylim=[14 30];
        elseif sleep==1 && ll==9 && trialkc==-1
          skipplot=1; % combine with alpha above
        elseif sleep==1 && ll==5 && trialkc==0
          cfg.ylim=[10 18];
        elseif sleep==1 && ll==9 && trialkc==0
          cfg.ylim=[10 18];
        else
          disp('get ylim right per ll beta')
          tmp.freq(find(mean(mean(tmp.mask,1),3)))
          keyboard
        end
        if skipplot==0
          %         if sleep==0 && ll==4 && iter==31
          %           masktime_tmp=find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),:),2),1)));
          %           masktime_tmp=masktime_tmp(26:end); % 1:25 is covered by alpha plot
          %           difftimes=diff(find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),:),2),1))));
          %         else
          masktime_tmp=find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),:),2),1)));
          difftimes=diff(find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),:),2),1))));
          %         end
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
          cfg.parameter='powspctrm';
          cfg.layout='elec1010.lay';
          cfg.maskalpha=0.5;
          %       cfg.zlim=[-1.4 1.4];
          cfg.highlight='on';
          cfg.baseline=baseline2;
          if commentson
            cfg.comment='auto';
          else
            cfg.comment='no';
          end
          for dd=1:length(masktime)
            cfg.xlim=[tmp.time(masktime{dd}(1)) tmp.time(masktime{dd}(end))];
            cfg.highlightchannel=tmp.label(find(ceil(mean(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2)) ),2),3))));
            cfg.zlim=[-1.1 1.1];
            cfg.highlightsize=25;
            cfg.highlightsymbol='.';
            figure(10*ll+5);
            ft_topoplotTFR(cfg,tmpuA);
            print(10*ll+5,[fdir 'tfrlo_topoU_beta_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '_time' num2str(dd) '.eps'],'-depsc2')
            figure(10*ll+6);
            ft_topoplotTFR(cfg,tmpsA);
            print(10*ll+6,[fdir 'tfrlo_topoM_beta_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '_time' num2str(dd) '.eps'],'-depsc2')
            cfg.zlim=[-1.1 1.1];
            figure(10*ll+7);
            ft_topoplotTFR(cfg,tmpA);
            print(10*ll+7,[fdir 'tfrlo_topoDiff_beta_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '_time' num2str(dd) '.eps'],'-depsc2')
          end
          
          
          if 0 %~isempty(setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2})))
            chansel=setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2}));
            zlim=[-2 2];
            figinds=10*ll+1000+2;
            figstrings=[];
            figstrings{1}=[fdir 'tfrlo_final_Xbeta_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
            baseline=[tmp.time(1) tmp.time(9)];
            tfr_subchannel_3cond_plot_pow(tmp,tmpu,tmps,chansel,zlim,figinds,figstrings,baseline)
          end
        end % skipplot
      end
      
    end % plottfrflag

%%        Useful, across all asynch regardless of significance

figure;imagesc(squeeze(mean(stattl_mc_comb1.stat(match_str(stattl_mc_comb1.label,chanplot{1}),:,:),1)));axis xy;caxis([-2.5 2.5])
figure;imagesc(squeeze(mean(stattl_mc_comb1plvabsDepT.stat(match_str(stattl_mc_comb1plvabsDepT.label,chanplot{1}),:,:),1)));axis xy;caxis([-2.5 2.5])
cfg=[];
cfg.ylim=[4 6.5];
cfg.parameter='stat';
cfg.layout='elec1010.lay';
cfg.maskalpha=0.5000;
cfg.highlight='on';
% cfg.baseline=[-0.6500 -0.5700];
cfg.comment='no';
cfg.xlim=[0.15 0.3];
cfg.zlim=[-3.5 3.5];
figure;
ft_topoplotTFR(cfg,stattl_mc_comb1)
cfg.zlim=[-1.5 1.5];
figure
ft_topoplotTFR(cfg,stattl_mc_comb1plvabsDepT)
% % High
% cfg=[];
% cfg.ylim=[30 70];
% cfg.parameter='stat';
% cfg.layout='elec1010.lay';
% cfg.maskalpha=0.5000;
% cfg.highlight='on';
% % cfg.baseline=[-0.6500 -0.5700];
% cfg.comment='no';
% cfg.xlim=[0.7400 1.3000];
% cfg.zlim=[-3.5 3.5];
% ft_topoplotTFR(cfg,statth_mc_comb1)
% cfg.zlim=[-1.5 1.5];
% ft_topoplotTFR(cfg,statth_mc_comb1plvabsDepT)

%%
    
    if plotplvflag
      %    %  % %%%%%%% PLV   %%%%%%%%
%       for adda=3
%         if combval=1
%           if adda==1
%             cfg=[];
%             cfg.latency=[stattl_mc_comb1plvadi.time(1) stattl_mc_comb1plvadi.time(end)];
%             tmp=ft_selectdata(cfg,gravelo_TPA_MSPN_comb1);
%             tmp.mask=stattl_mc_comb1plvadi.mask;
%             tmps=ft_selectdata(cfg,gravelo_tacMSpN_comb1);
%             tmps.mask=stattl_mc_comb1plvadi.mask;
%             tmpu=ft_selectdata(cfg,gravelo_tacPaud_comb1);
%             tmpu.mask=stattl_mc_comb1plvadi.mask;
%             
%             cfg=[];
%             tmpA=ft_selectdata(cfg,gravelo_TPA_MSPN_comb1);
%             tmpA.mask=logical(ones(size(tmpA.powspctrm)));
%             tmpA.mask(:,:,dsearchn(tmpA.time',stattl_mc_comb1plvadi.time(1)):dsearchn(tmpA.time',stattl_mc_comb1plvadi.time(end)))=stattl_mc_comb1plvadi.mask;
%             tmpsA=ft_selectdata(cfg,gravelo_tacMSpN_comb1);
%             tmpsA.mask=logical(ones(size(tmpsA.powspctrm)));
%             tmpsA.mask(:,:,dsearchn(tmpsA.time',stattl_mc_comb1plvadi.time(1)):dsearchn(tmpsA.time',stattl_mc_comb1plvadi.time(end)))=stattl_mc_comb1plvadi.mask;
%             tmpuA=ft_selectdata(cfg,gravelo_tacPaud_comb1);
%             tmpuA.mask=logical(ones(size(tmpuA.powspctrm)));
%             tmpuA.mask(:,:,dsearchn(tmpuA.time',stattl_mc_comb1plvadi.time(1)):dsearchn(tmpuA.time',stattl_mc_comb1plvadi.time(end)))=stattl_mc_comb1plvadi.mask;
%           elseif adda==2
%             cfg=[];
%             cfg.latency=[stattl_mc_comb1plvdai.time(1) stattl_mc_comb1plvdai.time(end)];
%             tmp=ft_selectdata(cfg,gravelo_TPA_MSPN_comb1);
%             tmp.mask=stattl_mc_comb1plvdai.mask;
%             tmps=ft_selectdata(cfg,gravelo_tacMSpN_comb1);
%             tmps.mask=stattl_mc_comb1plvdai.mask;
%             tmpu=ft_selectdata(cfg,gravelo_tacPaud_comb1);
%             tmpu.mask=stattl_mc_comb1plvdai.mask;
%             
%             cfg=[];
%             tmpA=ft_selectdata(cfg,gravelo_TPA_MSPN_comb1);
%             tmpA.mask=logical(ones(size(tmpA.powspctrm)));
%             tmpA.mask(:,:,dsearchn(tmpA.time',stattl_mc_comb1plvdai.time(1)):dsearchn(tmpA.time',stattl_mc_comb1plvdai.time(end)))=stattl_mc_comb1plvdai.mask;
%             tmpsA=ft_selectdata(cfg,gravelo_tacMSpN_comb1);
%             tmpsA.mask=logical(ones(size(tmpsA.powspctrm)));
%             tmpsA.mask(:,:,dsearchn(tmpsA.time',stattl_mc_comb1plvdai.time(1)):dsearchn(tmpsA.time',stattl_mc_comb1plvdai.time(end)))=stattl_mc_comb1plvdai.mask;
%             tmpuA=ft_selectdata(cfg,gravelo_tacPaud_comb1);
%             tmpuA.mask=logical(ones(size(tmpuA.powspctrm)));
%             tmpuA.mask(:,:,dsearchn(tmpuA.time',stattl_mc_comb1plvdai.time(1)):dsearchn(tmpuA.time',stattl_mc_comb1plvdai.time(end)))=stattl_mc_comb1plvdai.mask;
%           elseif adda==3
            cfg=[];
            cfg.latency=[stattl_mc_comb1plvabsDepT.time(1) stattl_mc_comb1plvabsDepT.time(end)];
            tmp=ft_selectdata(cfg,gravelo_TPA_MSPN_comb1);
            tmp.mask=stattl_mc_comb1plvabsDepT.mask;
            tmps=ft_selectdata(cfg,gravelo_tacMSpN_comb1);
            tmps.mask=stattl_mc_comb1plvabsDepT.mask;
            tmpu=ft_selectdata(cfg,gravelo_tacPaud_comb1);
            tmpu.mask=stattl_mc_comb1plvabsDepT.mask;
            
            cfg=[];
            tmpA=ft_selectdata(cfg,gravelo_TPA_MSPN_comb1);
            tmpA.mask=logical(ones(size(tmpA.powspctrm)));
            tmpA.mask(:,:,dsearchn(tmpA.time',stattl_mc_comb1plvabsDepT.time(1)):dsearchn(tmpA.time',stattl_mc_comb1plvabsDepT.time(end)))=stattl_mc_comb1plvabsDepT.mask;
            tmpsA=ft_selectdata(cfg,gravelo_tacMSpN_comb1);
            tmpsA.mask=logical(ones(size(tmpsA.powspctrm)));
            tmpsA.mask(:,:,dsearchn(tmpsA.time',stattl_mc_comb1plvabsDepT.time(1)):dsearchn(tmpsA.time',stattl_mc_comb1plvabsDepT.time(end)))=stattl_mc_comb1plvabsDepT.mask;
            tmpuA=ft_selectdata(cfg,gravelo_tacPaud_comb1);
            tmpuA.mask=logical(ones(size(tmpuA.powspctrm)));
            tmpuA.mask(:,:,dsearchn(tmpuA.time',stattl_mc_comb1plvabsDepT.time(1)):dsearchn(tmpuA.time',stattl_mc_comb1plvabsDepT.time(end)))=stattl_mc_comb1plvabsDepT.mask;
%           end
%         elseif combval==2
%           if adda==1
%             cfg=[];
%             cfg.latency=[stattl_mc_comb2plvadi.time(1) stattl_mc_comb2plvadi.time(end)];
%             tmp=ft_selectdata(cfg,gravelo_TPA_MSPN_comb2);
%             tmp.mask=stattl_mc_comb2plvadi.mask;
%             tmps=ft_selectdata(cfg,gravelo_tacMSpN_comb2);
%             tmps.mask=stattl_mc_comb2plvadi.mask;
%             tmpu=ft_selectdata(cfg,gravelo_tacPaud_comb2);
%             tmpu.mask=stattl_mc_comb2plvadi.mask;
%           elseif adda==2
%             cfg=[];
%             cfg.latency=[stattl_mc_comb2plvdai.time(1) stattl_mc_comb2plvdai.time(end)];
%             tmp=ft_selectdata(cfg,gravelo_TPA_MSPN_comb2);
%             tmp.mask=stattl_mc_comb2plvdai.mask;
%             tmps=ft_selectdata(cfg,gravelo_tacMSpN_comb2);
%             tmps.mask=stattl_mc_comb2plvdai.mask;
%             tmpu=ft_selectdata(cfg,gravelo_tacPaud_comb2);
%             tmpu.mask=stattl_mc_comb2plvdai.mask;
%           elseif adda==3
%             cfg=[];
%             cfg.latency=[stattl_mc_comb2plvabsDepT.time(1) stattl_mc_comb2plvabsDepT.time(end)];
%             tmp=ft_selectdata(cfg,gravelo_TPA_MSPN_comb2);
%             tmp.mask=stattl_mc_comb2plvabsDepT.mask;
%             tmps=ft_selectdata(cfg,gravelo_tacMSpN_comb2);
%             tmps.mask=stattl_mc_comb2plvabsDepT.mask;
%             tmpu=ft_selectdata(cfg,gravelo_tacPaud_comb2);
%             tmpu.mask=stattl_mc_comb2plvabsDepT.mask;
%           end
%         end
        
        tmp.plvavgangwrap=wrapToPi(tmp.plvavgang);
        tmps.plvavgangwrap=wrapToPi(tmps.plvavgang);
        tmpu.plvavgangwrap=wrapToPi(tmpu.plvavgang);
        tmpA.plvavgangwrap=wrapToPi(tmpA.plvavgang);
        tmpsA.plvavgangwrap=wrapToPi(tmpsA.plvavgang);
        tmpuA.plvavgangwrap=wrapToPi(tmpuA.plvavgang);
        %       gravelo_tNulAlone_comb.plvavgangwrap=wrapToPi(gravelo_tNulAlone_comb.plvavgang);
        %       gravelo_tTacAlone_comb.plvavgangwrap=wrapToPi(gravelo_tTacAlone_comb.plvavgang);
        %       gravelo_tAudAlone_comb.plvavgangwrap=wrapToPi(gravelo_tAudAlone_comb.plvavgang);
        %       gravelo_tMSAlone_comb.plvavgangwrap= wrapToPi(gravelo_tMSAlone_comb.plvavgang);
        
        %       if 0
        %         cfg=[];
        %         cfg.avgoverchan='yes';
        %         if ~isempty(tmp.label(find(ceil(mean(mean(tmp.mask(:,:,:),2),3)))))
        %           cfg.channel=tmp.label(find(ceil(mean(mean(tmp.mask(:,:,:),2),3))));
        %         else
        %           cfg.channel='all';
        %         end
        %         tmp1=ft_selectdata(cfg,tmp);
        %         tmp1.mask=logical(ceil(tmp1.mask));
        %
        %         figure(100*ll);
        %         cfg=[];
        %         cfg.parameter='plvavgabs';
        %         cfg.layout='elec1010.lay';
        %         cfg.maskparameter='mask';
        %         cfg.maskalpha=0.5;
        %         cfg.zlim='maxabs';
        %         ft_singleplotTFR(cfg,tmp1);
        %         print(100*ll,[fdir 'plvabslo_final_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
        %         figure(100*ll+10);
        %         cfg.parameter='plvavgang';
        %         ft_singleplotTFR(cfg,tmp1);
        %         print(100*ll+10,[fdir 'plvanglo_final_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
        %       end
        
        
        % plot TFR of plv abs and plv ang, for each suplot: u and m and difference.
        if pre30plot
          zlim=[0 0.6; -.25 .25; -20 400; -10 300];
        else
          if adda==3
            if sleep  % both trialkc =-1 and =0
              zlim=[0 0.34; -0.16 0.16;   -20 400; -10 300];
            else
              zlim=[0 0.61; -0.16 0.16;   -20 400; -10 300];
            end
          else
            zlim=[0 0.5; .1 .25;   -20 400; -10 300];
          end
        end
        
        if 0
          pow{1}=gravelo_tMSAlone_comb;
          pow{2}=tmpu;
          pow{3}=tmps;
          pow{4}=tmp;
          
          chansel=chanplot{1};
          figinds=[100*ll;  100*ll+10];
          figstrings{1}=[fdir 'plvabslo_crop_FC_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
          figstrings{2}=[fdir 'plvanglo_crop_FC_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
          tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
          
          chansel=chanplot{2};
          figinds=[100*ll+1;  100*ll+10+1];
          figstrings{1}=[fdir 'plvabslo_crop_OP_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
          figstrings{2}=[fdir 'plvanglo_crop_OP_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
          tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
          
          chansel=[chanplot{1} chanplot{2}];
          figinds=[100*ll+2;  100*ll+10+2];
          figstrings{1}=[fdir 'plvabslo_crop_FCOP_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
          figstrings{2}=[fdir 'plvanglo_crop_FCOP_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
          tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
        end
        
        clear pow
        if pre30plot
          pow{1}=gravelo_tMSAlone_comb;
          pow{2}=tmpuA;
          pow{3}=tmpsA;
          pow{4}=tmpA;
        else
          pow{1}=tmpuA;
          pow{2}=tmpsA;
          pow{3}=tmpA;
        end
        
        for pp=1:length(pow)
          
          %         cfg=[];
          %         cfg.baseline=baseline2;
          %         cfg.baselinetype='absolute';
          %         cfg.parameter={'powspctrm' 'mask'};
          %         pow{pp}=ft_freqbaseline(cfg, pow{pp})
          
          cfg=[];
          if ll==1 | ll==3 | ll==4 | ll==5
            cfg.latency=[-0.5 1.3];
          elseif ll==6
            cfg.latency=[-0.5+.02 1.3+.02];
          elseif ll==7
            cfg.latency=[-0.5+.07 1.3+.07];
          elseif ll==9
            cfg.latency=[-0.5+.50 1.3+.50];
          end
          pow{pp}=ft_selectdata(cfg,pow{pp});
        end
        
        %%
        
        %       disp('starting PLF plotting')
        %       keyboard
        chansel=chanplot{1};
        figinds=[100*ll;  100*ll+10];
        figstrings{1}=[fdir 'plvabslo_all_FC_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        figstrings{2}=[fdir 'plvanglo_all_FC_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings,adda)
        
        chansel=chanplot{2};
        figinds=[100*ll+1;  100*ll+10+1];
        figstrings{1}=[fdir 'plvabslo_all_OP_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        figstrings{2}=[fdir 'plvanglo_all_OP_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings,adda)
        
        chansel=[chanplot{1} chanplot{2}];
        figinds=[100*ll+2;  100*ll+10+2];
        figstrings{1}=[fdir 'plvabslo_all_FCOP_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        figstrings{2}=[fdir 'plvanglo_all_FCOP_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings,adda)
        
        masktime=find(squeeze(any(mean(tmp.mask(:,1:end,:),2),1)));
        if ~isempty(masktime)
          chansel=tmp.label(find(ceil(mean(mean(tmp.mask(:,:,dsearchn(tmp.time',tmp.time(masktime(1))):dsearchn(tmp.time',tmp.time(masktime(end)))),2),3))));
          figinds=[100*ll+9;  100*ll+10+9];
          figstrings{1}=[fdir 'plvabslo_all_allsig_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
          figstrings{2}=[fdir 'plvanglo_all_allsig_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
          tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings,adda)
        end
        
        % theta topo
        skipplot=0;
        if ~isempty(find(squeeze(any(mean(tmp.mask(:,1:2,:),2),1))))
          cfg=[];
          if sleep==0 && ll==3 && pre30plot==1
            cfg.ylim=[4 6.5];
          elseif sleep==0 && ll==6 && adda==1 && pre30plot==1
            cfg.ylim=[5 6.5];
          elseif sleep==0 && ll==7 && pre30plot==1
            cfg.ylim=[4 6.5];
          elseif sleep==0 && ll==3 && iter==31 && usetr==3 && adda==2
            cfg.ylim=[4 8.5];  %
          elseif sleep==0 && ll==3 && iter==31 && usetr==3 && adda==3
            cfg.ylim=[4 6.5];  % final
          elseif sleep==0 && ll==7 && iter==31 && usetr==3
            cfg.ylim=[4 8.5]; % final
          elseif sleep==0 && ll==7 && combval==1 && adda==2 && iter==32 && usetr==3
            cfg.ylim=[5.5 8.5];
          elseif sleep==0 && ll==3 && combval==1 && adda==2 && iter==32 && usetr==3
            cfg.ylim=[4 6.5];
          elseif sleep==0 && ll==3 && iter==27 && usetr==1
            cfg.ylim=[4 8.5];
          elseif sleep==0 && ll==7 && iter==27 && usetr==1
            cfg.ylim=[4 8.5];
          elseif sleep==1 && ll==9 && trialkc==-1
            cfg.ylim=[4 4.5];
            cfg.baseline=[-0.35 -0.27];  % NaN in original baseline window
          elseif sleep==1 && ll==5 && trialkc==0
            cfg.ylim=[4 8.5];
          elseif sleep==1 && ll==4 && trialkc==-1
            cfg.ylim=[4 6.5];
          elseif sleep==1 && ll==5 && trialkc==-1
            cfg.ylim=[4 6.5];
          else
            disp('get ylim right per ll theta plv')
            tmp.freq(find(mean(mean(tmp.mask,1),3)))
            keyboard
          end
          masktime=find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),:),2),1)));
          cfg.xlim=[tmp.time(masktime(1)) tmp.time(masktime(end))];
          %%
          if pptdemo
            cfg=[];
            cfg.ylim=[4 6.5];
            cfg.xlim=[.24 .49]+soades(ll);  % AT20 .24-.64;  AT0 0-.49;
          end
          if adda==2
            cfg.parameter='plvavgabs';
          else adda==3
            cfg.parameter='plvabs';
          end
          cfg.layout='elec1010.lay';
          cfg.maskalpha=0.5;
          cfg.zlim=[-.25 .25];
          cfg.highlight='on';
          if pptdemo
            load thetaITCchan.mat
            cfg.highlightchannel=thetaITCchan;
          else
            cfg.highlightchannel=tmp.label(find(ceil(mean(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2))),2),3))));
          end
          cfg.baseline=baseline2;
          if commentson
            cfg.comment='auto';
          else
            cfg.comment='no';
          end
          cfg.highlightsize=25;
          cfg.highlightsymbol='.';
          cfg.zlim=[-.25 .25];
          figure(100*ll+8);
          ft_topoplotTFR(cfg,tmpuA);
          print(100*ll+8,[fdir 'plvabslo_topoU_theta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
          figure(100*ll+3);
          ft_topoplotTFR(cfg,tmpsA);
          print(100*ll+3,[fdir 'plvabslo_topoM_theta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
          if sleep
            cfg.zlim=[-.06 .06];
          else
            cfg.zlim=[-.15 .15];
          end
          figure(100*ll+4);
          ft_topoplotTFR(cfg,tmpA);
          print(100*ll+4,[fdir 'plvabslo_topoDiff_theta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
          
          if ~isempty(setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2}))) || pptdemo
            if pptdemo
              load thetaITCchan.mat
              chansel=thetaITCchan;
            else
              chansel=cfg.highlightchannel;
            end
%             zlim=[-.3 .3; 0 0.35; -30 30; -60 60];
            zlim=[0 0.2; -.06 .06; -20 400; -10 300];
            figinds=[100*ll+20+1;  100*ll+30+1];
            figstrings{1}=[fdir 'plvabslo_all_Xtheta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
            figstrings{2}=[fdir 'plvanglo_all_Xtheta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
            tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings,adda)
          end
          
          %%
          cfg.parameter='plvavgangwrap';
          if sleep==0 && ll==3 && pre30plot==1
            cfg.xlim=[0.19 0.19];
          elseif sleep==0 && ll==6 && combval==1 && adda==1 && pre30plot==1
            cfg.xlim=[0.41 0.41];
          elseif sleep==0 && ll==6 && combval==2 && adda==1 && pre30plot==1
            cfg.xlim=[0.32 0.32];
          elseif sleep==0 && ll==7 && combval==1 && pre30plot==1
            cfg.xlim=[.21 .21];
          elseif sleep==0 && ll==7 && combval==2 && pre30plot==1
            cfg.xlim=[.23 .23];
          elseif sleep==0 && ll==3 && iter==31 && usetr==3
            cfg.xlim=[.11 .11]; % final
          elseif sleep==0 && ll==7 && iter==31 && usetr==3
            cfg.xlim=[.2 .2]; % final
          elseif sleep==0 && ll==3 && iter==32 && usetr==3
            cfg.xlim=[.03 .03];
          elseif sleep==0 && ll==7 && iter==32 && usetr==3
            cfg.xlim=[.11 .11];
          elseif sleep==0 && ll==3 && iter==27 && usetr==1
            cfg.xlim=[.02 .02];
          elseif sleep==0 && ll==7 && iter==27 && usetr==1
            cfg.xlim=[.11 .11];
            %         elseif sleep==1 && ll==9 && trialkc==-1
            %           cfg.xlim=[0.5 1.11];
            %         elseif sleep==1 && ll==5 && trialkc==0
            %           cfg.xlim=[0.01 0.62];
          else
            disp('get xlim right per ll theta plv')
            tmp.time(find(mean(mean(tmp.mask,1),2)))
            keyboard
          end
          cfg.zlim=[-4 4];
          cfg.highlightsize=25;
          cfg.highlightsymbol='.';
          figure(100*ll+8+10);
          ft_topoplotTFR(cfg,tmpuA);
          print(100*ll+8+10,[fdir 'plvanglo_topoU_theta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
          figure(100*ll+3+10);
          ft_topoplotTFR(cfg,tmpsA);
          print(100*ll+3+10,[fdir 'plvanglo_topoM_theta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
          figure(100*ll+4+10);
          ft_topoplotTFR(cfg,tmpA);
          print(100*ll+4+10,[fdir 'plvanglo_topoDiff_theta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
        end
        
        % alpha topo
        skipplot=0;
        if ~isempty(find(squeeze(any(mean(tmp.mask(:,3:5,:),2),1))))
          cfg=[];
          if sleep==0 && ll==6 && [adda==1 || (combval==2 && adda==2)] && pre30plot==1
            cfg.ylim=[8 11];
          elseif sleep==0 && ll==6 && combval==1 && adda==2 && pre30plot==1
            cfg.ylim=[8 14]; % slightly into beta but then no beta
          elseif sleep==0 && ll==3 && combval==1 && adda==1 && pre30plot==1
            cfg.ylim=[8 12];
          elseif sleep==0 && ll==3 && adda==2 && pre30plot==1
            cfg.ylim=[9 12];
          elseif sleep==0 && ll==3 && combval==2 && adda==1 && pre30plot==1
            cfg.ylim=[8 9];
          elseif sleep==0 && ll==3 && combval==1 && adda==2 && iter==31 && usetr==3
            skipplot=1;
          elseif sleep==0 && ll==5 && combval==1 && adda==2 && iter==31 && usetr==3
            cfg.ylim=[8 14];
          elseif sleep==0 && ll==6 && combval==1 && adda==2 && iter==31 && usetr==3
            cfg.ylim=[8 15]; % 2 clusters: early beta and later alpha
          elseif sleep==0 && ll==7 && combval==1 && adda==2 && iter==31 && usetr==3
            skipplot=1;
          elseif sleep==0 && ll==7 && combval==1 && adda==3 && iter==31 && usetr==3
            skipplot=1;
          elseif sleep==0 && ll==3 && combval==1 && adda==2 && iter==32 && usetr==3
            cfg.ylim=[8 11];
          elseif sleep==0 && ll==5 && combval==1 && adda==2 && iter==32 && usetr==3
            cfg.ylim=[8 14];
          elseif sleep==0 && ll==6 && combval==1 && adda==2 && iter==32 && usetr==3
            cfg.ylim=[8 14];
          elseif sleep==0 && ll==7 && combval==1 && adda==2 && iter==32 && usetr==3
            % only at 8Hz, so will tie in with theta
            skipplot=1;
          elseif sleep==0 && ll==3 && iter==27 && usetr==1
            skipplot=1;
          elseif sleep==0 && ll==6 && iter==27 && usetr==1
            cfg.ylim=[8 14]; % 2 clusters: early beta and later alpha
          elseif sleep==0 && ll==7 && iter==27 && usetr==1
            skipplot=1; % only at 8Hz, so will tie in with theta
          elseif sleep==1 && ll==3 && trialkc==-1
            cfg.ylim=[8 14];
          elseif sleep==1 && ll==5 && trialkc==0
            skipplot=1; % tie in with theta
          elseif sleep==1 && ll==6 && trialkc==0
            skipplot=1; % tie in with beta
          else
            disp('get ylim right per ll alpha plv')
            tmp.freq(find(mean(mean(tmp.mask,1),3)))
            keyboard
          end
          if skipplot==0
            masktime_tmp=find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),:),2),1)));
            difftimes=diff(find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),:),2),1))));
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
            for dd=1:length(masktime)
              
              if adda==2
                cfg.parameter='plvavgabs';
              else adda==3
                cfg.parameter='plvabs';
              end
              cfg.layout='elec1010.lay';
              cfg.maskalpha=0.5;
              %     cfg.zlim='maxabs';
              cfg.highlight='on';
              cfg.xlim=[tmp.time(masktime{dd}(1)) tmp.time(masktime{dd}(end))];
              cfg.highlightchannel=tmp.label(find(ceil(mean(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2)) ),2),3))));
              cfg.baseline=baseline2;
              cfg.zlim=[-.15 .15];
              if commentson
                cfg.comment='auto';
              else
                cfg.comment='no';
              end
              cfg.highlightsize=25;
              cfg.highlightsymbol='.';
              figure(100*ll+5);
              ft_topoplotTFR(cfg,tmpuA);
              print(100*ll+5,[fdir 'plvabslo_topoU_alpha_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '_time' num2str(dd) '.eps'],'-depsc2')
              figure(100*ll+6);
              ft_topoplotTFR(cfg,tmpsA);
              print(100*ll+6,[fdir 'plvabslo_topoM_alpha_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '_time' num2str(dd) '.eps'],'-depsc2')
              cfg.zlim=[-.1 .1];
              figure(100*ll+7);
              ft_topoplotTFR(cfg,tmpA);
              print(100*ll+7,[fdir 'plvabslo_topoDiff_alpha_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '_time' num2str(dd) '.eps'],'-depsc2')
            end
            
            
            if ~isempty(setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2})))
              chansel=cfg.highlightchannel;
              zlim=[-.3 .3; 0 0.35; -30 30; -60 60];
              zlim=[0 0.6; -.25 .25; -20 400; -10 300];
              figinds=[100*ll+40+1;  100*ll+50+1];
              figstrings{1}=[fdir 'plvabslo_all_Xalpha_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
              figstrings{2}=[fdir 'plvanglo_all_Xalpha_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
              tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings,adda)
            end
            
            cfg.parameter='plvavgangwrap';
            if sleep==0 && ll==6 && pre30plot==1
              cfg.xlim=[.34 .34];
              %           elseif sleep==0 && ll==3 && combval==1 && adda==1 && pre30plot==1
              %             cfg.xlim=[0 0];
              %           elseif sleep==0 && ll==3 && adda==2 && pre30plot==1
              %             cfg.xlim=[-.06 -.06];
              %           elseif sleep==0 && ll==3 && combval==2 && adda==1 && pre30plot==1
              %             cfg.xlim=[.21 .21];
              %           elseif sleep==0 && ll==5 && combval==1 && adda==2 && iter==31 && usetr==3
              %             cfg.xlim=[.2 .2];
              %           elseif sleep==0 && ll==6 && combval==1 && adda==2 && iter==31 && usetr==3
              %             cfg.xlim=[.2 .2];
              %           elseif sleep==0 && ll==3 && combval==1 && adda==2 && iter==32 && usetr==3
              %             cfg.xlim=[.17 .17];
              %           elseif sleep==0 && ll==5 && combval==1 && adda==2 && iter==32 && usetr==3
              %             cfg.xlim=[.11 .11];
              %           elseif sleep==0 && ll==6 iter==27 && usetr==1
              %             cfg.xlim=[.2 .2];
            elseif sleep==1 && ll==3 && trialkc==-1
              cfg.xlim=[1.03 1.03];
            else
              disp('get xlim right per ll alpha plv')
              tmp.time(find(mean(mean(tmp.mask,1),2)))
              keyboard
            end
            cfg.highlightsize=25;
            cfg.highlightsymbol='.';
            cfg.zlim=[-4 4];
            figure(100*ll+5+10);
            ft_topoplotTFR(cfg,tmpuA);
            print(100*ll+5+10,[fdir 'plvanglo_topoU_alpha_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
            figure(100*ll+6+10);
            ft_topoplotTFR(cfg,tmpsA);
            print(100*ll+6+10,[fdir 'plvanglo_topoM_alpha_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
            figure(100*ll+7+10);
            ft_topoplotTFR(cfg,tmpA);
            print(100*ll+7+10,[fdir 'plvanglo_topoDiff_alpha_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
          end % skipplot
        end
        
        % beta topo
        skipplot=0;
        if ~isempty(find(squeeze(any(mean(tmp.mask(:,6:10,:),2),1))))
          cfg=[];
          if sleep==0 && ll==6 && combval==1 && adda==2 && pre30plot==1
            cfg.ylim=[13 15];
          elseif sleep==0 && ll==3 && pre30plot==1
            cfg.ylim=[13 15];
          elseif sleep==1 && ll==1 && pre30plot==1
            cfg.ylim=[13 21];
          elseif sleep==0 && ll==5 && combval==1 && adda==2 && iter==31 && usetr==3
            % only 14 hz, tied on with alpha; ignore
            skipplot=1;
          elseif sleep==0 && ll==6 && combval==1 && adda==2 && iter==31 && usetr==3
            cfg.ylim=[11 27]; % there are 2 clusters: early beta and later alpha
          elseif sleep==0 && ll==7 && combval==1 && adda==2 && iter==31 && usetr==3
            cfg.ylim=[15 25]; % final
          elseif sleep==0 && ll==5 && combval==1 && adda==2 && iter==32 && usetr==3
            % only '6' (14 hz), tied on with alpha; ignore
            skipplot=1;
          elseif sleep==0 && ll==6 && combval==1 && adda==2 && iter==32 && usetr==3
            % only '6' (14 hz), tied on with alpha; ignore
            skipplot=1;
          elseif sleep==0 && ll==7 && combval==1 && adda==2 && iter==32 && usetr==3
            cfg.ylim=[15 25];
            %         elseif sleep==0 && ll==7 && combval==1 && adda==2 && iter==31 && usetr==2
            %           cfg.ylim=[13 19];
          elseif sleep==0 && ll==6 && iter==27 && usetr==1
            cfg.ylim=[14 18];
          elseif sleep==1 && ll==6 && trialkc==0
            cfg.ylim=[12 18];
          elseif sleep==1 && ll==7 && trialkc==0
            cfg.ylim=[12 18]; % p=0.056
          else
            disp('get ylim right per ll beta plv')
            tmp.freq(find(mean(mean(tmp.mask,1),3)))
            keyboard
          end
          if skipplot==0
            masktime=find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),:),2),1)));
            if sleep==0 && ll==6 && combval==1 && adda==2 && iter==31 && usetr==3
              masktime=masktime(1:10);   % early beta blob only
            end
            cfg.xlim=[tmp.time(masktime(1)) tmp.time(masktime(end))];
            %%
            if pptdemo
              cfg=[];
              cfg.ylim=[12 18];
              cfg.xlim=[.15 .3]+soades(ll);
            end
            if adda==2
              cfg.parameter='plvavgabs';
            else adda==3
              cfg.parameter='plvabs';
            end
            cfg.layout='elec1010.lay';
            cfg.maskalpha=0.5;
            %     cfg.zlim='maxabs';
            cfg.highlight='on';
            if pptdemo
              load betaITCtrialkc0.mat
              cfg.highlightchannel=betaITCtrialkc0_intersect;
              commentson=1;
            else
              cfg.highlightchannel=tmp.label(find(ceil(mean(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2)) ),2),3))));
            end
%             cfg.baseline=baseline2;
            cfg.zlim=[.1 .2];
            if commentson
              cfg.comment='auto';
            else
              cfg.comment='no';
            end
            cfg.highlightsize=25;
            cfg.highlightsymbol='.';
            figure(100*ll+5);
            ft_topoplotTFR(cfg,tmpuA);
            print(100*ll+5,[fdir 'plvabslo_topoU_beta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
            figure(100*ll+6);
            ft_topoplotTFR(cfg,tmpsA);
            print(100*ll+6,[fdir 'plvabslo_topoM_beta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
            cfg.zlim=[-.1 .1];
            figure(100*ll+7);
            ft_topoplotTFR(cfg,tmpA);
            print(100*ll+7,[fdir 'plvabslo_topoDiff_beta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
            
            if ~isempty(setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2}))) || pptdemo
              if pptdemo
                load betaITCtrialkc0.mat
                chansel=betaITCtrialkc0_intersect;
              else
                chansel=cfg.highlightchannel;
              end
%               zlim=[-.3 .3; 0 0.35; -30 30; -60 60];
              zlim=[.1 .2; -.1 .1; -20 400; -10 300];
              figinds=[100*ll+40+1;  100*ll+50+1];
              figstrings{1}=[fdir 'plvabslo_final_Xbeta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
              figstrings{2}=[fdir 'plvanglo_final_Xbeta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
              tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings,adda)
%               betaITCtrialkc0AT20=chansel;
              % save betaITCtrialkc0AT20.mat betaITCtrialkc0AT20
            end
            %%
            
            cfg.parameter='plvavgangwrap';
            if sleep==0 && ll==6 && combval==1 && adda==2 && pre30plot==1
              cfg.xlim=[.32 .32];
              %           elseif sleep==0 && ll==3 && combval==1 && adda==1 && pre30plot==1
              %             cfg.xlim=[0 0];
              %           elseif sleep==0 && ll==3 && adda==2 && pre30plot==1
              %             cfg.xlim=[-.06 -.06];
              %           elseif sleep==1 && ll==1 && adda==2 && pre30plot==1
              %             cfg.xlim=[0.05 0.61];
              %           elseif sleep==0 && ll==6 && combval==1 && adda==2 && iter==31 && usetr==3
              %             cfg.xlim=[.12 .12];
              %           elseif sleep==0 && ll==7 && combval==1 && adda==2 && iter==31 && usetr==3
              %             cfg.xlim=[.11 .11]; % this is beginning time for sign period; is that consistent with done above/previously?
              %           elseif sleep==0 && ll==7 && combval==1 && adda==2 && iter==32 && usetr==3
              %             cfg.xlim=[.11 .11]; % this is beginning time for sign period; is that consistent with done above/previously?
              %           elseif sleep==0 && ll==6 && iter==27 && usetr==1
              %             cfg.xlim=[.03 .03];
              %           elseif sleep==1 && ll==6 && trialkc==0
              %             cfg.xlim=[.05 .62];
            else
              disp('get xlim right per ll beta plv')
              tmp.time(find(mean(mean(tmp.mask,1),2)))
              keyboard
            end
            cfg.zlim=[-4 4];
            cfg.highlightsize=25;
            cfg.highlightsymbol='.';
            figure(100*ll+5+10);
            ft_topoplotTFR(cfg,tmpuA);
            print(100*ll+5+10,[fdir 'plvanglo_topoU_beta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
            figure(100*ll+6+10);
            ft_topoplotTFR(cfg,tmpsA);
            print(100*ll+6+10,[fdir 'plvanglo_topoM_beta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
            figure(100*ll+7+10);
            ft_topoplotTFR(cfg,tmpA);
            print(100*ll+7+10,[fdir 'plvanglo_topoDiff_beta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
          end % skippplot
        end
        
%       end % adda
      
    end % plotplvflag
    
%   end % combval
  
  % %
  
  %   if synchasynch && ll<5
  %     for combval=1
  %       close all
  %       clear adda
  %       % power
  %       if combval==1
  %         cfg=[];
  %         cfg.latency=[stattl_synch_comb1.time(1) stattl_synch_comb1.time(end)];
  %         tmp=ft_selectdata(cfg,gravelo_TMSs_TMSa_comb1);
  %         tmp.mask=stattl_synch_comb1.mask;
  %         tmps=ft_selectdata(cfg,gravelo_tMSasynch_comb1);
  %         tmps.mask=stattl_synch_comb1.mask;
  %         tmpu=ft_selectdata(cfg,gravelo_tMSsynch_comb1);
  %         tmpu.mask=stattl_synch_comb1.mask;
  %
  %         cfg=[];
  %         tmpA=ft_selectdata(cfg,gravelo_TMSs_TMSa_comb1);
  %         tmpA.mask=logical(ones(size(tmpA.powspctrm)));
  %         tmpA.mask(:,:,dsearchn(tmpA.time',stattl_synch_comb1.time(1)):dsearchn(tmpA.time',stattl_synch_comb1.time(end)))=stattl_synch_comb1.mask;
  %         tmpsA=ft_selectdata(cfg,gravelo_tMSasynch_comb1);
  %         tmpsA.mask=logical(ones(size(tmpsA.powspctrm)));
  %         tmpsA.mask(:,:,dsearchn(tmpsA.time',stattl_synch_comb1.time(1)):dsearchn(tmpsA.time',stattl_synch_comb1.time(end)))=stattl_synch_comb1.mask;
  %         tmpuA=ft_selectdata(cfg,gravelo_tMSsynch_comb1);
  %         tmpuA.mask=logical(ones(size(tmpuA.powspctrm)));
  %         tmpuA.mask(:,:,dsearchn(tmpuA.time',stattl_synch_comb1.time(1)):dsearchn(tmpuA.time',stattl_synch_comb1.time(end)))=stattl_synch_comb1.mask;
  %       elseif combval==2
  %         cfg=[];
  %         cfg.latency=[stattl_synch_comb2.time(1) stattl_synch_comb2.time(end)];
  %         tmp=ft_selectdata(cfg,gravelo_TMSs_TMSa_comb2);
  %         tmp.mask=stattl_synch_comb2.mask;
  %         tmps=ft_selectdata(cfg,gravelo_tMSasynch_comb2);
  %         tmps.mask=stattl_synch_comb2.mask;
  %         tmpu=ft_selectdata(cfg,gravelo_tMSsynch_comb2);
  %         tmpu.mask=stattl_synch_comb2.mask;
  %       end
  %       gravelo_tMSAlone_comb.mask= logical(ones(size(gravelo_tMSAlone_comb.powspctrm)));
  %
  %
  %       % plot TFR of powspctrm, for each suplot: u and m and difference.
  %       zlim=[-3 3];
  %       baseline3=[-.15 -.07];
  %
  %       if 0
  %         pow{1}=gravelo_tMSAlone_comb;
  %         pow{2}=tmpu;
  %         pow{3}=tmps;
  %         pow{4}=tmp;
  %
  %         chansel=chanplot{1};
  %         figinds=10*ll;
  %         figstrings=[];
  %         figstrings{1}=[fdir 'tfrlo_synch_crop_FC_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
  %         tfr_subchannel_3cond_plot_pow(pow,chansel,zlim,figinds,figstrings,baseline3)
  %
  %         chansel=chanplot{2};
  %         figinds=10*ll+1;
  %         figstrings=[];
  %         figstrings{1}=[fdir 'tfrlo_synch_crop_OP_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
  %         tfr_subchannel_3cond_plot_pow(pow,chansel,zlim,figinds,figstrings,baseline3)
  %
  %         chansel=[chanplot{1} chanplot{2}];
  %         figinds=10*ll+2;
  %         figstrings=[];
  %         figstrings{1}=[fdir 'tfrlo_synch_crop_FCOP_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
  %         tfr_subchannel_3cond_plot_pow(pow,chansel,zlim,figinds,figstrings,baseline3)
  %       end
  %
  %       pow{1}=gravelo_tMSAlone_comb;
  %       pow{2}=tmpuA;
  %       pow{3}=tmpsA;
  %       pow{4}=tmpA;
  %
  %       chansel=chanplot{1};
  %       figinds=10*ll;
  %       figstrings=[];
  %       figstrings{1}=[fdir 'tfrlo_synch_all_FC_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
  %       tfr_subchannel_3cond_plot_pow(pow,chansel,zlim,figinds,figstrings,baseline3)
  %
  %       chansel=chanplot{2};
  %       figinds=10*ll+1;
  %       figstrings=[];
  %       figstrings{1}=[fdir 'tfrlo_synch_all_OP_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
  %       tfr_subchannel_3cond_plot_pow(pow,chansel,zlim,figinds,figstrings,baseline3)
  %
  %       chansel=[chanplot{1} chanplot{2}];
  %       figinds=10*ll+2;
  %       figstrings=[];
  %       figstrings{1}=[fdir 'tfrlo_synch_all_FCOP_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
  %       tfr_subchannel_3cond_plot_pow(pow,chansel,zlim,figinds,figstrings,baseline3)
  %
  %       masktime=find(squeeze(any(mean(tmp.mask(:,1:end,:),2),1)));
  %       if ~isempty(masktime)
  %         chansel=tmp.label(find(ceil(mean(mean(tmp.mask(:,:,dsearchn(tmp.time',tmp.time(masktime(1))):dsearchn(tmp.time',tmp.time(masktime(end)))),2),3))));
  %         figinds=10*ll+9;
  %         figstrings=[];
  %         figstrings{1}=[fdir 'tfrlo_synch_all_allsig_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
  %         tfr_subchannel_3cond_plot_pow(pow,chansel,zlim,figinds,figstrings,baseline3)
  %       end
  %
  %       % theta topo
  %       if ~isempty(find(squeeze(any(mean(tmp.mask(:,1:2,:),2),1))))
  %         masktime=find(squeeze(any(mean(tmp.mask(:,1:2,:),2),1)));
  %         cfg=[];
  %         cfg.parameter='powspctrm';
  %         cfg.layout='elec1010.lay';
  %         cfg.maskalpha=0.5;
  %         if ll==4
  %           cfg.ylim=[4 8.5];
  %         else
  %           cfg.ylim=[4 6.5];
  %         end
  %         cfg.zlim=[-1.4 1.4];
  %         cfg.highlight='on';
  %         cfg.xlim=[tmp.time(masktime(1)) tmp.time(masktime(end))];
  %         %     if ll==1
  %         %       cfg.xlim=[.1 .35];
  %         %     elseif ll==3
  %         %       cfg.xlim=[.06 .42];
  %         %     elseif ll==5
  %         %       cfg.xlim=[-.04 .36];
  %         %     elseif ll==7
  %         %       cfg.xlim=[.1 .44];
  %         %     end
  %         cfg.highlightchannel=tmp.label(find(ceil(mean(mean(tmp.mask(:,1:2,dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2))),2),3))));
  %         cfg.baseline=baseline3;
  %         figure(10*ll+8);
  %         ft_topoplotTFR(cfg,tmpuA);
  %         print(10*ll+8,[fdir 'tfrlo_synch_topoU_theta_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %         figure(10*ll+3);
  %         ft_topoplotTFR(cfg,tmpsA);
  %         print(10*ll+3,[fdir 'tfrlo_synch_topoM_theta_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %         cfg.zlim=[-1.4 1.4];
  %         %     cfg.zlim='maxabs';
  %         %     cfg.baseline='no';
  %         figure(10*ll+4);
  %         ft_topoplotTFR(cfg,tmpA);
  %         print(10*ll+4,[fdir 'tfrlo_synch_topoDiff_theta_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %
  %         if 0 %~isempty(setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2})))
  %           chansel=setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2}));
  %           zlim=[-3 3];
  %           figinds=10*ll+1000;
  %           figstrings=[];
  %           figstrings{1}=[fdir 'tfrlo_synch_final_Xtheta_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
  %           baseline=[tmp.time(1) tmp.time(9)];
  %           tfr_subchannel_3cond_plot_pow(tmp,tmpu,tmps,chansel,zlim,figinds,figstrings,baseline)
  %         end
  %
  %       end
  %
  %       % alpha topo
  %       if ~isempty(find(squeeze(any(mean(tmp.mask(:,3:5,:),2),1))))
  %         cfg=[];
  %         if ll==3
  %           cfg.ylim=[8 11];
  %         elseif ll==4
  %           cfg.ylim=[10 14];
  %         else
  %           keyboard
  %         end
  %         disp('get ylim right per ll alpha synch')
  %         masktime=find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),:),2),1)));
  %         if ll==3 || ll==4
  %           cfg.xlim=[tmp.time(masktime(1)) .5];
  %         else
  %           cfg.xlim=[tmp.time(masktime(1)) tmp.time(masktime(end))];
  %         end
  %         cfg.parameter='powspctrm';
  %         cfg.layout='elec1010.lay';
  %         cfg.maskalpha=0.5;
  %         cfg.zlim=[-9 9];
  %         cfg.highlight='on';
  %         cfg.highlightchannel=tmp.label(find(ceil(mean(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2)) ),2),3))));
  %         cfg.baseline=baseline3;
  %         figure(10*ll+5);
  %         ft_topoplotTFR(cfg,tmpuA);
  %         print(10*ll+5,[fdir 'tfrlo_synch_topoU_alpha_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %         figure(10*ll+6);
  %         ft_topoplotTFR(cfg,tmpsA);
  %         print(10*ll+6,[fdir 'tfrlo_synch_topoM_alpha_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %         cfg.zlim=[-9 9];
  %         %     cfg.zlim='maxabs';
  %         %     cfg.baseline='no';
  %         figure(10*ll+7);
  %         ft_topoplotTFR(cfg,tmpA);
  %         print(10*ll+7,[fdir 'tfrlo_synch_topoDiff_alpha_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %         if ll==3 || ll==4
  %           cfg.xlim=[.5 tmp.time(masktime(end))];
  %           cfg.highlightchannel=tmp.label(find(ceil(mean(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2)) ),2),3))));
  %           figure(10*ll+5);
  %           ft_topoplotTFR(cfg,tmpuA);
  %           print(10*ll+5,[fdir 'tfrlo_synch_topoU_alpha2_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %           figure(10*ll+6);
  %           ft_topoplotTFR(cfg,tmpsA);
  %           print(10*ll+6,[fdir 'tfrlo_synch_topoM_alpha2_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %           cfg.zlim=[-9 9];
  %           %     cfg.zlim='maxabs';
  %           %     cfg.baseline='no';
  %           figure(10*ll+7);
  %           ft_topoplotTFR(cfg,tmpA);
  %           print(10*ll+7,[fdir 'tfrlo_synch_topoDiff_alpha2_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %         end
  %
  %         if 0 %~isempty(setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2})))
  %           chansel=setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2}));
  %           zlim=[-3 3];
  %           figinds=10*ll+1000+1;
  %           figstrings=[];
  %           figstrings{1}=[fdir 'tfrlo_synch_final_Xalpha_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
  %           baseline=[tmp.time(1) tmp.time(9)];
  %           tfr_subchannel_3cond_plot_pow(tmp,tmpu,tmps,chansel,zlim,figinds,figstrings,baseline)
  %         end
  %       end
  %
  %       % beta topo
  %       if ~isempty(find(squeeze(any(mean(tmp.mask(:,6:10,:),2),1))))
  %         cfg=[];
  %         if ll==3
  %           cfg.ylim=[14 28];
  %         elseif ll==4
  %           cfg.ylim=[16 30];
  %         else
  %           keyboard
  %         end
  %         disp('get ylim right per ll beta synch')
  %         masktime=find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),:),2),1)));
  %         if ll==3 || ll==4
  %           cfg.xlim=[tmp.time(masktime(1)) .5];
  %         else
  %           cfg.xlim=[tmp.time(masktime(1)) tmp.time(masktime(end))];
  %         end
  %         cfg.parameter='powspctrm';
  %         cfg.layout='elec1010.lay';
  %         cfg.maskalpha=0.5;
  %         cfg.zlim=[-1.1 1.1];
  %         cfg.highlight='on';
  %         cfg.highlightchannel=tmp.label(find(ceil(mean(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2)) ),2),3))));
  %         cfg.baseline=baseline3;
  %         figure(10*ll+5);
  %         ft_topoplotTFR(cfg,tmpuA);
  %         print(10*ll+5,[fdir 'tfrlo_synch_topoU_beta_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %         figure(10*ll+6);
  %         ft_topoplotTFR(cfg,tmpsA);
  %         print(10*ll+6,[fdir 'tfrlo_synch_topoM_beta_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %         cfg.zlim=[-1.1 1.1];
  %         %     cfg.baseline='no';
  %         %     cfg.zlim='maxabs';
  %         figure(10*ll+7);
  %         ft_topoplotTFR(cfg,tmpA);
  %         print(10*ll+7,[fdir 'tfrlo_synch_topoDiff_beta_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %         if ll==3 || ll==4
  %           cfg.xlim=[.5 tmp.time(masktime(end))];
  %           cfg.highlightchannel=tmp.label(find(ceil(mean(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2)) ),2),3))));
  %           figure(10*ll+5);
  %           ft_topoplotTFR(cfg,tmpuA);
  %           print(10*ll+5,[fdir 'tfrlo_synch_topoU_beta2_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %           figure(10*ll+6);
  %           ft_topoplotTFR(cfg,tmpsA);
  %           print(10*ll+6,[fdir 'tfrlo_synch_topoM_beta2_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %           cfg.zlim=[-1.1 1.1];
  %           %     cfg.baseline='no';
  %           %     cfg.zlim='maxabs';
  %           figure(10*ll+7);
  %           ft_topoplotTFR(cfg,tmpA);
  %           print(10*ll+7,[fdir 'tfrlo_synch_topoDiff_beta2_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %         end
  %
  %
  %         if 0 %~isempty(setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2})))
  %           chansel=setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2}));
  %           zlim=[-2 2];
  %           figinds=10*ll+1000+2;
  %           figstrings=[];
  %           figstrings{1}=[fdir 'tfrlo_synch_final_Xbeta_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
  %           baseline=[tmp.time(1) tmp.time(9)];
  %           tfr_subchannel_3cond_plot_pow(tmp,tmpu,tmps,chansel,zlim,figinds,figstrings,baseline)
  %         end
  %       end
  %
  %
  %
  %
  %
  %       %    %  % %%%%%%% PLV   %%%%%%%%
  %       for adda=2
  %         if combval==1
  %           if adda==1
  %             cfg=[];
  %             cfg.latency=[stattl_synch_comb1plvadi.time(1) stattl_synch_comb1plvadi.time(end)];
  %             tmp=ft_selectdata(cfg,gravelo_TMSs_TMSa_comb1);
  %             tmp.mask=stattl_synch_comb1plvadi.mask;
  %             tmps=ft_selectdata(cfg,gravelo_tMSasynch_comb1);
  %             tmps.mask=stattl_synch_comb1plvadi.mask;
  %             tmpu=ft_selectdata(cfg,gravelo_tMSsynch_comb1);
  %             tmpu.mask=stattl_synch_comb1plvadi.mask;
  %
  %             cfg=[];
  %             tmpA=ft_selectdata(cfg,gravelo_TMSs_TMSa_comb1);
  %             tmpA.mask=logical(ones(size(tmpA.powspctrm)));
  %             tmpA.mask(:,:,dsearchn(tmpA.time',stattl_synch_comb1plvadi.time(1)):dsearchn(tmpA.time',stattl_synch_comb1plvadi.time(end)))=stattl_synch_comb1plvadi.mask;
  %             tmpsA=ft_selectdata(cfg,gravelo_tMSasynch_comb1);
  %             tmpsA.mask=logical(ones(size(tmpsA.powspctrm)));
  %             tmpsA.mask(:,:,dsearchn(tmpsA.time',stattl_synch_comb1plvadi.time(1)):dsearchn(tmpsA.time',stattl_synch_comb1plvadi.time(end)))=stattl_synch_comb1plvadi.mask;
  %             tmpuA=ft_selectdata(cfg,gravelo_tMSsynch_comb1);
  %             tmpuA.mask=logical(ones(size(tmpuA.powspctrm)));
  %             tmpuA.mask(:,:,dsearchn(tmpuA.time',stattl_synch_comb1plvadi.time(1)):dsearchn(tmpuA.time',stattl_synch_comb1plvadi.time(end)))=stattl_synch_comb1plvadi.mask;
  %           elseif adda==2
  %             cfg=[];
  %             cfg.latency=[stattl_synch_comb1plvdai.time(1) stattl_synch_comb1plvdai.time(end)];
  %             tmp=ft_selectdata(cfg,gravelo_TMSs_TMSa_comb1);
  %             tmp.mask=stattl_synch_comb1plvdai.mask;
  %             tmps=ft_selectdata(cfg,gravelo_tMSasynch_comb1);
  %             tmps.mask=stattl_synch_comb1plvdai.mask;
  %             tmpu=ft_selectdata(cfg,gravelo_tMSsynch_comb1);
  %             tmpu.mask=stattl_synch_comb1plvdai.mask;
  %
  %             cfg=[];
  %             tmpA=ft_selectdata(cfg,gravelo_TMSs_TMSa_comb1);
  %             tmpA.mask=logical(ones(size(tmpA.powspctrm)));
  %             tmpA.mask(:,:,dsearchn(tmpA.time',stattl_synch_comb1plvdai.time(1)):dsearchn(tmpA.time',stattl_synch_comb1plvdai.time(end)))=stattl_synch_comb1plvdai.mask;
  %             tmpsA=ft_selectdata(cfg,gravelo_tMSasynch_comb1);
  %             tmpsA.mask=logical(ones(size(tmpsA.powspctrm)));
  %             tmpsA.mask(:,:,dsearchn(tmpsA.time',stattl_synch_comb1plvdai.time(1)):dsearchn(tmpsA.time',stattl_synch_comb1plvdai.time(end)))=stattl_synch_comb1plvdai.mask;
  %             tmpuA=ft_selectdata(cfg,gravelo_tMSsynch_comb1);
  %             tmpuA.mask=logical(ones(size(tmpuA.powspctrm)));
  %             tmpuA.mask(:,:,dsearchn(tmpuA.time',stattl_synch_comb1plvdai.time(1)):dsearchn(tmpuA.time',stattl_synch_comb1plvdai.time(end)))=stattl_synch_comb1plvdai.mask;
  %           end
  %         elseif combval==2
  %           if adda==1
  %             cfg=[];
  %             cfg.latency=[stattl_synch_comb2plvadi.time(1) stattl_synch_comb2plvadi.time(end)];
  %             tmp=ft_selectdata(cfg,gravelo_TMSs_TMSa_comb2);
  %             tmp.mask=stattl_synch_comb2plvadi.mask;
  %             tmps=ft_selectdata(cfg,gravelo_tMSasynch_comb2);
  %             tmps.mask=stattl_synch_comb2plvadi.mask;
  %             tmpu=ft_selectdata(cfg,gravelo_tMSsynch_comb2);
  %             tmpu.mask=stattl_synch_comb2plvadi.mask;
  %           elseif adda==2
  %             cfg=[];
  %             cfg.latency=[stattl_synch_comb2plvdai.time(1) stattl_synch_comb2plvdai.time(end)];
  %             tmp=ft_selectdata(cfg,gravelo_TMSs_TMSa_comb2);
  %             tmp.mask=stattl_synch_comb2plvdai.mask;
  %             tmps=ft_selectdata(cfg,gravelo_tMSasynch_comb2);
  %             tmps.mask=stattl_synch_comb2plvdai.mask;
  %             tmpu=ft_selectdata(cfg,gravelo_tMSsynch_comb2);
  %             tmpu.mask=stattl_synch_comb2plvdai.mask;
  %           end
  %         end
  %
  %         tmp.plvavgangwrap=wrapToPi(tmp.plvavgang);
  %         tmps.plvavgangwrap=wrapToPi(tmps.plvavgang);
  %         tmpu.plvavgangwrap=wrapToPi(tmpu.plvavgang);
  %         tmpA.plvavgangwrap=wrapToPi(tmpA.plvavgang);
  %         tmpsA.plvavgangwrap=wrapToPi(tmpsA.plvavgang);
  %         tmpuA.plvavgangwrap=wrapToPi(tmpuA.plvavgang);
  %
  %         % plot TFR of plv abs and plv ang, for each suplot: u and m and difference.
  %         zlim=[0 0.6; -.25 .25; -20 400; -10 300];
  %
  %
  %         if 0
  %           pow{1}=gravelo_tMSAlone_comb;
  %           pow{2}=tmpu;
  %           pow{3}=tmps;
  %           pow{4}=tmp;
  %
  %           chansel=chanplot{1};
  %           figinds=[100*ll;  100*ll+10];
  %           figstrings{1}=[fdir 'plvabslo_synch_crop_FC_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
  %           figstrings{2}=[fdir 'plvanglo_synch_crop_FC_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
  %           tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
  %
  %           chansel=chanplot{2};
  %           figinds=[100*ll+1;  100*ll+10+1];
  %           figstrings{1}=[fdir 'plvabslo_synch_crop_OP_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
  %           figstrings{2}=[fdir 'plvanglo_synch_crop_OP_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
  %           tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
  %
  %           chansel=[chanplot{1} chanplot{2}];
  %           figinds=[100*ll+2;  100*ll+10+2];
  %           figstrings{1}=[fdir 'plvabslo_synch_crop_FCOP_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
  %           figstrings{2}=[fdir 'plvanglo_synch_crop_FCOP_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
  %           tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
  %         end
  %
  %         pow{1}=gravelo_tMSAlone_comb;
  %         pow{2}=tmpuA;
  %         pow{3}=tmpsA;
  %         pow{4}=tmpA;
  %
  %         chansel=chanplot{1};
  %         figinds=[100*ll;  100*ll+10];
  %         figstrings{1}=[fdir 'plvabslo_synch_all_FC_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
  %         figstrings{2}=[fdir 'plvanglo_synch_all_FC_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
  %         tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
  %
  %         chansel=chanplot{2};
  %         figinds=[100*ll+1;  100*ll+10+1];
  %         figstrings{1}=[fdir 'plvabslo_synch_all_OP_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
  %         figstrings{2}=[fdir 'plvanglo_synch_all_OP_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
  %         tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
  %
  %         chansel=[chanplot{1} chanplot{2}];
  %         figinds=[100*ll+2;  100*ll+10+2];
  %         figstrings{1}=[fdir 'plvabslo_synch_all_FCOP_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
  %         figstrings{2}=[fdir 'plvanglo_synch_all_FCOP_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
  %         tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
  %
  %         masktime=find(squeeze(any(mean(tmp.mask(:,1:end,:),2),1)));
  %         if ~isempty(masktime)
  %           chansel=tmp.label(find(ceil(mean(mean(tmp.mask(:,:,dsearchn(tmp.time',tmp.time(masktime(1))):dsearchn(tmp.time',tmp.time(masktime(end)))),2),3))));
  %           figinds=[100*ll+9;  100*ll+10+9];
  %           figstrings{1}=[fdir 'plvabslo_synch_all_allsig_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
  %           figstrings{2}=[fdir 'plvanglo_synch_all_allsig_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
  %           tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
  %         end
  %
  %         % theta topo
  %         if ~isempty(find(squeeze(any(mean(tmp.mask(:,1:2,:),2),1))))
  %           cfg=[];
  %           if ll==3 && adda==1
  %             cfg.ylim=[4 6.5];
  %           elseif ll==3 && adda==2
  %             cfg.ylim=[5.5 6.5];
  %           else
  %             keyboard
  %           end
  %           disp('get ylim right per ll theta synch plv')
  %           masktime=find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),:),2),1)));
  %           cfg.xlim=[tmp.time(masktime(1)) tmp.time(masktime(end))];
  %           cfg.parameter='plvavgabs';
  %           cfg.layout='elec1010.lay';
  %           cfg.maskalpha=0.5;
  %           cfg.zlim=[-.25 .25];
  %           cfg.highlight='on';
  %           cfg.highlightchannel=tmp.label(find(ceil(mean(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2))),2),3))));
  %           cfg.baseline=baseline3;
  %           cfg.zlim=[-.25 .25];
  %           figure(100*ll+8);
  %           ft_topoplotTFR(cfg,tmpuA);
  %           print(100*ll+8,[fdir 'plvabslo_synch_topoU_theta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %           figure(100*ll+3);
  %           ft_topoplotTFR(cfg,tmpsA);
  %           print(100*ll+3,[fdir 'plvabslo_synch_topoM_theta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %           cfg.zlim=[-.15 .15];
  %           figure(100*ll+4);
  %           ft_topoplotTFR(cfg,tmpA);
  %           print(100*ll+4,[fdir 'plvabslo_synch_topoDiff_theta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %
  %           if ~isempty(setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2})))
  %             chansel=cfg.highlightchannel;
  %             zlim=[0 0.6; -.25 .25; -20 400; -10 300];
  %             figinds=[100*ll+20+1;  100*ll+30+1];
  %             figstrings{1}=[fdir 'plvabslo_synch_final_Xtheta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
  %             figstrings{2}=[fdir 'plvanglo_synch_final_Xtheta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
  %             tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
  %           end
  %
  %           cfg.parameter='plvavgangwrap';
  %           if ll==3 && adda==1
  %             cfg.xlim=[0.19 0.19];
  %           elseif ll==3 && adda==2
  %             cfg.xlim=[0 0];
  %           else
  %             keyboard
  %           end
  %           disp('get xlim right per ll theta synch plv')
  %           cfg.zlim=[-4 4];
  %           figure(100*ll+8+10);
  %           ft_topoplotTFR(cfg,tmpuA);
  %           print(100*ll+8+10,[fdir 'plvanglo_synch_topoU_theta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %           figure(100*ll+3+10);
  %           ft_topoplotTFR(cfg,tmpsA);
  %           print(100*ll+3+10,[fdir 'plvanglo_synch_topoM_theta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %           figure(100*ll+4+10);
  %           ft_topoplotTFR(cfg,tmpA);
  %           print(100*ll+4+10,[fdir 'plvanglo_synch_topoDiff_theta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %         end
  %
  %         % alpha topo
  %         if ~isempty(find(squeeze(any(mean(tmp.mask(:,3:5,:),2),1))))
  %           cfg=[];
  %           if ll==3 && adda==1
  %             cfg.ylim=[8 9];
  %           else
  %             keyboard
  %           end
  %           disp('get ylim right per ll alpha synch plv')
  %           masktime=find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),:),2),1)));
  %           cfg.xlim=[tmp.time(masktime(1)) tmp.time(masktime(end))];
  %           cfg.parameter='plvavgabs';
  %           cfg.layout='elec1010.lay';
  %           cfg.maskalpha=0.5;
  %           %     cfg.zlim='maxabs';
  %           cfg.highlight='on';
  %           cfg.highlightchannel=tmp.label(find(ceil(mean(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2)) ),2),3))));
  %           cfg.baseline=baseline3;
  %           cfg.zlim=[-.15 .15];
  %           figure(100*ll+5);
  %           ft_topoplotTFR(cfg,tmpuA);
  %           print(100*ll+5,[fdir 'plvabslo_synch_topoU_alpha_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %           figure(100*ll+6);
  %           ft_topoplotTFR(cfg,tmpsA);
  %           print(100*ll+6,[fdir 'plvabslo_synch_topoM_alpha_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %           cfg.zlim=[-.1 .1];
  %           figure(100*ll+7);
  %           ft_topoplotTFR(cfg,tmpA);
  %           print(100*ll+7,[fdir 'plvabslo_synch_topoDiff_alpha_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %
  %           if ~isempty(setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2})))
  %             chansel=cfg.highlightchannel;
  %             %             zlim=[-.3 .3; 0 0.35; -30 30; -60 60];
  %             zlim=[0 0.6; -.25 .25; -20 400; -10 300];
  %             figinds=[100*ll+40+1;  100*ll+50+1];
  %             figstrings{1}=[fdir 'plvabslo_synch_final_Xalpha_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
  %             figstrings{2}=[fdir 'plvanglo_synch_final_Xalpha_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
  %             tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
  %           end
  %
  %           cfg.parameter='plvavgangwrap';
  %           if ll==3 && adda==1
  %             cfg.xlim=[0.19 0.19];
  %           else
  %             keyboard
  %           end
  %           disp('get xlim right per ll alpha synch plv')
  %           cfg.zlim=[-4 4];
  %           figure(100*ll+5+10);
  %           ft_topoplotTFR(cfg,tmpuA);
  %           print(100*ll+5+10,[fdir 'plvanglo_synch_topoU_alpha_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %           figure(100*ll+6+10);
  %           ft_topoplotTFR(cfg,tmpsA);
  %           print(100*ll+6+10,[fdir 'plvanglo_synch_topoM_alpha_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %           figure(100*ll+7+10);
  %           ft_topoplotTFR(cfg,tmpA);
  %           print(100*ll+7+10,[fdir 'plvanglo_synch_topoDiff_alpha_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %         end
  %
  %         % beta topo
  %         if ~isempty(find(squeeze(any(mean(tmp.mask(:,6:10,:),2),1))))
  %           cfg=[];
  %           if ll==6 && combval==1 && adda==2
  %             cfg.ylim=[13 15];
  %             %         elseif ll==3
  %             %           cfg.ylim=[13 20];
  %             %           keyboard
  %           else
  %             keyboard
  %           end
  %           disp('get ylim right per ll beta synch plv')
  %           masktime=find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),:),2),1)));
  %           cfg.xlim=[tmp.time(masktime(1)) tmp.time(masktime(end))];
  %           cfg.parameter='plvavgabs';
  %           cfg.layout='elec1010.lay';
  %           cfg.maskalpha=0.5;
  %           %     cfg.zlim='maxabs';
  %           cfg.highlight='on';
  %           cfg.highlightchannel=tmp.label(find(ceil(mean(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2)) ),2),3))));
  %           cfg.baseline=baseline;
  %           cfg.zlim=[-.15 .15];
  %           figure(100*ll+5);
  %           ft_topoplotTFR(cfg,tmpuA);
  %           print(100*ll+5,[fdir 'plvabslo_synch_topoU_beta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %           figure(100*ll+6);
  %           ft_topoplotTFR(cfg,tmpsA);
  %           print(100*ll+6,[fdir 'plvabslo_synch_topoM_beta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %           cfg.zlim=[-.1 .1];
  %           figure(100*ll+7);
  %           ft_topoplotTFR(cfg,tmpA);
  %           print(100*ll+7,[fdir 'plvabslo_synch_topoDiff_beta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %
  %           if ~isempty(setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2})))
  %             chansel=cfg.highlightchannel;
  %             zlim=[-.3 .3; 0 0.35; -30 30; -60 60];
  %             zlim=[0 0.6; -.25 .25; -20 400; -10 300];
  %             figinds=[100*ll+40+1;  100*ll+50+1];
  %             figstrings{1}=[fdir 'plvabslo_synch_final_Xbeta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
  %             figstrings{2}=[fdir 'plvanglo_synch_final_Xbeta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
  %             tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
  %           end
  %
  %           cfg.parameter='plvavgangwrap';
  %           if ll==6 && combval==1 && adda==2
  %             cfg.xlim=[.32 .32];
  %           elseif ll==3 && combval==1 && adda==1
  %             cfg.xlim=[0 0];
  %           elseif ll==3 && adda==2
  %             cfg.xlim=[-.06 -.06];
  %           else
  %             keyboard
  %           end
  %           disp('get xlim right per ll beta synch plv')
  %           cfg.zlim=[-4 4];
  %           figure(100*ll+5+10);
  %           ft_topoplotTFR(cfg,tmpuA);
  %           print(100*ll+5+10,[fdir 'plvanglo_synch_topoU_beta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %           figure(100*ll+6+10);
  %           ft_topoplotTFR(cfg,tmpsA);
  %           print(100*ll+6+10,[fdir 'plvanglo_synch_topoM_beta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %           figure(100*ll+7+10);
  %           ft_topoplotTFR(cfg,tmpA);
  %           print(100*ll+7+10,[fdir 'plvanglo_synch_topoDiff_beta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
  %         end
  %
  %       end % adda
  %
  %     end % combval
  %   end  % if ll<5
  
end % ll
