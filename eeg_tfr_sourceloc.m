% source loc of TFR

clear all
close all
if ispc
  edir='D:\audtac\eeg_data\';
  esdir='D:\audtac\source_data\';
  ddir='D:\audtac\legomagic\diaries\';
  bdir='D:\audtac\behav_data\';
  sdir='D:\audtac\spss_stuff\';
  fdir='D:\audtac\figs\';
  mdir='D:\audtac\structural_MRI\';
  pdir='D:\audtac\polhemus\';
else
  [~,hostname]=system('hostname');
  if ~isempty(strfind(hostname,'les')) | ~isempty(strfind(hostname,'LES')) % either COLLES-151401 or LES-LINUX_FS3
    edir='/home/zumerj/audtac/eeg_data/';
    esdir='/home/zumerj/audtac/source_data/';
    ddir='/home/zumerj/audtac/legomagic/diaries/';
    bdir='/home/zumerj/audtac/behav_data/';
    sdir='/home/zumerj/audtac/spss_stuff/';
    fdir='/home/zumerj/audtac/figs/';
    mdir='/home/zumerj/audtac/structural_MRI/';
    pdir='/home/zumerj/audtac/polhemus/';
  else % assume on VM linux of psychl-132432
    edir='/mnt/hgfs/D/audtac/eeg_data/';
    esdir='/mnt/hgfs/D/audtac/source_data/';
    ddir='/mnt/hgfs/D/audtac/legomagic/diaries/';
    bdir='/mnt/hgfs/D/audtac/behav_data/';
    sdir='/mnt/hgfs/D/audtac/spss_stuff/';
    fdir='/mnt/hgfs/D/audtac/figs/';
    mdir='/mnt/hgfs/D/audtac/structural_MRI/';
    pdir='/mnt/hgfs/D/audtac/polhemus/';
  end
end

sub{5}='e05'; % ma.a. 25/06/14
sub{6}='e06'; % e.u.  01/07/14
sub{7}='e07'; % a.s.  07/07/14
sub{8}='e08'; % k.t.  09/07/14
sub{9}='e09';% d.a.  14/07/14
sub{10}='e10';% k.l.  15/07/14
sub{11}='e11';% ab.m.  16/07/14  % from here on, had EOG/EMG
sub{12}='e12';% b.s.  17/07/14
sub{13}='e13';% d.t.  21/07/14
sub{14}='e14';% f.g.  22/07/14
sub{15}='e15';% r.m.  23/07/14
sub{16}='e16';% t.p.  24/07/14 % from here on, had attempt Polhemus
sub{17}='e17';% t.t.  28/07/14
sub{18}='e18';% k.n.v.  29/07/14
sub{19}='e19';% j.b.  30/07/14
sub{20}='e20';% n.m.  31/07/14
sub{21}='e21';% l.c.  04/08/14
sub{22}='e22';% a.b.  05/08/14
sub{23}='e23';% r.c.
sub{24}='e24';% a.d.
sub{25}='e25';% j.c.
sub{26}='e26';% r.s.
sub{27}='e27';% a.p.
sub{28}='e28';% w.p.
sub{29}='e29';% i.r.
sub{30}='e30';% o.y.l.
sub{31}='e31';% r.b.
sub{32}='e32';% i.f.

if ispc
  rmpath(genpath('D:\matlab\spm8\external\fieldtrip\'))
  rmpath(genpath('D:\fieldtrip_svn\'))
  addpath('D:\fieldtrip_svn\')
else
  rmpath(genpath('/mnt/hgfs/D/matlab/spm8/external/fieldtrip/'))
  rmpath(genpath('/mnt/hgfs/D/fieldtrip_svn/'))
  addpath('/mnt/hgfs/D/fieldtrip_svn/')
end
which ft_defaults.m
ft_defaults;

fwd=fileparts(which('ft_defaults.m'))
if ispc
  fwd=[fwd '\'];
else
  fwd=[fwd '/'];
end

load([edir 'iikeep.mat'])

soalist=[1 3 4 5 6 7 9];
freqlim=[4 7; 8 12; 14 30; 30 45; 55 80];
soades=[-.5 nan -.07 -.02 0 .02 .07 nan .5];
tacbasemax=min(-.15,soades-.15);
audbasemax=min(-.15,fliplr(soades)-.15);
msstart=max(0,soades);

iiSuse=setdiff(iiSuse,11);

cd(esdir)

%%

aal=ft_read_atlas([fwd 'template/atlas/aal/ROI_MNI_V4.nii']);
[indx,indy,indz]=ndgrid(1:size(aal.tissue,1),1:size(aal.tissue,2),1:size(aal.tissue,3));
hglmni=aal.transform*[indx(find(aal.tissue==79)) indy(find(aal.tissue==79)) indz(find(aal.tissue==79)) ones(length(find(aal.tissue==79)),1)]';
hgrmni=aal.transform*[indx(find(aal.tissue==80)) indy(find(aal.tissue==80)) indz(find(aal.tissue==80)) ones(length(find(aal.tissue==80)),1)]';
stlmni=aal.transform*[indx(find(aal.tissue==81)) indy(find(aal.tissue==81)) indz(find(aal.tissue==81)) ones(length(find(aal.tissue==81)),1)]';
strmni=aal.transform*[indx(find(aal.tissue==82)) indy(find(aal.tissue==82)) indz(find(aal.tissue==82)) ones(length(find(aal.tissue==82)),1)]';
hglmni=hglmni(1:3,:)';
hgrmni=hgrmni(1:3,:)';
stlmni=stlmni(1:3,:)';
strmni=strmni(1:3,:)';

plotflag=0;
saveflag=1;
supflag=1;
dicsflag=0;
mneflag=0;

% for ii=setdiff(iiSuse,[8:11])
for ii=iiSuse
  % for ii=8
  clearvars -except ii* sub *dir sleep hostname soa* *max freqlim msstart *flag *mni
  load([esdir 'lf_' sub{ii} '.mat']);
%   load([esdir 'elec_' sub{ii} '.mat']);  % major bug to use this one instead of one below!
  load([mdir sub{ii} '/elec.mat']); 
  load([mdir sub{ii} '/vol.mat']);
  load([mdir sub{ii} '/mrinorm.mat']);
  cd([edir sub{ii} ])
  
  
  hglsub=ft_warp_apply(mrinormwarp.params, hglmni, 'sn2individual');
  hgrsub=ft_warp_apply(mrinormwarp.params, hgrmni, 'sn2individual');
  stlsub=ft_warp_apply(mrinormwarp.params, stlmni, 'sn2individual');
  strsub=ft_warp_apply(mrinormwarp.params, strmni, 'sn2individual');
  
  %   % assess if lf looks 'reasonable'; this is done already in eeg_mriheadmodels.m
  %   tlockcm=[];
  %   tlockcm.time=[1 2 3];
  %   tlockcm.avg=lf.leadfield{dsearchn(lf.pos,[30 0 40])};
  %   tlockcm.avg=lf.leadfield{dsearchn(lf.pos,[-20 0 50])};
  %   tlockcm.label=lf.label;
  %   tlockcm.dimord='chan_time';
  %   cfg=[];
  %   cfg.layout='elec1010.lay';
  %   cfg.xlim=[0.9 1.1];
  %   figure;
  %   ft_topoplotER(cfg,tlockcm);
  %   cfg.xlim=[1.9 2.1];
  %   figure;
  %   ft_topoplotER(cfg,tlockcm);
  %   cfg.xlim=[2.9 3.1];
  %   figure;
  %   ft_topoplotER(cfg,tlockcm);
  
  for sleep=[0]
    %     tmp0=load(['freq_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud0.mat']);
    tmp1=load(['freq_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud1.mat']);
    
    if sleep
      ssuse=[12 13 23];
    else
      ssuse=10;
    end
    for tt=2
      for ss=ssuse
        for ll=soalist
          for ff=1:size(freqlim,1)
            
            %             keyboard % use/test _comb as well
            if ff<4
              freqN=tmp1.freqlo_tNulAlone_fftadd{ll,tt,ss};
              freqT=tmp1.freqlo_tTacAlone_fftadd{ll,tt,ss};
              freqA=tmp1.freqlo_tAudAlone_fftadd{ll,tt,ss};
              freqS=tmp1.freqlo_tMSAlone_fftadd{ll,tt,ss};
              freqU=tmp1.freqlo_tacPaud_fftadd{ll,tt,ss};
              freqM=tmp1.freqlo_tacMSpN_fftadd{ll,tt,ss};
            else
              freqN=tmp1.freqhi_tNulAlone_fftadd{ll,tt,ss};
              freqT=tmp1.freqhi_tTacAlone_fftadd{ll,tt,ss};
              freqA=tmp1.freqhi_tAudAlone_fftadd{ll,tt,ss};
              freqS=tmp1.freqhi_tMSAlone_fftadd{ll,tt,ss};
              freqU=tmp1.freqhi_tacPaud_fftadd{ll,tt,ss};
              freqM=tmp1.freqhi_tacMSpN_fftadd{ll,tt,ss};
            end
            
            cfg=[];
            cfg.parameter='powspctrm';
            freqall=ft_freqgrandaverage(cfg,freqN,freqT,freqA,freqS,freqU,freqM);
            freqU.labelcmb=freqN.labelcmb;
            freqM.labelcmb=freqN.labelcmb;
            cfg=[];
            cfg.parameter='crsspctrm';
            tmp=ft_appendfreq(cfg,freqN,freqT,freqA,freqS,freqU,freqM);
            freqall.crsspctrm=squeeze(mean(tmp.crsspctrm,1));
            
            %           freq.time=[round(1000*freq.time(1)):round(1000*(mean(diff(freq.time)))):round(1000*freq.time(end))]/1000;
            %           cfg=[];
            %           cfg.layout='elec1010.lay';
            %           cfg.baseline='yes';
            %           cfg.baselinewindow=[-.5 -.9];
            %           cfg.baselinetype='relchange';
            %           %         cfg.ylim=[4 6];
            %           ft_multiplotTFR(cfg,freq)
            
            if supflag
              cfg=[];
              cfg.grid=lf;
              cfg.vol=voldipolimm;
              cfg.elec=elec;
              cfg.frequency=freqlim(ff,:);
              cfg.latency=[tacbasemax(ll)-.35; tacbasemax(ll)+1.15];
              cfg.keepfilter='yes';
              cfg.method='pcc';
              cfg.pcc.lambda='5%';
              cfg.pcc.realfilter='yes';
              cfg.pcc.feedback='none';
              cfg.pcc.supdip=lf.pos(unique(dsearchn(lf.pos,hglsub)),:);
              cfg.pcc.supdip=lf.pos(unique(dsearchn(lf.pos,mean(hglsub))),:);
              source_alltime{ll,tt,ff}=ft_sourceanalysis(cfg,freqall);
              
              for yy=1:length(source_alltime{ll,tt,ff}.filter)
                try
                  filter{ll,tt,ff}{yy}=source_alltime{ll,tt,ff}.filter{yy}(1:3,:);
                catch
                  filter{ll,tt,ff}{yy}=[];
                end
              end
              %               cfg.pcc.fixedori='yes';
              %               source_pcc_alltime1{ll,tt,ff}=ft_sourceanalysis(cfg,freqall);
              %               filter{ll,tt,ff}=source_alltime{ll,tt,ff}.filter;
            elseif dicsflag
              cfg=[];
              cfg.grid=lf;
              cfg.vol=voldipolimm;
              cfg.elec=elec;
              cfg.method='dics';
              cfg.dics.lambda='5%';
              cfg.dics.realfilter='yes';
              cfg.dics.feedback='none';
              cfg.frequency=freqlim(ff,:);
              cfg.latency=[tacbasemax(ll)-.35; tacbasemax(ll)+1.15];
              cfg.keepfilter='yes';
              source_alltime{ll,tt,ff}=ft_sourceanalysis(cfg,freqall);
              %               cfg.dics.fixedori='yes';
              %               source_dics_alltime1{ll,tt,ff}=ft_sourceanalysis(cfg,freqall);
              filter{ll,tt,ff}=source_alltime{ll,tt,ff}.filter;
            elseif mneflag
              cfg=[];
              cfg.grid=lf;
              cfg.vol=voldipolimm;
              cfg.elec=elec;
              cfg.method='mne';
              cfg.frequency=freqlim(ff,:);
              %               cfg.latency=[tacbasemax(ll)-.35; tacbasemax(ll)+1.15];
              cfg.keepfilter='yes';
              cfg.mne.snr=1;
              source_alltime{ll,tt,ff}=ft_sourceanalysis(cfg,freqall);
              %               cfg.dics.fixedori='yes';
              %               source_dics_alltime1{ll,tt,ff}=ft_sourceanalysis(cfg,freqall);
              filter{ll,tt,ff}=source_alltime{ll,tt,ff}.filter;
            end
            
            % assess if filter looks 'reasonable'
            if plotflag
              if ll==1
                tlockcm=[];
                tlockcm.time=[1 2 3];
                tlockcm.avg=source_alltime{ll,tt,ff}.filter{dsearchn(source_alltime{ll,tt,ff}.pos,mean(hgrsub))}(1:3,:)';
                %                 tlockcm.avg=source_dics_alltime{ll,tt,ff}.filter{dsearchn(source_alltime{ll,tt,ff}.pos,mean(hgrsub))}(1:3,:)';
                tlockcm.label=source_alltime{ll,tt,ff}.cfg.channel;
                %           tlockcm.label=lf.label;
                tlockcm.dimord='chan_time';
                cfg=[];
                cfg.layout='elec1010.lay';
                cfg.xlim=[0.9 1.1];
                figure;
                ft_topoplotER(cfg,tlockcm);
                cfg.xlim=[1.9 2.1];
                figure;
                ft_topoplotER(cfg,tlockcm);
                cfg.xlim=[2.9 3.1];
                figure;
                ft_topoplotER(cfg,tlockcm);
              end
            end
            
            
            %           if 0
            %             timesteps=-.6:.02:.7;
            %             cfg=[];
            %             cfg.grid=lf;
            %             cfg.grid.filter=source_dics_alltime.filter;
            %             cfg.vol=voldipolimm;
            %             cfg.elec=elec;
            %             cfg.method='dics';
            %             cfg.frequency=[4 6];
            %             cfg.keepfilter='no';
            %             for mm=1:length(timesteps)
            %               cfg.latency=timesteps(mm);
            %               source_dics_each=ft_sourceanalysis(cfg,freq);
            %               if mm==1
            %                 source_dics=source_dics_each;
            %               end
            %               source_dics.avg.pow(mm,:)=source_dics_each.avg.pow;
            %             end
            %             source_dics.time=timesteps;
            %             source_dics.pow=source_dics.avg.pow;
            %             source_dics.dimord='time_pos';
            %
            %             cfg=[];
            %             cfg.funparameter='pow';
            %             ft_sourceplot(cfg,source_dics);
            %
            %             post=mean(source_dics.avg.pow(dsearchn(timesteps',.05):dsearchn(timesteps',.15),:),1);
            %             pre=mean(source_dics.avg.pow(dsearchn(timesteps',-.6):dsearchn(timesteps',-.3),:),1);
            %             source_relative=source_dics;
            %             source_relative=rmfield(source_relative,'time');
            %             source_relative.avg.pow=post./pre-1;
            %
            %           else
            
            cfg=[];
            cfg.grid=lf;
            cfg.grid.filter=filter{ll,tt,ff};
            cfg.vol=voldipolimm;
            cfg.elec=elec;
            cfg.method='dics';
            cfg.dics.feedback='none';
            cfg.frequency=freqlim(ff,:);
            
            cfg.latency=[tacbasemax(ll)-.35; tacbasemax(ll)+.05]; % 400ms window before anything.
            source_Npre{ll,tt,ff}=ft_sourceanalysis(cfg,freqN);
            source_Tpre{ll,tt,ff}=ft_sourceanalysis(cfg,freqT);
            source_Apre{ll,tt,ff}=ft_sourceanalysis(cfg,freqA);
            source_Spre{ll,tt,ff}=ft_sourceanalysis(cfg,freqS);
            source_Upre{ll,tt,ff}=ft_sourceanalysis(cfg,freqU);
            source_Mpre{ll,tt,ff}=ft_sourceanalysis(cfg,freqM);
            
            cfg.latency=[0 .15]; % early poststim
            source_Npost1{ll,tt,ff}=ft_sourceanalysis(cfg,freqN);
            source_Tpost1{ll,tt,ff}=ft_sourceanalysis(cfg,freqT);
            cfg.latency=[soades(ll); soades(ll)+.15];
            source_Apost1{ll,tt,ff}=ft_sourceanalysis(cfg,freqA);
            cfg.latency=[msstart(ll); msstart(ll)+.15];
            source_Spost1{ll,tt,ff}=ft_sourceanalysis(cfg,freqS);
            source_Upost1{ll,tt,ff}=ft_sourceanalysis(cfg,freqU);
            source_Mpost1{ll,tt,ff}=ft_sourceanalysis(cfg,freqM);
            
            cfg.latency=[.15 .3]; % late poststim
            source_Npost2{ll,tt,ff}=ft_sourceanalysis(cfg,freqN);
            source_Tpost2{ll,tt,ff}=ft_sourceanalysis(cfg,freqT);
            cfg.latency=[soades(ll)+.15; soades(ll)+.3];
            source_Apost2{ll,tt,ff}=ft_sourceanalysis(cfg,freqA);
            cfg.latency=[msstart(ll)+.15; msstart(ll)+.3];
            source_Spost2{ll,tt,ff}=ft_sourceanalysis(cfg,freqS);
            source_Upost2{ll,tt,ff}=ft_sourceanalysis(cfg,freqU);
            source_Mpost2{ll,tt,ff}=ft_sourceanalysis(cfg,freqM);
            
            source_Npre{ll,tt,ff}=verwijder_onnodige_velden(source_Npre{ll,tt,ff});
            source_Tpre{ll,tt,ff}=verwijder_onnodige_velden(source_Tpre{ll,tt,ff});
            source_Apre{ll,tt,ff}=verwijder_onnodige_velden(source_Apre{ll,tt,ff});
            source_Spre{ll,tt,ff}=verwijder_onnodige_velden(source_Spre{ll,tt,ff});
            source_Upre{ll,tt,ff}=verwijder_onnodige_velden(source_Upre{ll,tt,ff});
            source_Mpre{ll,tt,ff}=verwijder_onnodige_velden(source_Mpre{ll,tt,ff});
            source_Npost1{ll,tt,ff}=verwijder_onnodige_velden(source_Npost1{ll,tt,ff});
            source_Tpost1{ll,tt,ff}=verwijder_onnodige_velden(source_Tpost1{ll,tt,ff});
            source_Apost1{ll,tt,ff}=verwijder_onnodige_velden(source_Apost1{ll,tt,ff});
            source_Spost1{ll,tt,ff}=verwijder_onnodige_velden(source_Spost1{ll,tt,ff});
            source_Upost1{ll,tt,ff}=verwijder_onnodige_velden(source_Upost1{ll,tt,ff});
            source_Mpost1{ll,tt,ff}=verwijder_onnodige_velden(source_Mpost1{ll,tt,ff});
            source_Npost2{ll,tt,ff}=verwijder_onnodige_velden(source_Npost2{ll,tt,ff});
            source_Tpost2{ll,tt,ff}=verwijder_onnodige_velden(source_Tpost2{ll,tt,ff});
            source_Apost2{ll,tt,ff}=verwijder_onnodige_velden(source_Apost2{ll,tt,ff});
            source_Spost2{ll,tt,ff}=verwijder_onnodige_velden(source_Spost2{ll,tt,ff});
            source_Upost2{ll,tt,ff}=verwijder_onnodige_velden(source_Upost2{ll,tt,ff});
            source_Mpost2{ll,tt,ff}=verwijder_onnodige_velden(source_Mpost2{ll,tt,ff});
            
            
            % % These are not per ll but do it here anyway for each (yes for ff, tt, ss)
            %             if ll==5
            % each condition's post versus it's own pre
            sourcerelchange_Npost1pre{ll,tt,ff}=source_Npost1{ll,tt,ff};
            sourcerelchange_Npost1pre{ll,tt,ff}.avg.pow=source_Npost1{ll,tt,ff}.avg.pow./source_Npre{ll,tt,ff}.avg.pow-1;
            sourcerelchange_Npost2pre{ll,tt,ff}=source_Npost1{ll,tt,ff};
            sourcerelchange_Npost2pre{ll,tt,ff}.avg.pow=source_Npost2{ll,tt,ff}.avg.pow./source_Npre{ll,tt,ff}.avg.pow-1;
            sourcerelchange_Tpost1pre{ll,tt,ff}=source_Npost1{ll,tt,ff};
            sourcerelchange_Tpost1pre{ll,tt,ff}.avg.pow=source_Tpost1{ll,tt,ff}.avg.pow./source_Tpre{ll,tt,ff}.avg.pow-1;
            sourcerelchange_Tpost2pre{ll,tt,ff}=source_Npost1{ll,tt,ff};
            sourcerelchange_Tpost2pre{ll,tt,ff}.avg.pow=source_Tpost2{ll,tt,ff}.avg.pow./source_Tpre{ll,tt,ff}.avg.pow-1;
            sourcerelchange_Apost1pre{ll,tt,ff}=source_Npost1{ll,tt,ff};
            sourcerelchange_Apost1pre{ll,tt,ff}.avg.pow=source_Apost1{ll,tt,ff}.avg.pow./source_Apre{ll,tt,ff}.avg.pow-1;
            sourcerelchange_Apost2pre{ll,tt,ff}=source_Npost1{ll,tt,ff};
            sourcerelchange_Apost2pre{ll,tt,ff}.avg.pow=source_Apost2{ll,tt,ff}.avg.pow./source_Apre{ll,tt,ff}.avg.pow-1;
            
            % each stimulus condition versus null
            sourcerelchange_Tpost1Npost1{ll,tt,ff}=source_Npost1{ll,tt,ff};
            sourcerelchange_Tpost1Npost1{ll,tt,ff}.avg.pow=source_Tpost1{ll,tt,ff}.avg.pow./source_Npost1{ll,tt,ff}.avg.pow-1;
            sourcerelchange_Tpost2Npost2{ll,tt,ff}=source_Npost1{ll,tt,ff};
            sourcerelchange_Tpost2Npost2{ll,tt,ff}.avg.pow=source_Tpost2{ll,tt,ff}.avg.pow./source_Npost2{ll,tt,ff}.avg.pow-1;
            sourcerelchange_Apost1Npost1{ll,tt,ff}=source_Npost1{ll,tt,ff};
            sourcerelchange_Apost1Npost1{ll,tt,ff}.avg.pow=source_Apost1{ll,tt,ff}.avg.pow./source_Npost1{ll,tt,ff}.avg.pow-1;
            sourcerelchange_Apost2Npost2{ll,tt,ff}=source_Npost1{ll,tt,ff};
            sourcerelchange_Apost2Npost2{ll,tt,ff}.avg.pow=source_Apost2{ll,tt,ff}.avg.pow./source_Npost2{ll,tt,ff}.avg.pow-1;
            %             end
            
            % % Now the rest are also dependent on ll
            % each condition's post versus it's own pre
            sourcerelchange_Spost1pre{ll,tt,ff}=source_Npost1{ll,tt,ff};
            sourcerelchange_Spost1pre{ll,tt,ff}.avg.pow=source_Spost1{ll,tt,ff}.avg.pow./source_Spre{ll,tt,ff}.avg.pow-1;
            sourcerelchange_Spost2pre{ll,tt,ff}=source_Npost1{ll,tt,ff};
            sourcerelchange_Spost2pre{ll,tt,ff}.avg.pow=source_Spost2{ll,tt,ff}.avg.pow./source_Spre{ll,tt,ff}.avg.pow-1;
            sourcerelchange_Mpost1pre{ll,tt,ff}=source_Npost1{ll,tt,ff};
            sourcerelchange_Mpost1pre{ll,tt,ff}.avg.pow=source_Mpost1{ll,tt,ff}.avg.pow./source_Mpre{ll,tt,ff}.avg.pow-1;
            sourcerelchange_Mpost2pre{ll,tt,ff}=source_Npost1{ll,tt,ff};
            sourcerelchange_Mpost2pre{ll,tt,ff}.avg.pow=source_Mpost2{ll,tt,ff}.avg.pow./source_Mpre{ll,tt,ff}.avg.pow-1;
            sourcerelchange_Upost1pre{ll,tt,ff}=source_Npost1{ll,tt,ff};
            sourcerelchange_Upost1pre{ll,tt,ff}.avg.pow=source_Upost1{ll,tt,ff}.avg.pow./source_Upre{ll,tt,ff}.avg.pow-1;
            sourcerelchange_Upost2pre{ll,tt,ff}=source_Npost1{ll,tt,ff};
            sourcerelchange_Upost2pre{ll,tt,ff}.avg.pow=source_Upost2{ll,tt,ff}.avg.pow./source_Upre{ll,tt,ff}.avg.pow-1;
            
            % each stimulus condition versus null
            sourcerelchange_Spost1Npost1{ll,tt,ff}=source_Npost1{ll,tt,ff};
            sourcerelchange_Spost1Npost1{ll,tt,ff}.avg.pow=source_Spost1{ll,tt,ff}.avg.pow./source_Npost1{ll,tt,ff}.avg.pow-1;
            sourcerelchange_Spost2Npost2{ll,tt,ff}=source_Npost1{ll,tt,ff};
            sourcerelchange_Spost2Npost2{ll,tt,ff}.avg.pow=source_Spost2{ll,tt,ff}.avg.pow./source_Npost2{ll,tt,ff}.avg.pow-1;
            
            % combination for ROI(and time) selection?
            sourcerelchange_Tpost1Spost1Npost1Apost1{ll,tt,ff}=source_Npost1{ll,tt,ff}; % localises T
            sourcerelchange_Tpost1Spost1Npost1Apost1{ll,tt,ff}.avg.pow=(source_Tpost1{ll,tt,ff}.avg.pow+source_Spost1{ll,tt,ff}.avg.pow)./(source_Npost1{ll,tt,ff}.avg.pow+source_Apost1{ll,tt,ff}.avg.pow)-1;
            sourcerelchange_Tpost2Spost2Npost2Apost2{ll,tt,ff}=source_Npost1{ll,tt,ff}; % localises T
            sourcerelchange_Tpost2Spost2Npost2Apost2{ll,tt,ff}.avg.pow=(source_Tpost2{ll,tt,ff}.avg.pow+source_Spost2{ll,tt,ff}.avg.pow)./(source_Npost2{ll,tt,ff}.avg.pow+source_Apost2{ll,tt,ff}.avg.pow)-1;
            sourcerelchange_Apost1Spost1Npost1Tpost1{ll,tt,ff}=source_Npost1{ll,tt,ff}; % localises A
            sourcerelchange_Apost1Spost1Npost1Tpost1{ll,tt,ff}.avg.pow=(source_Apost1{ll,tt,ff}.avg.pow+source_Spost1{ll,tt,ff}.avg.pow)./(source_Npost1{ll,tt,ff}.avg.pow+source_Tpost1{ll,tt,ff}.avg.pow)-1;
            sourcerelchange_Apost2Spost2Npost2Tpost2{ll,tt,ff}=source_Npost1{ll,tt,ff}; % localises A
            sourcerelchange_Apost2Spost2Npost2Tpost2{ll,tt,ff}.avg.pow=(source_Apost2{ll,tt,ff}.avg.pow+source_Spost2{ll,tt,ff}.avg.pow)./(source_Npost2{ll,tt,ff}.avg.pow+source_Tpost2{ll,tt,ff}.avg.pow)-1;
            
            % finally, the interesting contrast
            sourcerelchange_Mpost1Upost1{ll,tt,ff}=source_Npost1{ll,tt,ff};
            sourcerelchange_Mpost1Upost1{ll,tt,ff}.avg.pow=source_Mpost1{ll,tt,ff}.avg.pow./source_Upost1{ll,tt,ff}.avg.pow-1;
            sourcerelchange_Mpost2Upost2{ll,tt,ff}=source_Npost1{ll,tt,ff};
            sourcerelchange_Mpost2Upost2{ll,tt,ff}.avg.pow=source_Mpost2{ll,tt,ff}.avg.pow./source_Upost2{ll,tt,ff}.avg.pow-1;
            
            
            %           end
            
            if plotflag
              if ll==5
                close all
                cfg=[];
                cfg.funparameter='pow';
                ft_sourceplot(cfg,sourcerelchange_Npost1pre{ll,tt,ff});
                ft_sourceplot(cfg,sourcerelchange_Npost2pre{ll,tt,ff});
                ft_sourceplot(cfg,sourcerelchange_Tpost1pre{ll,tt,ff});
                ft_sourceplot(cfg,sourcerelchange_Tpost2pre{ll,tt,ff});
                ft_sourceplot(cfg,sourcerelchange_Apost1pre{ll,tt,ff});
                ft_sourceplot(cfg,sourcerelchange_Apost2pre{ll,tt,ff});
                ft_sourceplot(cfg,sourcerelchange_Tpost1Npost1{ll,tt,ff});
                ft_sourceplot(cfg,sourcerelchange_Tpost2Npost2{ll,tt,ff});
                ft_sourceplot(cfg,sourcerelchange_Apost1Npost1{ll,tt,ff});
                ft_sourceplot(cfg,sourcerelchange_Apost2Npost2{ll,tt,ff});
                
                ft_sourceplot(cfg,sourcerelchange_Spost1pre{ll,tt,ff}); % too top of head? correlated sources?
                ft_sourceplot(cfg,sourcerelchange_Spost2pre{ll,tt,ff});
                ft_sourceplot(cfg,sourcerelchange_Mpost1pre{ll,tt,ff}); % too top of head? correlated sources?
                ft_sourceplot(cfg,sourcerelchange_Mpost2pre{ll,tt,ff});
                ft_sourceplot(cfg,sourcerelchange_Upost1pre{ll,tt,ff});
                ft_sourceplot(cfg,sourcerelchange_Upost2pre{ll,tt,ff});
                ft_sourceplot(cfg,sourcerelchange_Spost1Npost1{ll,tt,ff});
                ft_sourceplot(cfg,sourcerelchange_Spost2Npost2{ll,tt,ff});
                
                ft_sourceplot(cfg,sourcerelchange_Tpost1Spost1Npost1Apost1{ll,tt,ff}); % works
                ft_sourceplot(cfg,sourcerelchange_Tpost2Spost2Npost2Apost2{ll,tt,ff});
                ft_sourceplot(cfg,sourcerelchange_Apost1Spost1Npost1Tpost1{ll,tt,ff}); % doesn't work (center head)
                ft_sourceplot(cfg,sourcerelchange_Apost2Spost2Npost2Tpost2{ll,tt,ff});
                
                ft_sourceplot(cfg,sourcerelchange_Mpost1Upost1{ll,tt,ff});
                ft_sourceplot(cfg,sourcerelchange_Mpost2Upost2{ll,tt,ff});
              end
            end
            
          end % ff
        end % ll
      end % ss
    end % tt
    
    clear source_alltime
    if saveflag
      save([esdir  'sourceTFR_' sub{ii} '_tacaud1_sleep0_supflag' num2str(supflag) '.mat'],'source_*','filter','-v7.3');
    end
    
  end % sleep
end % ii

return

%%   Group level

plotflag=1;
printflag=1;
statsflag=1;
supflag=1;
% audtacflag=0;

soalist=[1 3 4 5 6 7 9];
soadesc={'Aud first by 500ms' '' 'Aud first by 70ms' 'Aud first by 20ms' 'Simultaneous' 'Tac first by 20ms' 'Tac first by 70ms' '' 'Tac first by 500ms'};
load standard_sourcemodel3d8mm
load standard_mri

cd(esdir)
for sleep=[0]
  for tacaud=[1]
    for tt=2
      for ll=soalist
        %       for ll=9
        for ff=1:5
          %         for ff=3:5
          clearvars -except ff ll tt sub *dir ii* sleep *flag soa* ylim* neigh* *basemax tacaud sourcemodel mri
          
          if sleep
            subuseall=iiBuse;
          else
            subuseall=setdiff(iiSuse,[])
          end
          submin=subuseall(1)-1;
          subuseind=0;
          
          for ii=subuseall
            try
              tmp=load([esdir  'sourceTFR_' sub{ii} '_tacaud1_sleep0_supflag' num2str(supflag) '.mat'],'source_*');
            catch
              if supflag==0
                tmp=load(['sourceTFR_' sub{ii} '_tacaud' num2str(tacaud) '_sleep' num2str(sleep) '.mat'],'source*');
              end
            end
            svars=fieldnames(tmp);
            
            % THis is preliminary...how best to included all stages later on?
            if sleep==0
              ssuse=10; % awake
              sleepcond='Awake W';
            elseif sleep==1
              ssuse=12; % this is concatenation of N2 and N3
              sleepcond='Sleep N2';
              %             ssuse=23; % this is concatenation of N2 and N3
              %             sleepcond='Sleep (N2+N3)';
            end
            
            ss=ssuse;
            subuse=subuseall; % reset to all for each sleep stage
            subuseind=subuseind+1;
            
            for vv=1:length(svars)
              tmp.(svars{vv}){ll,tt,ff}.pos=sourcemodel.pos;
              sourceall.(svars{vv}){subuseind}=tmp.(svars{vv}){ll,tt,ff};
            end
            
            % % These are not per ll but do it here anyway for each (yes for ff, tt, ss)
            %             if ll==5
            % each condition's post versus it's own pre
            sourcerelchange.Npost1pre{subuseind}=tmp.source_Npost1{ll,tt,ff};
            sourcerelchange.Npost1pre{subuseind}.avg.pow=tmp.source_Npost1{ll,tt,ff}.avg.pow./tmp.source_Npre{ll,tt,ff}.avg.pow-1;
            sourcerelchange.Npost2pre{subuseind}=tmp.source_Npost1{ll,tt,ff};
            sourcerelchange.Npost2pre{subuseind}.avg.pow=tmp.source_Npost2{ll,tt,ff}.avg.pow./tmp.source_Npre{ll,tt,ff}.avg.pow-1;
            sourcerelchange.Tpost1pre{subuseind}=tmp.source_Npost1{ll,tt,ff};
            sourcerelchange.Tpost1pre{subuseind}.avg.pow=tmp.source_Tpost1{ll,tt,ff}.avg.pow./tmp.source_Tpre{ll,tt,ff}.avg.pow-1;
            sourcerelchange.Tpost2pre{subuseind}=tmp.source_Npost1{ll,tt,ff};
            sourcerelchange.Tpost2pre{subuseind}.avg.pow=tmp.source_Tpost2{ll,tt,ff}.avg.pow./tmp.source_Tpre{ll,tt,ff}.avg.pow-1;
            sourcerelchange.Apost1pre{subuseind}=tmp.source_Npost1{ll,tt,ff};
            sourcerelchange.Apost1pre{subuseind}.avg.pow=tmp.source_Apost1{ll,tt,ff}.avg.pow./tmp.source_Apre{ll,tt,ff}.avg.pow-1;
            sourcerelchange.Apost2pre{subuseind}=tmp.source_Npost1{ll,tt,ff};
            sourcerelchange.Apost2pre{subuseind}.avg.pow=tmp.source_Apost2{ll,tt,ff}.avg.pow./tmp.source_Apre{ll,tt,ff}.avg.pow-1;
            
            % each stimulus condition versus null
            sourcerelchange.Tpost1Npost1{subuseind}=tmp.source_Npost1{ll,tt,ff};
            sourcerelchange.Tpost1Npost1{subuseind}.avg.pow=tmp.source_Tpost1{ll,tt,ff}.avg.pow./tmp.source_Npost1{ll,tt,ff}.avg.pow-1;
            sourcerelchange.Tpost2Npost2{subuseind}=tmp.source_Npost1{ll,tt,ff};
            sourcerelchange.Tpost2Npost2{subuseind}.avg.pow=tmp.source_Tpost2{ll,tt,ff}.avg.pow./tmp.source_Npost2{ll,tt,ff}.avg.pow-1;
            sourcerelchange.Apost1Npost1{subuseind}=tmp.source_Npost1{ll,tt,ff};
            sourcerelchange.Apost1Npost1{subuseind}.avg.pow=tmp.source_Apost1{ll,tt,ff}.avg.pow./tmp.source_Npost1{ll,tt,ff}.avg.pow-1;
            sourcerelchange.Apost2Npost2{subuseind}=tmp.source_Npost1{ll,tt,ff};
            sourcerelchange.Apost2Npost2{subuseind}.avg.pow=tmp.source_Apost2{ll,tt,ff}.avg.pow./tmp.source_Npost2{ll,tt,ff}.avg.pow-1;
            %             end
            
            % % Now the rest are also dependent on ll
            % each condition's post versus it's own pre
            sourcerelchange.Spost1pre{subuseind}=tmp.source_Npost1{ll,tt,ff};
            sourcerelchange.Spost1pre{subuseind}.avg.pow=tmp.source_Spost1{ll,tt,ff}.avg.pow./tmp.source_Spre{ll,tt,ff}.avg.pow-1;
            sourcerelchange.Spost2pre{subuseind}=tmp.source_Npost1{ll,tt,ff};
            sourcerelchange.Spost2pre{subuseind}.avg.pow=tmp.source_Spost2{ll,tt,ff}.avg.pow./tmp.source_Spre{ll,tt,ff}.avg.pow-1;
            sourcerelchange.Mpost1pre{subuseind}=tmp.source_Npost1{ll,tt,ff};
            sourcerelchange.Mpost1pre{subuseind}.avg.pow=tmp.source_Mpost1{ll,tt,ff}.avg.pow./tmp.source_Mpre{ll,tt,ff}.avg.pow-1;
            sourcerelchange.Mpost2pre{subuseind}=tmp.source_Npost1{ll,tt,ff};
            sourcerelchange.Mpost2pre{subuseind}.avg.pow=tmp.source_Mpost2{ll,tt,ff}.avg.pow./tmp.source_Mpre{ll,tt,ff}.avg.pow-1;
            sourcerelchange.Upost1pre{subuseind}=tmp.source_Npost1{ll,tt,ff};
            sourcerelchange.Upost1pre{subuseind}.avg.pow=tmp.source_Upost1{ll,tt,ff}.avg.pow./tmp.source_Upre{ll,tt,ff}.avg.pow-1;
            sourcerelchange.Upost2pre{subuseind}=tmp.source_Npost1{ll,tt,ff};
            sourcerelchange.Upost2pre{subuseind}.avg.pow=tmp.source_Upost2{ll,tt,ff}.avg.pow./tmp.source_Upre{ll,tt,ff}.avg.pow-1;
            
            % each stimulus condition versus null
            sourcerelchange.Spost1Npost1{subuseind}=tmp.source_Npost1{ll,tt,ff};
            sourcerelchange.Spost1Npost1{subuseind}.avg.pow=tmp.source_Spost1{ll,tt,ff}.avg.pow./tmp.source_Npost1{ll,tt,ff}.avg.pow-1;
            sourcerelchange.Spost2Npost2{subuseind}=tmp.source_Npost1{ll,tt,ff};
            sourcerelchange.Spost2Npost2{subuseind}.avg.pow=tmp.source_Spost2{ll,tt,ff}.avg.pow./tmp.source_Npost2{ll,tt,ff}.avg.pow-1;
            
            % combination for ROI(and time) selection?
            sourcerelchange.Tpost1Spost1Npost1Apost1{subuseind}=tmp.source_Npost1{ll,tt,ff}; % localises T
            sourcerelchange.Tpost1Spost1Npost1Apost1{subuseind}.avg.pow=(tmp.source_Tpost1{ll,tt,ff}.avg.pow+tmp.source_Spost1{ll,tt,ff}.avg.pow)./(tmp.source_Npost1{ll,tt,ff}.avg.pow+tmp.source_Apost1{ll,tt,ff}.avg.pow)-1;
            sourcerelchange.Tpost2Spost2Npost2Apost2{subuseind}=tmp.source_Npost1{ll,tt,ff}; % localises T
            sourcerelchange.Tpost2Spost2Npost2Apost2{subuseind}.avg.pow=(tmp.source_Tpost2{ll,tt,ff}.avg.pow+tmp.source_Spost2{ll,tt,ff}.avg.pow)./(tmp.source_Npost2{ll,tt,ff}.avg.pow+tmp.source_Apost2{ll,tt,ff}.avg.pow)-1;
            sourcerelchange.Apost1Spost1Npost1Tpost1{subuseind}=tmp.source_Npost1{ll,tt,ff}; % localises A
            sourcerelchange.Apost1Spost1Npost1Tpost1{subuseind}.avg.pow=(tmp.source_Apost1{ll,tt,ff}.avg.pow+tmp.source_Spost1{ll,tt,ff}.avg.pow)./(tmp.source_Npost1{ll,tt,ff}.avg.pow+tmp.source_Tpost1{ll,tt,ff}.avg.pow)-1;
            sourcerelchange.Apost2Spost2Npost2Tpost2{subuseind}=tmp.source_Npost1{ll,tt,ff}; % localises A
            sourcerelchange.Apost2Spost2Npost2Tpost2{subuseind}.avg.pow=(tmp.source_Apost2{ll,tt,ff}.avg.pow+tmp.source_Spost2{ll,tt,ff}.avg.pow)./(tmp.source_Npost2{ll,tt,ff}.avg.pow+tmp.source_Tpost2{ll,tt,ff}.avg.pow)-1;
            
            % finally, the interesting contrast
            sourcerelchange.Mpost1Upost1{subuseind}=tmp.source_Npost1{ll,tt,ff};
            sourcerelchange.Mpost1Upost1{subuseind}.avg.pow=tmp.source_Mpost1{ll,tt,ff}.avg.pow./tmp.source_Upost1{ll,tt,ff}.avg.pow-1;
            sourcerelchange.Mpost2Upost2{subuseind}=tmp.source_Npost1{ll,tt,ff};
            sourcerelchange.Mpost2Upost2{subuseind}.avg.pow=tmp.source_Mpost2{ll,tt,ff}.avg.pow./tmp.source_Upost2{ll,tt,ff}.avg.pow-1;
          end % ii
          
          freq=sourcerelchange.Mpost2Upost2{subuseind}.freq;
          
          rvars=fieldnames(sourcerelchange);
          
          cfg=[];
          cfg.keepindividual='yes';
          for vv=1:length(svars)
            grind_orig.(svars{vv})=ft_sourcegrandaverage(cfg,sourceall.(svars{vv}){:});
            grind_orig.(svars{vv}).powdimord='rpt_pos';
          end
          for vv=1:length(rvars)
            grind_con.(rvars{vv})=ft_sourcegrandaverage(cfg,sourcerelchange.(rvars{vv}){:});
            grind_con.(rvars{vv}).powdimord='rpt_pos';
          end
          
          cfg=[];
          for vv=1:length(svars)
            grave_orig.(svars{vv})=ft_sourcegrandaverage(cfg,sourceall.(svars{vv}){:});
          end
          close all
          scfg=[];
          scfg.funparameter='pow';
          scfg.funcolorlim=[-.5 .5];
          for vv=1:length(rvars)
            grave_con.(rvars{vv})=ft_sourcegrandaverage(cfg,sourcerelchange.(rvars{vv}){:});
            if plotflag
              ft_sourceplot(scfg,grave_con.(rvars{vv}));
              if printflag
                print(vv,[fdir 'grave_sourcecon_' rvars{vv} '_cond' num2str(ll) num2str(tt) num2str(ss) num2str(tacaud) num2str(supflag) num2str(round(freq)) '.png'],'-dpng')
              end
            end
          end
          
          nsub=size(grind_orig.(svars{1}).pow,1);
          
          if statsflag
            cfg=[];
            cfg.parameter='pow';
            cfg.method='montecarlo';
            cfg.numrandomization=2000;
            % cfg.correctm='holm';
            cfg.correctm='cluster';
            cfg.clusteralpha = 0.05;
            cfg.clusterstatistic = 'maxsum';
            cfg.minnbchan = 2;
            cfg.statistic='depsamplesT';
            cfg.design=zeros(2,2*nsub);
            cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
            cfg.design(2,:)=[1:nsub 1:nsub];
            cfg.ivar=1;
            cfg.uvar=2;
            stat.Npost1pre=ft_sourcestatistics(cfg, grind_orig.source_Npost1, grind_orig.source_Npre);
            stat.Tpost1pre=ft_sourcestatistics(cfg, grind_orig.source_Tpost1, grind_orig.source_Tpre);
            stat.Apost1pre=ft_sourcestatistics(cfg, grind_orig.source_Apost1, grind_orig.source_Apre);
            stat.Spost1pre=ft_sourcestatistics(cfg, grind_orig.source_Spost1, grind_orig.source_Spre);
            stat.Mpost1pre=ft_sourcestatistics(cfg, grind_orig.source_Mpost1, grind_orig.source_Mpre);
            stat.Upost1pre=ft_sourcestatistics(cfg, grind_orig.source_Upost1, grind_orig.source_Upre);
            stat.Npost2pre=ft_sourcestatistics(cfg, grind_orig.source_Npost2, grind_orig.source_Npre);
            stat.Tpost2pre=ft_sourcestatistics(cfg, grind_orig.source_Tpost2, grind_orig.source_Tpre);
            stat.Apost2pre=ft_sourcestatistics(cfg, grind_orig.source_Apost2, grind_orig.source_Apre);
            stat.Spost2pre=ft_sourcestatistics(cfg, grind_orig.source_Spost2, grind_orig.source_Spre);
            stat.Mpost2pre=ft_sourcestatistics(cfg, grind_orig.source_Mpost2, grind_orig.source_Mpre);
            stat.Upost2pre=ft_sourcestatistics(cfg, grind_orig.source_Upost2, grind_orig.source_Upre);
            stat.Tpost1Npost1=ft_sourcestatistics(cfg, grind_orig.source_Tpost1, grind_orig.source_Npost1);
            stat.Tpost2Npost2=ft_sourcestatistics(cfg, grind_orig.source_Tpost2, grind_orig.source_Npost2);
            stat.Apost1Npost1=ft_sourcestatistics(cfg, grind_orig.source_Apost1, grind_orig.source_Npost1);
            stat.Apost2Npost2=ft_sourcestatistics(cfg, grind_orig.source_Apost2, grind_orig.source_Npost2);
            stat.Spost1Npost1=ft_sourcestatistics(cfg, grind_orig.source_Spost1, grind_orig.source_Npost1);
            stat.Spost2Npost2=ft_sourcestatistics(cfg, grind_orig.source_Spost2, grind_orig.source_Npost2);
            stat.Mpost1Upost1=ft_sourcestatistics(cfg, grind_orig.source_Mpost1, grind_orig.source_Upost1);
            stat.Mpost2Upost2=ft_sourcestatistics(cfg, grind_orig.source_Mpost2, grind_orig.source_Upost2);
            
            cfg.design=zeros(2,4*nsub);
            cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub) ones(1,nsub) 2*ones(1,nsub)];
            cfg.design(2,:)=[1:nsub 1:nsub 1:nsub 1:nsub];
            stat.Tpost1Spost1Npost1Apost1=ft_sourcestatistics(cfg, grind_orig.source_Tpost1, grind_orig.source_Npost1, grind_orig.source_Spost1, grind_orig.source_Apost1);
            stat.Tpost2Spost2Npost2Apost2=ft_sourcestatistics(cfg, grind_orig.source_Tpost2, grind_orig.source_Npost2, grind_orig.source_Spost2, grind_orig.source_Apost2);
            stat.Apost1Spost1Npost1Tpost1=ft_sourcestatistics(cfg, grind_orig.source_Apost1, grind_orig.source_Npost1, grind_orig.source_Spost1, grind_orig.source_Tpost1);
            stat.Apost2Spost2Npost2Tpost2=ft_sourcestatistics(cfg, grind_orig.source_Apost2, grind_orig.source_Npost2, grind_orig.source_Spost2, grind_orig.source_Tpost2);
            
            
            if plotflag
              scfg=[];
              scfg.funparameter='pow';
              scfg.maskparameter='mask';
              scfg.funcolorlim=[-.5 .5];
              
              close all
              for vv=1:length(rvars)
                grave_con.(rvars{vv}).mask=stat.(rvars{vv}).mask;
                %                 if any(grave_con.(rvars{vv}).mask)
                %                   ft_sourceplot(scfg,grave_con.(rvars{vv}));
                %                 else
                %                   figure
                %                 end
                %                 if printflag
                %                   print(vv,[fdir 'grave_sourceconMasked_' rvars{vv} '_cond' num2str(ll) num2str(tt) num2str(ss) num2str(tacaud) num2str(round(freq)) '.png'],'-dpng')
                %                 end
              end
              
              close all
              for vv=1:length(rvars)
                if any(grave_con.(rvars{vv}).mask)
                  cfg            = [];
                  cfg.parameter = {'pow' 'mask'};
                  cfg.downsample=2;
                  source_interp  = ft_sourceinterpolate(cfg, grave_con.(rvars{vv}) , mri);
                  cfg=[];
                  cfg.funparameter='pow';
                  cfg.maskparameter='mask';
                  ft_sourceplot(cfg,source_interp);
                else
                  figure
                end
                if printflag
                  print(vv,[fdir 'grave_sourceconMRImasked_' rvars{vv} '_cond' num2str(ll) num2str(tt) num2str(ss) num2str(tacaud) num2str(supflag) num2str(round(freq)) '.png'],'-dpng')
                end
              end
            end
            
          end % statsflag
          
          save([esdir 'statsgrave_source_cond' num2str(ll) num2str(tt) num2str(ss) num2str(tacaud) num2str(supflag) num2str(round(freq)) '.mat'],'stat*','grave_con*');
          
          clear stat gr*
          
        end % ff
      end % ll
    end % tt
  end % tacaud
end % sleep


