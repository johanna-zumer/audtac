% source loc of ERP

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
    root='/home/zumerj/';
  else % assume on VM linux of psychl-132432
    root='/mnt/hgfs/D/';
  end    
    edir=[root 'audtac/eeg_data/'];
    esdir=[root 'audtac/source_data/'];
    ddir=[root 'audtac/legomagic/diaries/'];
    bdir=[root 'audtac/behav_data/'];
    sdir=[root 'audtac/spss_stuff/'];
    fdir=[root 'audtac/figs/'];
    mdir=[root 'audtac/structural_MRI/'];
    pdir=[root 'audtac/polhemus/'];
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
  rmpath(genpath([root 'matlab/spm8/external/fieldtrip/']))
  rmpath(genpath([root 'fieldtrip_svn/']))
  addpath([root 'fieldtrip_svn/'])
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

iiSuse=setdiff(iiSuse,11); % no structural info for e11

cd(esdir)

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

%%
% Skip straight to group level, since saving out individual sourceERP too
% memory-harddrive-intensive.
% Could save out filter first and apply later, but not that huge of
% time-savings.


plotflag=1;
printflag=1;
statsflag=1;
supflag=0;
% audtacflag=0;
soalist=[1 3 4 5 6 7 9];
soadesc={'Aud first by 500ms' '' 'Aud first by 70ms' 'Aud first by 20ms' 'Simultaneous' 'Tac first by 20ms' 'Tac first by 70ms' '' 'Tac first by 500ms'};
load standard_sourcemodel3d8mm
load standard_mri

sourcemodel=ft_convert_units(sourcemodel,'mm');

cd(esdir)

for sleep=[0]
  if sleep
    ssuse=[12 13 23];
  else
    ssuse=10;
  end
  for tt=2
    for ss=ssuse
      for ll=soalist
        clearvars -except ll tt ss* sub *dir ii* sleep *flag soa* *basemax sourcemodel mri *mni filter
        
        if sleep
          subuseall=iiBuse;
        else
          subuseall=setdiff(iiSuse,[]);
        end
        submin=subuseall(1)-1;
        subuseind=0;
        
        for ii=subuseall
%           if ii==18 || ii==24 || ii==25 || ii==26 || ii==28 || || ii==29
%             disp('figure out why has only 62 channels in leadfield')
%             % not just ft_selectdata but also must rereference and recompute cov
%             continue
%           end
          subuseind=subuseind+1;
          cd([edir sub{ii} ])
          load([esdir 'lf_' sub{ii} '.mat']);
%           load([esdir 'elec_' sub{ii} '.mat']);
          load([mdir sub{ii} '/elec.mat']);
          load([mdir sub{ii} '/vol.mat']);
          load([mdir sub{ii} '/mrinorm.mat']);
          
          hglsub=ft_warp_apply(mrinormwarp.params, hglmni, 'sn2individual');
          hgrsub=ft_warp_apply(mrinormwarp.params, hgrmni, 'sn2individual');
          stlsub=ft_warp_apply(mrinormwarp.params, stlmni, 'sn2individual');
          strsub=ft_warp_apply(mrinormwarp.params, strmni, 'sn2individual');
          
          tmp=load(['tlock_diffs_averef_' sub{ii} '_sleep0']);
          
          
          %           cfg=[];
          %           cfg.grid=lf;
          %           cfg.vol=voldipolimm;
          %           cfg.elec=elec;
          %           cfg.method='mne';
          %           cfg.mne.snr=1;
          %           cfg.mne.projectmom=1;
          %           erp_mne{ll,tt,ss}=ft_sourceanalysis(cfg,tmp.tlock_tAudAlone{ll,tt,ss});
          %           erp_mne{ll,tt,ss}=mat_your_mom(erp_mne{ll,tt,ss});
          %
          %           cfg=[];
          %           cfg.funparameter='mommat';
          %           cfg.funcolorlim=[-1e-35 1e-35];
          %           ft_sourceplot(cfg,erp_mne{ll,tt,ss});
          
          
          tlockN=tmp.tlock_tNulAlone{ll,tt,ss};
          tlockT=tmp.tlock_tTacAlone{ll,tt,ss};
          tlockA=tmp.tlock_tAudAlone{ll,tt,ss};
          tlockS=tmp.tlock_tMSAlone{ll,tt,ss};
          tlockU=tmp.tlock_tacPaud{ll,tt,ss};
          tlockM=tmp.tlock_tacMSpN{ll,tt,ss};
          clear tmp
          
          cfg=[];
          cfg.latency=[-0.6 1.4];
          % % NOTE: the .cov is still from the timewindow specified in eeg_legomagic_erp_stats2_sepTacAud.m
          tlockN=ft_selectdata(cfg,tlockN);
          tlockT=ft_selectdata(cfg,tlockT);
          tlockA=ft_selectdata(cfg,tlockA);
          tlockS=ft_selectdata(cfg,tlockS);
          tlockU=ft_selectdata(cfg,tlockU);
          tlockM=ft_selectdata(cfg,tlockM);
          
          cfg=[];
          cfg.parameter='cov';
          tlockall=ft_timelockgrandaverage(cfg,tlockN,tlockT,tlockA,tlockS,tlockU,tlockM);
          tlockall.cov=tlockall.avg;
          tlockall=rmfield(tlockall,'avg');
          tlockall=rmfield(tlockall,'var');
          %             tlockall.dimord='chan_chan';
          tlockall.avg=zeros(length(tlockall.label),1);
          
          
          cfg=[];
          cfg.grid=lf;
          cfg.vol=voldipolimm;
          cfg.elec=elec_use;
          cfg.keepfilter='yes';
          cfg.inwardshift=10;
          
          cfg.method='lcmv';
          cfg.lcmv.lambda='5%';
          cfg.lcmv.projectmom=1;
          erpall_lcmv=ft_sourceanalysis(cfg,tlockall);
          filter_lcmv{ll,tt,ss}=erpall_lcmv.filter;
          
          cfg.method='mne';
          keyboard % tweak this snr!?
          cfg.mne.snr=1;
          cfg.mne.projectmom=1;
          erpall_mne=ft_sourceanalysis(cfg,tlockall);
          filter_mne{ll,tt,ss}=erpall_mne.filter;
          
          if plotflag
            if ll==1
              tlockcm=[];
              tlockcm.time=[1 2 3];
              tlockcm.avg=erpall_lcmv.filter{dsearchn(erpall_lcmv.pos,mean(hgrsub))}';
              tlockcm.label=erpall_lcmv.cfg.channel;
              tlockcm.dimord='chan_time';
              cfg=[];
              cfg.layout='elec1010.lay';
              cfg.xlim=[0.9 1.1];
              figure;
              ft_topoplotER(cfg,tlockcm);

              tlockcm=[];
              tlockcm.time=[1 2 3];
              tlockcm.avg=erpall_mne.filter{dsearchn(erpall_mne.pos,mean(hgrsub))}';
              tlockcm.label=erpall_mne.cfg.channel;
              tlockcm.dimord='chan_time';
              cfg=[];
              cfg.layout='elec1010.lay';
              cfg.xlim=[0.9 1.1];
              figure;
              ft_topoplotER(cfg,tlockcm);
              %               cfg.xlim=[1.9 2.1];
              %               figure;
              %               ft_topoplotER(cfg,tlockcm);
              %               cfg.xlim=[2.9 3.1];
              %               figure;
              %               ft_topoplotER(cfg,tlockcm);
            end
          end
          
          
          
          %           % original covariance, no regularisation
          %           erp_lcmvc=ft_sourceanalysis(cfg,tmp.tlock_tAudAlone{ll,tt,ss});
          % %           mineigUnave=min(eig(tmp.tlock_tAudAlone{ll,tt,ss}.cov));
          % %           tmpnocov=tmp.tlock_tAudAlone{ll,tt,ss};
          % %           tmpnocov=rmfield(tmpnocov,'cov');
          % %           tmpavgcov=tmpnocov;
          % %           tmpavgcov.cov=cov(tmpavgcov.avg');
          % %           cfg.lambda='15%';
          % %           % original covariance, 15% regularisation
          % %           erp_lcmvcr=ft_sourceanalysis(cfg,tmp.tlock_tAudAlone{ll,tt,ss});
          % %           cfg.lambda=mineigUnave;
          % %           % original covariance, mineigUnave regularisation
          % %           erp_lcmvcme=ft_sourceanalysis(cfg,tmp.tlock_tAudAlone{ll,tt,ss});
          % % %           cfg.lambda=0;
          % % %           % no covariance (identity)
          % % %           erp_lcmvnc=ft_sourceanalysis(cfg,tmpnocov);
          % % %           % covariance of the average, no reg
          % % %           erp_lcmvac=ft_sourceanalysis(cfg,tmpavgcov);
          % %           cfg.lambda='15%';
          % %           % covariance of the average, 15% reg
          % %           erp_lcmvacr=ft_sourceanalysis(cfg,tmpavgcov);
          % % %           cfg.lambda=mineigUnave;
          % % %           % covariance of the average, mineigUnave reg
          % % %           erp_lcmvacme=ft_sourceanalysis(cfg,tmpavgcov);
          %
          %           erp_lcmvc=mat_your_mom(erp_lcmvc);
          % %           erp_lcmvacr=mat_your_mom(erp_lcmvacr);
          % %           erp_lcmvcr=mat_your_mom(erp_lcmvcr);
          % %           erp_lcmvcme=mat_your_mom(erp_lcmvcme);
          % %           erp_lcmvnc=mat_your_mom(erp_lcmvnc);
          % %           erp_lcmvac=mat_your_mom(erp_lcmvac);
          % %           erp_lcmvacme=mat_your_mom(erp_lcmvacme);
          %
          %           cfg=[];
          %           cfg.funparameter='mommat';
          %           cfg.funcolorlim=[-3e4 3e4];
          %           ft_sourceplot(cfg,erp_lcmvc) % best in terms of space/timecourse tradeoff.
          % %           ft_sourceplot(cfg,erp_lcmvcme) % same as lcmvc
          % %           ft_sourceplot(cfg,erp_lcmvacr) % more widespread
          % %           ft_sourceplot(cfg,erp_lcmvcr) % even more widespread
          % %           ft_sourceplot(cfg,erp_lcmvnc)
          % %           ft_sourceplot(cfg,erp_lcmvac)
          % %           ft_sourceplot(cfg,erp_lcmvacme)
          
          %           cfg=[];
          %           cfg.funparameter='pow';
          %           ft_sourcemovie(cfg,erp_lcmvc{ll,tt,ss})
          
          %         cfg=[];
          %         cfg.grid=lf;
          %         cfg.vol=voldipolimm;
          %         cfg.elec=elec;
          %         cfg.gridsearch='no';
          %         cfg.dip.pos=[ -61.1271  -90.8194  -40.8452]
          % %         cfg.numdipoles=2;
          % %         cfg.symmetry = 'y';
          %         erp_dip{ll,tt,ss}=ft_dipolefitting(cfg,tmp.tlock_tAudAlone{ll,tt,ss});
          
          cfg=[];
          cfg.grid=lf;
          cfg.grid.filter=filter_lcmv{ll,tt,ss};
          cfg.vol=voldipolimm;
          cfg.elec=elec_use;
          
          cfg.method='lcmv';
          cfg.lcmv.projectmom=1;
          cfg.lcmv.feedback='none';          
          if 0
          
          source_N{subuseind}=ft_sourceanalysis(cfg,tlockN);
          source_T{subuseind}=ft_sourceanalysis(cfg,tlockT);
          source_A{subuseind}=ft_sourceanalysis(cfg,tlockA);
          source_S{subuseind}=ft_sourceanalysis(cfg,tlockS);
          source_U{subuseind}=ft_sourceanalysis(cfg,tlockU);
          source_M{subuseind}=ft_sourceanalysis(cfg,tlockM);
          source_N{subuseind}=verwijder_onnodige_velden(source_N{subuseind},0);
          source_T{subuseind}=verwijder_onnodige_velden(source_T{subuseind},0);
          source_A{subuseind}=verwijder_onnodige_velden(source_A{subuseind},0);
          source_S{subuseind}=verwijder_onnodige_velden(source_S{subuseind},0);
          source_U{subuseind}=verwijder_onnodige_velden(source_U{subuseind},0);
          source_M{subuseind}=verwijder_onnodige_velden(source_M{subuseind},0);
          source_N{subuseind}=mat_your_mom(source_N{subuseind});
          source_T{subuseind}=mat_your_mom(source_T{subuseind});
          source_A{subuseind}=mat_your_mom(source_A{subuseind});
          source_S{subuseind}=mat_your_mom(source_S{subuseind});
          source_U{subuseind}=mat_your_mom(source_U{subuseind});
          source_M{subuseind}=mat_your_mom(source_M{subuseind});
          source_N{subuseind}.pos=sourcemodel.pos;
          source_T{subuseind}.pos=sourcemodel.pos;
          source_A{subuseind}.pos=sourcemodel.pos;
          source_S{subuseind}.pos=sourcemodel.pos;
          source_U{subuseind}.pos=sourcemodel.pos;
          source_M{subuseind}.pos=sourcemodel.pos;          
%           source_N{subuseind}=tidy_up_grind(source_N{subuseind});
%           source_T{subuseind}=tidy_up_grind(source_T{subuseind});
%           source_A{subuseind}=tidy_up_grind(source_A{subuseind});
%           source_S{subuseind}=tidy_up_grind(source_S{subuseind});
%           source_U{subuseind}=tidy_up_grind(source_U{subuseind});
%           source_M{subuseind}=tidy_up_grind(source_M{subuseind});
          end
          
          cfg.method='mne';
          cfg.mne.projectmom=1;
          cfg.grid.filter=filter_mne{ll,tt,ss};
          source_Nm{subuseind}=ft_sourceanalysis(cfg,tlockN);
          source_Tm{subuseind}=ft_sourceanalysis(cfg,tlockT);
          source_Am{subuseind}=ft_sourceanalysis(cfg,tlockA);
          source_Sm{subuseind}=ft_sourceanalysis(cfg,tlockS);
          source_Um{subuseind}=ft_sourceanalysis(cfg,tlockU);
          source_Mm{subuseind}=ft_sourceanalysis(cfg,tlockM);
          source_Nm{subuseind}=verwijder_onnodige_velden(source_Nm{subuseind},0);
          source_Tm{subuseind}=verwijder_onnodige_velden(source_Tm{subuseind},0);
          source_Am{subuseind}=verwijder_onnodige_velden(source_Am{subuseind},0);
          source_Sm{subuseind}=verwijder_onnodige_velden(source_Sm{subuseind},0);
          source_Um{subuseind}=verwijder_onnodige_velden(source_Um{subuseind},0);
          source_Mm{subuseind}=verwijder_onnodige_velden(source_Mm{subuseind},0);
          source_Nm{subuseind}=mat_your_3mom(source_Nm{subuseind});
          source_Tm{subuseind}=mat_your_3mom(source_Tm{subuseind});
          source_Am{subuseind}=mat_your_3mom(source_Am{subuseind});
          source_Sm{subuseind}=mat_your_3mom(source_Sm{subuseind});
          source_Um{subuseind}=mat_your_3mom(source_Um{subuseind});
          source_Mm{subuseind}=mat_your_3mom(source_Mm{subuseind});
          source_Nm{subuseind}.pos=sourcemodel.pos;
          source_Tm{subuseind}.pos=sourcemodel.pos;
          source_Am{subuseind}.pos=sourcemodel.pos;
          source_Sm{subuseind}.pos=sourcemodel.pos;
          source_Um{subuseind}.pos=sourcemodel.pos;
          source_Mm{subuseind}.pos=sourcemodel.pos;          
%           source_N{subuseind}=tidy_up_grind(source_N{subuseind});
%           source_T{subuseind}=tidy_up_grind(source_T{subuseind});
%           source_A{subuseind}=tidy_up_grind(source_A{subuseind});
%           source_S{subuseind}=tidy_up_grind(source_S{subuseind});
%           source_U{subuseind}=tidy_up_grind(source_U{subuseind});
%           source_M{subuseind}=tidy_up_grind(source_M{subuseind});
          
        end % ii
        
        if 0
        cfg=[];
        cfg.parameter='mommat';
        grave_all_N=ft_sourcegrandaverage(cfg,source_N{:});
        grave_all_T=ft_sourcegrandaverage(cfg,source_T{:});
        grave_all_A=ft_sourcegrandaverage(cfg,source_A{:});
        grave_all_S=ft_sourcegrandaverage(cfg,source_S{:});
        grave_all_U=ft_sourcegrandaverage(cfg,source_U{:});
        grave_all_M=ft_sourcegrandaverage(cfg,source_M{:});
        end
        
        cfg=[];
        cfg.parameter='mommat';
        grave_all_Nm=ft_sourcegrandaverage(cfg,source_Nm{:});
        grave_all_Tm=ft_sourcegrandaverage(cfg,source_Tm{:});
        grave_all_Am=ft_sourcegrandaverage(cfg,source_Am{:});
        grave_all_Sm=ft_sourcegrandaverage(cfg,source_Sm{:});
        grave_all_Um=ft_sourcegrandaverage(cfg,source_Um{:});
        grave_all_Mm=ft_sourcegrandaverage(cfg,source_Mm{:});

        cfg=[];
        cfg.parameter='mommat';
        cfg.keepindividual='yes';
        if 0
        grind_all_N=ft_sourcegrandaverage(cfg,source_N{:});
        clear source_N
        grind_all_T=ft_sourcegrandaverage(cfg,source_T{:});
        clear source_T
        grind_all_A=ft_sourcegrandaverage(cfg,source_A{:});
        clear source_A
        grind_all_S=ft_sourcegrandaverage(cfg,source_S{:});
        clear source_S
        grind_all_U=ft_sourcegrandaverage(cfg,source_U{:});
        clear source_U
        grind_all_M=ft_sourcegrandaverage(cfg,source_M{:});
        clear source_M
        end
        grind_all_N=ft_sourcegrandaverage(cfg,source_Nm{:});
        clear source_Nm
        grind_all_T=ft_sourcegrandaverage(cfg,source_Tm{:});
        clear source_Tm
        grind_all_A=ft_sourcegrandaverage(cfg,source_Am{:});
        clear source_Am
        grind_all_S=ft_sourcegrandaverage(cfg,source_Sm{:});
        clear source_Sm
        grind_all_U=ft_sourcegrandaverage(cfg,source_Um{:});
        clear source_Um
        grind_all_M=ft_sourcegrandaverage(cfg,source_Mm{:});
        clear source_Mm
        
        if 0
        grind_all_N.mommatdimord='rpt_pos_time';
        grind_all_T.mommatdimord='rpt_pos_time';
        grind_all_A.mommatdimord='rpt_pos_time';
        grind_all_S.mommatdimord='rpt_pos_time';
        grind_all_U.mommatdimord='rpt_pos_time';
        grind_all_M.mommatdimord='rpt_pos_time';
        end
        grind_all_Nm.mommatdimord='rpt_pos_time';
        grind_all_Tm.mommatdimord='rpt_pos_time';
        grind_all_Am.mommatdimord='rpt_pos_time';
        grind_all_Sm.mommatdimord='rpt_pos_time';
        grind_all_Um.mommatdimord='rpt_pos_time';
        grind_all_Mm.mommatdimord='rpt_pos_time';
        
%         grind_all_N=tidy_up_grind(grind_all_N);
%         grind_all_T=tidy_up_grind(grind_all_T);
%         grind_all_A=tidy_up_grind(grind_all_A);
%         grind_all_S=tidy_up_grind(grind_all_S);
%         grind_all_U=tidy_up_grind(grind_all_U);
%         grind_all_M=tidy_up_grind(grind_all_M);
        
        
        if plotflag
          cfg=[];
          cfg.funparameter='mommat';
          if 0
          ft_sourceplot(cfg,grave_all_N);
          ft_sourceplot(cfg,grave_all_T);
          ft_sourceplot(cfg,grave_all_A);
          ft_sourceplot(cfg,grave_all_S);
          ft_sourceplot(cfg,grave_all_U);
          ft_sourceplot(cfg,grave_all_M);
          end
          ft_sourceplot(cfg,grave_all_Nm);
          ft_sourceplot(cfg,grave_all_Tm);
          ft_sourceplot(cfg,grave_all_Am);
          ft_sourceplot(cfg,grave_all_Sm);
          ft_sourceplot(cfg,grave_all_Um);
          ft_sourceplot(cfg,grave_all_Mm);
        end % plotflag
        
        
        keyboard
        nsub=size(grind_all_N.mommat,1);
        if statsflag
          cfg=[];
          cfg.parameter='mommat';
          cfg.method='analytic';
          %           cfg.correctm='holm';
          cfg.minnbchan = 2;
          cfg.alpha=.005; % loose on purpose just to create vague mask
          cfg.statistic='depsamplesT';
          cfg.design=zeros(2,2*nsub);
          cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
          cfg.design(2,:)=[1:nsub 1:nsub];
          cfg.ivar=1;
          cfg.uvar=2;
          stat.TvsN=ft_sourcestatistics(cfg, grind_all_T, grind_all_N);
          stat.AvsN=ft_sourcestatistics(cfg, grind_all_A, grind_all_N);
          stat.SvsN=ft_sourcestatistics(cfg, grind_all_S, grind_all_N);
          cfg.design=zeros(2,4*nsub);
          cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub) ones(1,nsub) 2*ones(1,nsub)];
          cfg.design(2,:)=[1:nsub 1:nsub 1:nsub 1:nsub];
          stat.TSvsNA=ft_sourcestatistics(cfg, grind_all_T, grind_all_N, grind_all_S, grind_all_A);
          stat.ASvsNT=ft_sourcestatistics(cfg, grind_all_A, grind_all_N, grind_all_S, grind_all_T);
          
          % get tac ROIs based on unisensory
%           [mx,tactimeind]=max(mean(stat.TvsN.mask));
          tactimeind=dsearchn(grave_all_Nm.time',0.08);
          meanpostac=mean(stat.TvsN.pos(find(all([stat.TvsN.mask(:,tactimeind-13:tactimeind+13)]')),:));
          roitac=find(nut_rownorm(nut_coord_diff(grind_all_N.pos,meanpostac))<10);
          
          audtimeind=dsearchn(grave_all_N.time',0.08);
          meanposaud=mean(stat.AvsN.pos(find(all([stat.AvsN.mask(:,audtimeind-13:audtimeind+13)]')),:));
          roitac=find(nut_rownorm(nut_coord_diff(grind_all_N.pos,meanpostac))<10);

          cfg=[];
          cfg.coordinate=grind_all_N.pos(roitac,:)+.1;
          %         cfg.avgoverpos='yes';
          grind_all_N_roitac=ft_selectdata(cfg,grind_all_N);
          grind_all_T_roitac=ft_selectdata(cfg,grind_all_T);
          grind_all_A_roitac=ft_selectdata(cfg,grind_all_A);
          grind_all_S_roitac=ft_selectdata(cfg,grind_all_S);
          grind_all_U_roitac=ft_selectdata(cfg,grind_all_U);
          grind_all_M_roitac=ft_selectdata(cfg,grind_all_M);
          
          % this doesn't work. don't know why.
          %           figure;plot(stat.TvsN.time,mean(stat.TvsN.mask & stat.SvsN.mask))
          %           figure;plot(stat.AvsN.time,mean(stat.AvsN.mask & stat.SvsN.mask))
          %           [mx,tind]=max(mean(stat.TvsN.mask & stat.SvsN.mask))
          %           [mx,aind]=max(mean(stat.AvsN.mask & stat.SvsN.mask))
          %           clear tmp
          %           tmp=stat.TvsN;
          %           tmp.mommat=double(repmat(stat.TvsN.mask(:,tind),[1 length(stat.TvsN.time)]));
          %           cfg=[];
          %           cfg.funparameter='mommat';
          %           ft_sourceplot(cfg,tmp);
          
          
          
          cfg=[];
          cfg.parameter='mommat';
          cfg.method='montecarlo';
          cfg.numrandomization=100;
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
          statmc.TvsN=ft_sourcestatistics(cfg, grind_all_T_roitac, grind_all_N_roitac);
          statmc.TvsN=ft_sourcestatistics(cfg, grind_all_T, grind_all_N);
          statmc.AvsN=ft_sourcestatistics(cfg, grind_all_A, grind_all_N);
          statmc.SvsN=ft_sourcestatistics(cfg, grind_all_S, grind_all_N);
          statmc.UvsM=ft_sourcestatistics(cfg, grind_all_U, grind_all_M);
          cfg.design=zeros(2,4*nsub);
          cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub) ones(1,nsub) 2*ones(1,nsub)];
          cfg.design(2,:)=[1:nsub 1:nsub 1:nsub 1:nsub];
          statmc.TSvsNA=ft_sourcestatistics(cfg, grind_all_T, grind_all_N, grind_all_S, grind_all_A);
          statmc.ASvsNT=ft_sourcestatistics(cfg, grind_all_A, grind_all_N, grind_all_S, grind_all_T);
          
          
        end % statsflag
        
      end % ll
    end % ss
  end % tt
  
  if saveflag
    %       save([esdir  'sourceERP_' sub{ii} '_tacaud1_sleep' num2str(sleep) '_supflag' num2str(supflag) '.mat'],'source_*','filter','-v7.3');
    save([esdir  'filterERP_' sub{ii} '_tacaud1_sleep' num2str(sleep) '_supflag' num2str(supflag) '.mat'],'filter','-v7.3');
  end
  
end % sleep

return

%% Group level




cfg            = [];
cfg.downsample = 2;
cfg.parameter = 'avg.pow';
erp_source  = ft_sourceinterpolate(cfg, erp_mne{ll,tt,ss} , mrinormlin);

cfg=[];
cfg.funparameter='avg.pow';
ft_sourceplot(cfg,erp_source);

