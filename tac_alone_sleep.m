
%% Group level: awake and asleep testing tactile alone response
% This is not for the main finding, but rather testing strange
% finding of no seeming Tactile response for N23 bed data.

plotflag=1;
printflag=1;
statsflag=0;
tacaud=1; % tacaud=1 means triggered on tactile; tacaud=0 means triggered on auditory

tt=3;
sleepss_tacAll=zeros(13,11,max([iiSuse iiBuse]));
sleepss_tac19T=zeros(13,11,max([iiSuse iiBuse]));
sleepss_tacAlone=zeros(13,11,max([iiSuse iiBuse]));
for sleep=[1]
  clearvars -except tt sub edir ddir ii*use *flag sleep* tacaud
  if sleep==1
    subuse=setdiff(iiBuse,3:7);
    iter=11;
    trialkc=-1; % -1 all, 0 no Kc, 1 only Kc
  elseif sleep==0
    subuse=iiSuse;
    iter=27;
  end
  subind=1;
  for ii=subuse
    cd([edir sub{ii} ])
    %   load(['tlock_diffs_' sub{ii} '.mat']);
    %       load(['tlock_diffs_averef_' sub{ii} '.mat']);
    %     load(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '.mat'])
    if sleep
      load(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_trialkc' num2str(trialkc) '.mat'])
      tk0=load(['trialkept_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(0) '_trialkc' num2str(trialkc) '.mat']);
      tk1=load(['trialkept_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '_trialkc' num2str(trialkc) '.mat']);
      ssuse=tk1.tr.stageuse;
    else
      load(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud1' '_iter' num2str(iter) '.mat'])
      tk1=load(['trialkept_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '_iter' num2str(iter) '.mat']);
      ssuse=tk1.tr.stageuse;
    end
    %         load(['trialkept_tt' num2str(tt) '_sleep' num2str(sleep) '.mat'],'tr')
    
    for ss=ssuse+10
      try
        sleepss_tacAll(ss,sleep+10,subind)=tlock_tacAll{tt,ss}.dof(1,dsearchn(tlock_tacAll{tt,ss}.time',0));
        %         if sleepss_tacAll(ss,sleep+10,subind)>=20
        tlock_tacAll_each{ss,sleep+10,subind}=tlock_tacAll{tt,ss};
        %         end
      end
      try
        sleepss_tac19T(ss,sleep+10,subind)=tlock_tac19T{tt,ss}.dof(1,dsearchn(tlock_tac19T{tt,ss}.time',0));
        %         if sleepss_tac19T(ss,sleep+10,subind)>=20
        tlock_tac19T_each{ss,sleep+10,subind}=tlock_tac19T{tt,ss};
        %         end
      end
      try
        sleepss_tacAlone(ss,sleep+10,subind)=tlock_tTacAlone{5,tt,ss}.dof(1,dsearchn(tlock_tTacAlone{5,tt,ss}.time',0));
        %         if sleepss_tacAlone(ss,sleep+10,subind)>=20
        tlock_tacAlone_each{ss,sleep+10,subind}=tlock_tTacAlone{5,tt,ss}; % the 'll' is arbitrary here
        %         end
      end
    end % ss
    
    clear tlock_tacAll tlock_tac19T tlock_tTacAlone tlock_a* tlock_tacMSpN tlock_tacPaud tlock_t*Alone
  end % sleep
  % if any(any(sleepss_tacAll(:,:,subind))) || any(any(sleepss_tac19T(:,:,subind))) || any(any(sleepss_Alone(:,:,subind)))
  subind=subind+1;
end % ii
subind=subind-1;

thresh=10;
for ii=1:subind
  sleep=0;
  ss=10;
  if sleepss_tacAll(ss,sleep+10,ii)>=thresh
    tacAll_s0{ii}=tlock_tacAll_each{ss,sleep+10,ii};
  else
    tacAll_s0{ii}=[];
  end
  if sleepss_tac19T(ss,sleep+10,ii)>=thresh
    tac19T_s0{ii}=tlock_tac19T_each{ss,sleep+10,ii};
  else
    tac19T_s0{ii}=[];
  end
  if sleepss_tacAlone(ss,sleep+10,ii)>=thresh
    tacAlone_s0{ii}=tlock_tacAlone_each{ss,sleep+10,ii};
  else
    tacAlone_s0{ii}=[]
  end
  
  sleep=0;
  ss=11;
  if sleepss_tacAll(ss,sleep+10,ii)>=thresh
    tacAll_s1{ii}=tlock_tacAll_each{ss,sleep+10,ii};
  else
    tacAll_s1{ii}=[];
  end
  if sleepss_tac19T(ss,sleep+10,ii)>=thresh
    tac19T_s1{ii}=tlock_tac19T_each{ss,sleep+10,ii};
  else
    tac19T_s1{ii}=[];
  end
  if sleepss_tacAlone(ss,sleep+10,ii)>=thresh
    tacAlone_s1{ii}=tlock_tacAlone_each{ss,sleep+10,ii};
  else
    tacAlone_s1{ii}=[];
  end
  
  sleep=0;
  ss=12;
  if sleepss_tacAll(ss,sleep+10,ii)>=thresh
    tacAll_s2{ii}=tlock_tacAll_each{ss,sleep+10,ii};
  else
    tacAll_s2{ii}=[];
  end
  if sleepss_tac19T(ss,sleep+10,ii)>=thresh
    tac19T_s2{ii}=tlock_tac19T_each{ss,sleep+10,ii};
  else
    tac19T_s2{ii}=[];
  end
  if sleepss_tacAlone(ss,sleep+10,ii)>=thresh
    tacAlone_s2{ii}=tlock_tacAlone_each{ss,sleep+10,ii};
  else
    tacAlone_s2{ii}=[];
  end
  
  sleep=1;
  ss=10;
  if sleepss_tacAll(ss,sleep+10,ii)>=thresh
    tacAll_b0{ii}=tlock_tacAll_each{ss,sleep+10,ii};
  else
    tacAll_b0{ii}=[];
  end
  if sleepss_tac19T(ss,sleep+10,ii)>=thresh
    tac19T_b0{ii}=tlock_tac19T_each{ss,sleep+10,ii};
  else
    tac19T_b0{ii}=[];
  end
  if sleepss_tacAlone(ss,sleep+10,ii)>=thresh
    tacAlone_b0{ii}=tlock_tacAlone_each{ss,sleep+10,ii};
  else
    tacAlone_b0{ii}=[];
  end
  
  sleep=1;
  ss=11;
  if sleepss_tacAll(ss,sleep+10,ii)>=thresh
    tacAll_b1{ii}=tlock_tacAll_each{ss,sleep+10,ii};
  else
    tacAll_b1{ii}=[];
  end
  if sleepss_tac19T(ss,sleep+10,ii)>=thresh
    tac19T_b1{ii}=tlock_tac19T_each{ss,sleep+10,ii};
  else
    tac19T_b1{ii}=[];
  end
  if sleepss_tacAlone(ss,sleep+10,ii)>=thresh
    tacAlone_b1{ii}=tlock_tacAlone_each{ss,sleep+10,ii};
  else
    tacAlone_b1{ii}=[];
  end
  
  sleep=1;
  ss=12;
  if sleepss_tacAll(ss,sleep+10,ii)>=thresh
    tacAll_b2{ii}=tlock_tacAll_each{ss,sleep+10,ii};
  else
    tacAll_b2{ii}=[];
  end
  if sleepss_tac19T(ss,sleep+10,ii)>=thresh
    tac19T_b2{ii}=tlock_tac19T_each{ss,sleep+10,ii};
  else
    tac19T_b2{ii}=[];
  end
  if sleepss_tacAlone(ss,sleep+10,ii)>=thresh
    tacAlone_b2{ii}=tlock_tacAlone_each{ss,sleep+10,ii};
  else
    tacAlone_b2{ii}=[];
  end
  
  sleep=1;
  ss=13;
  if sleepss_tacAll(ss,sleep+10,ii)>=thresh
    tacAll_b3{ii}=tlock_tacAll_each{ss,sleep+10,ii};
  else
    tacAll_b3{ii}=[];
  end
  if sleepss_tac19T(ss,sleep+10,ii)>=thresh
    tac19T_b3{ii}=tlock_tac19T_each{ss,sleep+10,ii};
  else
    tac19T_b3{ii}=[];
  end
  if sleepss_tacAlone(ss,sleep+10,ii)>=thresh
    tacAlone_b3{ii}=tlock_tacAlone_each{ss,sleep+10,ii};
  else
    tacAlone_b3{ii}=[];
  end
end % ii


% insert these into lines below.
figure(100);
figcfg=[];
figcfg.xlim=[-0.7 1.1];
figcfg.channel={'Fz' 'Cz' 'F1' 'F2' 'FC1' 'FC2' 'C1' 'C2'};
figcfg.ylim=[-7 7];

cfg=[];
tmp={tacAll_s0{squeeze(sleepss_tacAll(10,10,:))>=thresh}};
if ~isempty(tmp)
  grave_tacAll_s0=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,1);ft_singleplotER(figcfg,grave_tacAll_s0);
  ylabel('TacAll')
  title(num2str(sum(squeeze(sleepss_tacAll(10,10,:)))))
else
end
cfg=[];
tmp={tacAll_s1{squeeze(sleepss_tacAll(11,10,:))>=thresh}};
if ~isempty(tmp)
  grave_tacAll_s1=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,3);ft_singleplotER(figcfg,grave_tacAll_s1);
  title(num2str(sum(squeeze(sleepss_tacAll(11,10,:)))))
else
end
cfg=[];
tmp={tacAll_s2{squeeze(sleepss_tacAll(12,10,:))>=thresh}};
if ~isempty(tmp)
  grave_tacAll_s2=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,5);ft_singleplotER(figcfg,grave_tacAll_s2);
  title(num2str(sum(squeeze(sleepss_tacAll(12,10,:)))))
else
end
% cfg=[];
% tmp={tacAll_s3{squeeze(sleepss_tacAll(13,10,:))>=thresh}};
% grave_tacAll_s3=ft_timelockgrandaverage(cfg,tmp{:});
% subplot(3,8,7);ft_singleplotER(figcfg,grave_tacAll_s3);

cfg=[];
tmp={tacAll_b0{squeeze(sleepss_tacAll(10,11,:))>=thresh}};
if ~isempty(tmp)
  grave_tacAll_b0=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,2);ft_singleplotER(figcfg,grave_tacAll_b0);
  title(num2str(sum(squeeze(sleepss_tacAll(10,11,:)))))
else
end
cfg=[];
tmp={tacAll_b1{squeeze(sleepss_tacAll(11,11,:))>=thresh}};
if ~isempty(tmp)
  grave_tacAll_b1=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,4);ft_singleplotER(figcfg,grave_tacAll_b1);
  title(num2str(sum(squeeze(sleepss_tacAll(11,11,:)))))
else
end
cfg=[];
tmp={tacAll_b2{squeeze(sleepss_tacAll(12,11,:))>=thresh}};
if ~isempty(tmp)
  grave_tacAll_b2=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,6);ft_singleplotER(figcfg,grave_tacAll_b2);
  title(num2str(sum(squeeze(sleepss_tacAll(12,11,:)))))
else
end
cfg=[];
tmp={tacAll_b3{squeeze(sleepss_tacAll(13,11,:))>=thresh}};
if ~isempty(tmp)
  grave_tacAll_b3=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,8);ft_singleplotER(figcfg,grave_tacAll_b3);
  title(num2str(sum(squeeze(sleepss_tacAll(13,11,:)))))
else
end

cfg=[];
tmp={tac19T_s0{squeeze(sleepss_tac19T(10,10,:))>=thresh}};
if ~isempty(tmp)
  grave_tac19T_s0=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,9); ft_singleplotER(figcfg,grave_tac19T_s0);
  title(num2str(sum(squeeze(sleepss_tac19T(10,10,:)))))
  ylabel('Tac PlusMinus500 and Alone')
else
end
cfg=[];
tmp={tac19T_s1{squeeze(sleepss_tac19T(11,10,:))>=thresh}};
if ~isempty(tmp)
  grave_tac19T_s1=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,11);ft_singleplotER(figcfg,grave_tac19T_s1);
  title(num2str(sum(squeeze(sleepss_tac19T(11,10,:)))))
else
end
cfg=[];
tmp={tac19T_s2{squeeze(sleepss_tac19T(12,10,:))>=thresh}};
if ~isempty(tmp)
  grave_tac19T_s2=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,13);ft_singleplotER(figcfg,grave_tac19T_s2);
  title(num2str(sum(squeeze(sleepss_tac19T(12,10,:)))))
else
end
% cfg=[];
% tmp={tac19T_s3{squeeze(sleepss_tac19T(13,10,:))>=thresh}};
% grave_tac19T_s3=ft_timelockgrandaverage(cfg,tmp{:});
% subplot(3,8,15);ft_singleplotER(figcfg,grave_tac19T_s3);

cfg=[];
tmp={tac19T_b0{squeeze(sleepss_tac19T(10,11,:))>=thresh}};
if ~isempty(tmp)
  grave_tac19T_b0=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,10);ft_singleplotER(figcfg,grave_tac19T_b0);
  title(num2str(sum(squeeze(sleepss_tac19T(10,11,:)))))
else
end
cfg=[];
tmp={tac19T_b1{squeeze(sleepss_tac19T(11,11,:))>=thresh}};
if ~isempty(tmp)
  grave_tac19T_b1=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,12);ft_singleplotER(figcfg,grave_tac19T_b1);
  title(num2str(sum(squeeze(sleepss_tac19T(11,11,:)))))
else
end
cfg=[];
tmp={tac19T_b2{squeeze(sleepss_tac19T(12,11,:))>=thresh}};
if ~isempty(tmp)
  grave_tac19T_b2=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,14);ft_singleplotER(figcfg,grave_tac19T_b2);
  title(num2str(sum(squeeze(sleepss_tac19T(12,11,:)))))
else
end
cfg=[];
tmp={tac19T_b3{squeeze(sleepss_tac19T(13,11,:))>=thresh}};
if ~isempty(tmp)
  grave_tac19T_b3=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,16);ft_singleplotER(figcfg,grave_tac19T_b3);
  title(num2str(sum(squeeze(sleepss_tac19T(13,11,:)))))
else
end

cfg=[];
tmp={tacAlone_s0{squeeze(sleepss_tacAlone(10,10,:))>=thresh}};
if ~isempty(tmp)
  grave_tacAlone_s0=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,17);ft_singleplotER(figcfg,grave_tacAlone_s0);
  title(num2str(sum(squeeze(sleepss_tacAlone(10,10,:)))))
  ylabel('Tac Alone')
else
end
cfg=[];
tmp={tacAlone_s1{squeeze(sleepss_tacAlone(11,10,:))>=thresh}};
if ~isempty(tmp)
  grave_tacAlone_s1=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,19);ft_singleplotER(figcfg,grave_tacAlone_s1);
  title(num2str(sum(squeeze(sleepss_tacAlone(11,10,:)))))
else
end
cfg=[];
tmp={tacAlone_s2{squeeze(sleepss_tacAlone(12,10,:))>=thresh}};
if ~isempty(tmp)
  grave_tacAlone_s2=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,21);ft_singleplotER(figcfg,grave_tacAlone_s2);
  title(num2str(sum(squeeze(sleepss_tacAlone(12,10,:)))))
else
end
% cfg=[];
% tmp={tacAlone_s3{squeeze(sleepss_tacAlone(13,10,:))>=thresh}};
% grave_tacAlone_s3=ft_timelockgrandaverage(cfg,tmp{:});
% subplot(3,8,23);ft_singleplotER(figcfg,grave_tacAlone_s3);

cfg=[];
tmp={tacAlone_b0{squeeze(sleepss_tacAlone(10,11,:))>=thresh}};
if ~isempty(tmp)
  grave_tacAlone_b0=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,18);ft_singleplotER(figcfg,grave_tacAlone_b0);
  title(num2str(sum(squeeze(sleepss_tacAlone(10,11,:)))))
else
end
cfg=[];
tmp={tacAlone_b1{squeeze(sleepss_tacAlone(11,11,:))>=thresh}};
if ~isempty(tmp)
  grave_tacAlone_b1=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,20);ft_singleplotER(figcfg,grave_tacAlone_b1);
  title(num2str(sum(squeeze(sleepss_tacAlone(11,11,:)))))
else
end
cfg=[];
tmp={tacAlone_b2{squeeze(sleepss_tacAlone(12,11,:))>=thresh}};
if ~isempty(tmp)
  grave_tacAlone_b2=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,22);ft_singleplotER(figcfg,grave_tacAlone_b2);
  title(num2str(sum(squeeze(sleepss_tacAlone(12,11,:)))))
else
end
cfg=[];
tmp={tacAlone_b3{squeeze(sleepss_tacAlone(13,11,:))>=thresh}};
if ~isempty(tmp)
  grave_tacAlone_b3=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,24);ft_singleplotER(figcfg,grave_tacAlone_b3);
  title(num2str(sum(squeeze(sleepss_tacAlone(13,11,:)))))
else
end

%% More testing of tactile-alone response, for sleep=1 only, using tacaloneproc=1 flag
% sorting ERP amplitudes

% Documentation:
% tk flag refers to tk=1 is trialkc-1 and tk=2 is trialkc0
% 'F' at end of name refers to selection of Fz/Fc (5 total) channels averaged over

tt=3;
sleep=1;
iiuse=setdiff(iiBuse,[3:7]);
tacaud=1;
iter=11;
trialkc=-1;  % redo for -1 also

clearvars -except tt sub* *dir ii*use *flag sleep* tacaud iter trialkc

load eeg1010_neighb
submin=iiuse(1)-1;
subuseind=0;
subuse=iiuse;

fsample=1000;
latuse=[-0.5 1];
time=latuse(1):1/fsample:latuse(end);
coloruse=varycolor(length(subuse));
cmap=colormap('parula');

% numtr=10; % for BACN
numtr=20;

for ii=iiuse
  cd([edir sub{ii} ])
  subuseind=subuseind+1;
  for tk=1:3
    try
      load(['tlock_tacaloneERPamp_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_trialkc' num2str(tk-2) '.mat'])
    catch
      tlock_tacT_tlockall(:,:,10:12,tk,subuseind)= nan(63,length(time),3);
      tlock_tac9T_tlockall(:,:,10:12,tk,subuseind)= nan(63,length(time),3);
      tlock_audA_tlockall(:,:,10:12,tk,subuseind)= nan(63,length(time),3);
      tlock_aud1A_tlockall(:,:,10:12,tk,subuseind)= nan(63,length(time),3);
      tlock_nul_tlockall(:,:,10:12,tk,subuseind)= nan(63,length(time),3);
      continue
    end
    for ss=10:12
      clear *tmp
      cfg=[];
      cfg.latency=latuse;
      if ~isempty(tlock_tacT_tlock{tt,ss}) && tlock_tacT_tlock{tt,ss}.dof(17,dsearchn(tlock_tacT_tlock{tt,ss}.time',0))>numtr
        tlock_tacT_tlocktmp=ft_selectdata(cfg,tlock_tacT_tlock{tt,ss});
        numtr_tacT(ss,tk,subuseind)=tlock_tacT_tlock{tt,ss}.dof(17,dsearchn(tlock_tacT_tlock{tt,ss}.time',0));
      elseif ~isempty(tlock_tacT_tlock{tt,ss})
        tlock_tacT_tlocktmp.avg=nan(63,length(time));
        numtr_tacT(ss,tk,subuseind)=tlock_tacT_tlock{tt,ss}.dof(17,dsearchn(tlock_tacT_tlock{tt,ss}.time',0));
      else
        tlock_tacT_tlocktmp.avg=nan(63,length(time));
        numtr_tacT(ss,tk,subuseind)=0;
      end
      if ~isempty(tlock_tac9T_tlock{tt,ss}) && tlock_tac9T_tlock{tt,ss}.dof(17,dsearchn(tlock_tac9T_tlock{tt,ss}.time',0))>numtr
        tlock_tac9T_tlocktmp=ft_selectdata(cfg,tlock_tac9T_tlock{tt,ss});
        numtr_tac9T(ss,tk,subuseind)=tlock_tac9T_tlock{tt,ss}.dof(17,dsearchn(tlock_tac9T_tlock{tt,ss}.time',0));
      elseif ~isempty(tlock_tac9T_tlock{tt,ss})
        tlock_tac9T_tlocktmp.avg=nan(63,length(time));
        numtr_tac9T(ss,tk,subuseind)=tlock_tac9T_tlock{tt,ss}.dof(17,dsearchn(tlock_tac9T_tlock{tt,ss}.time',0));
      else
        tlock_tac9T_tlocktmp.avg=nan(63,length(time));
        numtr_tac9T(ss,tk,subuseind)=0;
      end
      if ~isempty(tlock_audA_tlock{tt,ss}) && tlock_audA_tlock{tt,ss}.dof(17,dsearchn(tlock_audA_tlock{tt,ss}.time',0))>numtr
        tlock_audA_tlocktmp=ft_selectdata(cfg,tlock_audA_tlock{tt,ss});
        numtr_audA(ss,tk,subuseind)=tlock_audA_tlock{tt,ss}.dof(17,dsearchn(tlock_audA_tlock{tt,ss}.time',0));
      elseif ~isempty(tlock_audA_tlock{tt,ss})
        tlock_audA_tlocktmp.avg=nan(63,length(time));
        numtr_audA(ss,tk,subuseind)=tlock_audA_tlock{tt,ss}.dof(17,dsearchn(tlock_audA_tlock{tt,ss}.time',0));
      else
        tlock_audA_tlocktmp.avg=nan(63,length(time));
        numtr_audA(ss,tk,subuseind)=0;
      end
      if ~isempty(tlock_aud1A_tlock{tt,ss}) && tlock_aud1A_tlock{tt,ss}.dof(17,dsearchn(tlock_aud1A_tlock{tt,ss}.time',0))>numtr
        tlock_aud1A_tlocktmp=ft_selectdata(cfg,tlock_aud1A_tlock{tt,ss});
        numtr_aud1A(ss,tk,subuseind)=tlock_aud1A_tlock{tt,ss}.dof(17,dsearchn(tlock_aud1A_tlock{tt,ss}.time',0));
      elseif ~isempty(tlock_aud1A_tlock{tt,ss})
        tlock_aud1A_tlocktmp.avg=nan(63,length(time));
        numtr_aud1A(ss,tk,subuseind)=tlock_aud1A_tlock{tt,ss}.dof(17,dsearchn(tlock_aud1A_tlock{tt,ss}.time',0));
      else
        tlock_aud1A_tlocktmp.avg=nan(63,length(time));
        numtr_aud1A(ss,tk,subuseind)=0;
      end
      if ~isempty(tlock_nul_tlock{tt,ss}) && tlock_nul_tlock{tt,ss}.dof(17,dsearchn(tlock_nul_tlock{tt,ss}.time',0))>numtr
        tlock_nul_tlocktmp=ft_selectdata(cfg,tlock_nul_tlock{tt,ss});
        numtr_nul(ss,tk,subuseind)=tlock_nul_tlock{tt,ss}.dof(17,dsearchn(tlock_nul_tlock{tt,ss}.time',0));
      elseif ~isempty(tlock_nul_tlock{tt,ss})
        tlock_nul_tlocktmp.avg=nan(63,length(time));
        numtr_nul(ss,tk,subuseind)=tlock_nul_tlock{tt,ss}.dof(17,dsearchn(tlock_nul_tlock{tt,ss}.time',0));
      else
        tlock_nul_tlocktmp.avg=nan(63,length(time));
        numtr_nul(ss,tk,subuseind)=0;
      end
      
      tlock_tacT_tlockall(:,:,ss,tk,subuseind)= tlock_tacT_tlocktmp.avg;
      tlock_tac9T_tlockall(:,:,ss,tk,subuseind)= tlock_tac9T_tlocktmp.avg;
      tlock_audA_tlockall(:,:,ss,tk,subuseind)= tlock_audA_tlocktmp.avg;
      tlock_aud1A_tlockall(:,:,ss,tk,subuseind)= tlock_aud1A_tlocktmp.avg;
      tlock_nul_tlockall(:,:,ss,tk,subuseind)= tlock_nul_tlocktmp.avg;
    end % ss
  end % tk
end % ii

% figure;plot(time,squeeze(tlock_tacT_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,10,:)))
% figure;plot(time,squeeze(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,10,:)))
% figure;plot(time,squeeze(tlock_audA_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,10,:)))
% figure;plot(time,squeeze(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,10,:)))
% figure;plot(time,squeeze(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,10,:)))
% figure;imagesc(time,1:subuseind,squeeze(tlock_tacT_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,10,:))');caxis([-5 5])
% figure;imagesc(time,1:subuseind,squeeze(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,10,:))');caxis([-5 5])
% figure;imagesc(time,1:subuseind,squeeze(tlock_audA_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,10,:))');caxis([-5 5])
% figure;imagesc(time,1:subuseind,squeeze(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,10,:))');caxis([-5 5])
% figure;imagesc(time,1:subuseind,squeeze(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,10,:))');caxis([-5 5])

% Get max/min of peaks from awake data
[maxval,maxind]=max(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,10,1,:),5));
[minval,minind]=min(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,10,1,:),5));

% max versus min (peak P200 vs peak N100)
for ss=10:12
  [sortval(:,ss),sortind(:,ss)]=sort(squeeze(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),maxind,ss,1,:)-tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),minind,ss,1,:)));
end

% wideband P200 versus prestim
for ss=10:12
  [sortwideval(:,ss),sortwideind(:,ss)]=sort(squeeze(mean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),maxind-30:maxind+30,ss,1,:),2)-mean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),1:500,ss,1,:),2)));
end

% wideband P200 versus wideband N100
for ss=10:12
  [sortwideP2N1val(:,ss),sortwideP2N1ind(:,ss)]=sort(squeeze(mean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),maxind-30:maxind+30,ss,1,:),2)-mean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),minind-30:minind+30,ss,1,:),2)));
end

load([edir 'spindles.mat']);
for ss=10:12
  tmpuse=length(find(~isnan(sortval(:,ss))));
  tmp=spindles(subuse(sortind(:,ss)),:);        [cc,pp]=corr(sortval(1:tmpuse,ss),tmp(1:tmpuse,:))
  tmp=spindles(subuse(sortwideind(:,ss)),:);    [cc,pp]=corr(sortwideval(1:tmpuse,ss),tmp(1:tmpuse,:))
  tmp=spindles(subuse(sortwideP2N1ind(:,ss)),:);[cc,pp]=corr(sortwideP2N1val(1:tmpuse,ss),tmp(1:tmpuse,:))
end


% for tk=1:3
%   figure(tk)
%   for ss=10:12
%     subplot(3,3,1+(ss-10)*3);plot(time,squeeze(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,:)));axis([-inf inf -5 5])
%     hold on;plot(time,squeeze(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,:),5)),'k');
%     subplot(3,3,2+(ss-10)*3);plot(time,squeeze(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,:)));axis([-inf inf -5 5])
%     hold on;plot(time,squeeze(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,:),5)),'k');
%     subplot(3,3,3+(ss-10)*3);plot(time,squeeze(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,:)));axis([-inf inf -5 5])
%     hold on;plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,:),5)),'k');
%   end
% end

sortTacN2=subuse(sortind(:,12))';
% save([edir 'sortTacN2.mat'],'sort*');
save([edir 'sortTacN2_numtr' num2str(numtr) '.mat'],'sort*');

% tophalf=sortind(end-8:end,12);
% bothalf=sortind(1:9,12);
% figure(10);
% for ss=10:12
%   subplot(3,3,1+(ss-10)*3);axis([-inf inf -5 5])
%   hold on;plot(time,squeeze(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,tophalf),5)),'Color',cmapuse(end,:));
%   hold on;plot(time,squeeze(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,bothalf),5)),'Color',cmapuse(1,:));
%   subplot(3,3,2+(ss-10)*3);axis([-inf inf -5 5])
%   hold on;plot(time,squeeze(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,tophalf),5)),'Color',cmapuse(end,:));
%   hold on;plot(time,squeeze(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,bothalf),5)),'Color',cmapuse(1,:));
%   subplot(3,3,3+(ss-10)*3);axis([-inf inf -5 5])
%   hold on;plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,tophalf),5)),'Color',cmapuse(end,:));
%   hold on;plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,bothalf),5)),'Color',cmapuse(1,:));
% end
%
% tophalfwide=sortwideind(end-8:end,12);
% bothalfwide=sortwideind(1:9,12);
% figure(11);
% for ss=10:12
%   subplot(3,3,1+(ss-10)*3);axis([-inf inf -5 5])
%   hold on;plot(time,squeeze(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,tophalfwide),5)),'Color',cmapuse(end,:));
%   hold on;plot(time,squeeze(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,bothalfwide),5)),'Color',cmapuse(1,:));
%   subplot(3,3,2+(ss-10)*3);axis([-inf inf -5 5])
%   hold on;plot(time,squeeze(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,tophalfwide),5)),'Color',cmapuse(end,:));
%   hold on;plot(time,squeeze(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,bothalfwide),5)),'Color',cmapuse(1,:));
%   subplot(3,3,3+(ss-10)*3);axis([-inf inf -5 5])
%   hold on;plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,tophalfwide),5)),'Color',cmapuse(end,:));
%   hold on;plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,bothalfwide),5)),'Color',cmapuse(1,:));
% end
%
% tophalfwideP2N1=sortwideP2N1ind(end-8:end,12);
% bothalfwideP2N1=sortwideP2N1ind(1:9,12);
% figure(12);
% for ss=10:12
%   subplot(3,3,1+(ss-10)*3);axis([-inf inf -5 5])
%   hold on;plot(time,squeeze(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,tophalfwideP2N1),5)),'Color',cmapuse(end,:));
%   hold on;plot(time,squeeze(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,bothalfwideP2N1),5)),'Color',cmapuse(1,:));
%   subplot(3,3,2+(ss-10)*3);axis([-inf inf -5 5])
%   hold on;plot(time,squeeze(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,tophalfwideP2N1),5)),'Color',cmapuse(end,:));
%   hold on;plot(time,squeeze(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,bothalfwideP2N1),5)),'Color',cmapuse(1,:));
%   subplot(3,3,3+(ss-10)*3);axis([-inf inf -5 5])
%   hold on;plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,tophalfwideP2N1),5)),'Color',cmapuse(end,:));
%   hold on;plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,bothalfwideP2N1),5)),'Color',cmapuse(1,:));
% end

tk=1;
% chanuse={'Fz' 'Cz'};
chanuse={'Fz' 'C1' 'C2' 'FC1' 'FC2'};

tophalfwide=sortwideind(end-8:end,12);
bothalfwide=sortwideind(1:9,12);
sortplot=[bothalfwide  tophalfwide];
cmapuse=cmap(1:4:4*length(sortplot),:);
figure(10);
for ss=10:12
  subplot(3,3,1+(ss-10)*3);axis([-inf inf -5 5])
  %   hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,chanuse),:,ss,tk,tophalfwide),1),5)),'Color',cmapuse(end,:));set(lh(2),'linestyle',':','linewidth',3);
  %   hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,chanuse),:,ss,tk,bothalfwide),1),5)),'Color',cmapuse(1,:));set(lh(2),'linestyle',':','linewidth',3);
  %   subplot(3,3,2+(ss-10)*3);axis([-inf inf -5 5])
  %   hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,chanuse),:,ss,tk,tophalfwide),1),5)),'Color',cmapuse(end,:));set(lh(2),'linestyle',':','linewidth',3);
  %   hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,chanuse),:,ss,tk,bothalfwide),1),5)),'Color',cmapuse(1,:));set(lh(2),'linestyle',':','linewidth',3);
  %   subplot(3,3,3+(ss-10)*3);axis([-inf inf -5 5])
  %   hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,chanuse),:,ss,tk,tophalfwide),1),5)),'Color',cmapuse(end,:));set(lh(2),'linestyle',':','linewidth',3);
  %   hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,chanuse),:,ss,tk,bothalfwide),1),5)),'Color',cmapuse(1,:));set(lh(2),'linestyle',':','linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,chanuse),:,ss,tk,tophalfwide),1),5)),'Color',cmapuse(end,:));set(lh(1),'linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,chanuse),:,ss,tk,bothalfwide),1),5)),'Color',cmapuse(1,:));set(lh(1),'linewidth',3);
  subplot(3,3,2+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,chanuse),:,ss,tk,tophalfwide),1),5)),'Color',cmapuse(end,:));set(lh(1),'linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,chanuse),:,ss,tk,bothalfwide),1),5)),'Color',cmapuse(1,:));set(lh(1),'linewidth',3);
  subplot(3,3,3+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,chanuse),:,ss,tk,tophalfwide),1),5)),'Color',cmapuse(end,:));set(lh(1),'linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,chanuse),:,ss,tk,bothalfwide),1),5)),'Color',cmapuse(1,:));set(lh(1),'linewidth',3);
end

tophalfwideP2N1=sortwideP2N1ind(end-8:end,12);
bothalfwideP2N1=sortwideP2N1ind(1:9,12);
sortplot=[bothalfwideP2N1  tophalfwideP2N1];
cmapuse=cmap(1:4:4*length(sortplot),:);
figure(11);
for ss=10:12
  subplot(3,3,1+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,chanuse),:,ss,tk,tophalfwideP2N1),1),5)),'Color',cmapuse(end,:));set(lh(1),'linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,chanuse),:,ss,tk,bothalfwideP2N1),1),5)),'Color',cmapuse(1,:));set(lh(1),'linewidth',3);
  subplot(3,3,2+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,chanuse),:,ss,tk,tophalfwideP2N1),1),5)),'Color',cmapuse(end,:));set(lh(1),'linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,chanuse),:,ss,tk,bothalfwideP2N1),1),5)),'Color',cmapuse(1,:));set(lh(1),'linewidth',3);
  subplot(3,3,3+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,chanuse),:,ss,tk,tophalfwideP2N1),1),5)),'Color',cmapuse(end,:));set(lh(1),'linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,chanuse),:,ss,tk,bothalfwideP2N1),1),5)),'Color',cmapuse(1,:));set(lh(1),'linewidth',3);
end

tophalf=sortind(end-8:end,12);
bothalf=sortind(1:9,12);
sortplot=[bothalf  tophalf];
cmapuse=cmap(1:4:4*length(sortplot),:);
figure(12);
for ss=10:12
  subplot(3,3,1+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,chanuse),:,ss,tk,tophalf),1),5)),'Color',cmapuse(end,:));set(lh(1),'linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,chanuse),:,ss,tk,bothalf),1),5)),'Color',cmapuse(1,:));set(lh(1),'linewidth',3);
  subplot(3,3,2+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,chanuse),:,ss,tk,tophalf),1),5)),'Color',cmapuse(end,:));set(lh(1),'linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,chanuse),:,ss,tk,bothalf),1),5)),'Color',cmapuse(1,:));set(lh(1),'linewidth',3);
  subplot(3,3,3+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,chanuse),:,ss,tk,tophalf),1),5)),'Color',cmapuse(end,:));set(lh(1),'linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,chanuse),:,ss,tk,bothalf),1),5)),'Color',cmapuse(1,:));set(lh(1),'linewidth',3);
end

% Is bottom half of auditory-alone the *opposite*?


% contrast tophalf bothalf for conditions not used for sorting.
for ss=10:12
  for tk=1:2
    tlock_tac9T_gndavgtop{ss,tk}=[];
    tlock_tac9T_gndavgtop{ss,tk}.label=tlock_nul_tlocktmp.label;
    tlock_tac9T_gndavgtop{ss,tk}.dimord='subj_chan_time';
    tlock_tac9T_gndavgtop{ss,tk}.time=tlock_nul_tlocktmp.time;
    
    tlock_tac9T_gndavgbot{ss,tk}=tlock_tac9T_gndavgtop{ss,tk};
    tlock_aud1A_gndavgtop{ss,tk}=tlock_tac9T_gndavgtop{ss,tk};
    tlock_aud1A_gndavgbot{ss,tk}=tlock_tac9T_gndavgtop{ss,tk};
    tlock_nul_gndavgtop{ss,tk}=tlock_tac9T_gndavgtop{ss,tk};
    tlock_nul_gndavgbot{ss,tk}=tlock_tac9T_gndavgtop{ss,tk};
    
    %     tlock_tac9T_gndavgtop{ss,tk}.individual=permute(squeeze(tlock_tac9T_tlockall(:,:,10,1,tophalf)),[3 1 2]);    % bug that this was '10' always?
    %     tlock_tac9T_gndavgbot{ss,tk}.individual=permute(squeeze(tlock_tac9T_tlockall(:,:,10,1,bothalf)),[3 1 2]);
    %     tlock_aud1A_gndavgtop{ss,tk}.individual=permute(squeeze(tlock_aud1A_tlockall(:,:,10,1,tophalf)),[3 1 2]);
    %     tlock_aud1A_gndavgbot{ss,tk}.individual=permute(squeeze(tlock_aud1A_tlockall(:,:,10,1,bothalf)),[3 1 2]);
    %     tlock_nul_gndavgtop{ss,tk}.individual  =permute(squeeze(tlock_nul_tlockall(:,:,10,1,tophalf)),[3 1 2]);
    %     tlock_nul_gndavgbot{ss,tk}.individual  =permute(squeeze(tlock_nul_tlockall(:,:,10,1,bothalf)),[3 1 2]);
    tlock_tac9T_gndavgtop{ss,tk}.individual=permute(squeeze(tlock_tac9T_tlockall(:,:,ss,tk,tophalf)),[3 1 2]);
    tlock_tac9T_gndavgbot{ss,tk}.individual=permute(squeeze(tlock_tac9T_tlockall(:,:,ss,tk,bothalf)),[3 1 2]);
    tlock_aud1A_gndavgtop{ss,tk}.individual=permute(squeeze(tlock_aud1A_tlockall(:,:,ss,tk,tophalf)),[3 1 2]);
    tlock_aud1A_gndavgbot{ss,tk}.individual=permute(squeeze(tlock_aud1A_tlockall(:,:,ss,tk,bothalf)),[3 1 2]);
    tlock_nul_gndavgtop{ss,tk}.individual  =permute(squeeze(tlock_nul_tlockall(:,:,ss,tk,tophalf)),[3 1 2]);
    tlock_nul_gndavgbot{ss,tk}.individual  =permute(squeeze(tlock_nul_tlockall(:,:,ss,tk,bothalf)),[3 1 2]);
  end
end
for ss=10:12
  for tk=1:2
    tlock_tac9T_gndavgtopwide{ss,tk}=[];
    tlock_tac9T_gndavgtopwide{ss,tk}.label=tlock_nul_tlocktmp.label;
    tlock_tac9T_gndavgtopwide{ss,tk}.dimord='subj_chan_time';
    tlock_tac9T_gndavgtopwide{ss,tk}.time=tlock_nul_tlocktmp.time;
    
    tlock_tac9T_gndavgbotwide{ss,tk}=tlock_tac9T_gndavgtopwide{ss,tk};
    tlock_aud1A_gndavgtopwide{ss,tk}=tlock_tac9T_gndavgtopwide{ss,tk};
    tlock_aud1A_gndavgbotwide{ss,tk}=tlock_tac9T_gndavgtopwide{ss,tk};
    tlock_nul_gndavgtopwide{ss,tk}=tlock_tac9T_gndavgtopwide{ss,tk};
    tlock_nul_gndavgbotwide{ss,tk}=tlock_tac9T_gndavgtopwide{ss,tk};
    
    tlock_tac9T_gndavgtopwide{ss,tk}.individual=permute(squeeze(tlock_tac9T_tlockall(:,:,ss,tk,tophalfwide )),[3 1 2]);
    tlock_tac9T_gndavgbotwide{ss,tk}.individual=permute(squeeze(tlock_tac9T_tlockall(:,:,ss,tk,bothalfwide)),[3 1 2]);
    tlock_aud1A_gndavgtopwide{ss,tk}.individual=permute(squeeze(tlock_aud1A_tlockall(:,:,ss,tk,tophalfwide)),[3 1 2]);
    tlock_aud1A_gndavgbotwide{ss,tk}.individual=permute(squeeze(tlock_aud1A_tlockall(:,:,ss,tk,bothalfwide)),[3 1 2]);
    tlock_nul_gndavgtopwide{ss,tk}.individual  =permute(squeeze(tlock_nul_tlockall(:,:,ss,tk,tophalfwide)),[3 1 2]);
    tlock_nul_gndavgbotwide{ss,tk}.individual  =permute(squeeze(tlock_nul_tlockall(:,:,ss,tk,bothalfwide)),[3 1 2]);
  end
end
for ss=10:12
  for tk=1:2
    tlock_tac9T_gndavgtopwideP2N1{ss,tk}=[];
    tlock_tac9T_gndavgtopwideP2N1{ss,tk}.label=tlock_nul_tlocktmp.label;
    tlock_tac9T_gndavgtopwideP2N1{ss,tk}.dimord='subj_chan_time';
    tlock_tac9T_gndavgtopwideP2N1{ss,tk}.time=tlock_nul_tlocktmp.time;
    
    tlock_tac9T_gndavgbotwideP2N1{ss,tk}=tlock_tac9T_gndavgtopwideP2N1{ss,tk};
    tlock_aud1A_gndavgtopwideP2N1{ss,tk}=tlock_tac9T_gndavgtopwideP2N1{ss,tk};
    tlock_aud1A_gndavgbotwideP2N1{ss,tk}=tlock_tac9T_gndavgtopwideP2N1{ss,tk};
    tlock_nul_gndavgtopwideP2N1{ss,tk}=tlock_tac9T_gndavgtopwideP2N1{ss,tk};
    tlock_nul_gndavgbotwideP2N1{ss,tk}=tlock_tac9T_gndavgtopwideP2N1{ss,tk};
    
    tlock_tac9T_gndavgtopwideP2N1{ss,tk}.individual=permute(squeeze(tlock_tac9T_tlockall(:,:,ss,tk,tophalfwideP2N1 )),[3 1 2]);
    tlock_tac9T_gndavgbotwideP2N1{ss,tk}.individual=permute(squeeze(tlock_tac9T_tlockall(:,:,ss,tk,bothalfwideP2N1)),[3 1 2]);
    tlock_aud1A_gndavgtopwideP2N1{ss,tk}.individual=permute(squeeze(tlock_aud1A_tlockall(:,:,ss,tk,tophalfwideP2N1)),[3 1 2]);
    tlock_aud1A_gndavgbotwideP2N1{ss,tk}.individual=permute(squeeze(tlock_aud1A_tlockall(:,:,ss,tk,bothalfwideP2N1)),[3 1 2]);
    tlock_nul_gndavgtopwideP2N1{ss,tk}.individual  =permute(squeeze(tlock_nul_tlockall(:,:,ss,tk,tophalfwideP2N1)),[3 1 2]);
    tlock_nul_gndavgbotwideP2N1{ss,tk}.individual  =permute(squeeze(tlock_nul_tlockall(:,:,ss,tk,bothalfwideP2N1)),[3 1 2]);
  end
end

nsub=length(tophalf);

cfg=[];
cfg.latency=[.03 .43];
cfg.neighbours=neighbours;
cfg.parameter='individual';
cfg.method='montecarlo';
cfg.numrandomization=2000;
cfg.correctm='cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2;
cfg.statistic='indepsamplesT';
cfg.design=zeros(2,2*nsub);
cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
% cfg.design(2,:)=[1:nsub 1:nsub];
cfg.ivar=1;
% cfg.uvar=2;
for ss=10:12
  for tk=1:2
    statt_tactopbot{ss,tk}=ft_timelockstatistics(cfg, tlock_tac9T_gndavgtop{ss,tk}, tlock_tac9T_gndavgbot{ss,tk});
    statt_audtopbot{ss,tk}=ft_timelockstatistics(cfg, tlock_aud1A_gndavgtop{ss,tk}, tlock_aud1A_gndavgbot{ss,tk});
    statt_nultopbot{ss,tk}=ft_timelockstatistics(cfg, tlock_nul_gndavgtop{ss,tk},   tlock_nul_gndavgbot{ss,tk});
  end
end
% numtr=20; not signif for  tac or aud; but yes for nul in 12,1 (p=0.047)
for ss=10:12
  for tk=1:2
    statt_tactopbotwide{ss,tk}=ft_timelockstatistics(cfg, tlock_tac9T_gndavgtopwide{ss,tk}, tlock_tac9T_gndavgbotwide{ss,tk});
    statt_audtopbotwide{ss,tk}=ft_timelockstatistics(cfg, tlock_aud1A_gndavgtopwide{ss,tk}, tlock_aud1A_gndavgbotwide{ss,tk});
    statt_nultopbotwide{ss,tk}=ft_timelockstatistics(cfg, tlock_nul_gndavgtopwide{ss,tk},   tlock_nul_gndavgbotwide{ss,tk});
  end
end
% numtr=20;  tac: 10,2 p=0.07;    aud 11,1 p=0.045  11,2 p=0.02  12,2 p=0.04;    nul  none
for ss=10:12
  for tk=1:2
    statt_tactopbotwideP2N1{ss,tk}=ft_timelockstatistics(cfg, tlock_tac9T_gndavgtopwideP2N1{ss,tk}, tlock_tac9T_gndavgbotwideP2N1{ss,tk});
    statt_audtopbotwideP2N1{ss,tk}=ft_timelockstatistics(cfg, tlock_aud1A_gndavgtopwideP2N1{ss,tk}, tlock_aud1A_gndavgbotwideP2N1{ss,tk});
    statt_nultopbotwideP2N1{ss,tk}=ft_timelockstatistics(cfg, tlock_nul_gndavgtopwideP2N1{ss,tk},   tlock_nul_gndavgbotwideP2N1{ss,tk});
  end
end
% numtr=20;   tac: none;   aud:  none;   nul: 12,1 p=0.036

% try for just one channel of main effect.
cfg=[];
% cfg.channel='Fz';
cfg.channel=chanuse;
cfg.avgoverchan='yes';
cfg.latency=[.03 .43];
cfg.parameter='individual';
cfg.method='montecarlo';
cfg.numrandomization=2000;
cfg.correctm='cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.statistic='indepsamplesT';
cfg.design=zeros(2,2*nsub);
cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
cfg.ivar=1;
for ss=10:12
  for tk=1:2
    statt_tactopbotF{ss,tk}=ft_timelockstatistics(cfg, tlock_tac9T_gndavgtop{ss,tk}, tlock_tac9T_gndavgbot{ss,tk});
    statt_audtopbotF{ss,tk}=ft_timelockstatistics(cfg, tlock_aud1A_gndavgtop{ss,tk}, tlock_aud1A_gndavgbot{ss,tk});
    statt_nultopbotF{ss,tk}=ft_timelockstatistics(cfg, tlock_nul_gndavgtop{ss,tk},   tlock_nul_gndavgbot{ss,tk});
  end
end
% numtr=20;  tac 12,1 p=0.02;   aud 10,1 p=0.04  10,2 p=0.06;    nul 12,2 p=0.06
for ss=10:12
  for tk=1:2
    statt_tactopbotwideF{ss,tk}=ft_timelockstatistics(cfg, tlock_tac9T_gndavgtopwide{ss,tk}, tlock_tac9T_gndavgbotwide{ss,tk});
    statt_audtopbotwideF{ss,tk}=ft_timelockstatistics(cfg, tlock_aud1A_gndavgtopwide{ss,tk}, tlock_aud1A_gndavgbotwide{ss,tk});
    statt_nultopbotwideF{ss,tk}=ft_timelockstatistics(cfg, tlock_nul_gndavgtopwide{ss,tk},   tlock_nul_gndavgbotwide{ss,tk});
  end
end
% numtr=20;   tac 12,1 p=0.052;    aud 10,1 p=0.04  11,1 p=0.04  12,1 p=0.04  10,2 p=0.058  11,2 p=0.02  12,2 p=0.061;    nul  10,1  p=0.0005  10,2 p=0.006   12,2 p=0.04
for ss=10:12
  for tk=1:2
    statt_tactopbotwideP2N1F{ss,tk}=ft_timelockstatistics(cfg, tlock_tac9T_gndavgtopwideP2N1{ss,tk}, tlock_tac9T_gndavgbotwideP2N1{ss,tk});
    statt_audtopbotwideP2N1F{ss,tk}=ft_timelockstatistics(cfg, tlock_aud1A_gndavgtopwideP2N1{ss,tk}, tlock_aud1A_gndavgbotwideP2N1{ss,tk});
    statt_nultopbotwideP2N1F{ss,tk}=ft_timelockstatistics(cfg, tlock_nul_gndavgtopwideP2N1{ss,tk},   tlock_nul_gndavgbotwideP2N1{ss,tk});
  end
end
% numtr=20;    tac 12,1 p=0.051     aud 10,2 p=0.074    nul None

% % try for just one channel of main effect and around 100ms and 200ms peaks
% cfg=[];
% cfg.channel='Fz';
% cfg.latency=[.08 .23];
% cfg.parameter='individual';
% cfg.method='montecarlo';
% cfg.numrandomization=2000;
% cfg.correctm='cluster';
% cfg.clusteralpha = 0.05;
% cfg.clusterstatistic = 'maxsum';
% cfg.statistic='indepsamplesT';
% cfg.design=zeros(2,2*nsub);
% cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
% cfg.ivar=1;
% for ss=10:12
%   for tk=1:2
%     statt_tactopbotF12{ss,tk}=ft_timelockstatistics(cfg, tlock_tac9T_gndavgtop{ss,tk}, tlock_tac9T_gndavgbot{ss,tk});
%     statt_audtopbotF12{ss,tk}=ft_timelockstatistics(cfg, tlock_aud1A_gndavgtop{ss,tk}, tlock_aud1A_gndavgbot{ss,tk});
%     statt_nultopbotF12{ss,tk}=ft_timelockstatistics(cfg, tlock_nul_gndavgtop{ss,tk},   tlock_nul_gndavgbot{ss,tk});
%   end
% end
% % numtr:10 tac is now signif, aud trend, but nul highly signif!  ???



% %%%%%%%%%%%%%%%%% this code was run with 'sitting' data loaded on one machine and 'lying' data loaded on other machine.

% Uta suggested sorting based on P200 tactile awake data peak versus prestim
tt=3;
sleep=0;
iiuse=intersect(setdiff(iiBuse,[3:7]),iiSuse);
tacaud=1;
iter=27;
trialkc=-1;  % redo for -1 also

submin=iiuse(1)-1;
subuseind=0;
subuse=iiuse;

clearvars -except tt sub* edir ddir ii*use *flag sleep* tacaud iter trialkc

% numtr=10;
for ii=iiuse
  cd([edir sub{ii} ])
  subuseind=subuseind+1;
  try
    load(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'])
  catch
    load(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '.mat'])
  end
  cfg=[];
  cfg.latency=[-.5 0.5];
  tlock_tacall{subuseind}=ft_selectdata(cfg,tlock_tTacAlone{5,3,10});
  tlock_tacavg(:,:,subuseind)=tlock_tacall{subuseind}.avg;
end % ii
time=tlock_tacall{subuseind}.time;

% Get max/min of peaks from awake data
[maxval,maxind]=max(nanmean(tlock_tacavg(match_str(tlock_tacall{subuseind}.label,'Fz'),:,:),3));
[minval,minind]=min(nanmean(tlock_tacavg(match_str(tlock_tacall{subuseind}.label,'Fz'),:,:),3));

% wide P200 vs prestim
[sortwideWval,sortwideWind]=sort(squeeze(mean(tlock_tacavg(match_str(tlock_tacall{subuseind}.label,'Fz'),maxind-30:maxind+30,:),2)-mean(tlock_tacavg(match_str(tlock_tacall{subuseind}.label,'Fz'),1:500,:),2)));

% peak P200 vs peak N100
[sortP2N1Wval,sortP2N1Wind]=sort(squeeze(tlock_tacavg(match_str(tlock_tacall{subuseind}.label,'Fz'),maxind,:)-tlock_tacavg(match_str(tlock_tacall{subuseind}.label,'Fz'),minind,:)));

% wide P200 vs wide N100
[sortwideP2N1Wval,sortwideP2N1Wind]=sort(squeeze(mean(tlock_tacavg(match_str(tlock_tacall{subuseind}.label,'Fz'),maxind-30:maxind+30,:),2)-mean(tlock_tacavg(match_str(tlock_tacall{subuseind}.label,'Fz'),minind-30:minind+30,:),2)));

load([edir 'spindles.mat'])
[cc,pp]=corr(sortwideP2N1Wval,spindles(subuse(sortwideP2N1Wind),:))
% none significant


tophalfWwide=sortwideWind(end-7:end);
bothalfWwide=sortwideWind(1:8);

tophalfWP2N1=sortP2N1Wind(end-7:end);
bothalfWP2N1=sortP2N1Wind(1:8);

tophalfWwideP2N1=sortwideP2N1Wind(end-7:end);
bothalfWwideP2N1=sortwideP2N1Wind(1:8);

iiuse_tophalfWwide=iiuse(tophalfWwide);
iiuse_bothalfWwide=iiuse(bothalfWwide);
iiuse_tophalfWP2N1=iiuse(tophalfWP2N1);
iiuse_bothalfWP2N1=iiuse(bothalfWP2N1);
iiuse_tophalfWwideP2N1=iiuse(tophalfWwideP2N1);
iiuse_bothalfWwideP2N1=iiuse(bothalfWwideP2N1);
save([edir 'sortTacWP200.mat'],'iiuse_*');

% % % % % %

load([edir 'sortTacWP200.mat']);
iiBfinal=iiBuse(5:end);

coloruse=varycolor(length(subuse));
sortplot=[iiuse_bothalfWwideP2N1  iiuse_tophalfWwideP2N1]
cmapuse=cmap(1:4:4*length(sortplot),:);

figure(20);
for ss=10:12
  subplot(3,3,1+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_tophalfWwide')),5)),'Color',cmapuse(end,:));set(lh(2),'linestyle',':','linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_bothalfWwide')),5)),'Color',cmapuse(1,:));set(lh(2),'linestyle',':','linewidth',3);
  subplot(3,3,2+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_tophalfWwide')),5)),'Color',cmapuse(end,:));set(lh(2),'linestyle',':','linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_bothalfWwide')),5)),'Color',cmapuse(1,:));set(lh(2),'linestyle',':','linewidth',3);
  subplot(3,3,3+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_tophalfWwide')),5)),'Color',cmapuse(end,:));set(lh(2),'linestyle',':','linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_bothalfWwide')),5)),'Color',cmapuse(1,:));set(lh(2),'linestyle',':','linewidth',3);
end
% Conclusion:  P200 in W sitting data partly explains P200 in W lying data, but not N2 lying data


% FINAL in BACN 2015 poster.
tk=1;
figure(42);
for ss=10:12
  subplot(3,3,1+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_tacT_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'C1' 'C2' 'FC1' 'FC2'}),:,ss,tk,dsearchn(iiBfinal',iiuse_tophalfWwideP2N1')),5),1)),'Color',cmapuse(end,:));set(lh(1),'linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_tacT_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'C1' 'C2' 'FC1' 'FC2'}),:,ss,tk,dsearchn(iiBfinal',iiuse_bothalfWwideP2N1')),5),1)),'Color',cmapuse(1,:));set(lh(1),'linewidth',3);
  if ss==12,xlabel('Time (s)');end
  if ss==10,title('Tactile alone');end
  if ss==10,ylabel('Awake'),elseif ss==11,ylabel('N1');elseif ss==12,ylabel('N2');end
  set(gca,'FontSize',20)
  subplot(3,3,2+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_audA_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'C1' 'C2' 'FC1' 'FC2'}),:,ss,tk,dsearchn(iiBfinal',iiuse_tophalfWwideP2N1')),5),1)),'Color',cmapuse(end,:));set(lh(1),'linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_audA_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'C1' 'C2' 'FC1' 'FC2'}),:,ss,tk,dsearchn(iiBfinal',iiuse_bothalfWwideP2N1')),5),1)),'Color',cmapuse(1,:));set(lh(1),'linewidth',3);
  if ss==12,xlabel('Time (s)');end
  if ss==10,title('Auditory alone');end
  set(gca,'FontSize',20)
  subplot(3,3,3+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'C1' 'C2' 'FC1' 'FC2'}),:,ss,tk,dsearchn(iiBfinal',iiuse_tophalfWwideP2N1')),5),1)),'Color',cmapuse(end,:));set(lh(1),'linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'C1' 'C2' 'FC1' 'FC2'}),:,ss,tk,dsearchn(iiBfinal',iiuse_bothalfWwideP2N1')),5),1)),'Color',cmapuse(1,:));set(lh(1),'linewidth',3);
  if ss==12,xlabel('Time (s)');end
  if ss==10,title('Null');end
  if ss==12,lg=legend({'Top half (Part 2, P200-N100)' 'Bottom half (Part 2, P200-N100)'});set(lg,'fontsize',14);end
  set(gca,'FontSize',20)
end
print(42,[fdir 'unisensory_mediansplitWideP2N1_sleep1_numtr' num2str(numtr) '.png'],'-dpng');
print(42,[fdir 'unisensory_mediansplitWideP2N1_sleep1_numtr' num2str(numtr) '.eps'],'-depsc');

figure(43);
for ss=10:12
  subplot(3,3,1+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'C1' 'C2' 'FC1' 'FC2'}),:,ss,tk,dsearchn(iiBfinal',iiuse_tophalfWwideP2N1')),5),1)),'Color',cmapuse(end,:));set(lh(1),'linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'C1' 'C2' 'FC1' 'FC2'}),:,ss,tk,dsearchn(iiBfinal',iiuse_bothalfWwideP2N1')),5),1)),'Color',cmapuse(1,:));set(lh(1),'linewidth',3);
  if ss==12,xlabel('Time (s)');end
  if ss==10,title('Tactile alone');end
  if ss==10,ylabel('Awake'),elseif ss==11,ylabel('N1');elseif ss==12,ylabel('N2');end
  set(gca,'FontSize',20)
  subplot(3,3,2+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'C1' 'C2' 'FC1' 'FC2'}),:,ss,tk,dsearchn(iiBfinal',iiuse_tophalfWwideP2N1')),5),1)),'Color',cmapuse(end,:));set(lh(1),'linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'C1' 'C2' 'FC1' 'FC2'}),:,ss,tk,dsearchn(iiBfinal',iiuse_bothalfWwideP2N1')),5),1)),'Color',cmapuse(1,:));set(lh(1),'linewidth',3);
  if ss==12,xlabel('Time (s)');end
  if ss==10,title('Auditory alone');end
  set(gca,'FontSize',20)
  subplot(3,3,3+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'C1' 'C2' 'FC1' 'FC2'}),:,ss,tk,dsearchn(iiBfinal',iiuse_tophalfWwideP2N1')),5),1)),'Color',cmapuse(end,:));set(lh(1),'linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'C1' 'C2' 'FC1' 'FC2'}),:,ss,tk,dsearchn(iiBfinal',iiuse_bothalfWwideP2N1')),5),1)),'Color',cmapuse(1,:));set(lh(1),'linewidth',3);
  if ss==12,xlabel('Time (s)');end
  if ss==10,title('Null');end
  if ss==12,lg=legend({'Top half (Part 2, P200-N100)' 'Bottom half (Part 2, P200-N100)'});set(lg,'fontsize',14);end
  set(gca,'FontSize',20)
end
print(43,[fdir 'unisensory_mediansplitWideP2N1_sleep1_numtr' num2str(numtr) '.png'],'-dpng');
print(43,[fdir 'unisensory_mediansplitWideP2N1_sleep1_numtr' num2str(numtr) '.eps'],'-depsc');

tk=1;
figure(44);
timeinduse=1:1001;
for ss=10:12
  subplot(3,3,1+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time(timeinduse),squeeze(nanmean(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'C1' 'C2' 'FC1' 'FC2'}),(timeinduse),ss,tk,dsearchn(iiBfinal',iiuse_tophalfWwideP2N1')),5),1)),'Color',cmapuse(end,:));set(lh(1),'linewidth',3);
  hold on;lh=plot(time(timeinduse),squeeze(nanmean(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'C1' 'C2' 'FC1' 'FC2'}),(timeinduse),ss,tk,dsearchn(iiBfinal',iiuse_bothalfWwideP2N1')),5),1)),'Color',cmapuse(1,:));set(lh(1),'linewidth',3);
  if ss==12,xlabel('Time (s)');end
  if ss==10,title('Tactile alone');end
  if ss==10,ylabel('Awake'),elseif ss==11,ylabel('N1');elseif ss==12,ylabel('N2');end
  set(gca,'FontSize',20)
  subplot(3,3,2+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time(timeinduse),squeeze(nanmean(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'C1' 'C2' 'FC1' 'FC2'}),(timeinduse),ss,tk,dsearchn(iiBfinal',iiuse_tophalfWwideP2N1')),5),1)),'Color',cmapuse(end,:));set(lh(1),'linewidth',3);
  hold on;lh=plot(time(timeinduse),squeeze(nanmean(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'C1' 'C2' 'FC1' 'FC2'}),(timeinduse),ss,tk,dsearchn(iiBfinal',iiuse_bothalfWwideP2N1')),5),1)),'Color',cmapuse(1,:));set(lh(1),'linewidth',3);
  if ss==12,xlabel('Time (s)');end
  if ss==10,title('Auditory alone');end
  set(gca,'FontSize',20)
  subplot(3,3,3+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time(timeinduse),squeeze(nanmean(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'C1' 'C2' 'FC1' 'FC2'}),(timeinduse),ss,tk,dsearchn(iiBfinal',iiuse_tophalfWwideP2N1')),5),1)),'Color',cmapuse(end,:));set(lh(1),'linewidth',3);
  hold on;lh=plot(time(timeinduse),squeeze(nanmean(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'C1' 'C2' 'FC1' 'FC2'}),(timeinduse),ss,tk,dsearchn(iiBfinal',iiuse_bothalfWwideP2N1')),5),1)),'Color',cmapuse(1,:));set(lh(1),'linewidth',3);
  if ss==12,xlabel('Time (s)');end
  if ss==10,title('Null');end
  if ss==12,lg=legend({'Top half (Part 2, P200-N100)' 'Bottom half (Part 2, P200-N100)'});set(lg,'fontsize',14);end
  set(gca,'FontSize',20)
end
print(44,[fdir 'unisensory_mediansplitWideP2N1_sleep1_numtr' num2str(numtr) '.png'],'-dpng');
print(44,[fdir 'unisensory_mediansplitWideP2N1_sleep1_numtr' num2str(numtr) '.eps'],'-depsc');

tk=1;
figure(22);
for ss=10:12
  subplot(3,3,1+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(tlock_tacT_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_tophalfWP2N1')),5)),'Color',cmapuse(end,:));set(lh(2),'linestyle',':','linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(tlock_tacT_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_bothalfWP2N1')),5)),'Color',cmapuse(1,:));set(lh(2),'linestyle',':','linewidth',3);
  subplot(3,3,2+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(tlock_audA_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_tophalfWP2N1')),5)),'Color',cmapuse(end,:));set(lh(2),'linestyle',':','linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(tlock_audA_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_bothalfWP2N1')),5)),'Color',cmapuse(1,:));set(lh(2),'linestyle',':','linewidth',3);
  subplot(3,3,3+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_tophalfWP2N1')),5)),'Color',cmapuse(end,:));set(lh(2),'linestyle',':','linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_bothalfWP2N1')),5)),'Color',cmapuse(1,:));set(lh(2),'linestyle',':','linewidth',3);
end

tk=1;
figure(23);
for ss=10:12
  subplot(3,3,1+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(tlock_tacT_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_tophalfWP2N1')),5)),'Color',cmapuse(end,:));set(lh(2),'linestyle',':','linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(tlock_tacT_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_bothalfWP2N1')),5)),'Color',cmapuse(1,:));set(lh(2),'linestyle',':','linewidth',3);
  subplot(3,3,2+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(tlock_audA_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_tophalfWP2N1')),5)),'Color',cmapuse(end,:));set(lh(2),'linestyle',':','linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(tlock_audA_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_bothalfWP2N1')),5)),'Color',cmapuse(1,:));set(lh(2),'linestyle',':','linewidth',3);
  subplot(3,3,3+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_tophalfWP2N1')),5)),'Color',cmapuse(end,:));set(lh(2),'linestyle',':','linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_bothalfWP2N1')),5)),'Color',cmapuse(1,:));set(lh(2),'linestyle',':','linewidth',3);
end

% choose which 'iiuse_*' to use
% decision: iiuse_tophalfWP2N1 iiuse_bothalfWP2N1

% test if difference can be explained in trial numbers
for ss=10:12
  [hh,pp]=ttest2(squeeze(numtr_tac9T(ss,1,dsearchn(iiuse',iiuse_bothalfWwideP2N1'))),squeeze(numtr_tac9T(ss,1,dsearchn(iiuse',iiuse_tophalfWwideP2N1'))))
end

sortplot=[iiuse_bothalfWwideP2N1  iiuse_tophalfWwideP2N1]
cmap=colormap('parula');
cmapuse=cmap(1:4:4*length(sortplot),:);

% this also in final BACN poster
tk=1;figure(1)
for ss=10:12
  subplot(3,3,1+(ss-10)*3);
  for ii=1:length(sortplot),plot(time,squeeze(tlock_tacT_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,dsearchn(subuse',sortplot(ii)))),'Color',cmapuse(ii,:));hold on;end;axis([-inf inf -7 7])
  hold on;
  plot(time,squeeze(nanmean(tlock_tacT_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,:),5)),'k','LineWidth',3');
  if ss==12,xlabel('Time (s)');end
  if ss==10,title('Tactile alone');end
  if ss==10,ylabel('Awake'),elseif ss==11,ylabel('N1');elseif ss==12,ylabel('N2');end
  set(gca,'FontSize',20)
  subplot(3,3,2+(ss-10)*3);
  for ii=1:length(sortplot),plot(time,squeeze(tlock_audA_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,dsearchn(subuse',sortplot(ii)))),'Color',cmapuse(ii,:));hold on;end;axis([-inf inf -7 7])
  hold on;
  plot(time,squeeze(nanmean(tlock_audA_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,:),5)),'k','LineWidth',3');
  if ss==12,xlabel('Time (s)');end
  if ss==10,title('Auditory alone');end
  set(gca,'FontSize',20)
  subplot(3,3,3+(ss-10)*3);
  for ii=1:length(sortplot),plot(time,squeeze(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,dsearchn(subuse',sortplot(ii)))),'Color',cmapuse(ii,:));hold on;end;axis([-inf inf -7 7])
  hold on;
  plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,:),5)),'k','LineWidth',3');
  if ss==12,xlabel('Time (s)');end
  if ss==10,title('Null');end
  set(gca,'FontSize',20)
end
print(1,[fdir 'unisensory_sortedWideP2N1_sleep1_numtr' num2str(numtr) '.png'],'-dpng');
print(1,[fdir 'unisensory_sortedWideP2N1_sleep1_numtr' num2str(numtr) '.eps'],'-depsc');

tk=1;figure(2)
for ss=10:12
  subplot(3,3,1+(ss-10)*3);
  for ii=1:length(sortplot),plot(time,squeeze(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,dsearchn(subuse',sortplot(ii)))),'Color',cmapuse(ii,:));hold on;end;axis([-inf inf -7 7])
  hold on;
  plot(time,squeeze(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,:),5)),'k','LineWidth',3');
  if ss==12,xlabel('Time (s)');end
  if ss==10,title('Tactile alone');end
  if ss==10,ylabel('Awake'),elseif ss==11,ylabel('N1');elseif ss==12,ylabel('N2');end
  set(gca,'FontSize',20)
  subplot(3,3,2+(ss-10)*3);
  for ii=1:length(sortplot),plot(time,squeeze(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,dsearchn(subuse',sortplot(ii)))),'Color',cmapuse(ii,:));hold on;end;axis([-inf inf -7 7])
  hold on;
  plot(time,squeeze(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,:),5)),'k','LineWidth',3');
  if ss==12,xlabel('Time (s)');end
  if ss==10,title('Auditory alone');end
  set(gca,'FontSize',20)
  subplot(3,3,3+(ss-10)*3);
  for ii=1:length(sortplot),plot(time,squeeze(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,dsearchn(subuse',sortplot(ii)))),'Color',cmapuse(ii,:));hold on;end;axis([-inf inf -7 7])
  hold on;
  plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,:),5)),'k','LineWidth',3');
  if ss==12,xlabel('Time (s)');end
  if ss==10,title('Null');end
  set(gca,'FontSize',20)
end
print(2,[fdir 'unisensory_sortedWideP2N1_sleep1_numtr' num2str(numtr) '.png'],'-dpng');
print(2,[fdir 'unisensory_sortedWideP2N1_sleep1_numtr' num2str(numtr) '.eps'],'-depsc');

tk=1;figure(3)
for ss=10:12
  subplot(3,3,1+(ss-10)*3);
  for ii=1:length(sortplot),plot(time(timeinduse),squeeze(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),timeinduse,ss,tk,dsearchn(subuse',sortplot(ii)))),'Color',cmapuse(ii,:));hold on;end;axis([-inf inf -7 7])
  hold on;
  plot(time(timeinduse),squeeze(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),timeinduse,ss,tk,:),5)),'k','LineWidth',3');
  if ss==12,xlabel('Time (s)');end
  if ss==10,title('Tactile alone');end
  if ss==10,ylabel('Awake'),elseif ss==11,ylabel('N1');elseif ss==12,ylabel('N2');end
  set(gca,'FontSize',20)
  subplot(3,3,2+(ss-10)*3);
  for ii=1:length(sortplot),plot(time(timeinduse),squeeze(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),timeinduse,ss,tk,dsearchn(subuse',sortplot(ii)))),'Color',cmapuse(ii,:));hold on;end;axis([-inf inf -7 7])
  hold on;
  plot(time(timeinduse),squeeze(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),timeinduse,ss,tk,:),5)),'k','LineWidth',3');
  if ss==12,xlabel('Time (s)');end
  if ss==10,title('Auditory alone');end
  set(gca,'FontSize',20)
  subplot(3,3,3+(ss-10)*3);
  for ii=1:length(sortplot),plot(time(timeinduse),squeeze(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),timeinduse,ss,tk,dsearchn(subuse',sortplot(ii)))),'Color',cmapuse(ii,:));hold on;end;axis([-inf inf -7 7])
  hold on;
  plot(time(timeinduse),squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),timeinduse,ss,tk,:),5)),'k','LineWidth',3');
  if ss==12,xlabel('Time (s)');end
  if ss==10,title('Null');end
  set(gca,'FontSize',20)
end
print(3,[fdir 'unisensory_sortedWideP2N1_sleep1_numtr' num2str(numtr) '.png'],'-dpng');
print(3,[fdir 'unisensory_sortedWideP2N1_sleep1_numtr' num2str(numtr) '.eps'],'-depsc');

% tk=1;figure(2)
% for ss=10:12
%   subplot(3,3,1+(ss-10)*3);
%   for ii=1:length(subuse),plot(time,squeeze(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,sortwideind(ii,12))),'Color',cmapuse(ii,:));hold on;end;axis([-inf inf -7 7])
%   hold on;
%   plot(time,squeeze(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,:),5)),'k','LineWidth',3');
%   if ss==12,xlabel('Time (s)');end
%   if ss==10,title('Tactile alone');end
%   if ss==10,ylabel('Awake'),elseif ss==11,ylabel('N1');elseif ss==12,ylabel('N2');end
%   subplot(3,3,2+(ss-10)*3);
%   for ii=1:length(subuse),plot(time,squeeze(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,sortwideind(ii,12))),'Color',cmapuse(ii,:));hold on;end;axis([-inf inf -7 7])
%   hold on;
%   plot(time,squeeze(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,:),5)),'k','LineWidth',3');
%   if ss==12,xlabel('Time (s)');end
%   if ss==10,title('Auditory alone');end
%   subplot(3,3,3+(ss-10)*3);
%   for ii=1:length(subuse),plot(time,squeeze(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,sortwideind(ii,12))),'Color',cmapuse(ii,:));hold on;end;axis([-inf inf -7 7])
%   hold on;
%   plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,:),5)),'k','LineWidth',3');
%   if ss==12,xlabel('Time (s)');end
%   if ss==10,title('Null');end
% end

% contrast tophalf bothalf for conditions not used for sorting.

for ss=10:12
  for tk=1:2
    tlock_tac9T_gndavgtopWwideP2N1{ss,tk}=[];
    tlock_tac9T_gndavgtopWwideP2N1{ss,tk}.label=tlock_nul_tlocktmp.label;
    tlock_tac9T_gndavgtopWwideP2N1{ss,tk}.dimord='subj_chan_time';
    tlock_tac9T_gndavgtopWwideP2N1{ss,tk}.time=tlock_nul_tlocktmp.time;
    
    tlock_tac9T_gndavgbotWwideP2N1{ss,tk}=tlock_tac9T_gndavgtopWwideP2N1{ss,tk};
    tlock_aud1A_gndavgtopWwideP2N1{ss,tk}=tlock_tac9T_gndavgtopWwideP2N1{ss,tk};
    tlock_aud1A_gndavgbotWwideP2N1{ss,tk}=tlock_tac9T_gndavgtopWwideP2N1{ss,tk};
    tlock_nul_gndavgtopWwideP2N1{ss,tk}=tlock_tac9T_gndavgtopWwideP2N1{ss,tk};
    tlock_nul_gndavgbotWwideP2N1{ss,tk}=tlock_tac9T_gndavgtopWwideP2N1{ss,tk};
    
    tlock_tac9T_gndavgtopWwideP2N1{ss,tk}.individual=permute(squeeze(tlock_tac9T_tlockall(:,:,ss,tk,dsearchn(subuse',iiuse_tophalfWwideP2N1'))),[3 1 2]);
    tlock_tac9T_gndavgbotWwideP2N1{ss,tk}.individual=permute(squeeze(tlock_tac9T_tlockall(:,:,ss,tk,dsearchn(subuse',iiuse_bothalfWwideP2N1'))),[3 1 2]);
    tlock_aud1A_gndavgtopWwideP2N1{ss,tk}.individual=permute(squeeze(tlock_aud1A_tlockall(:,:,ss,tk,dsearchn(subuse',iiuse_tophalfWwideP2N1'))),[3 1 2]);
    tlock_aud1A_gndavgbotWwideP2N1{ss,tk}.individual=permute(squeeze(tlock_aud1A_tlockall(:,:,ss,tk,dsearchn(subuse',iiuse_bothalfWwideP2N1'))),[3 1 2]);
    tlock_nul_gndavgtopWwideP2N1{ss,tk}.individual  =permute(squeeze(tlock_nul_tlockall(:,:,ss,tk,dsearchn(subuse',iiuse_tophalfWwideP2N1'))),[3 1 2]);
    tlock_nul_gndavgbotWwideP2N1{ss,tk}.individual  =permute(squeeze(tlock_nul_tlockall(:,:,ss,tk,dsearchn(subuse',iiuse_bothalfWwideP2N1'))),[3 1 2]);
  end
end
for ss=10:12
  for tk=1:2
    tlock_tacT_gndavgtopWwideP2N1{ss,tk}=[];
    tlock_tacT_gndavgtopWwideP2N1{ss,tk}.label=tlock_nul_tlocktmp.label;
    tlock_tacT_gndavgtopWwideP2N1{ss,tk}.dimord='subj_chan_time';
    tlock_tacT_gndavgtopWwideP2N1{ss,tk}.time=tlock_nul_tlocktmp.time;
    
    tlock_tacT_gndavgbotWwideP2N1{ss,tk}=tlock_tacT_gndavgtopWwideP2N1{ss,tk};
    tlock_audA_gndavgtopWwideP2N1{ss,tk}=tlock_tacT_gndavgtopWwideP2N1{ss,tk};
    tlock_audA_gndavgbotWwideP2N1{ss,tk}=tlock_tacT_gndavgtopWwideP2N1{ss,tk};
    tlock_nul_gndavgtopWwideP2N1{ss,tk}=tlock_tacT_gndavgtopWwideP2N1{ss,tk};
    tlock_nul_gndavgbotWwideP2N1{ss,tk}=tlock_tacT_gndavgtopWwideP2N1{ss,tk};
    
    tlock_tacT_gndavgtopWwideP2N1{ss,tk}.individual=permute(squeeze(tlock_tacT_tlockall(:,:,ss,tk,dsearchn(subuse',iiuse_tophalfWwideP2N1'))),[3 1 2]);
    tlock_tacT_gndavgbotWwideP2N1{ss,tk}.individual=permute(squeeze(tlock_tacT_tlockall(:,:,ss,tk,dsearchn(subuse',iiuse_bothalfWwideP2N1'))),[3 1 2]);
    tlock_audA_gndavgtopWwideP2N1{ss,tk}.individual=permute(squeeze(tlock_audA_tlockall(:,:,ss,tk,dsearchn(subuse',iiuse_tophalfWwideP2N1'))),[3 1 2]);
    tlock_audA_gndavgbotWwideP2N1{ss,tk}.individual=permute(squeeze(tlock_audA_tlockall(:,:,ss,tk,dsearchn(subuse',iiuse_bothalfWwideP2N1'))),[3 1 2]);
    tlock_nul_gndavgtopWwideP2N1{ss,tk}.individual  =permute(squeeze(tlock_nul_tlockall(:,:,ss,tk,dsearchn(subuse',iiuse_tophalfWwideP2N1'))),[3 1 2]);
    tlock_nul_gndavgbotWwideP2N1{ss,tk}.individual  =permute(squeeze(tlock_nul_tlockall(:,:,ss,tk,dsearchn(subuse',iiuse_bothalfWwideP2N1'))),[3 1 2]);
  end
end

% % contrast tophalf bothalf for conditions not used for sorting.
% for ss=10:12
%   for tk=1:2
%     tlock_tac9T_gndavgtopWwide{ss,tk}=[];
%     tlock_tac9T_gndavgtopWwide{ss,tk}.label=tlock_nul_tlocktmp.label;
%     tlock_tac9T_gndavgtopWwide{ss,tk}.dimord='subj_chan_time';
%     tlock_tac9T_gndavgtopWwide{ss,tk}.time=tlock_nul_tlocktmp.time;
%
%     tlock_tac9T_gndavgbotWwide{ss,tk}=tlock_tac9T_gndavgtopWwide{ss,tk};
%     tlock_aud1A_gndavgtopWwide{ss,tk}=tlock_tac9T_gndavgtopWwide{ss,tk};
%     tlock_aud1A_gndavgbotWwide{ss,tk}=tlock_tac9T_gndavgtopWwide{ss,tk};
%     tlock_nul_gndavgtopWwide{ss,tk}=tlock_tac9T_gndavgtopWwide{ss,tk};
%     tlock_nul_gndavgbotWwide{ss,tk}=tlock_tac9T_gndavgtopWwide{ss,tk};
%
%     tlock_tac9T_gndavgtopWwide{ss,tk}.individual=permute(squeeze(tlock_tac9T_tlockall(:,:,ss,tk,dsearchn(subuse',iiuse_tophalfWwide'))),[3 1 2]);
%     tlock_tac9T_gndavgbotWwide{ss,tk}.individual=permute(squeeze(tlock_tac9T_tlockall(:,:,ss,tk,dsearchn(subuse',iiuse_bothalfWwide'))),[3 1 2]);
%     tlock_aud1A_gndavgtopWwide{ss,tk}.individual=permute(squeeze(tlock_aud1A_tlockall(:,:,ss,tk,dsearchn(subuse',iiuse_tophalfWwide'))),[3 1 2]);
%     tlock_aud1A_gndavgbotWwide{ss,tk}.individual=permute(squeeze(tlock_aud1A_tlockall(:,:,ss,tk,dsearchn(subuse',iiuse_bothalfWwide'))),[3 1 2]);
%     tlock_nul_gndavgtopWwide{ss,tk}.individual  =permute(squeeze(tlock_nul_tlockall(:,:,ss,tk,dsearchn(subuse',iiuse_tophalfWwide'))),[3 1 2]);
%     tlock_nul_gndavgbotWwide{ss,tk}.individual  =permute(squeeze(tlock_nul_tlockall(:,:,ss,tk,dsearchn(subuse',iiuse_bothalfWwide'))),[3 1 2]);
%   end
% end

nsub=size(tlock_tac9T_gndavgtopWwideP2N1{ss,tk}.individual,1);
% try for just one channel of main effect.
cfg=[];
cfg.channel={'Fz' 'C1' 'C2' 'FC1' 'FC2'};
cfg.avgoverchan='yes';
cfg.latency=[.05 .25];
cfg.parameter='individual';
cfg.method='montecarlo';
cfg.numrandomization=2000;
% cfg.correctm='fdr';
cfg.correctm='cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.statistic='indepsamplesT';
cfg.design=zeros(2,2*nsub);
cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
cfg.ivar=1;
% for ss=10:12
%   for tk=1:2
%     statt_tactopbotWP2N1F{ss,tk}=ft_timelockstatistics(cfg, tlock_tac9T_gndavgtopWP2N1{ss,tk}, tlock_tac9T_gndavgbotWP2N1{ss,tk});
%     statt_audtopbotWP2N1F{ss,tk}=ft_timelockstatistics(cfg, tlock_aud1A_gndavgtopWP2N1{ss,tk}, tlock_aud1A_gndavgbotWP2N1{ss,tk});
%     statt_nultopbotWP2N1F{ss,tk}=ft_timelockstatistics(cfg, tlock_nul_gndavgtopWP2N1{ss,tk},   tlock_nul_gndavgbotWP2N1{ss,tk});
%   end
% end

for ss=10:12
  for tk=1:2
    statt_tactopbotWwideP2N1F{ss,tk}=ft_timelockstatistics(cfg, tlock_tac9T_gndavgtopWwideP2N1{ss,tk}, tlock_tac9T_gndavgbotWwideP2N1{ss,tk});
    statt_audtopbotWwideP2N1F{ss,tk}=ft_timelockstatistics(cfg, tlock_aud1A_gndavgtopWwideP2N1{ss,tk}, tlock_aud1A_gndavgbotWwideP2N1{ss,tk});
    statt_nultopbotWwideP2N1F{ss,tk}=ft_timelockstatistics(cfg, tlock_nul_gndavgtopWwideP2N1{ss,tk},   tlock_nul_gndavgbotWwideP2N1{ss,tk});
  end
end

% numtr=10:
% numtr=20:  tac 10,1  p=0.03   12,1 p=0.02   10,2 p=0.051  12,2 p=0.003     aud:  11,1 p=0.09  11,2 p=0.07     nul:  None
% numtr=20: N  tac 10 7v6  11 8v8 12,1 8v8      aud:  10 6v4   11 8v8  12 8v8     nul:  10 1v5   11 6v7  12  8v8
for ss=10:12
  for tk=1:2
    statt_tacTtopbotWwideP2N1F{ss,tk}=ft_timelockstatistics(cfg, tlock_tacT_gndavgtopWwideP2N1{ss,tk}, tlock_tacT_gndavgbotWwideP2N1{ss,tk});
    statt_audAtopbotWwideP2N1F{ss,tk}=ft_timelockstatistics(cfg, tlock_audA_gndavgtopWwideP2N1{ss,tk}, tlock_audA_gndavgbotWwideP2N1{ss,tk});
  end
end
% numtr=20;  tac 11,1 p=0.055  12,1 p=0.01   11,2 p=0.054  12,2 p=0.049     aud:  11,1 p=0.057
% numtr=20:  N  tac 10  6v1  11  7v6   8v8;     aud:  10  2v6  11 7v7  12 8v8
save([edir 'stat_medsplit_sleep1_numtr' num2str(numtr) '.mat'],'stat*')

load([edir 'stat_medsplit_sleep1_numtr' num2str(numtr) '.mat'],'stat*')
cfg=[];
cfg.latency=[.05 .25];
cfg.parameter='individual';
cfg.method='montecarlo';
cfg.numrandomization=2000;
cfg.neighbours=neighbours;
% cfg.correctm='fdr';
cfg.correctm='cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.statistic='indepsamplesT';
cfg.design=zeros(2,2*nsub);
cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
cfg.ivar=1;
for ss=10:12
  for tk=1:2
    statt_tactopbotWwideP2N1{ss,tk}=ft_timelockstatistics(cfg, tlock_tac9T_gndavgtopWwideP2N1{ss,tk}, tlock_tac9T_gndavgbotWwideP2N1{ss,tk});
    statt_audtopbotWwideP2N1{ss,tk}=ft_timelockstatistics(cfg, tlock_aud1A_gndavgtopWwideP2N1{ss,tk}, tlock_aud1A_gndavgbotWwideP2N1{ss,tk});
    statt_nultopbotWwideP2N1{ss,tk}=ft_timelockstatistics(cfg, tlock_nul_gndavgtopWwideP2N1{ss,tk},   tlock_nul_gndavgbotWwideP2N1{ss,tk});
  end
end
% T not 9T
for ss=10:12
  for tk=1:2
    statt_tacTtopbotWwideP2N1{ss,tk}=ft_timelockstatistics(cfg, tlock_tacT_gndavgtopWwideP2N1{ss,tk}, tlock_tacT_gndavgbotWwideP2N1{ss,tk});
    statt_audAtopbotWwideP2N1{ss,tk}=ft_timelockstatistics(cfg, tlock_audA_gndavgtopWwideP2N1{ss,tk}, tlock_audA_gndavgbotWwideP2N1{ss,tk});
  end
end

cfg.latency=[.05 .5];  % full possible time window if including 9T
for ss=10:12
  for tk=1:2
    statt_tactopbotWwideP2N15{ss,tk}=ft_timelockstatistics(cfg, tlock_tac9T_gndavgtopWwideP2N1{ss,tk}, tlock_tac9T_gndavgbotWwideP2N1{ss,tk});
    statt_audtopbotWwideP2N15{ss,tk}=ft_timelockstatistics(cfg, tlock_aud1A_gndavgtopWwideP2N1{ss,tk}, tlock_aud1A_gndavgbotWwideP2N1{ss,tk});
    statt_nultopbotWwideP2N15{ss,tk}=ft_timelockstatistics(cfg, tlock_nul_gndavgtopWwideP2N1{ss,tk},   tlock_nul_gndavgbotWwideP2N1{ss,tk});
  end
end
for ss=10:12
  for tk=1:2
    statt_tacTtopbotWwideP2N15{ss,tk}=ft_timelockstatistics(cfg, tlock_tacT_gndavgtopWwideP2N1{ss,tk}, tlock_tacT_gndavgbotWwideP2N1{ss,tk});
    statt_audAtopbotWwideP2N15{ss,tk}=ft_timelockstatistics(cfg, tlock_audA_gndavgtopWwideP2N1{ss,tk}, tlock_audA_gndavgbotWwideP2N1{ss,tk});
  end
end
cfg.latency=[.05 .8];  % hint of effect around 500-700ms from plot; can't test this with 9T
for ss=10:12
  for tk=1:2
    statt_tacTtopbotWwideP2N18{ss,tk}=ft_timelockstatistics(cfg, tlock_tacT_gndavgtopWwideP2N1{ss,tk}, tlock_tacT_gndavgbotWwideP2N1{ss,tk});
    statt_audAtopbotWwideP2N18{ss,tk}=ft_timelockstatistics(cfg, tlock_audA_gndavgtopWwideP2N1{ss,tk}, tlock_audA_gndavgbotWwideP2N1{ss,tk});
  end
end


save([edir 'stat_medsplit_sleep1_numtr' num2str(numtr) '.mat'],'stat*')
% see also sleep_specific/file_record_keeping.xls (tab topBot)


cfg=[];cfg.xlim=[.2 .3];
cfg.layout='EEG1020';cfg.parameter='stat'; cfg.zlim='maxabs';
ft_topoplotER(cfg,statt_tacTtopbotWwideP2N15{12,2})


% for ss=10:12
%   for tk=1:2
%     statt_tactopbotWwideF{ss,tk}=ft_timelockstatistics(cfg, tlock_tac9T_gndavgtopWwide{ss,tk}, tlock_tac9T_gndavgbotWwide{ss,tk});
%     statt_audtopbotWwideF{ss,tk}=ft_timelockstatistics(cfg, tlock_aud1A_gndavgtopWwide{ss,tk}, tlock_aud1A_gndavgbotWwide{ss,tk});
%     statt_nultopbotWwideF{ss,tk}=ft_timelockstatistics(cfg, tlock_nul_gndavgtopWwide{ss,tk},   tlock_nul_gndavgbotWwide{ss,tk});
%   end
% end

% examine peaks across participants at times when statistically significant
ss=10;
timeind=dsearchn(tlock_tac9T_gndavgbotWwideP2N1{ss,1}.time',statt_tactopbotWwideP2N1F{ss,1}.time((statt_tactopbotWwideP2N1F{ss,1}.prob<.05))');
chanind=match_str(tlock_tac9T_gndavgbotWwideP2N1{ss,1}.label,cfg.channel);

[mx,mxind]=max(nanmean(nanmean(tlock_aud1A_tlockall(chanind,:,ss,1,:),1),5))
[mn,mnind]=min(nanmean(nanmean(tlock_aud1A_tlockall(chanind,:,ss,1,:),1),5))

sub_aud1A_tactime=squeeze(mean(mean(tlock_aud1A_tlockall(chanind,timeind,ss,1,:),2),1));
sub_aud1A_amintime=squeeze(mean(mean(tlock_aud1A_tlockall(chanind,mnind-length(timeind)/2:mnind+length(timeind)/2,ss,1,:),2),1));
sub_aud1A_amaxtime=squeeze(mean(mean(tlock_aud1A_tlockall(chanind,mxind-length(timeind)/2:mxind+length(timeind)/2,ss,1,:),2),1));

[cc,pp]=corr(sub_aud1A_tactime(~isnan(sub_aud1A_amaxtime)),sub_tac9T(~isnan(sub_aud1A_amaxtime)))
[cc,pp]=corr(sub_aud1A_amintime(~isnan(sub_aud1A_amaxtime)),sub_tac9T(~isnan(sub_aud1A_amaxtime)))
[cc,pp]=corr(sub_aud1A_amaxtime(~isnan(sub_aud1A_amaxtime)),sub_tac9T(~isnan(sub_aud1A_amaxtime)))

sub_aud1A_tactime=squeeze(mean(mean(tlock_aud1A_tlockall(chanind,timeind,ss,1,:),2),1));
sub_aud1A_amintime=squeeze(mean(mean(tlock_aud1A_tlockall(chanind,mnind-length(timeind)/2:mnind+length(timeind)/2,ss,1,:),2),1));
sub_aud1A_amaxtime=squeeze(mean(mean(tlock_aud1A_tlockall(chanind,mxind-length(timeind)/2:mxind+length(timeind)/2,ss,1,:),2),1));

[mx,mxind]=max(nanmean(nanmean(tlock_tac9T_tlockall(chanind,501:1000,ss,1,:),1),5))
[mn,mnind]=min(nanmean(nanmean(tlock_tac9T_tlockall(chanind,501:1000,ss,1,:),1),5))

sub_tac9T_tactime=squeeze(mean(mean(tlock_tac9T_tlockall(chanind,timeind,ss,1,:),2),1));
sub_tac9T_amintime=squeeze(mean(mean(tlock_tac9T_tlockall(chanind,mnind-length(timeind)/2:mnind+length(timeind)/2,ss,1,:),2),1));
sub_tac9T_amaxtime=squeeze(mean(mean(tlock_tac9T_tlockall(chanind,mxind-length(timeind)/2:mxind+length(timeind)/2,ss,1,:),2),1));

[cc,pp]=corr(sub_aud1A_tactime(~isnan(sub_aud1A_amaxtime)),sub_tac9T_tactime(~isnan(sub_aud1A_amaxtime)))
[cc,pp]=corr(sub_aud1A_amintime(~isnan(sub_aud1A_amaxtime)),sub_tac9T_tactime(~isnan(sub_aud1A_amaxtime)))
[cc,pp]=corr(sub_aud1A_amaxtime(~isnan(sub_aud1A_amaxtime)),sub_tac9T_tactime(~isnan(sub_aud1A_amaxtime)))

[cc,pp]=corr(sub_aud1A_tactime(~isnan(sub_aud1A_amaxtime)),sub_tac9T_amintime(~isnan(sub_aud1A_amaxtime)))
[cc,pp]=corr(sub_aud1A_amintime(~isnan(sub_aud1A_amaxtime)),sub_tac9T_amintime(~isnan(sub_aud1A_amaxtime)))
[cc,pp]=corr(sub_aud1A_amaxtime(~isnan(sub_aud1A_amaxtime)),sub_tac9T_amintime(~isnan(sub_aud1A_amaxtime)))

[cc,pp]=corr(sub_aud1A_tactime(~isnan(sub_aud1A_amaxtime)),sub_tac9T_amaxtime(~isnan(sub_aud1A_amaxtime))) % significant
[cc,pp]=corr(sub_aud1A_amintime(~isnan(sub_aud1A_amaxtime)),sub_tac9T_amaxtime(~isnan(sub_aud1A_amaxtime))) % significant
[cc,pp]=corr(sub_aud1A_amaxtime(~isnan(sub_aud1A_amaxtime)),sub_tac9T_amaxtime(~isnan(sub_aud1A_amaxtime)))

figure;plot(sub_aud1A_amintime(~isnan(sub_aud1A_amintime)),sub_tac9T_amaxtime(~isnan(sub_aud1A_amaxtime)),'o')

% compare across stages


% compare to spindles
load([edir 'spindles.mat'])
sortplot=[iiuse_bothalfWwideP2N1  iiuse_tophalfWwideP2N1];

%% More testing of tactile-alone response, using tacaloneproc=1 flag, phase relationship

tt=3;
sleep=1;
iiuse=setdiff(iiBuse,[3:7]);
tacaud=1;
iter=11;
trialkc=-1;  % redo for -1 also

submin=iiuse(1)-1;
subuseind=0;
subuse=iiuse;
erpamp_phasefft_aud1A_all=nan(4,13,2,12,length(iiuse));
erpamp_phasefft_audA_all=nan(4,13,2,12,length(iiuse));
erpamp_phasefft_tac9T_all=nan(4,13,2,12,length(iiuse));
erpamp_phasefft_tacT_all=nan(4,13,2,12,length(iiuse));
erpamp_phasefft_nul_all=nan(4,13,2,12,length(iiuse));
erpamp_phasehilb_aud1A_all=nan(4,3,2,12,length(iiuse));
erpamp_phasehilb_audA_all=nan(4,3,2,12,length(iiuse));
erpamp_phasehilb_tac9T_all=nan(4,3,2,12,length(iiuse));
erpamp_phasehilb_tacT_all=nan(4,3,2,12,length(iiuse));
erpamp_phasehilb_nul_all=nan(4,3,2,12,length(iiuse));
erpamp_powrat_audA_all=nan(26,3,2,length(iiuse));
erpamp_powrat_aud1A_all=nan(26,3,2,length(iiuse));
erpamp_powrat_tac9T_all=nan(26,3,2,length(iiuse));
erpamp_powrat_tacT_all=nan(26,3,2,length(iiuse));
erpamp_powrat_nul_all=nan(26,3,2,length(iiuse));
for ii=iiuse
  cd([edir sub{ii} ])
  load(['tlock_tacaloneERPamp_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'])
  subuseind=subuseind+1;
  for ss=10:12
    if ~isempty(erpamp_phasefft_aud1A{tt,ss}), erpamp_phasefft_aud1A_all(:,:,:,ss,subuseind)=erpamp_phasefft_aud1A{tt,ss}; end
    if ~isempty(erpamp_phasefft_audA{tt,ss}),  erpamp_phasefft_audA_all(:,:,:,ss,subuseind) =erpamp_phasefft_audA{tt,ss};  end
    if ~isempty(erpamp_phasefft_tac9T{tt,ss}),  erpamp_phasefft_tac9T_all(:,:,:,ss,subuseind)=erpamp_phasefft_tac9T{tt,ss};   end
    if ~isempty(erpamp_phasefft_tacT{tt,ss}),  erpamp_phasefft_tacT_all(:,:,:,ss,subuseind) =erpamp_phasefft_tacT{tt,ss};  end
    if ~isempty(erpamp_phasefft_nul{tt,ss}),   erpamp_phasefft_nul_all(:,:,:,ss,subuseind)  =erpamp_phasefft_nul{tt,ss};   end
    if ~isempty(erpamp_phasehilb_aud1A{tt,ss}), erpamp_phasehilb_aud1A_all(:,:,:,ss,subuseind)=erpamp_phasehilb_aud1A{tt,ss}; end
    if ~isempty(erpamp_phasehilb_audA{tt,ss}),  erpamp_phasehilb_audA_all(:,:,:,ss,subuseind) =erpamp_phasehilb_audA{tt,ss};  end
    if ~isempty(erpamp_phasehilb_tac9T{tt,ss}), erpamp_phasehilb_tac9T_all(:,:,:,ss,subuseind)=erpamp_phasehilb_tac9T{tt,ss};   end
    if ~isempty(erpamp_phasehilb_tacT{tt,ss}),  erpamp_phasehilb_tacT_all(:,:,:,ss,subuseind) =erpamp_phasehilb_tacT{tt,ss};  end
    if ~isempty(erpamp_phasehilb_nul{tt,ss}),   erpamp_phasehilb_nul_all(:,:,:,ss,subuseind)  =erpamp_phasehilb_nul{tt,ss};   end
  end
  erpamp_powrat_audA_all(:,:,:,subuseind) =erpamp_powrat_audA;
  erpamp_powrat_aud1A_all(:,:,:,subuseind)=erpamp_powrat_aud1A;
  erpamp_powrat_tac9T_all(:,:,:,subuseind)=erpamp_powrat_tac9T;
  erpamp_powrat_tacT_all(:,:,:,subuseind) =erpamp_powrat_tacT;
  erpamp_powrat_nul_all(:,:,:,subuseind)  =erpamp_powrat_nul;
end

for ss=[10:12 24:26]
  ss
  for pk=1:2
    pk
    [hh,pp]=ttest(squeeze(erpamp_powrat_tac9T_all(ss,1,pk,:)),squeeze(erpamp_powrat_tac9T_all(ss,3,pk,:)))
    [hh,pp]=ttest(squeeze(erpamp_powrat_tacT_all(ss,1,pk,:)), squeeze(erpamp_powrat_tacT_all(ss,3,pk,:)))
    [hh,pp]=ttest(squeeze(erpamp_powrat_aud1A_all(ss,1,pk,:)),squeeze(erpamp_powrat_aud1A_all(ss,3,pk,:)))
    [hh,pp]=ttest(squeeze(erpamp_powrat_audA_all(ss,1,pk,:)), squeeze(erpamp_powrat_audA_all(ss,3,pk,:)))
    [hh,pp]=ttest(squeeze(erpamp_powrat_nul_all(ss,1,pk,:)),  squeeze(erpamp_powrat_nul_all(ss,3,pk,:)))
  end
end
% results:
% for trialkc=0
% for ss=24, pk=1, then audA is sign
% for ss=26, pk=2, then tac9T is sign (tacT borders sign):  Does this say anything more than for N2 the P2 goes away?
%
% for trialkc=-1
% for ss=10, pk=1, tac9T
% for ss=24, pk=1, then tacT, and pk=2 tacT and aud1A
% for ss=26, pk=1, then tac9T (tacT borders) and aud1A
%
% for trialkc=-1 (only tac9T and aud1A and nul)
% for ss=10, pk=1, tac9T
% for ss=24, pk=2, aud1A
% for ss=26, pk=1, then tac9T (tacT borders) and aud1A


close all
alphaval=.005;
for ss=[10:12]
  for ff=1:13
    for pk=1:2
      ppfa1(ff,pk,ss)=anova1(squeeze(erpamp_phasefft_aud1A_all(:,ff,pk,ss,:))',[],'off');
      %       ppfa(ff,pk,ss)=anova1(squeeze(erpamp_phasefft_audA_all(:,ff,pk,ss,:))',[],'off');
      ppft9(ff,pk,ss)=anova1(squeeze(erpamp_phasefft_tac9T_all(:,ff,pk,ss,:))',[],'off');
      %       ppft(ff,pk,ss)=anova1(squeeze(erpamp_phasefft_tacT_all(:,ff,pk,ss,:))',[],'off');
      ppfn(ff,pk,ss)=anova1(squeeze(erpamp_phasefft_nul_all(:,ff,pk,ss,:))',[],'off');
      if ff<4
        ppha1(ff,pk,ss)=anova1(squeeze(erpamp_phasehilb_aud1A_all(:,ff,pk,ss,:))',[],'off');
        %         ppha(ff,pk,ss)=anova1(squeeze(erpamp_phasehilb_audA_all(:,ff,pk,ss,:))',[],'off');
        ppht9(ff,pk,ss)=anova1(squeeze(erpamp_phasehilb_tac9T_all(:,ff,pk,ss,:))',[],'off');
        %         ppht(ff,pk,ss)=anova1(squeeze(erpamp_phasehilb_tacT_all(:,ff,pk,ss,:))',[],'off');
        pphn(ff,pk,ss)=anova1(squeeze(erpamp_phasehilb_nul_all(:,ff,pk,ss,:))',[],'off');
      end
    end
    figure(1);subplot(3,13,(ss-10)*13+ff);bar(squeeze(nanmean(erpamp_phasefft_tac9T_all(:,ff,:,ss,:),5)));axis([0 5 -3 3])
    if ppft9(ff,1,ss)<alphaval, hold on; subplot(3,13,(ss-10)*13+ff);plot(2,-2.5,'b*');end
    if ppft9(ff,2,ss)<alphaval, hold on; subplot(3,13,(ss-10)*13+ff);plot(3,-2.5,'r*');end
    figure(2);subplot(3,13,(ss-10)*13+ff);bar(squeeze(nanmean(erpamp_phasefft_aud1A_all(:,ff,:,ss,:),5)));axis([0 5 -3 3])
    if ppfa1(ff,1,ss)<alphaval, hold on; subplot(3,13,(ss-10)*13+ff);plot(2,-2.5,'b*');end
    if ppfa1(ff,2,ss)<alphaval, hold on; subplot(3,13,(ss-10)*13+ff);plot(3,-2.5,'r*');end
    figure(3);subplot(3,13,(ss-10)*13+ff);bar(squeeze(nanmean(erpamp_phasefft_nul_all(:,ff,:,ss,:),5)));axis([0 5 -3 3])
    if ppfn(ff,1,ss)<alphaval, hold on; subplot(3,13,(ss-10)*13+ff);plot(2,-2.5,'b*');end
    if ppfn(ff,2,ss)<alphaval, hold on; subplot(3,13,(ss-10)*13+ff);plot(3,-2.5,'r*');end
    if ff<4
      figure(4);subplot(3,3,(ss-10)*3+ff);bar(squeeze(nanmean(erpamp_phasehilb_tac9T_all(:,ff,:,ss,:),5)));axis([0 5 -3 3])
      if ppht9(ff,1,ss)<alphaval, hold on; subplot(3,3,(ss-10)*3+ff);plot(2,-2.5,'b*');end
      if ppht9(ff,2,ss)<alphaval, hold on; subplot(3,3,(ss-10)*3+ff);plot(3,-2.5,'r*');end
      figure(5);subplot(3,3,(ss-10)*3+ff);bar(squeeze(nanmean(erpamp_phasehilb_aud1A_all(:,ff,:,ss,:),5)));axis([0 5 -3 3])
      if ppha1(ff,1,ss)<alphaval, hold on; subplot(3,3,(ss-10)*3+ff);plot(2,-2.5,'b*');end
      if ppha1(ff,2,ss)<alphaval, hold on; subplot(3,3,(ss-10)*3+ff);plot(3,-2.5,'r*');end
      figure(6);subplot(3,3,(ss-10)*3+ff);bar(squeeze(nanmean(erpamp_phasehilb_nul_all(:,ff,:,ss,:),5)));axis([0 5 -3 3])
      if pphn(ff,1,ss)<alphaval, hold on; subplot(3,3,(ss-10)*3+ff);plot(2,-2.5,'b*');end
      if pphn(ff,2,ss)<alphaval, hold on; subplot(3,3,(ss-10)*3+ff);plot(3,-2.5,'r*');end
    end
  end
end
% interesting!  some N2 delta phase dependence and W alpha phase dependence;
% however, how much of it is above and beyond what is seen in nul trials?

