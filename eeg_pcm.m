%% start

audtac_startup

%% Load data

soalist=[1 3 4 5 6 7 9];

ftver=0;  % 0 means the svn version on the local machine that is already on path
mcseed=13;
chanuse_sleep0={'all' '-F4'};
chanuse_sleep1={'all' '-AF7' '-AF3' '-Fp1'};
sleep=0;
plotflag=1;

if sleep
  chanuse=chanuse_sleep1;
  iteruse=11;
  trialkcflag=1; % 0 for ignore trialkc (old results; all trials), or 1 for use trialkc value
  if trialkcflag
    trialkc=0;  % 0 (no Kc) or 1 (only Kc) or -1 (all trials)
  else
    trialkc=-1;
  end
  ss=12;
  
else
  chanuse=chanuse_sleep0;
  %   warning('change me back to 27!!')
  trialkcflag=0; % 0 for ignore trialkc (old results; all trials), or 1 for use trialkc value
  iteruse=27;
  %   iteruse=32;
  trialkc=-1;
  ss=10;
end
iter=iteruse;
tt=3;

close all

if sleep
  if tophalfflag
    load sortTacN2.mat
    subuseall=sort(sortTacN2(end-8:end)');
  else
    subuseall=setdiff(iiBuse,[3:7]);
  end
else
  subuseall=iiSuse;
end

% numtimesteps=25;
% numtimesteps=numel(50:30:440);  % 14
numtimesteps=numel(50:25:450);  % 17
subcnt=0;
for ii=subuseall
  subcnt=subcnt+1;
  clear tlockdiff_slide data_condmean
  %   cd([edir sub{ii} ])
  
  % Load individual unaveraged data
  load([edir sub{ii} '/tlock_pcm_' sub{ii} '_sleep' num2str(sleep) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat']);
  
  % Extract each difference across time windows of ERP, normalise.
  llcnt=0;
  trlcnt=0;
  for ll=soalist
    llcnt=llcnt+1;
    switch numtimesteps
      case 25
        for jj=1:numtimesteps; % every 20 ms
          tlockdiff_slide{llcnt}(:,:,jj)=mean(pcm_TPAmMSPN{ll,tt,ss}.trial(:,:,jj*20-19:jj*20+1),3);
        end
      case 14
        for jj=1:numtimesteps;
          tlockdiff_slide{llcnt}(:,:,jj)=mean(pcm_TPAmMSPN{ll,tt,ss}.trial(:,:,20+jj*30-49:20+jj*30+50),3);  % every 30 ms over a window of 100ms duration
        end
      case 17
        for jj=1:numtimesteps;
          [~,ssi]=sort(rand(1,size(pcm_TPAmMSPN{ll,tt,ss}.trial,1)));
          tlockdiff_slide{llcnt}(:,:,jj)=mean(pcm_TPAmMSPN{ll,tt,ss}.trial(ssi,:,25+jj*25-24:25+jj*25+25),3);  % every 25 ms over a window of 50ms duration
        end
    end
    designX{subcnt}(1+trlcnt:size(tlockdiff_slide{llcnt},1)+trlcnt,llcnt)=1;
    trlcnt=trlcnt+size(tlockdiff_slide{llcnt},1);
    data_condmean_squeeze{llcnt}=mean(tlockdiff_slide{llcnt},1);
    data_condmean{llcnt}=repmat(mean(tlockdiff_slide{llcnt},1),[size(tlockdiff_slide{llcnt},1) 1 1]);
  end
  datacat{subcnt}=cat(1,tlockdiff_slide{:}); % non-equal number of trials per condition; okay?
  partVec2{subcnt}=[rem(1:size(datacat{subcnt},1),2)+1]';  % this is nonequal number of trials per partition & condition; okay?
  partVec5{subcnt}=[rem(1:size(datacat{subcnt},1),5)+1]';  % this is nonequal number of trials per partition & condition; okay?
  datameancat=cat(1,data_condmean{:});
  datacat_meansubt{subcnt}=datacat{subcnt}-datameancat;
  datameansqueezecat=cat(1,data_condmean_squeeze{:});
  if size(designX{subcnt},1)~=size(datacat{subcnt},1)
    error('something went wrong trlcnt');
  end
  
  % Mahalanobis distance (from Mate)
  sigma = compCovEEG(datacat{subcnt},designX{subcnt});
  %   pdistvec =  NaN(7*(7-1)/2,numtimesteps);
  %   for iTimePoint = 1:numtimesteps
  % %     actTimePointFeats = squeeze(feat_avg(:,iTimePoint,:))';
  % %     actTimePointFeats = squeeze(datameancat(:,:,iTimePoint));
  %     actTimePointFeats = squeeze(datameansqueezecat(:,:,iTimePoint));
  %     pdistvec(:,iTimePoint) = pdist(actTimePointFeats,'mahalanobis',sigma(:,:,iTimePoint));
  %   end
  %   if plotflag
  %     for jj=1:25,subplot(5,5,jj);imagesc(squareform(pdistvec(:,jj)));caxis([0 1.5]);end
  %   end
  
  
  % Multivariate noise normalisation
  for jj=1:numtimesteps
    datacat_meansubt_noisenorm{subcnt}(:,:,jj)=datacat_meansubt{subcnt}(:,:,jj)*(sigma(:,:,jj)^-0.5);
  end
  for jj=1:numtimesteps
    datacat_noisenorm{subcnt}(:,:,jj)=datacat{subcnt}(:,:,jj)*(sigma(:,:,jj)^-0.5);
  end
  
end % ii
switch numtimesteps
  case 25
    save([edir 'pcm_datacat.mat'],'datacat*','designX','partVec*');
  case 14
    save([edir 'pcm_datacat_NumTime14.mat'],'datacat*','designX','partVec*');
  case 17
    save([edir 'pcm_datacat_NumTime17.mat'],'datacat*','designX','partVec*');
end

return

%% create models

ncond=7;
plotflag=0;
% create_sim1models
% create_sim2_featAllModels % with Identity as 7 features, each diagonal element
% create_sim3_featAllModels  % with Identity as 1 feature = eye
% create_sim4_featAllModels  % with Identity as 4 features (in pairs) and also background all conditions
create_sim5_featAllModels  % with Identity as 4 features (in pairs) and also background all conditions

%% create 'beta value' average of trials

% load([edir 'pcm_datacat_NumTime14.mat']);
% numtimesteps=14;
load([edir 'pcm_datacat_NumTime17.mat']);
numtimesteps=17;

% for numtruse=[5 10 11 12 13 14 15],  % to test parameter space
% for numPart=2:5, % to test parameter space
for numtruse=9,
  for numPart=4,
    clear datacat_chunk* Xchunk partVecch numchunk iiNo iiUse
    for ii=1:22,iiNo(ii)=any(sum(designX{ii})<numPart*numtruse);end
    iiUse=setdiff(1:22,find(iiNo));
    iiInd=0;
    for ii=iiUse
      iiInd=iiInd+1;
      clear numchunk
      for nc=1:7
        induse=find(designX{ii}(:,nc));
        %     maxuse=floor(length(induse)/numtruse)*numtruse;
        numchunk(nc)=floor(length(induse)/numtruse);
        for ch=1:numchunk(nc)
          if nc==1
            datacat_chunk{iiInd}(ch,:,:)=mean(datacat_noisenorm{ii}(induse(((ch-1)*numtruse+1):((ch-1)*numtruse+numtruse)),:,:),1);
            Xchunk{iiInd}(ch,nc)=1;
          else
            %         if nc==2,keyboard,end;
            datacat_chunk{iiInd}(sum(numchunk(1:nc-1))+ch,:,:)=mean(datacat_noisenorm{ii}(induse(((ch-1)*numtruse+1):((ch-1)*numtruse+numtruse)),:,:),1);
            Xchunk{iiInd}(sum(numchunk(1:nc-1))+ch,nc)=1;
          end
        end
        
      end
      
      % multivariate noise normalisation
      sigma = compCovEEG(datacat_chunk{iiInd},Xchunk{iiInd});
      for jj=1:numtimesteps
        datacat_chunk_nn{iiInd}(:,:,jj)=datacat_chunk{iiInd}(:,:,jj)*(sigma(:,:,jj)^-0.5);
      end
      
      partVecch{iiInd}=[rem(1:size(datacat_chunk_nn{iiInd},1),numPart)+1]';  % this is nonequal number of trials per partition & condition; okay?
      for jj=1:numtimesteps
        Gemp_ch(:,:,jj,iiInd)=cov([Xchunk{iiInd}'*datacat_chunk_nn{iiInd}(:,:,jj)]',1);
        GempCVch_p4(:,:,jj,iiInd)=pcm_estGCrossval(datacat_chunk_nn{iiInd}(:,:,jj),partVecch{iiInd},Xchunk{iiInd});
      end
      
      
    end
    
    for nc=1:7,for jj=1:numtimesteps,for ii=1:numel(iiUse),gcvdiag(:,jj,ii)=GempCVch_p4(nc,nc,jj,ii);end;end;end
    lngcvd(numtruse,numPart)=length(find(gcvdiag(:)<0));
    
    % (this is for old 14 timepoints)
    % % % % % [5,4]: 16  [7,4]: 2  [8,4]:0  [10,4]: 0
    % % % % % [8,2]:0  [8,3]:0  [8,4]:0 [8,5]:0 [8,6]:0  [8,7]:0
    % % % % % [20,4]:0
    for jj=1:numtimesteps,for ii=1:numel(iiUse),
        eigGCV=eig(GempCVch_p4(:,:,jj,ii));
        negeig(numtruse,numPart,jj,ii)=any(eigGCV(abs(eigGCV)>0.05*sum(abs(eigGCV)))<0);
      end
    end
    % (this is for old 14 timepoints)
    % % % % % [8,3]: 307   [8,4]: 308  [8,5]: 307  [8,6]: 306  [8,7]: 308
    % % % % % [10,3]: 281  [10,4]:269  [10,5]: 262 [10,6]:268
    % % % % % [12,3]:183  [15,4]: 67  [18,4]:9  [20,4]:3  [20,3]:23  [20,5]:22
    % % % % % [19,4]:17
  end
end
% lngcvd([5 10:15],2:5)
%     57    51    51    48
%     12     9    11    10
%      7     4     5     5
%      3     0     0     1
%      2     1     2     1
%      1     1     1     1
%      0     0     0     0
%      0
%      0
%      0
% sum(sum(negeig([5 10:15],2:5,:,:),4),3)
%    374   374   373   373
%    270   250   247   241
%    233   208   201   200
%    162   155   140   138
%    120   110   100    91
%     88    76    78    72
%     35    26    29    27
%     26
%     24
%      3


% % Conclusion: use numtruse=15, and numPart=3;
% BUT...subjects 10, 20, 21, 22 have only 2 partitions in some conditions
% % Final Conclusion: use numtruse=18, and numPart=2;

% save([edir 'pcm_datacat_nt14_nc' num2str(numtruse) '_np' num2str(numPart) '.mat'],'datacat_chunk_nn','Xchunk','partVecch');
save([edir 'pcm_datacat_nt17_nc' num2str(numtruse) '_np' num2str(numPart) '.mat'],'datacat_chunk_nn','Xchunk','partVecch');



% timesteps=50:30:440;
timesteps=50:25:450;
nfig=2;
set(0,'DefaultFigurePosition',[15 450 1300 350+100*(nfig-3)])
jjtest=1:numtimesteps;
% jjtest=[1 6 12];
jjcnt=0;
for jj=jjtest
  jjcnt=jjcnt+1;
  figure(300);
  subplot(nfig,numel(jjtest),jjcnt);imagesc(mean(Gemp_ch(:,:,jj,:),4));caxis([0.5*min(Gemp_ch(:)) 0.4*max(Gemp_ch(:))]);title([num2str(timesteps(jj)) ' ms'])
  if jj==1,ylabel(['Gemp ' num2str(numtruse)]);end
  subplot(nfig,numel(jjtest),numel(jjtest)+jjcnt);imagesc(mean(GempCVch_p4(:,:,jj,:),4));caxis([1*min(GempCVch_p4(:)) .3*max(GempCVch_p4(:))]);
  if jj==1,ylabel(['GempCV ' num2str(numPart)]);end
end

% Plotting individual G matrices
set(0,'DefaultFigurePosition',[15 50 1300 900])
figure(1);nn=0;for jj=1:numtimesteps,for ii=1:11, nn=nn+1;subplot(numtimesteps,11,nn);imagesc(Gemp_ch(:,:,jj,ii));caxis([0 0.45*max(Gemp_ch(:))]);if ii==1 && jj==1,title(['chunk ' num2str(numtruse)]);end; end;end
figure(2);nn=0;for jj=1:numtimesteps,for ii=12:22,nn=nn+1;subplot(numtimesteps,11,nn);imagesc(Gemp_ch(:,:,jj,ii));caxis([0 0.45*max(Gemp_ch(:))]);if ii==12 && jj==1,title(['chunk ' num2str(numtruse)]);end;end;end
figure(3);nn=0;for jj=1:numtimesteps,for ii=1:11, nn=nn+1;subplot(numtimesteps,11,nn);imagesc(GempCVch_p4(:,:,jj,ii));caxis([0.45*min(GempCVch_p4(:)) 0.45*max(GempCVch_p4(:))]);if ii==1 && jj==1,title(['chunk ' num2str(numtruse)]);end;end;end
figure(4);nn=0;for jj=1:numtimesteps,for ii=12:22,nn=nn+1;subplot(numtimesteps,11,nn);imagesc(GempCVch_p4(:,:,jj,ii));caxis([0.45*min(GempCVch_p4(:)) 0.45*max(GempCVch_p4(:))]);if ii==12 && jj==1,title(['chunk ' num2str(numtruse)]);end;end;end



%% call Model Group Fit

nC=18;nP=2;
% nC=9;nP=4;

% load([edir 'pcm_datacat.mat']);
% load([edir 'pcm_datacat_NumTime14.mat']);
% load([edir 'pcm_datacat_nt14_nc' num2str(nC) '_np' num2str(nP) '.mat'],'datacat_chunk_nn','Xchunk','partVecch');
load([edir 'pcm_datacat_nt17_nc' num2str(nC) '_np' num2str(nP) '.mat'],'datacat_chunk_nn','Xchunk','partVecch');
% load([edir 'pcm_models.mat']);
% load([edir 'sim2_featmodels.mat']);
% load([edir 'sim3_featmodels.mat']);
% load([edir 'sim4_featmodels.mat']);
% load([edir 'sim4b_featmodels.mat']);
% load([edir 'sim4bNR_featmodels.mat']);
% load([edir 'sim4dNR_featmodels.mat']);
load([edir 'sim4eNR_featmodels.mat']);
% load([edir 'sim4fNR_featmodels.mat']);
% load([edir 'sim4gNR_featmodels.mat']);
% load([edir 'sim5_featmodels.mat']);
% load([edir 'sim6_featmodels.mat']);
% load([edir 'sim7_featmodels.mat']);
% load([edir 'sim7NR_featmodels.mat']);
numtimesteps=size(datacat_chunk_nn{1},3);
H = eye(7)-ones(7)/7;

% % new partVec for proper cross-val
% for subcnt=1:length(datacat)
%   partVec{subcnt}=repmat([1 2],[1 floor(size(datacat{subcnt},1)/2)])';
%   datacat{subcnt}=datacat{subcnt}(1:2*floor(size(datacat{subcnt},1)/2),:,:);
%   datacat_meansubt{subcnt}=datacat_meansubt{subcnt}(1:2*floor(size(datacat{subcnt},1)/2),:,:);
%   datacat_meansubt_noisenorm{subcnt}=datacat_meansubt_noisenorm{subcnt}(1:2*floor(size(datacat{subcnt},1)/2),:,:);
%   designX{subcnt}=designX{subcnt}(1:2*floor(size(datacat{subcnt},1)/2),:);
%   numobs(subcnt)=size(designX{subcnt},1);
% end
% numobsall=sum(numobs);
for subcnt=1:length(datacat_chunk_nn)
  try
    numobs(subcnt)=size(Xchunk{subcnt},1);
  catch
    numobs(subcnt)=size(designX{subcnt},1);
  end
end
numobsall=sum(numobs);

% compare methods for obtaining G (empirical, CV, fit)
Mind{1}=M{1};
Mind=M;
nM=numel(M);

% if 0
%   for jj=1:numtimesteps
%     for ii=1:length(datacat_chunk_nn)
%       Y{ii}=datacat_chunk_nn{ii}(:,:,jj); % keep mean in
%     end
%     [pcmind{jj}.T,pcmind{jj}.theta_hat,pcmind{jj}.G_pred]=pcm_fitModelIndivid(Y,Mind,partVecch,Xchunk,'runEffect','fixed');
%     %     [pcmind{jj}.D,pcmind{jj}.Tcross,pcmind{jj}.thetaCV]=pcm_fitModelIndividCrossval(Y,Mind,partVecch,Xchunk,'runEffect','fixed');
%   end
%   save([pcdir 'real_mean_kept_feature/pcmind.mat'],'pcmind','Mind')
%   load([pcdir 'real_mean_kept_feature/pcmF_ind5.mat'])
% end

for jj=1:numtimesteps
  %     load([pcdir 'real_mean_kept_feature_14TS_eye1/pcmF_datafit_p5_TS14_Meye1_nc' num2str(nC) 'np' num2str(nP) '_jj' num2str(jj) '.mat'])
  
  %     load([pcdir 'real_mean_kept_feature_17TS_eye1/pcmF_datafit_TS17_Meye1_nc' num2str(nC) 'np' num2str(nP) '_jj' num2str(jj) '.mat'])
%   load([pcdir 'real_mean_kept_feature_17TS_eye1/pcmF_datafit_TS17_Msim4bNR_nc' num2str(nC) 'np' num2str(nP) '_jj' num2str(jj) '.mat'])
%   load([pcdir 'real_mean_kept_feature_17TS_eye1/pcmF_datafit_TS17_Msim4bNRrand_nc' num2str(nC) 'np' num2str(nP) '_jj' num2str(jj) '.mat'])
%   load([pcdir 'real_mean_kept_feature_17TS_eye1/pcmF_datafit_TS17_Msim4dNRrand_nc' num2str(nC) 'np' num2str(nP) '_jj' num2str(jj) '.mat'])
  load([pcdir 'real_mean_kept_feature_17TS_eye1/pcmF_datafit_TS17_Msim4eNRrand_nc' num2str(nC) 'np' num2str(nP) '_jj' num2str(jj) '.mat'])
%   load([pcdir 'real_mean_kept_feature_17TS_eye1/pcmF_datafit_TS17_Msim4fNRrand_nc' num2str(nC) 'np' num2str(nP) '_jj' num2str(jj) '.mat'])
%   load([pcdir 'real_mean_kept_feature_17TS_eye1/pcmF_datafit_TS17_Msim4gNRrand_nc' num2str(nC) 'np' num2str(nP) '_jj' num2str(jj) '.mat'])
%         load([pcdir 'real_mean_kept_feature_17TS_eye1/pcmF_datafit_TS17_Msim4b_nc' num2str(nC) 'np' num2str(nP) '_jj' num2str(jj) '.mat'])
  %     load([pcdir 'real_mean_kept_feature_17TS_eye1/pcmF_datafit_TS17_Msim7_nc' num2str(nC) 'np' num2str(nP) '_jj' num2str(jj) '.mat'])
%     load([pcdir 'real_mean_kept_feature_17TS_eye1/pcmF_datafit_TS17_Msim7NR_nc' num2str(nC) 'np' num2str(nP) '_jj' num2str(jj) '.mat'])
%     load([pcdir 'real_mean_kept_feature_17TS_eye1/pcmF_datafit_TS17_Msim7NRrand_nc' num2str(nC) 'np' num2str(nP) '_jj' num2str(jj) '.mat'])
  pcmfit{nC,nP}{jj}=A;    clear A
  M=Mind;
  
  %     load([pcdir 'real_mean_kept_feature_14TS/pcmF_datafit_p5_TS14_Meye_nc5np4_jj' num2str(jj) '.mat'])
  %     pcmfit{nC,nP}{jj}=A;    clear A
  %     load([pcdir 'real_mean_kept_feature_14TS/pcmF_datafit_p5_TS14_Meye_jj' num2str(jj) '.mat'])
  %     pcmfite5{jj}=A;    clear A
  
  %     load([pcdir 'real_mean_kept_feature_14TS/pcmF_ind2_jj' num2str(jj) '.mat'])
  %     pcmind2{jj}=A{jj};
  %     load([pcdir 'real_mean_kept_feature_14TS/pcmF_ind3_jj' num2str(jj) '.mat'])
  %     pcmind3{jj}=A{jj};
  %     load([pcdir 'real_mean_kept_feature_14TS/pcmF_ind5_jj' num2str(jj) '.mat'])
  %     pcmind5{jj}=A{jj};
  
  %     load([pcdir 'real_mean_kept_feature_14TS_eye1/pcmF_ind5_Meye1_nc5np4_jj' num2str(jj) '.mat'])
  
  %     load([pcdir 'real_mean_kept_feature_17TS_eye1/pcmF_ind17_Meye1_nc18np2_jj' num2str(jj) '.mat'])
  %     load([pcdir 'real_mean_kept_feature_17TS_eye1/pcmF_ind17_Msim4b_nc18np2_jj' num2str(jj) '.mat'])
  %     load([pcdir 'real_mean_kept_feature_17TS_eye1/pcmF_ind17_Msim4bCV_nc18np2_jj' num2str(jj) '.mat'])
  %     load([pcdir 'real_mean_kept_feature_17TS_eye1/pcmF_ind17_Msim4bNR_nc18np2_jj' num2str(jj) '.mat'])
  %     load([pcdir 'real_mean_kept_feature_17TS_eye1/pcmF_ind17_Msim7_nc18np2_jj' num2str(jj) '.mat'])
  %     load([pcdir 'real_mean_kept_feature_17TS_eye1/pcmF_ind17_Msim7NR_nc18np2_jj' num2str(jj) '.mat'])
  %     pcmind{nC,nP}{jj}=A{jj};    clear A
  
  %     load([pcdir 'real_mean_kept_feature_14TS/pcmF_ind5_Meye_nc5np4_jj' num2str(jj) '.mat'])
  %     pcmind{nC,nP}{jj}=A{jj};    clear A
  
  %     load([pcdir 'real_mean_kept_feature_14TS/pcmF_ind5_Meye_jj' num2str(jj) '.mat'])
  %     pcminde5{jj}=A{jj};    clear A
  
  % M=Mind;
  %     for mm=1:nM
  %       clear GpredCV
  %       for nn=1:nP
  %         thetaCV=pcmind{nC,nP}{jj}.thetaCV{mm}(nn:nP:end,:);  % subj x numfeat
  %         for ss=1:size(thetaCV,1)
  %           GpredCV(:,:,ss,nn)=pcm_calculateG(M{mm},thetaCV(ss,:)');
  %         end
  %       end
  %       pcmind{nC,nP}{jj}.G_predCV{mm}=mean(GpredCV,4);
  %     end
  
end

designX=Xchunk;
datacat_noisenorm=datacat_chunk_nn;


if 0
  for jj=1:numtimesteps
    for ii=1:length(datacat_noisenorm)
      Y{ii}=datacat_noisenorm{ii}(:,:,jj); % keep mean in
    end
    [pcmind{jj}.D,pcmind{jj}.Tcross,pcmind{jj}.thetaCV]=pcm_fitModelIndividCrossval(Y,Mind,partVec5,designX,'runEffect','fixed');
  end
end


% sim4
% mtest=[3:10];
% mtest=[11:17];
% mtest=[3 11:17];

% sim4b
% mtest=[3:18];

% sim4bNR
% mtest=[3:18];
% mtest=[3:10];

% % sim4dNR
% mtest=[3:18];
% mtest=[3:10];

% % sim4eNR
% mtest=[3:34];
mtest=[3 5 6 9 28 31 32 34];

% % sim4fNR
% mtest=[3 7 10 14 18 22 23 24 25 26 27 28];

% % sim4gNR
% mtest=[3:18];
% mtest=[3 5 6 9 12 15 16 18];

% sim5
% mtest=[20 23:30];
% mtest=[23:25 27 28 30];
% mtest=30;

% % sim6
% mtest=[4:18];

% % sim7
% mtest=[3:257];
% for mm=3:257,mmExclude(mm)=any(M{mm}.featincl==8);end
% mmNo8=find(~mmExclude);
% mtest=mmNo8(3:end);

% % Plotting G
for jj=1:numtimesteps
  for iireal=1:length(datacat_noisenorm)
    G_emp(:,:,jj,iireal)=cov([designX{iireal}'*datacat_noisenorm{iireal}(:,:,jj)]',1);
    %     G_empCV2(:,:,jj,iireal)=pcm_estGCrossval(datacat_noisenorm{iireal}(:,:,jj),partVec2{iireal},designX{iireal});
    %     G_empCV3(:,:,jj,iireal)=pcm_estGCrossval(datacat_noisenorm{iireal}(:,:,jj),partVec3{iireal},designX{iireal});
    G_empCV5(:,:,jj,iireal)=pcm_estGCrossval(datacat_noisenorm{iireal}(:,:,jj),partVecch{iireal},designX{iireal});
  end
end
% Gmodel
for mm=3:nM;Gmodel(:,:,mm)=pcm_calculateG(M{mm},ones(1,M{mm}.numGparams)');end

%%  T-test of GempCV
for jj=1:17,
for cc=1:7,
  GH(:,cc,jj)=ttest(squeeze(G_empCV5(:,cc,jj,:))',zeros(size(squeeze(G_empCV5(:,cc,jj,:))')),'alpha',0.025,'tail','right')';
end
end



%%
% Plotting individual G matrices
set(0,'DefaultFigurePosition',[15 50 1300 1300])
figure;nn=0;for jj=1:numtimesteps,for ii=1:22,nn=nn+1;subplot(numtimesteps,22,nn);imagesc(G_emp(:,:,jj,ii));caxis([-2 10]);end;end
figure;nn=0;for jj=1:numtimesteps,for ii=1:22,nn=nn+1;subplot(numtimesteps,22,nn);imagesc(G_empCV5(:,:,jj,ii));caxis([-0.3 1.5]);end;end

figure;nn=0;for jj=1:numtimesteps,for ii=1:22,nn=nn+1;subplot(numtimesteps,22,nn);imagesc(G_empCV5(:,:,jj,ii));end;end

jj=3;figure;
for ss=1:22,subplot(2,22,ss);imagesc(G_empCV5(:,:,jj,ss));end
for ss=1:22,subplot(2,22,ss+22);imagesc(pcmfit{nC,nP}{jj}.G_predcv{1}(:,:,ss));end
% for ss=1:22,subplot(2,22,ss+22);imagesc(pcmind{nC,nP}{jj}.G_predCV{1}(:,:,ss));end

figure;nn=0;for jj=[1 6 12],for ii=1:22,nn=nn+1;subplot(3,22,nn);imagesc(G_emp(:,:,jj,ii));caxis([-10 60]);end;end

figure;
for mm=2:17,
  if mm==2,subplot(2,8,mm-1);imagesc(Gmodel(:,:,mm));end
  if mm==3,subplot(2,8,9);imagesc(Gmodel(:,:,mm));end
  if mm<11 && mm>3,subplot(2,8,mm-2);imagesc(Gmodel(:,:,mm));end
  if mm>10 && mm>3,subplot(2,8,8+mm-9);imagesc(Gmodel(:,:,mm));end
end

%% Plotting models for sim4bNR

% figure;for ff=1:size(M{10}.Ac,3),subplot(4,4,ff);imagesc(M{10}.Ac(:,:,ff)*M{10}.Ac(:,:,ff)');end

figure(53);
ff=1;subplot(4,4,1);imagesc(M{10}.Ac(:,ff,ff)*M{10}.Ac(:,ff,ff)');ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
  if ff==1,ylabel('Identity');end
ff=2;subplot(4,4,2);imagesc(M{10}.Ac(:,ff,ff)*M{10}.Ac(:,ff,ff)');ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
ff=7;subplot(4,4,4);imagesc(M{10}.Ac(:,ff,ff)*M{10}.Ac(:,ff,ff)');ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
ffind=4;for ff=[15:16 ];ffind=ffind+1;subplot(4,4,ffind);imagesc(M{10}.Ac(:,ff,ff)*M{10}.Ac(:,ff,ff)');ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
  if ffind==5,ylabel('TIW');end
end;
ffind=8;for ff=[11:14 ];ffind=ffind+1;subplot(4,4,ffind);imagesc(M{10}.Ac(:,ff,ff)*M{10}.Ac(:,ff,ff)');ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
  if ffind==9,ylabel('Asymmetry');end
end;
ffind=12;for ff=[8:10];ffind=ffind+1;subplot(4,4,ffind);imagesc(M{10}.Ac(:,ff,ff)*M{10}.Ac(:,ff,ff)');ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
  if ffind==13,ylabel('Sym. Pairs');end
end;

print(53,[fdir 'pcmmodels4bNR.eps'],'-depsc')
print(53,[fdir 'pcmmodels4bNR.png'],'-dpng')


%% Plotting models for sim4eNR

% figure;for ff=1:size(M{10}.Ac,3),subplot(4,4,ff);imagesc(M{10}.Ac(:,:,ff)*M{10}.Ac(:,:,ff)');end

figure(54);
ff=1;subplot(4,4,1);imagesc(M{10}.Ac(:,ff,ff)*M{10}.Ac(:,ff,ff)');ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
  if ff==1,ylabel('Identity');end
ff=2;subplot(4,4,2);imagesc(M{10}.Ac(:,ff,ff)*M{10}.Ac(:,ff,ff)');ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
ff=7;subplot(4,4,4);imagesc(M{10}.Ac(:,ff,ff)*M{10}.Ac(:,ff,ff)');ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
ffind=4;for ff=[15:16 ];ffind=ffind+1;subplot(4,4,ffind);
  imagesc(M{10}.Ac(:,ff,ff)*M{10}.Ac(:,ff,ff)');ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
  if ffind==5,ylabel('TIW');end
end;
subplot(4,4,7);imagesc(ones(7,7));caxis([0 1]);ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
ffind=8;for ff=[11:14 ];ffind=ffind+1;subplot(4,4,ffind);imagesc(M{10}.Ac(:,ff,ff)*M{10}.Ac(:,ff,ff)');ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
  if ffind==9,ylabel('Asymmetry');end
end;
ffind=12;for ff=[8:10];ffind=ffind+1;subplot(4,4,ffind);imagesc(M{10}.Ac(:,ff,ff)*M{10}.Ac(:,ff,ff)');ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
  if ffind==13,ylabel('Sym. Pairs');end
end;

print(54,[fdir 'pcmmodels4eNR.eps'],'-depsc')
print(54,[fdir 'pcmmodels4eNR.png'],'-dpng')

%% Plotting models for sim5

figure;for ff=8:21,subplot(7,2,ff-7);imagesc(M{20}.Ac(:,ff,ff)*M{20}.Ac(:,ff,ff)');end;

figure;
ffind=0;for ff=[8 10:13];ffind=ffind+1;subplot(4,5,ffind);imagesc(M{20}.Ac(:,ff,ff)*M{20}.Ac(:,ff,ff)');end;
ffind=6;for ff=[14:17 ];ffind=ffind+1;subplot(4,5,ffind);imagesc(M{20}.Ac(:,ff,ff)*M{20}.Ac(:,ff,ff)');end;
ffind=10;for ff=[9 18 19 ];ffind=ffind+1;subplot(4,5,ffind);imagesc(M{20}.Ac(:,ff,ff)*M{20}.Ac(:,ff,ff)');end;
ffind=16;for ff=[20 21 ];ffind=ffind+1;subplot(4,5,ffind);imagesc(M{20}.Ac(:,ff,ff)*M{20}.Ac(:,ff,ff)');end;





%%  AIC BIC Groupfit
clear aic*
G_predG_mix=zeros(7,7,numtimesteps);
for jj=1:numtimesteps
  for mm=1:nM
    for subcnt=1:length(datacat_noisenorm)
      [aicg(mm,jj,subcnt),bicg(mm,jj,subcnt)] = aicbic(pcmfit{nC,nP}{jj}.Tgroup.likelihood(subcnt,mm),M{mm}.numGparams,numobs(subcnt));
    end
  end
  maicg=mean(aicg(mtest,jj,:),3);
  [mmaicg]=min(maicg);
  [tmp,mnindg(jj,:)]=min(squeeze(aicg(mtest,jj,:)),[],1);
  aicgdif(:,jj)=maicg-mmaicg; %
  aicgrat(:,jj)=exp(.5*(mmaicg-maicg)); %
  
  aicdenomG=sum(exp(-0.5*aicgdif(:,jj)));
  mmuseG{jj}=[find(aicgdif(:,jj)<3)];
  for mmu=1:numel(mmuseG{jj})
    G_predG_mix(:,:,jj)=G_predG_mix(:,:,jj)+  pcmfit{nC,nP}{jj}.Gpred{mtest(mmuseG{jj}(mmu))};
    G_predG_mixw(:,:,jj)=G_predG_mix(:,:,jj)+ [exp(-0.5*aicgdif(mmuseG{jj}(mmu),jj))/aicdenomG]* mean(pcmfit{nC,nP}{jj}.Gpred{mtest(mmuseG{jj}(mmu))},3);
  end
end

%% AIC individual
clear aic*
G_predI_mixw=zeros(7,7,numtimesteps);
G_predI_mix=zeros(7,7,numtimesteps);
for jj=1:numtimesteps
  for mm=1:nM
    for subcnt=1:length(datacat_noisenorm)
      [aici(mm,jj,subcnt),bici(mm,jj,subcnt)] = aicbic(pcmind{nC,nP}{jj}.Tgroup.likelihood(subcnt,mm),M{mm}.numGparams,numobs(subcnt));
    end
  end
  maici=mean(aici(mtest,jj,:),3);
  [mmaici]=min(maici);
  [tmp,mnindi(jj,:)]=min(squeeze(aici(mtest,jj,:)),[],1);
  aicidif(:,jj)=maici-mmaici; %
  aicirat(:,jj)=exp(.5*(mmaici-maici)); %
  
  aicdenomI=sum(exp(-0.5*aicidif(:,jj)));
  mmuseI{jj}=[find(aicidif(:,jj)<3)];   % 1.9
  for mmu=1:numel(mmuseI{jj})
    G_predI_mixw(:,:,jj)=G_predI_mix(:,:,jj)+ [exp(-0.5*aicidif(mmuseI{jj}(mmu),jj))/aicdenomI]* mean(pcmind{nC,nP}{jj}.G_pred{mtest(mmuseI{jj}(mmu))},3);
    G_predI_mix(:,:,jj)=G_predI_mix(:,:,jj)+  mean(pcmind{nC,nP}{jj}.G_pred{mtest(mmuseI{jj}(mmu))},3);
  end
end

%% plot likelihood

figure(123);
for jj=1:17,
  subplot(3,6,jj);pcm_plotModelLikelihood(pcmfit{nC,nP}{jj}.Tcross,{M{[1 2 mtest]}},'normalize',1,'Nceil',1,'Nnull',2,'upperceil',pcmfit{nC,nP}{jj}.Tgroup.likelihood(:,1));
end

jj=7;pcm_plotModelLikelihood(pcmfit{nC,nP}{jj}.Tcross,{M{[1 2 mtest]}},'normalize',1,'Nceil',1,'Nnull',2,'upperceil',pcmfit{nC,nP}{jj}.Tgroup.likelihood(:,1));

%% plot likelihood  for paper
figure(124);  
for jj=1:17,
  subplot(6,3,jj);plotlike{jj}=pcm_plotModelLikelihood(pcmfit{nC,nP}{jj}.Tcross,{M{[1 2 mtest]}},'normalize',0,'Nceil',1,'Nnull',2,'upperceil',pcmfit{nC,nP}{jj}.Tgroup.likelihood(:,1));
  ylim([min(mean(plotlike{jj}.likelihood_norm(:,mtest),1)) max(mean(plotlike{jj}.likelihood_norm(:,mtest),1))+3])
  title([num2str(timesteps(jj)) ' ms'])
  [mx,mind]=max(mean(plotlike{jj}.likelihood_norm(:,mtest),1));
  hold on;plot(0:9,(mean(plotlike{jj}.likelihood_norm(:,mtest(mind)),1)-3)*ones(1,10));
end

figure(125);  
jj=7;
plotlike{jj}=pcm_plotModelLikelihood(pcmfit{nC,nP}{jj}.Tcross,{M{[1 2 mtest]}},'normalize',0,'Nceil',1,'Nnull',2,'upperceil',pcmfit{nC,nP}{jj}.Tgroup.likelihood(:,1));
  ylim([min(mean(plotlike{jj}.likelihood_norm(:,mtest),1))-3 max(mean(plotlike{jj}.likelihood_norm(:,mtest),1))+25])
  title([num2str(timesteps(jj)) ' ms'])
  ylabel('Relative log-likelihood')
  [mx,mind]=max(mean(plotlike{jj}.likelihood_norm(:,mtest),1));
  hold on;plot(0:9,(mean(plotlike{jj}.likelihood_norm(:,mtest(mind)),1)-3)*ones(1,10));

% mtest
  set(0,'DefaultFigurePosition',[15 150 500 400])
  for jj=1:numtimesteps
    figure(300+jj);
    Tcrosstmp.likelihood=pcmfit{nC,nP}{jj}.Tcross.likelihood(:,[1 mtest]);
    Tcrosstmp.SN=pcmfit{nC,nP}{jj}.Tcross.SN;
    plotlike{jj}=pcm_plotModelLikelihood(Tcrosstmp,{M{[1 mtest]}},'normalize',0,'Nceil',1,'Nnull',2,'upperceil',pcmfit{nC,nP}{jj}.Tgroup.likelihood(:,1));
%     close(300+jj);
%     figure(300+jj);    bar(mean(plotlike{jj}.likelihood_norm(:,3:end),1));
    ylim([min(mean(plotlike{jj}.likelihood_norm(:,2:end),1))-3 max(mean(plotlike{jj}.likelihood_norm(:,2:end),1))+10])
    title([num2str(timesteps(jj)) ' ms'])
    ylabel('Relative log-likelihood')
    [mx,mind]=max(mean(plotlike{jj}.likelihood_norm(:,2:end),1));
    hold on;plot(0:[length(mtest)+1],(mean(plotlike{jj}.likelihood_norm(:,mind+1),1)-3)*ones(1,length(mtest)+2),'LineWidth',2);
%     hold on;plot(0:[length(mtest)+1],(mean(plotlike{jj}.likelihood_norm(:,9),1)-3)*ones(1,length(mtest)+2),'LineWidth',2);
    %   ax=get(124,'Children');set(ax,'FontSize',16)
  end

set(0,'DefaultFigurePosition',[15 150 500 400])
figure(126);  
jj=7;
Tcrosstmp.likelihood=pcmfit{nC,nP}{jj}.Tcross.likelihood(:,[1 mtest]);
Tcrosstmp.SN=pcmfit{nC,nP}{jj}.Tcross.SN;
plotlike{jj}=pcm_plotModelLikelihood(Tcrosstmp,{M{[1 mtest]}},'normalize',0,'Nceil',1,'Nnull',2,'upperceil',pcmfit{nC,nP}{jj}.Tgroup.likelihood(:,1));
  ylim([min(mean(plotlike{jj}.likelihood_norm(:,2:end),1))-3 max(mean(plotlike{jj}.likelihood_norm(:,2:end),1))+10])
  title([num2str(timesteps(jj)) ' ms'])
  ylabel('Relative log-likelihood')
  [mx,mind]=max(mean(plotlike{jj}.likelihood_norm(:,2:end),1));
  hold on;plot(0:8,(mean(plotlike{jj}.likelihood_norm(:,mind+1),1)-3)*ones(1,9),'LineWidth',2);
  ax=get(126,'Children');set(ax,'FontSize',16)
  
  print(126,[fdir 'pcm_relloglike.eps'],'-depsc')
  print(126,[fdir 'pcm_relloglike.png'],'-dpng')

set(0,'DefaultFigurePosition',[15 150 500 400])
figure(127);  
jj=11;
Tcrosstmp.likelihood=pcmfit{nC,nP}{jj}.Tcross.likelihood(:,[1 11:18]);
Tcrosstmp.SN=pcmfit{nC,nP}{jj}.Tcross.SN;
plotlike{jj}=pcm_plotModelLikelihood(Tcrosstmp,{M{[1 11:18]}},'normalize',0,'Nceil',1,'Nnull',2,'upperceil',pcmfit{nC,nP}{jj}.Tgroup.likelihood(:,1));
  ylim([min(mean(plotlike{jj}.likelihood_norm(:,2:8),1))-3 max(mean(plotlike{jj}.likelihood_norm(:,2:8),1))+10])
  title([num2str(timesteps(jj)) ' ms'])
  ylabel('Relative log-likelihood')
  [mx,mind]=max(mean(plotlike{jj}.likelihood_norm(:,2:end),1));
  hold on;plot(0:8,(mean(plotlike{jj}.likelihood_norm(:,mind+1),1)-3)*ones(1,9),'LineWidth',2);
%   ax=get(127,'Children');set(ax,'FontSize',16)
  
  print(127,[fdir 'pcm_relloglike.eps'],'-depsc')
  print(127,[fdir 'pcm_relloglike.png'],'-dpng')


  
figure;for mm=3:10,
  subplot(2,4,mm-2);plot(pcmfit{nC,nP}{jj}.Tcross.likelihood(:,mm)-pcmfit{nC,nP}{jj}.Tcross.likelihood(:,2),'o');
  ylim([0 500]);
end
figure;for mm=3:10,
  subplot(2,4,mm-2);plot(pcmfit{nC,nP}{jj}.Tcross.likelihood(:,mm)-pcmfit{nC,nP}{jj}.Tcross.likelihood(:,3),'o');
  ylim([0 120]);
end
figure;for mm=3:10,
  subplot(2,4,mm-2);plot(pcmfit{nC,nP}{jj}.Tcross.likelihood(:,mm)-pcmfit{nC,nP}{jj}.Tcross.likelihood(:,6),'o');
  ylim([0 60]);
end


%%
timesteps=50:25:450;
  nfig=3;
set(0,'DefaultFigurePosition',[15 150 1300 300+100*(nfig-3)])
for jj=1:numtimesteps,
  mxg1(jj)=max(max(mean(pcmfit{nC,nP}{jj}.G_predcv{1},3)));
  try
    [mx,mxindG(jj)]=max(mean(pcmfit{nC,nP}{jj}.Tcross.likelihood(:,mtest)-repmat(pcmfit{nC,nP}{jj}.Tcross.likelihood(:,2),[1 numel(mtest)])));
  end
  try
    [mx,mxindGnonCV(jj)]=max(mean(pcmfit{nC,nP}{jj}.Tgroup.likelihood(:,mtest)-repmat(pcmfit{nC,nP}{jj}.Tgroup.likelihood(:,2),[1 numel(mtest)])));
  end
  for ss=1:22
    [mx,mxindGeach(jj,ss)]=max(pcmfit{nC,nP}{jj}.Tcross.likelihood(ss,mtest)-pcmfit{nC,nP}{jj}.Tcross.likelihood(ss,2));
  end
  [mxsortG{jj},sortindG{jj}]=sort(mean(pcmfit{nC,nP}{jj}.Tcross.likelihood(:,mtest)-repmat(pcmfit{nC,nP}{jj}.Tcross.likelihood(:,2),[1 numel(mtest)])));
  mxi3(jj)=mxsortG{jj}(end);
  mxg2(jj)=max(max(mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindG(jj))},3)));
  mxg2ncv(jj)=max(max(mean(pcmfit{nC,nP}{jj}.Gpred{mtest(mxindGnonCV(jj))},3)));
end
for jj=1:numtimesteps,
  figure(150);
  %   subplot(nfig,numtimesteps,jj);imagesc(mean(G_emp(:,:,jj,:),4));
  %   if jj==1,ylabel(['Gemp']);end
  subplot(nfig,numtimesteps,jj);imagesc(mean(G_empCV5(:,:,jj,:),4));  
  caxis([min(min(min(mean(G_empCV5,4)))) max(max(max(mean(G_empCV5,4))))])
  title([num2str(timesteps(jj)) ' ms'])
  ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
  if jj==1,ylabel(['GempCV']);end
  %   subplot(nfig,numtimesteps,numtimesteps+jj);imagesc(0.5*mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(sortindG{jj}(end-2))},3)+0.8*mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(sortindG{jj}(end-1))},3)+mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindG(jj))},3));  caxis([0 1.4*max(mxg2)])
  %   if jj==1,ylabel(['Gfit3 CV']);end
  %   subplot(nfig,numtimesteps,2*numtimesteps+jj);imagesc(0.8*mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(sortindG{jj}(end-1))},3)+mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindG(jj))},3));  caxis([0 1.4*max(mxg2)])
  %   if jj==1,ylabel(['Gfit2CV']);end
  %   subplot(nfig,numtimesteps,1*numtimesteps+jj);imagesc(mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindG(jj))},3));  caxis([min(mxg2) 1*max(mxg2)])
  subplot(nfig,numtimesteps,1*numtimesteps+jj);imagesc(mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindG(jj))},3));
%   if numel(M)==257
%     caxis([-0.05+min(min(mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindG(jj))},3))) 1.1*max(max(mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindG(jj))},3)))])
%   else
%     caxis([0+min(min(mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindG(jj))},3))) 0.6*max(max(mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindG(jj))},3)))])
%   end
  caxis([0+min(min(mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindG(jj))},3))) 1*max(max(mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindG(jj))},3)))])
  ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
  if jj==1,ylabel(['GmodelfitCV']);end

%   subplot(nfig,numtimesteps,2*numtimesteps+jj);imagesc(mean(pcmfit{nC,nP}{jj}.Gpred{mtest(mxindG(jj))},3));
%   caxis([0+min(min(mean(pcmfit{nC,nP}{jj}.Gpred{mtest(mxindG(jj))},3))) 1*max(max(mean(pcmfit{nC,nP}{jj}.Gpred{mtest(mxindG(jj))},3)))])
%   ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
%   if jj==1,ylabel(['GfitSame']);end

%   for ss=1:22,
%     Gpredcvplot(:,:,ss)=pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindGeach(jj,ss))}(:,:,ss);
%   end
%   %   subplot(nfig,numtimesteps,2*numtimesteps+jj);imagesc(mean(Gpredcvplot,3));  caxis([0 0.7*max(mxg2)])
%   subplot(nfig,numtimesteps,2*numtimesteps+jj);imagesc(mean(Gpredcvplot,3)); 
% %   if numel(M)==257
% %     caxis([-0.05+min(min(mean(Gpredcvplot,3))) 1.2*max(max(mean(Gpredcvplot,3)))])
% %   else
% %     caxis([0+min(min(mean(Gpredcvplot,3))) 0.6*max(max(mean(Gpredcvplot,3)))])
% %   end
%     caxis([0+min(min(mean(Gpredcvplot,3))) 1*max(max(mean(Gpredcvplot,3)))])
%   ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
%   if jj==1,ylabel(['GfitCVeach']);end

  %   subplot(nfig,numtimesteps,28+jj);imagesc(Gmodel(:,:,mtest(mxindG(jj))));
  %   if jj==1,ylabel(['Gmodel win']);end
  %   subplot(nfig,numtimesteps,42+jj);imagesc(mean(pcmfit{nC,nP}{jj}.Gpred{mode(mnindg(jj,:))},3));
  %   if jj==1,ylabel(['Gfit AICmin']);end
%   subplot(nfig,numtimesteps,3*numtimesteps+jj);imagesc(G_predG_mixw(:,:,jj));  caxis([min(G_predG_mixw(:)) 0.3*max(G_predG_mixw(:))])
%   ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
%   if jj==1,ylabel(['Gfit AICmix']);end
  subplot(nfig,numtimesteps,2*numtimesteps+jj);imagesc(mean(pcmfit{nC,nP}{jj}.G_predcv{1},3)); 
%   caxis([min(mxg1) 1*max(mxg1)])
  ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
  if jj==1,ylabel(['GFreeCV']);end
%   if numel(M)==257
% %     subplot(nfig,numtimesteps,4*numtimesteps+jj);imagesc(mean(pcmfit{nC,nP}{jj}.G_predcv{257},3)); % caxis([0 max(mxg1)])
% %     ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
% %     if jj==1,ylabel(['GfitCV257']);end
%     subplot(nfig,numtimesteps,4*numtimesteps+jj);imagesc(mean(pcmfit{nC,nP}{jj}.G_predcv{249},3)); % caxis([0 max(mxg1)])
%     ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
%     if jj==1,ylabel(['GfitCV249']);end
%   end
  %   subplot(nfig,numtimesteps,84+jj);imagesc(mean(pcmfit{nC,nP}{jj}.Gpred{17},3));
  %   if jj==1,ylabel(['Gfit All']);end
end

%% For paper (200ms G matrix)

nfig=3;
set(0,'DefaultFigurePosition',[15 150 100 300+100*(nfig-3)])
for jj=7
  figure(160);
  subplot(nfig,1,1);imagesc(mean(G_empCV5(:,:,jj,:),4));  
  caxis([min(min(min(mean(G_empCV5,4)))) max(max(max(mean(G_empCV5,4))))])
  title([num2str(timesteps(jj)) ' ms'])
  ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
  if jj==7,ylabel(['GempCV']);end
  subplot(nfig,1,2);imagesc(mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindG(jj))},3));
  caxis([0+min(min(mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindG(jj))},3))) 1*max(max(mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindG(jj))},3)))])
  ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
  if jj==7,ylabel(['GmodelfitCV']);end
  subplot(nfig,1,3);imagesc(mean(pcmfit{nC,nP}{jj}.G_predcv{1},3)); 
  ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
  if jj==7,ylabel(['GFreeCV']);end
end


print(160,[fdir 'pcm_Gemp_Gfit_freechol.eps'],'-depsc')
print(160,[fdir 'pcm_Gemp_Gfit_freechol.png'],'-dpng')


nfig=2;
set(0,'DefaultFigurePosition',[15 150 100 300+100*(nfig-3)])
for jj=7
  figure(161);
  subplot(nfig,1,1);imagesc(mean(G_empCV5(:,:,jj,:),4));  
  caxis([min(min(min(mean(G_empCV5,4)))) max(max(max(mean(G_empCV5,4))))])
  title([num2str(timesteps(jj)) ' ms'])
  ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
  if jj==7,ylabel(['GempCV']);end
  subplot(nfig,1,2);imagesc(mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindG(jj))},3));
  caxis([0+min(min(mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindG(jj))},3))) 1*max(max(mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindG(jj))},3)))])
  ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
  if jj==7,ylabel(['GmodelfitCV']);end
end


print(161,[fdir 'pcm_Gemp_Gfit.eps'],'-depsc')
print(161,[fdir 'pcm_Gemp_Gfit.png'],'-dpng')


%% Difference of G-matrices
figure;
caxlim=[-0.7 0.7];
for jj=1:numtimesteps
  subplot(3,numtimesteps,0*numtimesteps+jj);imagesc(mean(G_empCV5(:,:,jj,:),4)-mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindG(jj))},3));caxis(caxlim)
  if jj==1,ylabel('ECV-Fit');end
  title([num2str(timesteps(jj)) ' ms'])
  subplot(3,numtimesteps,1*numtimesteps+jj);imagesc(mean(G_empCV5(:,:,jj,:),4)-mean(pcmfit{nC,nP}{jj}.G_predcv{1},3));caxis(caxlim)
  if jj==1,ylabel('ECV-FC');end
  subplot(3,numtimesteps,2*numtimesteps+jj);imagesc(mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindG(jj))},3)-mean(pcmfit{nC,nP}{jj}.G_predcv{1},3));caxis(caxlim)
  if jj==1,ylabel('Fit-FC');end
  if jj==17,colorbar;end
end




%%  Sim4bNR for Joern
nfig=6;
set(0,'DefaultFigurePosition',[15 150 1900 350+100*(nfig-3)])
jj=find(timesteps==200);figure(jj);
subplot(nfig,24,1);imagesc(mean(G_empCV5(:,:,jj,:),4));ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};title('Group average')
%  caxis([min(min(min(mean(G_empCV5,4)))) 1*max(max(max(mean(G_empCV5,4))))])
subplot(nfig,24,1+1*24);imagesc(mean(pcmfit{nC,nP}{jj}.G_predcv{1},3)); ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
%  caxis([min(min(mean(pcmfit{nC,nP}{jj}.G_predcv{1},3))) 0.7*max(mxg1)])
subplot(nfig,24,1+2*24);imagesc(mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindG(jj))},3)); title(strrep(M{mtest(mxindG(jj))}.name(1:end-4), '_', ''));ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
%  caxis([0+min(min(mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindG(jj))},3))) 0.6*max(max(mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindG(jj))},3)))])
subplot(nfig,24,1+3*24);imagesc(mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindGeach(jj,ss))},3)); ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
%  caxis([0+min(min(mean(Gpredcvplot,3))) 0.5*max(max(mean(Gpredcvplot,3)))])


for ss=1:22,subplot(nfig,24,ss+2);imagesc(G_empCV5(:,:,jj,ss));  if ss==1,ylabel(['GempCV']);end;ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};title(['Sub ' num2str(ss)]);
%   caxis([min(min(min(mean(G_empCV5,4)))) 1*max(max(max(mean(G_empCV5,4))))])
end
for ss=1:22,subplot(nfig,24,ss+1*24+2);imagesc(pcmfit{nC,nP}{jj}.G_predcv{1}(:,:,ss)); if ss==1,ylabel(['GFreeCV']);end;ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
%   caxis([min(min(mean(pcmfit{nC,nP}{jj}.G_predcv{1},3))) 0.7*max(mxg1)])
end
for ss=1:22,subplot(nfig,24,ss+2*24+2);imagesc(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindG(jj))}(:,:,ss)); if ss==1,ylabel(['GfitCVsame']);end;ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
%   caxis([0+min(min(mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindG(jj))},3))) 0.6*max(max(mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindG(jj))},3)))]);
end
for ss=1:22,subplot(nfig,24,ss+3*24+2);imagesc(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindGeach(jj,ss))}(:,:,ss)); title(strrep(M{mtest(mxindGeach(jj,ss))}.name(1:end-4), '_', ''));if ss==1,ylabel(['GfitCVeach']);end;ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
%   caxis([0+min(min(mean(Gpredcvplot,3))) 0.5*max(max(mean(Gpredcvplot,3)))]);
end


%% Sim7NR for Joern
nfig=6;
set(0,'DefaultFigurePosition',[15 150 1900 350+100*(nfig-3)])
jj=find(timesteps==200);figure(jj);
subplot(nfig,24,1);imagesc(mean(G_empCV5(:,:,jj,:),4));ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};title('Group average')
%  caxis([min(min(min(mean(G_empCV5,4)))) 1*max(max(max(mean(G_empCV5,4))))])
subplot(nfig,24,1+1*24);imagesc(mean(pcmfit{nC,nP}{jj}.G_predcv{1},3)); ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
%  caxis([min(min(mean(pcmfit{nC,nP}{jj}.G_predcv{1},3))) 1*max(mxg1)])
subplot(nfig,24,1+2*24);imagesc(mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindG(jj))},3)); title(strrep(M{mtest(mxindG(jj))}.name(7:end), '_', ''));ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
%  caxis([0+min(min(mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindG(jj))},3))) 1*max(max(mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindG(jj))},3)))])
subplot(nfig,24,1+3*24);imagesc(mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindGeach(jj,ss))},3)); ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
%  caxis([0+min(min(mean(Gpredcvplot,3))) 1*max(max(mean(Gpredcvplot,3)))])
subplot(nfig,24,1+4*24);imagesc(mean(pcmfit{nC,nP}{jj}.G_predcv{257},3)); ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
%  caxis([min(min(mean(pcmfit{nC,nP}{jj}.G_predcv{257},3))) max(max(mean(pcmfit{nC,nP}{jj}.G_predcv{257},3)))])
subplot(nfig,24,1+5*24);imagesc(mean(pcmfit{nC,nP}{jj}.G_predcv{249},3)); ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
%  caxis([min(min(mean(pcmfit{nC,nP}{jj}.G_predcv{249},3))) max(max(mean(pcmfit{nC,nP}{jj}.G_predcv{249},3)))])

for ss=1:22,subplot(nfig,24,ss+2);imagesc(G_empCV5(:,:,jj,ss));  if ss==1,ylabel(['GempCV']);end;ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};title(['Sub ' num2str(ss)]); 
%   caxis([min(min(min(mean(G_empCV5,4)))) 1*max(max(max(mean(G_empCV5,4))))])
end
for ss=1:22,subplot(nfig,24,ss+1*24+2);imagesc(pcmfit{nC,nP}{jj}.G_predcv{1}(:,:,ss)); if ss==1,ylabel(['GFreeCV']);end;ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
%  caxis([min(min(mean(pcmfit{nC,nP}{jj}.G_predcv{1},3))) 1*max(mxg1)])
end
for ss=1:22,subplot(nfig,24,ss+2*24+2);imagesc(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindG(jj))}(:,:,ss)); if ss==1,ylabel(['GfitCVsame']);end;ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
%  caxis([0+min(min(mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindG(jj))},3))) 1*max(max(mean(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindG(jj))},3)))])
end
for ss=1:22,subplot(nfig,24,ss+3*24+2);imagesc(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindGeach(jj,ss))}(:,:,ss)); title(strrep(M{mtest(mxindGeach(jj,ss))}.name(7:end), '_', '')); if ss==1,ylabel(['GfitCVeach']);end;ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
%  caxis([0+min(min(mean(Gpredcvplot,3))) 1*max(max(mean(Gpredcvplot,3)))])
end
for ss=1:22,subplot(nfig,24,ss+4*24+2);imagesc(pcmfit{nC,nP}{jj}.G_predcv{257}(:,:,ss)); if ss==1,ylabel(['GfitCV257']);end;ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
%  caxis([min(min(mean(pcmfit{nC,nP}{jj}.G_predcv{257},3))) max(max(mean(pcmfit{nC,nP}{jj}.G_predcv{257},3)))])
end
for ss=1:22,subplot(nfig,24,ss+5*24+2);imagesc(pcmfit{nC,nP}{jj}.G_predcv{249}(:,:,ss)); if ss==1,ylabel(['GfitCV249']);end;ax=gca;ax.XTickLabel={'' '' '' '' '' '' ''};ax.YTickLabel={'' '' '' '' '' '' ''};
%  caxis([min(min(mean(pcmfit{nC,nP}{jj}.G_predcv{249},3))) max(max(mean(pcmfit{nC,nP}{jj}.G_predcv{249},3)))])
end

%%
% Sim7(NR) only
nfig=6;
set(0,'DefaultFigurePosition',[15 150 1300 350+100*(nfig-3)])
jj=find(timesteps==375);figure(jj);
for ss=1:22,subplot(nfig,22,ss);imagesc(G_empCV5(:,:,jj,ss));  if ss==1,ylabel(['GempCV']);end;end
for ss=1:22,subplot(nfig,22,ss+22);imagesc(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindG(jj))}(:,:,ss)); if ss==1,ylabel(['GfitCVsame']);end;end
for ss=1:22,subplot(nfig,22,ss+44);imagesc(pcmfit{nC,nP}{jj}.G_predcv{mtest(mxindGeach(jj,ss))}(:,:,ss)); if ss==1,ylabel(['GfitCVeach']);end;end
for ss=1:22,subplot(nfig,22,ss+66);imagesc(pcmfit{nC,nP}{jj}.G_predcv{257}(:,:,ss)); if ss==1,ylabel(['GfitCV257']);end;end
for ss=1:22,subplot(nfig,22,ss+88);imagesc(pcmfit{nC,nP}{jj}.G_predcv{249}(:,:,ss)); if ss==1,ylabel(['GfitCV249']);end;end
for ss=1:22,subplot(nfig,22,ss+110);imagesc(pcmfit{nC,nP}{jj}.G_predcv{1}(:,:,ss)); if ss==1,ylabel(['GFreeCV']);end;end


%%  Testing theta fit and likelihood output
% jj==3
% for ii=1:length(datacat_chunk_nn)
%       Y{ii}=datacat_chunk_nn{ii}(:,:,jj); % keep mean in
% end
% clear Mtheta
% Mtheta{1}=M{mtest(mxindGeach(jj,ss))};
% [A1{jj}.T,A1{jj}.theta_hat,A1{jj}.G_pred]=pcm_fitModelGroupCrossval(Y,Mtheta,partVecch,Xchunk,'runEffect','fixed');


%% for sim4eNR

addpath('C:\toolbox\spm12')
clear family* model*
likeuse='likCVind';
fmtest=[3 5 6 9 28 31 32 34];
for jj=1:numtimesteps
  switch likeuse
    case 'likCVind'
      lme=pcmfit{nC,nP}{jj}.Tcross.likelihood(:,fmtest);
  end
  familyT{jj}.infer='RFX';
  familyT{jj}.partition=[1 1 1 1 2 2 2 2];
  familyT{jj}.names={'without TIW','with TIW'};
  [familyT{jj},modelT{jj}] = spm_compare_families (lme,familyT{jj});
end

for jj=1:numtimesteps
  switch likeuse
    case 'likCVind'
      lme=pcmfit{nC,nP}{jj}.Tcross.likelihood(:,fmtest);
  end
  familyA{jj}.infer='RFX';
  familyA{jj}.partition=[1 2 1 2 1 2 1 2];
  familyA{jj}.names={'without Asym','with Asym'};
  [familyA{jj},modelA{jj}] = spm_compare_families (lme,familyA{jj});
end

for jj=1:numtimesteps
  switch likeuse
    case 'likCVind'
      lme=pcmfit{nC,nP}{jj}.Tcross.likelihood(:,fmtest);
  end
  familyS{jj}.infer='RFX';
  familyS{jj}.partition=[1 1 2 2 1 1 2 2];
  familyS{jj}.names={'without SymPair','with SymPair'};
  [familyS{jj},modelS{jj}] = spm_compare_families (lme,familyS{jj});
end

figind=34;
ymin=0;
yline=0.95;
ybot=0.05;
try close(figind);end
for jj=1:numtimesteps,xpt(jj)=familyT{jj}.xp(2);end;
for jj=1:numtimesteps,xpa(jj)=familyA{jj}.xp(2);end;
for jj=1:numtimesteps,xps(jj)=familyS{jj}.xp(2);end;
figure(figind);subplot(1,3,1);bar(xpt);ylim([ymin 1]);hold on;plot(0:18,yline*ones(1,19),'LineWidth',2);hold on;plot(0:18,ybot*ones(1,19),'g','LineWidth',2);title('Exc-Prob for TIW')
figure(figind);subplot(1,3,2);bar(xpa);ylim([ymin 1]);hold on;plot(0:18,yline*ones(1,19),'LineWidth',2);hold on;plot(0:18,ybot*ones(1,19),'g','LineWidth',2);title('Exc-Prob for Asym')
figure(figind);subplot(1,3,3);bar(xps);ylim([ymin 1]);hold on;plot(0:18,yline*ones(1,19),'LineWidth',2);hold on;plot(0:18,ybot*ones(1,19),'g','LineWidth',2);title('Exc-Prob for SymPair')
a=get(figind,'Children');
timechar=num2str([50:25:450]');
for pp=1:numel(a)
  a(pp).XTick=[3 7 11 15];%0:18;
  a(pp).XTickLabel={timechar(3,:) timechar(7,:) timechar(11,:) timechar(15,:) };
%   a(pp).XTickLabel={'' '' '' timechar(3,:) '' '' '' timechar(7,:) '' '' '' timechar(11,:) '' '' '' timechar(15,:) '' '' '' };
  a(pp).FontSize=18;
end

print(figind,[fdir 'pcm_familyXP_sim4eNR.eps'],'-depsc')
print(figind,[fdir 'pcm_familyXP_sim4eNR.png'],'-dpng')





%% for sim4b and sim4bNR and sim4bNRrand groupfit

addpath('C:\toolbox\spm12')
clear family* model*
likeuse='likCVind';
fmtest=mtest;
for jj=1:numtimesteps
  switch likeuse
    case 'aicind'
      lme=squeeze(-aicg(fmtest,jj,:))';
    case 'likind'
      lme=pcmfit{nC,nP}{jj}.Tcross.likelihood(:,fmtest);
    case 'likCVind'
      lme=pcmfit{nC,nP}{jj}.Tcross.likelihood(:,fmtest);
  end
  familyT{jj}.infer='RFX';
  familyT{jj}.partition=[1 2 1 1 2 2 1 2 1 2 1 1 2 2 1 2];
  familyT{jj}.names={'without TIW','with TIW'};
  [familyT{jj},modelT{jj}] = spm_compare_families (lme,familyT{jj});
end

for jj=1:numtimesteps
  switch likeuse
    case 'aicind'
      lme=squeeze(-aicg(fmtest,jj,:))';
    case 'likind'
      lme=pcmfit{nC,nP}{jj}.Tcross.likelihood(:,fmtest);
    case 'likCVind'
      lme=pcmfit{nC,nP}{jj}.Tcross.likelihood(:,fmtest);
  end
  familyA{jj}.infer='RFX';
  familyA{jj}.partition=[1 1 2 1 2 1 2 2 1 1 2 1 2 1 2 2];
  familyA{jj}.names={'without Asym','with Asym'};
  [familyA{jj},modelA{jj}] = spm_compare_families (lme,familyA{jj});
end


for jj=1:numtimesteps
  switch likeuse
    case 'aicind'
      lme=squeeze(-aicg(fmtest,jj,:))';
    case 'likind'
      lme=pcmfit{nC,nP}{jj}.Tcross.likelihood(:,fmtest);
    case 'likCVind'
      lme=pcmfit{nC,nP}{jj}.Tcross.likelihood(:,fmtest);
  end
  familyS{jj}.infer='RFX';
  familyS{jj}.partition=[1 1 1 2 1 2 2 2 1 1 1 2 1 2 2 2];
  familyS{jj}.names={'without SymPair','with SymPair'};
  [familyS{jj},modelS{jj}] = spm_compare_families (lme,familyS{jj});
end

figind=33;
ymin=0;
try close(figind);end
for jj=1:numtimesteps,xpt(jj)=familyT{jj}.xp(2);end;
for jj=1:numtimesteps,xpa(jj)=familyA{jj}.xp(2);end;
for jj=1:numtimesteps,xps(jj)=familyS{jj}.xp(2);end;
figure(figind);subplot(1,3,1);bar(xpt);ylim([ymin 1]);hold on;plot(0:18,0.75*ones(1,19));title('Exc-Prob for TIW')
figure(figind);subplot(1,3,2);bar(xpa);ylim([ymin 1]);hold on;plot(0:18,0.75*ones(1,19));title('Exc-Prob for Asym')
figure(figind);subplot(1,3,3);bar(xps);ylim([ymin 1]);hold on;plot(0:18,0.75*ones(1,19));title('Exc-Prob for SymPair')
a=get(figind,'Children');
timechar=num2str([50:25:450]');
for pp=1:numel(a)
  a(pp).XTick=0:18;
  a(pp).XTickLabel={'' '' timechar(2,:) '' '' timechar(5,:) '' '' timechar(8,:) '' '' timechar(11,:) '' '' timechar(14,:) '' '' timechar(17,:) '' };
end

%% sim7NR: specifc likelihoods

jj=7;
for mm=mtest
  for ff=1:8
      partition(mm-2,ff)=1+any(M{mm}.featincl==ff);
  end
end

lplot=mtest(find(partition(:,8)==2 & partition(:,2)==1 & partition(:,4)==1 & partition(:,6)==1));

set(0,'DefaultFigurePosition',[15 150 500 400])
figure(300+jj);  
Tcrosstmp.likelihood=pcmfit{nC,nP}{jj}.Tcross.likelihood(:,[1 lplot]);
Tcrosstmp.SN=pcmfit{nC,nP}{jj}.Tcross.SN;
plotlike{jj}=pcm_plotModelLikelihood(Tcrosstmp,{M{[1 lplot]}},'normalize',0,'Nceil',1,'Nnull',2,'upperceil',pcmfit{nC,nP}{jj}.Tgroup.likelihood(:,1));
  ylim([min(mean(plotlike{jj}.likelihood_norm(:,2:end),1))-3 max(mean(plotlike{jj}.likelihood_norm(:,2:end),1))+10])
  title([num2str(timesteps(jj)) ' ms'])
  ylabel('Relative log-likelihood')
  [mx,mind]=max(mean(plotlike{jj}.likelihood_norm(:,2:end),1));
  hold on;plot(0:size(plotlike{jj}.likelihood_norm,2),(mean(plotlike{jj}.likelihood_norm(:,mind+1),1)-3)*ones(1,size(plotlike{jj}.likelihood_norm,2)+1),'LineWidth',2);
%   ax=get(127,'Children');set(ax,'FontSize',16)


%% Family model comparison for sim7 group (all features)

clear family* model*
likeuse='likCVind';
addpath('C:\toolbox\spm12')
fmtest=3:257;

for jj=1:numtimesteps
  lme=pcmfit{nC,nP}{jj}.Tcross.likelihood(:,fmtest);
  for ff=1:8
    family{ff}{jj}.infer='RFX';
    family{ff}{jj}.partition=nan(1,length(fmtest));
    for mm=fmtest
      family{ff}{jj}.partition(mm-2)=1+any(M{mm}.featincl==ff);
    end
    family{ff}{jj}.names={['without ' num2str(ff)],['with ' num2str(ff)]};
    [family{ff}{jj},model{ff}{jj}] = spm_compare_families (lme,family{ff}{jj});
  end
end


for jj=1:numtimesteps,for ff=1:8,xpt(ff,jj)=family{ff}{jj}.xp(2);end;end;

figind=10;
ymin=0;
for ff=1:8
  figure(figind);subplot(4,2,ff);bar(xpt(ff,:));ylim([ymin 1]);hold on;plot(0:18,0.75*ones(1,19));title(['Exc-Prob for ' num2str(ff)])
end
a=get(figind,'Children');
timechar=num2str([50:25:450]');
for pp=1:numel(a)
  a(pp).XTick=0:18;
  a(pp).XTickLabel={'' '' timechar(2,:) '' '' timechar(5,:) '' '' timechar(8,:) '' '' timechar(11,:) '' '' timechar(14,:) '' '' timechar(17,:) '' };
end

%% Family model comparison for sim7 group: excluding feature8 (6x6)

clear family* model*
likeuse='likCVind';
fmtest=mtest;  % excluding 8
addpath('C:\toolbox\spm12')

for jj=1:numtimesteps
  lme=pcmfit{nC,nP}{jj}.Tcross.likelihood(:,fmtest);
  for ff=1:7
    family{ff}{jj}.infer='RFX';
    family{ff}{jj}.partition=nan(1,length(fmtest));
    mmind=0;
    for mm=fmtest
      mmind=mmind+1;
      family{ff}{jj}.partition(mmind)=1+any(M{mm}.featincl==ff);
    end
    family{ff}{jj}.names={['without ' num2str(ff)],['with ' num2str(ff)]};
    [family{ff}{jj},model{ff}{jj}] = spm_compare_families (lme,family{ff}{jj});
  end
end


for jj=1:numtimesteps,for ff=1:7,xpt(ff,jj)=family{ff}{jj}.xp(2);end;end;

figind=10;
ymin=0;
for ff=1:7
  figure(figind);subplot(4,2,ff);bar(xpt(ff,:));ylim([ymin 1]);hold on;plot(0:18,0.75*ones(1,19));title(['Exc-Prob for ' num2str(ff)])
end
a=get(figind,'Children');
timechar=num2str([50:25:450]');
for pp=1:numel(a)
  a(pp).XTick=0:18;
  a(pp).XTickLabel={'' '' timechar(2,:) '' '' timechar(5,:) '' '' timechar(8,:) '' '' timechar(11,:) '' '' timechar(14,:) '' '' timechar(17,:) '' };
end

%% Family model comparison for sim7 group: excluding feature8 (6x6); consolidation of features

clear family* model*
fmtest=mtest; 
addpath('C:\toolbox\spm12')

for jj=1:numtimesteps
  lme=pcmfit{nC,nP}{jj}.Tcross.likelihood(:,fmtest);

  % (Cross0) 1+2, 3+4+5, 6+7
  for ff=1:3,
    familyCross0{ff}{jj}.infer='RFX';
    familyCross0{ff}{jj}.partition=nan(1,length(fmtest));
    mmind=0;
    for mm=fmtest
      mmind=mmind+1;
      if ff==1
        familyCross0{ff}{jj}.partition(mmind)=1+any(M{mm}.featincl==1 | M{mm}.featincl==2);
        ffstr{ff}='1 2';
      elseif ff==2
        familyCross0{ff}{jj}.partition(mmind)=1+any(M{mm}.featincl==3 | M{mm}.featincl==4 | M{mm}.featincl==5);
        ffstr{ff}='3 4 5';
      elseif ff==3
        familyCross0{ff}{jj}.partition(mmind)=1+any(M{mm}.featincl==6 | M{mm}.featincl==7);
        ffstr{ff}='6 7';
      end
    end
    familyCross0{ff}{jj}.names={['without ' ffstr{ff}],['with ' ffstr{ff}]};
    [familyCross0{ff}{jj}] = spm_compare_families (lme,familyCross0{ff}{jj});
  end
  
  
  % (Edges) 1+2, 3+5+7, 4+6
  for ff=1:3,
    familyEdges{ff}{jj}.infer='RFX';
    familyEdges{ff}{jj}.partition=nan(1,length(fmtest));
    mmind=0;
    for mm=fmtest
      mmind=mmind+1;
      if ff==1
        familyEdges{ff}{jj}.partition(mmind)=1+any(M{mm}.featincl==1 | M{mm}.featincl==2);
        ffstr{ff}='1 2';
      elseif ff==2
        familyEdges{ff}{jj}.partition(mmind)=1+any(M{mm}.featincl==3 | M{mm}.featincl==5 | M{mm}.featincl==7);
        ffstr{ff}='3 5 7';
      elseif ff==3
        familyEdges{ff}{jj}.partition(mmind)=1+any(M{mm}.featincl==4 | M{mm}.featincl==6);
        ffstr{ff}='4 6';
      end
    end
    familyEdges{ff}{jj}.names={['without ' ffstr{ff}],['with ' ffstr{ff}]};
    [familyEdges{ff}{jj}] = spm_compare_families (lme,familyEdges{ff}{jj});
  end
  
  % (Size) 1+3+4, 5+6, 2+7
  for ff=1:3,
    familySize{ff}{jj}.infer='RFX';
    familySize{ff}{jj}.partition=nan(1,length(fmtest));
    mmind=0;
    for mm=fmtest
      mmind=mmind+1;
      if ff==1
        familySize{ff}{jj}.partition(mmind)=1+any(M{mm}.featincl==1 | M{mm}.featincl==3 | M{mm}.featincl==4);
        ffstr{ff}='1 3 4';
      elseif ff==2
        familySize{ff}{jj}.partition(mmind)=1+any(M{mm}.featincl==5 | M{mm}.featincl==6);
        ffstr{ff}='5 6';
      elseif ff==3
        familySize{ff}{jj}.partition(mmind)=1+any(M{mm}.featincl==2 | M{mm}.featincl==7);
        ffstr{ff}='2 7';
      end
    end
    familySize{ff}{jj}.names={['without ' ffstr{ff}],['with ' ffstr{ff}]};
    [familySize{ff}{jj}] = spm_compare_families (lme,familySize{ff}{jj});
  end
  
  % (TIWAsym) 1+2, 3-7
  for ff=1:2,
    familyTIWAsym{ff}{jj}.infer='RFX';
    familyTIWAsym{ff}{jj}.partition=nan(1,length(fmtest));
    mmind=0;
    for mm=fmtest
      mmind=mmind+1;
      if ff==1
        familyTIWAsym{ff}{jj}.partition(mmind)=1+any(M{mm}.featincl==1 | M{mm}.featincl==2);
        ffstr{ff}='1 2';
      elseif ff==2
        familyTIWAsym{ff}{jj}.partition(mmind)=1+any(M{mm}.featincl==3 | M{mm}.featincl==4 | M{mm}.featincl==5 | M{mm}.featincl==6 | M{mm}.featincl==7);
        ffstr{ff}='3 4 5 6 7';
      end
    end
    familyTIWAsym{ff}{jj}.names={['without ' ffstr{ff}],['with ' ffstr{ff}]};
    [familyTIWAsym{ff}{jj}] = spm_compare_families (lme,familyTIWAsym{ff}{jj});
  end
  
end
  
for jj=1:numtimesteps,
  for ff=1:3,xpC(ff,jj)=familyCross0{ff}{jj}.xp(2);end;
  for ff=1:3,xpE(ff,jj)=familyEdges{ff}{jj}.xp(2);end;
  for ff=1:3,xpS(ff,jj)=familySize{ff}{jj}.xp(2);end;
  for ff=1:2,xpT(ff,jj)=familyTIWAsym{ff}{jj}.xp(2);end;
end;

ymin=0;
figind=10; for ff=1:3
  figure(figind);subplot(1,3,ff);bar(xpC(ff,:));ylim([ymin 1]);hold on;plot(0:18,0.75*ones(1,19));title(['Exc-Prob ' familyCross0{ff}{jj}.names{2}])
end
figind=11; for ff=1:3
  figure(figind);subplot(1,3,ff);bar(xpE(ff,:));ylim([ymin 1]);hold on;plot(0:18,0.75*ones(1,19));title(['Exc-Prob ' familyEdges{ff}{jj}.names{2}])
end
figind=12; for ff=1:3
  figure(figind);subplot(1,3,ff);bar(xpS(ff,:));ylim([ymin 1]);hold on;plot(0:18,0.75*ones(1,19));title(['Exc-Prob ' familySize{ff}{jj}.names{2}])
end
figind=13; for ff=1:2
  figure(figind);subplot(1,3,ff);bar(xpT(ff,:));ylim([ymin 1]);hold on;plot(0:18,0.75*ones(1,19));title(['Exc-Prob ' familyTIWAsym{ff}{jj}.names{2}])
end

for figind=10:13
  a=get(figind,'Children');
  timechar=num2str([50:25:450]');
  for pp=1:numel(a)
    a(pp).XTick=0:18;
    a(pp).XTickLabel={'' '' timechar(2,:) '' '' timechar(5,:) '' '' timechar(8,:) '' '' timechar(11,:) '' '' timechar(14,:) '' '' timechar(17,:) '' };
  end
end



%%  For sim4

% % --->  spm_compare_families

clear family* model*
likeuse='aicind';
fmtest=3:32;
for jj=1:numtimesteps
  switch likeuse
    case 'aicind'
      lme=squeeze(-aici(3:32,jj,:))';
    case 'likind'
      lme=pcmind{nC,nP}{jj}.Tgroup.likelihood(:,3:32);
    case 'likCVind'
      lme=pcmind{nC,nP}{jj}.D.likelihood(:,3:32);
  end
  familyC{jj}.infer='RFX';
  familyC{jj}.partition=[ones(1,15) 2*ones(1,15)];
  familyC{jj}.names={'without allC','with allC'};
  [familyC{jj},modelC{jj}] = spm_compare_families (lme,familyC{jj})
end

for jj=1:numtimesteps  % super strong support for 'with eye'
  switch likeuse
    case 'aicind'
      lme=squeeze(-aici(3:32,jj,:))';
    case 'likind'
      lme=pcmind{nC,nP}{jj}.Tgroup.likelihood(:,3:32);
    case 'likCVind'
      lme=pcmind{nC,nP}{jj}.D.likelihood(:,3:32);
  end
  familyI{jj}.infer='RFX';
  familyI{jj}.partition=[2 ones(1,7) 2*ones(1,7) 2 ones(1,7) 2*ones(1,7)];
  familyI{jj}.names={'without eye','with eye'};
  [familyI{jj},modelI{jj}] = spm_compare_families (lme,familyI{jj})
end
for jj=1:numtimesteps
  switch likeuse
    case 'aicind'
      lme=squeeze(-aici(fmtest(familyI{jj}.partition==2),jj,:))';
    case 'likind'
      lme=pcmind{nC,nP}{jj}.Tgroup.likelihood(:,fmtest(familyI{jj}.partition==2));
    case 'likCVind'
      lme=pcmind{nC,nP}{jj}.D.likelihood(:,fmtest(familyI{jj}.partition==2));
  end
  familyC_I{jj}.infer='RFX';
  familyC_I{jj}.partition=[ones(1,8) 2*ones(1,8)];
  familyC_I{jj}.names={'without allC','with allC'};
  [familyC_I{jj},modelC_I{jj}] = spm_compare_families (lme,familyC_I{jj})
end

for jj=1:numtimesteps
  switch likeuse
    case 'aicind'
      lme=squeeze(-aici(fmtest(familyI{jj}.partition==2),jj,:))';
    case 'likind'
      lme=pcmind{nC,nP}{jj}.Tgroup.likelihood(:,fmtest(familyI{jj}.partition==2));
    case 'likCVind'
      lme=pcmind{nC,nP}{jj}.D.likelihood(:,fmtest(familyI{jj}.partition==2));
  end
  familyT_I{jj}.infer='RFX';
  familyT_I{jj}.partition=[1 2 1 1 2 2 1 2  1 2 1 1 2 2 1 2];
  familyT_I{jj}.names={'without TIW','with TIW'};
  [familyT_I{jj},modelT_I{jj}] = spm_compare_families (lme,familyT_I{jj})
end

for jj=1:numtimesteps
  switch likeuse
    case 'aicind'
      lme=squeeze(-aici(fmtest(familyI{jj}.partition==2),jj,:))';
    case 'likind'
      lme=pcmind{nC,nP}{jj}.Tgroup.likelihood(:,fmtest(familyI{jj}.partition==2));
    case 'likCVind'
      lme=pcmind{nC,nP}{jj}.D.likelihood(:,fmtest(familyI{jj}.partition==2));
  end
  familyA_I{jj}.infer='RFX';
  familyA_I{jj}.partition=[1 1 2 1 2 1 2 2 1 1 2 1 2 1 2 2];
  familyA_I{jj}.names={'without Asym','with Asym'};
  [familyA_I{jj},modelA_I{jj}] = spm_compare_families (lme,familyA_I{jj})
end
for jj=1:numtimesteps
  switch likeuse
    case 'aicind'
      lme=squeeze(-aici(fmtest(familyI{jj}.partition==2),jj,:))';
    case 'likind'
      lme=pcmind{nC,nP}{jj}.Tgroup.likelihood(:,fmtest(familyI{jj}.partition==2));
    case 'likCVind'
      lme=pcmind{nC,nP}{jj}.D.likelihood(:,fmtest(familyI{jj}.partition==2));
  end
  familyP_I{jj}.infer='RFX';
  familyP_I{jj}.partition=[1 1 1 2 1 2 2 2 1 1 1 2 1 2 2 2];
  familyP_I{jj}.names={'without SPair','with SPair'};
  [familyP_I{jj},modelP_I{jj}] = spm_compare_families (lme,familyP_I{jj})
end

fmatest=fmtest(familyI{jj}.partition==2);
for jj=1:numtimesteps
  switch likeuse
    case 'aicind'
      lme=squeeze(-aici(fmatest(familyA_I{jj}.partition==2),jj,:))';
    case 'likind'
      lme=pcmind{nC,nP}{jj}.Tgroup.likelihood(:,fmatest(familyA_I{jj}.partition==2));
    case 'likCVind'
      lme=pcmind{nC,nP}{jj}.D.likelihood(:,fmatest(familyA_I{jj}.partition==2));
  end
  familyT_AI{jj}.infer='RFX';
  familyT_AI{jj}.partition=[1 2 1 2 1 2 1 2];
  familyT_AI{jj}.names={'without TIW','with TIW'};
  [familyT_AI{jj},modelT_AI{jj}] = spm_compare_families (lme,familyT_AI{jj});
end

for jj=1:numtimesteps
  switch likeuse
    case 'aicind'
      lme=squeeze(-aici(fmatest(familyA_I{jj}.partition==2),jj,:))';
    case 'likind'
      lme=pcmind{nC,nP}{jj}.Tgroup.likelihood(:,fmatest(familyA_I{jj}.partition==2));
    case 'likCVind'
      lme=pcmind{nC,nP}{jj}.D.likelihood(:,fmatest(familyA_I{jj}.partition==2));
  end
  familyP_AI{jj}.infer='RFX';
  familyP_AI{jj}.partition=[1 1 2 2 1 1 2 2];
  familyP_AI{jj}.names={'without SPair','with SPair'};
  [familyP_AI{jj},modelP_AI{jj}] = spm_compare_families (lme,familyP_AI{jj});
end

save([pcdir 'familycompare_' likeuse],'family*')

likeuse='likCVind';
load([pcdir 'familycompare_' likeuse],'family*')
ymin=0;  % 0.4
switch likeuse
  case 'aicind'
    figind=12;
  case 'likind'
    figind=10;
  case 'likCVind'
    figind=11;
end
try close(figind);end

for jj=1:numtimesteps,xpi(jj)=familyI{jj}.xp(2);end;
figure(figind);subplot(2,4,1);bar(xpi);ylim([ymin 1]);hold on;plot(0:18,0.75*ones(1,19));title('Exc-Prob for Identity')
for jj=1:numtimesteps,xpci(jj)=familyC_I{jj}.xp(2);end;
figure(figind);subplot(2,4,5);bar(xpci);ylim([ymin 1]);hold on;plot(0:18,0.75*ones(1,19));title('Exc-Prob for OnesCol | I')
for jj=1:numtimesteps,xpti(jj)=familyT_I{jj}.xp(2);end;
figure(figind);subplot(2,4,2);bar(xpti);ylim([ymin 1]);hold on;plot(0:18,0.75*ones(1,19));title('Exc-Prob for TIW | I')
for jj=1:numtimesteps,xpai(jj)=familyA_I{jj}.xp(2);end;
figure(figind);subplot(2,4,4);bar(xpai);ylim([ymin 1]);hold on;plot(0:18,0.75*ones(1,19));title('Exc-Prob for Asymmetry | I')
for jj=1:numtimesteps,xppi(jj)=familyP_I{jj}.xp(2);end;
figure(figind);subplot(2,4,3);bar(xppi);ylim([ymin 1]);hold on;plot(0:18,0.75*ones(1,19));title('Exc-Prob for SymPairs | I')
for jj=1:numtimesteps,xptai(jj)=familyT_AI{jj}.xp(2);end;
figure(figind);subplot(2,4,6);bar(xptai);ylim([ymin 1]);hold on;plot(0:18,0.75*ones(1,19));title('Exc-Prob for TIW | I & A')
for jj=1:numtimesteps,xppai(jj)=familyP_AI{jj}.xp(2);end;
figure(figind);subplot(2,4,7);bar(xppai);ylim([ymin 1]);hold on;plot(0:18,0.75*ones(1,19));title('Exc-Prob for SymPairs | I & A')
a=get(figind,'Children');
timechar=num2str([50:25:450]');
for pp=1:numel(a)
  a(pp).XTick=0:18;
  a(pp).XTickLabel={'' '' timechar(2,:) '' '' timechar(5,:) '' '' timechar(8,:) '' '' timechar(11,:) '' '' timechar(14,:) '' '' timechar(17,:) '' };
end


%% Family model comparison for sim5

clear family* model*
likeuse='likCVind';
fmtest=18:31;

for jj=1:numtimesteps
  lme=pcmind{nC,nP}{jj}.D.likelihood(:,23:31);
  family3{jj}.infer='RFX';
  family3{jj}.partition=[2 1 1 1 2 1 1 2 1];
  family3{jj}.names={'without all3','with all3'};
  [family3{jj},model3{jj}] = spm_compare_families (lme,family3{jj});
end

for jj=1:numtimesteps
  lme=pcmind{nC,nP}{jj}.D.likelihood(:,23:31);
  family4{jj}.infer='RFX';
  family4{jj}.partition=[1 2 1 1 2 2 1 2 2];
  family4{jj}.names={'without all4','with all4'};
  [family4{jj},model4{jj}] = spm_compare_families (lme,family4{jj});
end

for jj=1:numtimesteps
  lme=pcmind{nC,nP}{jj}.D.likelihood(:,23:31);
  family5{jj}.infer='RFX';
  family5{jj}.partition=[1 1 2 1 1 2 2 2 2];
  family5{jj}.names={'without all5','with all5'};
  [family5{jj},model5{jj}] = spm_compare_families (lme,family5{jj});
end

for jj=1:numtimesteps
  lme=pcmind{nC,nP}{jj}.D.likelihood(:,23:31);
  family6{jj}.infer='RFX';
  family6{jj}.partition=[1 1 1 2 1 1 2 1 2];
  family6{jj}.names={'without all6','with all6'};
  [family6{jj},model6{jj}] = spm_compare_families (lme,family6{jj});
end


for jj=1:numtimesteps,xp3(jj)=family3{jj}.xp(2);end;
for jj=1:numtimesteps,xp4(jj)=family4{jj}.xp(2);end;
for jj=1:numtimesteps,xp5(jj)=family5{jj}.xp(2);end;
for jj=1:numtimesteps,xp6(jj)=family6{jj}.xp(2);end;

figind=10;
ymin=0;
figure(figind);subplot(4,1,1);bar(xp3);ylim([ymin 1]);hold on;plot(0:18,0.75*ones(1,19));title('Exc-Prob for 3x3')
figure(figind);subplot(4,1,2);bar(xp4);ylim([ymin 1]);hold on;plot(0:18,0.75*ones(1,19));title('Exc-Prob for 4x4')
figure(figind);subplot(4,1,3);bar(xp5);ylim([ymin 1]);hold on;plot(0:18,0.75*ones(1,19));title('Exc-Prob for 5x5')
figure(figind);subplot(4,1,4);bar(xp6);ylim([ymin 1]);hold on;plot(0:18,0.75*ones(1,19));title('Exc-Prob for 6x6')
a=get(figind,'Children');
timechar=num2str([50:25:450]');
for pp=1:numel(a)
  a(pp).XTick=0:18;
  a(pp).XTickLabel={'' '' timechar(2,:) '' '' timechar(5,:) '' '' timechar(8,:) '' '' timechar(11,:) '' '' timechar(14,:) '' '' timechar(17,:) '' };
end


%% Family model comparison for sim6

clear family* model*
likeuse='likCVind';
fmtest=4:18;

for jj=1:numtimesteps
  lme=pcmind{nC,nP}{jj}.D.likelihood(:,fmtest);
  familyT{jj}.infer='RFX';
  familyT{jj}.partition=[2 2 2 2 1 1 1 1 1 1 2 2 2 1 2];
  familyT{jj}.names={'without onTIW','with onTIW'};
  [familyT{jj},modelT{jj}] = spm_compare_families (lme,familyT{jj});
end

for jj=1:numtimesteps
  lme=pcmind{nC,nP}{jj}.D.likelihood(:,fmtest);
  familyA{jj}.infer='RFX';
  familyA{jj}.partition=[1 2 2 2 2 2 2 1 1 1 1 1 1 2 2];
  familyA{jj}.names={'without Asym','with Asym'};
  [familyA{jj},modelA{jj}] = spm_compare_families (lme,familyA{jj});
end

for jj=1:numtimesteps
  lme=pcmind{nC,nP}{jj}.D.likelihood(:,fmtest);
  familyG{jj}.infer='RFX';
  familyG{jj}.partition=[1 1 2 2 1 2 2 2 2 1 2 1 2 1 1];
  familyG{jj}.names={'without Green','with Green'};
  [familyG{jj},modelG{jj}] = spm_compare_families (lme,familyG{jj});
end

for jj=1:numtimesteps
  lme=pcmind{nC,nP}{jj}.D.likelihood(:,fmtest);
  familyC{jj}.infer='RFX';
  familyC{jj}.partition=[1 1 1 2 1 1 2 1 2 2 1 2 2 2 2];
  familyC{jj}.names={'without Cyan','with Cyan'};
  [familyC{jj},modelC{jj}] = spm_compare_families (lme,familyC{jj});
end

for jj=1:numtimesteps,xpt(jj)=familyT{jj}.xp(2);end;
for jj=1:numtimesteps,xpa(jj)=familyA{jj}.xp(2);end;
for jj=1:numtimesteps,xpg(jj)=familyG{jj}.xp(2);end;
for jj=1:numtimesteps,xpc(jj)=familyC{jj}.xp(2);end;

figind=10;
ymin=0;
figure(figind);subplot(4,1,1);bar(xpt);ylim([ymin 1]);hold on;plot(0:18,0.75*ones(1,19));title('Exc-Prob for TIW')
figure(figind);subplot(4,1,2);bar(xpa);ylim([ymin 1]);hold on;plot(0:18,0.75*ones(1,19));title('Exc-Prob for Asym')
figure(figind);subplot(4,1,3);bar(xpg);ylim([ymin 1]);hold on;plot(0:18,0.75*ones(1,19));title('Exc-Prob for Green')
figure(figind);subplot(4,1,4);bar(xpc);ylim([ymin 1]);hold on;plot(0:18,0.75*ones(1,19));title('Exc-Prob for Cyan')
a=get(figind,'Children');
timechar=num2str([50:25:450]');
for pp=1:numel(a)
  a(pp).XTick=0:18;
  a(pp).XTickLabel={'' '' timechar(2,:) '' '' timechar(5,:) '' '' timechar(8,:) '' '' timechar(11,:) '' '' timechar(14,:) '' '' timechar(17,:) '' };
end

for jj=1:numtimesteps,xpt(jj)=familyT{jj}.exp_r(2);end;
for jj=1:numtimesteps,xpa(jj)=familyA{jj}.exp_r(2);end;
for jj=1:numtimesteps,xpg(jj)=familyG{jj}.exp_r(2);end;
for jj=1:numtimesteps,xpc(jj)=familyC{jj}.exp_r(2);end;

figind=10;
ymin=0;
figure(figind);subplot(4,1,1);bar(xpt);ylim([ymin 1]);hold on;plot(0:18,0.75*ones(1,19));title('Exc-Prob for TIW')
figure(figind);subplot(4,1,2);bar(xpa);ylim([ymin 1]);hold on;plot(0:18,0.75*ones(1,19));title('Exc-Prob for Asym')
figure(figind);subplot(4,1,3);bar(xpg);ylim([ymin 1]);hold on;plot(0:18,0.75*ones(1,19));title('Exc-Prob for Green')
figure(figind);subplot(4,1,4);bar(xpc);ylim([ymin 1]);hold on;plot(0:18,0.75*ones(1,19));title('Exc-Prob for Cyan')
a=get(figind,'Children');
timechar=num2str([50:25:450]');
for pp=1:numel(a)
  a(pp).XTick=0:18;
  a(pp).XTickLabel={'' '' timechar(2,:) '' '' timechar(5,:) '' '' timechar(8,:) '' '' timechar(11,:) '' '' timechar(14,:) '' '' timechar(17,:) '' };
end

%% Family model comparison for sim7

clear family* model*
likeuse='likCVind';
fmtest=3:257;

for jj=1:numtimesteps
  lme=pcmind{nC,nP}{jj}.D.likelihood(:,fmtest);
  for ff=1:8
    family{ff}{jj}.infer='RFX';
    family{ff}{jj}.partition=nan(1,length(fmtest));
    for mm=fmtest
      family{ff}{jj}.partition(mm-2)=1+any(M{mm}.featincl==ff);
    end
    family{ff}{jj}.names={['without ' num2str(ff)],['with ' num2str(ff)]};
    [family{ff}{jj},model{ff}{jj}] = spm_compare_families (lme,family{ff}{jj});
  end
end


for jj=1:numtimesteps,for ff=1:8,xpt(ff,jj)=family{ff}{jj}.xp(2);end;end;

figind=10;
ymin=0;
for ff=1:8
  figure(figind);subplot(4,2,ff);bar(xpt(ff,:));ylim([ymin 1]);hold on;plot(0:18,0.75*ones(1,19));title(['Exc-Prob for ' num2str(ff)])
end
a=get(figind,'Children');
timechar=num2str([50:25:450]');
for pp=1:numel(a)
  a(pp).XTick=0:18;
  a(pp).XTickLabel={'' '' timechar(2,:) '' '' timechar(5,:) '' '' timechar(8,:) '' '' timechar(11,:) '' '' timechar(14,:) '' '' timechar(17,:) '' };
end



%%

% Ignore this old way of family comparison
% Plot AICi and AICg per subject
if 0
  m_eye=[3 11:17];
  m_tiw=[4 7 8 10 11 14 15 17];
  m_asm=[5 7 9 10 12 14 16 17];
  m_spr=[6 8 9 10 13 15 16 17];
  mt_yeye=dsearchn(mtest',m_eye');
  mt_neye=dsearchn(mtest',setdiff(3:17,m_eye)');
  mt_ytiw=dsearchn(mtest',m_tiw');
  mt_ntiw=dsearchn(mtest',setdiff(3:17,m_tiw)');
  mt_yasm=dsearchn(mtest',m_asm');
  mt_nasm=dsearchn(mtest',setdiff(3:17,m_asm)');
  mt_yspr=dsearchn(mtest',m_spr');
  mt_nspr=dsearchn(mtest',setdiff(3:17,m_spr)');
  mtest2=mtest(mt_nasm);
  for ss=1:size(aici,3)
    aiciss=aici(mtest,:,ss);
    aicissdif(:,:,ss)=aiciss-repmat(min(aiciss),[numel(mtest) 1]);
    figure(45);subplot(4,6,ss);imagesc(aicissdif(:,:,ss));caxis([0 15])
    aicgss=aicg(mtest,:,ss);
    aicgssdif(:,:,ss)=aicgss-repmat(min(aicgss),[numel(mtest) 1]);
    figure(46);subplot(4,6,ss);imagesc(aicgssdif(:,:,ss));caxis([0 15])
    figure(47);subplot(4,6,ss);imagesc(aicgssdif(:,:,ss)-aicissdif(:,:,ss));caxis([0 15])
    
    %   figure(50);subplot(4,6,ss);bar([mean(aicissdif(mt_yeye,:,ss),1); mean(aicissdif(mt_neye,:,ss),1)]');
    %   figure(51);subplot(4,6,ss);bar([mean(aicissdif(mt_ytiw,:,ss),1); mean(aicissdif(mt_ntiw,:,ss),1)]');
    %   figure(52);subplot(4,6,ss);bar([mean(aicissdif(mt_yasm,:,ss),1); mean(aicissdif(mt_nasm,:,ss),1)]');
    %   figure(53);subplot(4,6,ss);bar([mean(aicissdif(mt_yspr,:,ss),1); mean(aicissdif(mt_nspr,:,ss),1)]');
    %   figure(60);subplot(4,6,ss);bar([mean(aicgssdif(mt_yeye,:,ss),1); mean(aicgssdif(mt_neye,:,ss),1)]');
    %   figure(61);subplot(4,6,ss);bar([mean(aicgssdif(mt_ytiw,:,ss),1); mean(aicgssdif(mt_ntiw,:,ss),1)]');
    %   figure(62);subplot(4,6,ss);bar([mean(aicgssdif(mt_yasm,:,ss),1); mean(aicgssdif(mt_nasm,:,ss),1)]');
    %   figure(63);subplot(4,6,ss);bar([mean(aicgssdif(mt_yspr,:,ss),1); mean(aicgssdif(mt_nspr,:,ss),1)]');
  end
  figure(45);subplot(4,6,24);imagesc(mean(aicissdif,3));caxis([0 15])
  figure(46);subplot(4,6,24);imagesc(mean(aicgssdif,3));caxis([0 15])
  figure(47);subplot(4,6,24);imagesc(mean(aicgssdif,3)-mean(aicissdif,3));caxis([0 15])
  
  for jj=1:14,
    figure(90);subplot(4,14,00+jj);bar(squeeze(mean(aici(mtest(mt_yeye),jj,:),1)-mean(aici(mtest(mt_neye),jj,:),1))');xlim([0 23]);ylim([-5 5]);title([num2str(round(mean(squeeze(mean(aici(mtest(mt_yeye),jj,:),1)-mean(aici(mtest(mt_neye),jj,:),1))'),1))]);
    if jj==1,ylabel('Identity');end;
    figure(90);subplot(4,14,14+jj);bar(squeeze(mean(aici(mtest(mt_ytiw),jj,:),1)-mean(aici(mtest(mt_ntiw),jj,:),1))');xlim([0 23]);ylim([-5 5]);title([num2str(round(mean(squeeze(mean(aici(mtest(mt_ytiw),jj,:),1)-mean(aici(mtest(mt_ntiw),jj,:),1))'),1))]);
    if jj==1,ylabel('TIW');end;
    figure(90);subplot(4,14,28+jj);bar(squeeze(mean(aici(mtest(mt_yasm),jj,:),1)-mean(aici(mtest(mt_nasm),jj,:),1))');xlim([0 23]);ylim([-5 5]);title([num2str(round(mean(squeeze(mean(aici(mtest(mt_yasm),jj,:),1)-mean(aici(mtest(mt_nasm),jj,:),1))'),1))]);
    if jj==1,ylabel('Asym');end;
    figure(90);subplot(4,14,42+jj);bar(squeeze(mean(aici(mtest(mt_yspr),jj,:),1)-mean(aici(mtest(mt_nspr),jj,:),1))');xlim([0 23]);ylim([-5 5]);title([num2str(round(mean(squeeze(mean(aici(mtest(mt_yspr),jj,:),1)-mean(aici(mtest(mt_nspr),jj,:),1))'),1))]);
    if jj==1,ylabel('SymPair');end;
  end
  
  for jj=1:14,
    figure(70);subplot(4,4,jj);
    bar(squeeze([mean(aicissdif(mt_yeye,jj,:),1); mean(aicissdif(mt_neye,jj,:),1);])');xlim([0 23]);
    figure(71);subplot(4,4,jj);
    bar(squeeze([mean(aicissdif(mt_ytiw,jj,:),1); mean(aicissdif(mt_ntiw,jj,:),1);])');xlim([0 23]);
    figure(72);subplot(4,4,jj);
    bar(squeeze([mean(aicissdif(mt_yasm,jj,:),1); mean(aicissdif(mt_nasm,jj,:),1);])');xlim([0 23]);
    figure(73);subplot(4,4,jj);
    bar(squeeze([mean(aicissdif(mt_yspr,jj,:),1); mean(aicissdif(mt_nspr,jj,:),1);])');xlim([0 23]);
    figure(80);subplot(4,4,jj);
    bar(squeeze([mean(aicgssdif(mt_yeye,jj,:),1); mean(aicgssdif(mt_neye,jj,:),1);])');xlim([0 23]);
    figure(81);subplot(4,4,jj);
    bar(squeeze([mean(aicgssdif(mt_ytiw,jj,:),1); mean(aicgssdif(mt_ntiw,jj,:),1);])');xlim([0 23]);
    figure(82);subplot(4,4,jj);
    bar(squeeze([mean(aicgssdif(mt_yasm,jj,:),1); mean(aicgssdif(mt_nasm,jj,:),1);])');xlim([0 23]);
    figure(83);subplot(4,4,jj);
    bar(squeeze([mean(aicgssdif(mt_yspr,jj,:),1); mean(aicgssdif(mt_nspr,jj,:),1);])');xlim([0 23]);
  end
end

%%
% Winning models
timesteps=50:25:450;
nfig=4;  % Individual level fit6
set(0,'DefaultFigurePosition',[15 150 1300 300+100*(nfig-3)])
for jj=1:numtimesteps,
  mxi1(jj)=max(max(mean(pcmind{nC,nP}{jj}.G_predCV{1},3)));
  [mx,mxindI(jj)]=max(mean(pcmind{nC,nP}{jj}.D.likelihood(:,mtest)-repmat(pcmind{nC,nP}{jj}.D.likelihood(:,2),[1 numel(mtest)])));
  [mxsortI{jj},sortindI{jj}]=sort(mean(pcmind{nC,nP}{jj}.D.likelihood(:,mtest)-repmat(pcmind{nC,nP}{jj}.D.likelihood(:,2),[1 numel(mtest)])));
  mxi3(jj)=mxsortI{jj}(end);
  mxi2(jj)=max(max(mean(pcmind{nC,nP}{jj}.G_predCV{mtest(mxindI(jj))},3)));
end
for jj=1:numtimesteps,
  [mx,mxindIeach(jj,:)]=max(pcmind{nC,nP}{jj}.D.likelihood(:,mtest)-repmat(pcmind{nC,nP}{jj}.D.likelihood(:,2),[1 numel(mtest)]),[],2);
end


for jj=1:numtimesteps,
  figure(100);
  %   subplot(nfig,numtimesteps,jj);imagesc(mean(G_emp(:,:,jj,:),4)); % caxis([0 0.9*max(max(max(mean(G_emp,4))))])
  %   if jj==1,ylabel(['Gemp']);end
  subplot(nfig,numtimesteps,jj);imagesc(mean(G_empCV5(:,:,jj,:),4));   % caxis([0 0.9*max(max(max(mean(G_empCV5,4))))])
  if jj==1,ylabel(['GempCV']);end
  title([num2str(timesteps(jj)) ' ms'])
  %   subplot(nfig,numtimesteps,numtimesteps+jj);imagesc(H*mean(G_empCV5(:,:,jj,:),4)*H');   caxis([0 0.1*max(max(max(mean(G_empCV5,4))))])
  %   if jj==1,ylabel(['H*GeCV*H']);end
  %   subplot(nfig,numtimesteps,28+jj);imagesc(mean(pcminde5{jj}.G_predCV{mtest(mxindI(jj))},3));caxis([0 max(mxi2)])
  %  subplot(nfig,numtimesteps,28+jj);imagesc(0.8*mean(pcminde5{jj}.G_predCV{mtest(sortindI{jj}(end-1))},3)+mean(pcminde5{jj}.G_predCV{mtest(mxindI(jj))},3)); % caxis([0 max(mxi2)])
  %   subplot(nfig,numtimesteps,numtimesteps+jj);imagesc(0.5*mean(pcmind{nC,nP}{jj}.G_predCV{mtest(sortindI{jj}(end-2))},3)+0.8*mean(pcmind{nC,nP}{jj}.G_predCV{mtest(sortindI{jj}(end-1))},3)+mean(pcmind{nC,nP}{jj}.G_predCV{mtest(mxindI(jj))},3));  caxis([0 1.4*max(mxi2)])
  %   if jj==1,ylabel(['Gfit3 CV']);end
  %   subplot(nfig,numtimesteps,2*numtimesteps+jj);imagesc(0.8*mean(pcmind{nC,nP}{jj}.G_predCV{mtest(sortindI{jj}(end-1))},3)+mean(pcmind{nC,nP}{jj}.G_predCV{mtest(mxindI(jj))},3));  caxis([0 1.4*max(mxi2)])
  %   if jj==1,ylabel(['Gfit2 CV']);end
  
  for ss=1:22,
    GpredCVplot(:,:,ss)=pcmind{nC,nP}{jj}.G_predCV{mtest(mxindIeach(jj,ss))}(:,:,ss);
  end
  subplot(nfig,numtimesteps,1*numtimesteps+jj);imagesc(mean(GpredCVplot,3));  % caxis([0 0.4*max(mxi2)])
  
  %     subplot(nfig,numtimesteps,1*numtimesteps+jj);imagesc(mean(pcmind{nC,nP}{jj}.G_predCV{mtest(mxindI(jj))},3)); % caxis([0 0.7*max(mxi2)])
  if jj==1,ylabel(['Gfit1 CV']);end
  %   subplot(nfig,numtimesteps,28+jj);imagesc(Gmodel(:,:,mtest(mxindI(jj))));
  %   if jj==1,ylabel(['Gmodel win']);end
  %   subplot(nfig,numtimesteps,numtimesteps+jj);imagesc(mean(Gmodel(:,:,mtest(sortindI{jj}(end-0:end))),3));   caxis([0 3.5*max(mxi3)]) %caxis([0 3*max(mxi3)])
  %   if jj==1,ylabel(['Gmodel win']);end
  %   subplot(nfig,numtimesteps,42+jj);imagesc(mean(pcminde5{jj}.G_pred{mode(mnindi(jj,:))},3));
  %   if jj==1,ylabel(['Gfit AICmin']);end
  
  % subplot(nfig,numtimesteps,2*numtimesteps+jj);imagesc(G_predI_mixw(:,:,jj));  caxis([0 1*max(max(G_predI_mixw(:,:,jj)))])
  %   if jj==1,ylabel(['Gfit AICmix']);end
  %   subplot(nfig,numtimesteps,3*numtimesteps+jj);imagesc(mean(pcmind{nC,nP}{jj}.G_pred{1},3)); % caxis([0 max(mxi1)])
  subplot(nfig,numtimesteps,3*numtimesteps+jj);imagesc(mean(pcmind{nC,nP}{jj}.G_predCV{1},3)); % caxis([0 max(mxi1)])
  %   subplot(nfig,numtimesteps,70+jj);imagesc(H*mean(pcmind{nC,nP}{jj}.G_pred{1},3)*H'); caxis([0 max(mxi1)])
  if jj==1,ylabel(['Gfit Free']);end
  %   subplot(nfig,numtimesteps,84+jj);imagesc(mean(pcminde5{jj}.G_pred{17},3));
  %   if jj==1,ylabel(['Gfit All']);end
end


jj=5;figure;
for ss=1:22,
  GpredCVplot(:,:,ss)=pcmind{nC,nP}{jj}.G_predCV{mtest(mxindI(jj))}(:,:,ss);
  %   GpredCVplot(:,:,ss)=pcmind{nC,nP}{jj}.G_predCV{mtest(mxindIeach(jj,ss))}(:,:,ss);
  %   GpredCVplot(:,:,ss)=pcmind{nC,nP}{jj}.G_predCV{257}(:,:,ss);
end
for ss=1:22,subplot(3,22,ss);imagesc(G_empCV5(:,:,jj,ss));end
for ss=1:22,subplot(3,22,ss+22);imagesc(pcmind{nC,nP}{jj}.G_predCV{1}(:,:,ss));end
for ss=1:22,subplot(3,22,ss+44);imagesc(GpredCVplot(:,:,ss));end


%%
for jj=1:numtimesteps,
  GpredCVind(:,:,:,jj)=pcmind{nC,nP}{jj}.G_predCV{mtest(mxindI(jj))};
end
set(0,'DefaultFigurePosition',[15 50 1300 1300])
figure(101);
nn=0;for jj=1:numtimesteps,for ii=1:22,nn=nn+1;subplot(numtimesteps,22,nn);imagesc(GpredCVind(:,:,ii,jj));caxis([0 3]);end;end
figure(102);
nn=0;for jj=1:numtimesteps,for ii=1:22,nn=nn+1;subplot(numtimesteps,22,nn);imagesc(GpredCVind(:,:,ii,jj));end;end

for jj=1:numtimesteps,
  GpredCVindeach(:,:,:,jj)=pcmind{nC,nP}{jj}.G_predCV{mtest(mxindIeach(jj))};
end
set(0,'DefaultFigurePosition',[15 50 1300 1300])
figure(101);
nn=0;for jj=1:numtimesteps,for ii=1:22,nn=nn+1;subplot(numtimesteps,22,nn);imagesc(GpredCVindeach(:,:,ii,jj));caxis([0 2]);end;end
figure(102);
nn=0;for jj=1:numtimesteps,for ii=1:22,nn=nn+1;subplot(numtimesteps,22,nn);imagesc(GpredCVindeach(:,:,ii,jj));end;end


%% Plotting ThetaCV for sim7


% First, using model with all features
clear thetaCV*
mm=257;
for jj=1:numtimesteps,
  clear thetaCVtmp
  for nn=1:nP
    thetaCVtmp(:,:,nn)=pcmind{nC,nP}{jj}.thetaCV{mm}(nn:nP:end,:);  % subj x numfeat
    %       thetaCVrmstmp(:,1,nn)=rms(thetaCVtmp(:,12:13,nn),2); % TIW
  end
  thetaCV{mm}(:,:,jj)=mean(thetaCVtmp,3);
end

set(0,'DefaultFigurePosition',[15 50 1300 900])
figure(20);
induse=1;
for ff=1:8
  if ff<3
    indinc=0;
  else
    indinc=1;
  end
  subplot(4,2,ff);plot(50:25:450,  squeeze(mean(thetaCV{mm}(:,7+induse:7+induse+indinc,:),1)));
  ylim([-0.2 .2])
  induse=induse+1+indinc;
end


% Model with only 2 and 7, but also plot diagonals
clear thetaCV*
mm=22;
for jj=1:numtimesteps,
  clear thetaCVtmp
  for nn=1:nP
    thetaCVtmp(:,:,nn)=pcmind{nC,nP}{jj}.thetaCV{mm}(nn:nP:end,:);  % subj x numfeat
    %       thetaCVrmstmp(:,1,nn)=rms(thetaCVtmp(:,12:13,nn),2); % TIW
  end
  thetaCV{mm}(:,:,jj)=mean(thetaCVtmp,3);
end

set(0,'DefaultFigurePosition',[15 50 1300 900])
figure(21);
for ff=1:10
  subplot(5,2,ff);plot(50:25:450,  squeeze(mean(thetaCV{mm}(:,ff,:),1)));
  ylim([-0.2 .2])
end

% % Next, using best model for each time point. thetaCV may be NAN if feature not in that model.
% clear thetaCV*
% for jj=1:numtimesteps,
%     clear thetaCVtmp
%     for nn=1:nP
%       thetaCVtmp{jj}(:,:,nn)=pcmind{nC,nP}{jj}.thetaCV{mtest(mxindI(jj))}(nn:nP:end,:);  % subj x numfeat
% %       thetaCVrmstmp(:,1,nn)=rms(thetaCVtmp(:,12:13,nn),2); % TIW
%     end
%     thetaCV{jj}=mean(thetaCVtmp{jj},3);
% end





%%
% % Plotting ThetaCV for winning model (17 for sim4)
clear thetaCV*
mm=17;
for jj=1:numtimesteps,
  % %   for mm=unique(mtest(mxindI))
  clear thetaCVtmp
  for nn=1:nP
    thetaCVtmp(:,:,nn)=pcmind{nC,nP}{jj}.thetaCV{mm}(nn:nP:end,:);  % subj x numfeat
    thetaCVrmstmp(:,1,nn)=rms(thetaCVtmp(:,12:13,nn),2); % TIW
    thetaCVrmstmp(:,2,nn)=rms(thetaCVtmp(:,8:11,nn),2); % Asymmetry
    thetaCVrmstmp(:,3,nn)=rms(thetaCVtmp(:,5:7,nn),2); % SymPair
    thetaCVrmstmp(:,4,nn)=rms(thetaCVtmp(:,1:4,nn),2); % Identity
  end
  thetaCV{mm}(:,:,jj)=mean(thetaCVtmp,3);
  thetaCVrms{mm}(:,:,jj)=mean(thetaCVrmstmp,3);
  %   end
end

figure(20);  % thetaCV
subplot(4,1,1);plot(50:25:450,squeeze(mean(thetaCVrms{mm}(:,1,:),1)));
subplot(4,1,2);plot(50:25:450,squeeze(mean(thetaCVrms{mm}(:,2,:),1)));
subplot(4,1,3);plot(50:25:450,squeeze(mean(thetaCVrms{mm}(:,3,:),1)));
subplot(4,1,4);plot(50:25:450,squeeze(mean(thetaCVrms{mm}(:,4,:),1)));

figure(21);  % thetaCV
subplot(4,1,1);plot(50:25:450,squeeze(mean(rms(thetaCV{17}(:,12:13,:),2),1)));
subplot(4,1,2);plot(50:25:450,squeeze(mean(rms(thetaCV{17}(:,8:11,:),2),1)));
subplot(4,1,3);plot(50:25:450,squeeze(mean(rms(thetaCV{17}(:,5:7,:),2),1)));
subplot(4,1,4);plot(50:25:450,squeeze(mean(rms(thetaCV{17}(:,1:4,:),2),1)));

figure(22);  % thetaCV
subplot(4,1,1);plot(50:25:450,(mean([squeeze(mean(thetaCV{17}(:,12:13,:),1))].^2).^0.5));
subplot(4,1,2);plot(50:25:450,(mean([squeeze(mean(thetaCV{17}(:,8:11,:),1))].^2).^0.5));
subplot(4,1,3);plot(50:25:450,(mean([squeeze(mean(thetaCV{17}(:,5:7,:),1))].^2).^0.5));
subplot(4,1,4);plot(50:25:450,(mean([squeeze(mean(thetaCV{17}(:,1:4,:),1))].^2).^0.5));

figure(23);  % thetaCV
subplot(4,1,1);plot(50:25:450,(([squeeze(mean(thetaCV{17}(:,12:13,:),1))].^2).^0.5));
subplot(4,1,2);plot(50:25:450,(([squeeze(mean(thetaCV{17}(:,8:11,:),1))].^2).^0.5));
subplot(4,1,3);plot(50:25:450,(([squeeze(mean(thetaCV{17}(:,5:7,:),1))].^2).^0.5));
subplot(4,1,4);plot(50:25:450,(([squeeze(mean(thetaCV{17}(:,1:4,:),1))].^2).^0.5));

figure(24);  % thetaCV
subplot(4,1,1);plot(50:25:450,squeeze(mean(abs(thetaCV{17}(:,12:13,:)),1)))
subplot(4,1,2);plot(50:25:450,squeeze(mean(abs(thetaCV{17}(:,8:11,:)),1)))
subplot(4,1,3);plot(50:25:450,squeeze(mean(abs(thetaCV{17}(:,5:7,:)),1)))
subplot(4,1,4);plot(50:25:450,squeeze(mean(abs(thetaCV{17}(:,1:4,:)),1)))

figure(25);  % thetaCV
subplot(4,1,1);plot(50:25:450,squeeze(mean((thetaCV{17}(:,12:13,:)),1)))
subplot(4,1,2);plot(50:25:450,squeeze(mean((thetaCV{17}(:,8:11,:)),1)))
subplot(4,1,3);plot(50:25:450,squeeze(mean((thetaCV{17}(:,5:7,:)),1)))
subplot(4,1,4);plot(50:25:450,squeeze(mean((thetaCV{17}(:,1:4,:)),1)))

%%

% Winning models
timesteps=50:30:440;
nfig=2;
set(0,'DefaultFigurePosition',[15 450 1300 300+100*(nfig-3)])
jjtest=1:numtimesteps;
jjtest=[1 6 12];
jjcnt=0;
for jj=jjtest
  jjcnt=jjcnt+1;
  figure(300);
  subplot(nfig,numel(jjtest),jjcnt);imagesc(mean(G_emp(:,:,jj,:),4));caxis([-10 60]);title([num2str(timesteps(jj)) ' ms'])
  if jj==1,ylabel(['Gemp']);end
  subplot(nfig,numel(jjtest),numel(jjtest)+jjcnt);imagesc(mean(G_empCV5(:,:,jj,:),4));caxis([-.001 .01]);
  if jj==1,ylabel(['GempCV']);end
  %   [mx,mxind(jj)]=max(mean(pcminde5{jj}.D.likelihood-repmat(pcminde5{jj}.D.likelihood(:,2),[1 17])));
  %   subplot(nfig,14,numel(jjtest)*2+jj);imagesc(mean(pcminde5{jj}.G_predCV{mxind(jj)},3));
  %   if jj==1,ylabel(['GfitCV ML']);end
  %   subplot(nfig,14,numel(jjtest)*3+jj);imagesc(mean(pcminde5{jj}.G_predCV{1},3));
  %   if jj==1,ylabel(['GfitCV Free']);end
  
  %   subplot(nfig,14,numel(jjtest)*4+jj);imagesc(mean(pcminde5{jj}.G_predCV{17},3));
  %   if jj==1,ylabel(['GfitCV All']);end
  
end

%%
nfig=5;
set(0,'DefaultFigurePosition',[15 450 1300 300+100*(nfig-3)])
for jj=1:14,
  [mx,mxind(jj)]=max(mean(pcmind5{jj}.D.likelihood-repmat(pcmind5{jj}.D.likelihood(:,2),[1 nM])));
  figure(400);
  subplot(nfig,14,jj);imagesc(H*mean(G_emp(:,:,jj,:),4)*H');title([num2str(timesteps(jj)) ' ms'])
  if jj==1,ylabel(['Gemp']);end
  subplot(nfig,14,14+jj);imagesc(H*mean(G_empCV5(:,:,jj,:),4)*H');
  if jj==1,ylabel(['GempCV']);end
  subplot(nfig,14,28+jj);imagesc(H*mean(pcmind5{jj}.G_predCV{mxind(jj)},3)*H');
  if jj==1,ylabel(['GfitCV ML']);end
  subplot(nfig,14,42+jj);imagesc(H*mean(pcmind5{jj}.G_predCV{1},3)*H');
  if jj==1,ylabel(['GfitCV Free']);end
  subplot(nfig,14,56+jj);imagesc(H*mean(pcmind5{jj}.G_predCV{17},3)*H');
  if jj==1,ylabel(['GfitCV All']);end
  
end



nfig=4;
nspf=11; % number of subjects per figure
for jj=1:numtimesteps
  for iireal=1:length(datacat)
    %   for ii=1:nspf
    if iireal/nspf>1
      set(0,'DefaultFigurePosition',[750 150 700 850])
      figure(jj+3100);
      ii=iireal-11;
    else
      set(0,'DefaultFigurePosition',[50 150 700 850])
      figure(jj+3000);
      ii=iireal;
    end
    subplot(nspf,nfig,nfig*(ii-1)+1);imagesc(G_emp(:,:,jj,iireal));colorbar;if ii==1,title(['G emp TP' num2str(jj)]);end
    ylabel(['Subj ' num2str(iireal)])
    subplot(nspf,nfig,nfig*(ii-1)+2);imagesc(G_empCV3(:,:,jj,iireal));colorbar;if ii==1,title('G empCV');end
    subplot(nspf,nfig,nfig*(ii-1)+3);imagesc(pcmind3{jj}.G_pred{1}(:,:,iireal));colorbar;if ii==1,title('G fit ind');end
    subplot(nspf,nfig,nfig*(ii-1)+4);imagesc(pcmind3{jj}.G_predCV{1}(:,:,iireal));colorbar;if ii==1,title('G fitCV ind');end
    %     subplot(nspf,nfig,nfig*(ii-1)+5);imagesc(pcmfit5{jj}.Gpred{1});colorbar;if ii==1,title('G fit grp');end
  end
  for iireal=1:length(datacat)
    %   for ii=1:nspf
    if iireal/nspf>1
      set(0,'DefaultFigurePosition',[750 150 700 850])
      figure(jj+5100);
      ii=iireal-11;
    else
      set(0,'DefaultFigurePosition',[50 150 700 850])
      figure(jj+5000);
      ii=iireal;
    end
    subplot(nspf,nfig,nfig*(ii-1)+1);imagesc(G_emp(:,:,jj,iireal));colorbar;if ii==1,title(['G emp TP' num2str(jj)]);end
    ylabel(['Subj ' num2str(iireal)])
    subplot(nspf,nfig,nfig*(ii-1)+2);imagesc(G_empCV5(:,:,jj,iireal));colorbar;if ii==1,title('G empCV');end
    subplot(nspf,nfig,nfig*(ii-1)+3);imagesc(pcmind5{jj}.G_pred{1}(:,:,iireal));colorbar;if ii==1,title('G fit ind');end
    subplot(nspf,nfig,nfig*(ii-1)+4);imagesc(pcmind5{jj}.G_predCV{1}(:,:,iireal));colorbar;if ii==1,title('G fitCV ind');end
    %     subplot(nspf,nfig,nfig*(ii-1)+5);imagesc(pcmfit5{jj}.Gpred{1});colorbar;if ii==1,title('G fit grp');end
  end
end
for jj=1:numtimesteps
  nfig=7;
  for iireal=1:length(datacat)
    %   for ii=1:nspf
    if iireal/nspf>1
      set(0,'DefaultFigurePosition',[750 150 700 850])
      figure(jj+5100);
      ii=iireal-11;
    else
      set(0,'DefaultFigurePosition',[50 150 700 850])
      figure(jj+5000);
      ii=iireal;
    end
    subplot(nspf,nfig,nfig*(ii-1)+1);imagesc(G_emp(:,:,jj,iireal));colorbar;if ii==1,title(['G emp TP' num2str(jj)]);end
    ylabel(['Subj ' num2str(iireal)])
    subplot(nspf,nfig,nfig*(ii-1)+2);imagesc(G_empCV3(:,:,jj,iireal));colorbar;if ii,title('GempCV3P');end
    subplot(nspf,nfig,nfig*(ii-1)+3);imagesc(pcmind3{jj}.G_pred{1}(:,:,iireal));colorbar;if ii,title('GfitI3');end
    subplot(nspf,nfig,nfig*(ii-1)+4);imagesc(pcmind3{jj}.G_predCV{1}(:,:,iireal));colorbar;if ii,title('GfitICV3');end
    subplot(nspf,nfig,nfig*(ii-1)+5);imagesc(G_empCV5(:,:,jj,iireal));colorbar;if ii,title('GempCV5P');end
    subplot(nspf,nfig,nfig*(ii-1)+6);imagesc(pcmind5{jj}.G_pred{1}(:,:,iireal));colorbar;if ii,title('GfitI5');end
    subplot(nspf,nfig,nfig*(ii-1)+7);imagesc(pcmind5{jj}.G_predCV{1}(:,:,iireal));colorbar;if ii,title('GfitICV5');end
  end
  
  for iireal=1:length(datacat)
    %   for ii=1:nspf
    if iireal/nspf>1
      set(0,'DefaultFigurePosition',[750 150 700 850])
      figure(jj+1100);
      ii=iireal-11;
    else
      set(0,'DefaultFigurePosition',[50 150 700 850])
      figure(jj+1000);
      ii=iireal;
    end
    subplot(nspf,nfig,nfig*(ii-1)+1);imagesc(H*G_emp(:,:,jj,iireal)*H');colorbar;if ii==1,title(['G emp TP' num2str(jj)]);end
    ylabel(['Subj ' num2str(iireal)])
    subplot(nspf,nfig,nfig*(ii-1)+2);imagesc(H*G_empCV3(:,:,jj,iireal)*H');colorbar;if ii,title('GempCV3P');end
    subplot(nspf,nfig,nfig*(ii-1)+3);imagesc(H*pcmind3{jj}.G_pred{1}(:,:,iireal)*H');colorbar;if ii,title('GfitI3');end
    subplot(nspf,nfig,nfig*(ii-1)+4);imagesc(H*pcmind3{jj}.G_predCV{1}(:,:,iireal)*H');colorbar;if ii,title('GfitICV3');end
    subplot(nspf,nfig,nfig*(ii-1)+5);imagesc(H*G_empCV5(:,:,jj,iireal)*H');colorbar;if ii,title('GempCV5P');end
    subplot(nspf,nfig,nfig*(ii-1)+6);imagesc(H*pcmind5{jj}.G_pred{1}(:,:,iireal)*H');colorbar;if ii,title('GfitI5');end
    subplot(nspf,nfig,nfig*(ii-1)+7);imagesc(H*pcmind5{jj}.G_predCV{1}(:,:,iireal)*H');colorbar;if ii,title('GfitICV5');end
  end
  
end

% with H*G*H'
nfig=3;
nspf=11; % number of subjects per figure
for jj=1:numtimesteps
  for iireal=1:length(datacat)
    %   for ii=1:nspf
    if iireal/nspf>1
      set(0,'DefaultFigurePosition',[1200 50 700 950])
      figure(1000+jj+100);
      ii=iireal-11;
    else
      set(0,'DefaultFigurePosition',[500 50 700 950])
      figure(1000+jj);
      ii=iireal;
    end
    G_emp(:,:,jj,iireal)=cov([designX{iireal}'*datacat_noisenorm{iireal}(:,:,jj)]',1);
    G_empCV(:,:,jj,iireal)=pcm_estGCrossval(datacat_noisenorm{iireal}(:,:,jj),partVec5{iireal},designX{iireal});
    subplot(nspf,nfig,nfig*(ii-1)+1);imagesc(H*G_emp(:,:,jj,iireal)*H');colorbar;if ii==1,title(['G emp TP' num2str(jj)]);end
    ylabel(['Subj ' num2str(iireal)])
    subplot(nspf,nfig,nfig*(ii-1)+2);imagesc(H*G_empCV(:,:,jj,iireal)*H');colorbar;if ii==1,title('G empCV');end
    subplot(nspf,nfig,nfig*(ii-1)+3);imagesc(H*pcmind5{jj}.G_pred{1}(:,:,iireal)*H');colorbar;if ii==1,title('G fit ind');end
    %     subplot(nspf,nfig,nfig*(ii-1)+4);imagesc(H*pcmfit{jj}.Gpred{1}*H');colorbar;if ii==1,title('G fit grp');end
  end
end

if 0
  for jj=1:numtimesteps
    % Y=datacat;
    % Y=datacat_meansubt;
    for ii=1:length(datacat_meansubt_noisenorm)
      %       Y{ii}=datacat_meansubt_noisenorm{ii}(:,:,jj);
      Y{ii}=datacat_noisenorm{ii}(:,:,jj); % keep mean in
    end
    [pcmfit{jj}.Tgroup,pcmfit{jj}.theta,pcmfit{jj}.Gpred] = pcm_fitModelGroup(Y,M,partVec,designX,'runEffect','fixed','fitScale',1);
    [pcmfit{jj}.Tcross,pcmfit{jj}.thetaCr,pcmfit{jj}.G_predcv] = pcm_fitModelGroupCrossval(Y,M,partVec,designX,'runEffect','fixed','groupFit',pcmfit{jj}.theta,'fitScale',1);
    save([edir 'pcm_datafit.mat'],'pcmfit');
  end
else % bluebear
  % audtac_pcm_run
  for jj=1:numtimesteps
    load([pcdir '/real_mean_kept_feature/pcmF_datafit_jj' num2str(jj) '.mat']); % using real partVec
    pcmfit{jj}=A;
  end
end

p=[700 400 1200 550];
set(0,'DefaultFigurePosition',p)
for jj=1:25;figure(jj);mmind=0;
  for mm=[1 3:17],
    mmind=mmind+1;subplot(4,4,mmind);imagesc(pcmfit{jj}.Gpred{mm});
    colorbar
  end;
end

for jj=1:numtimesteps
  for mm=1:length(M)
    for subcnt=1:length(datacat)
      [aicg(mm,jj,subcnt),bicg(mm,jj,subcnt)] = aicbic(pcmfit{jj}.Tgroup.likelihood(subcnt,mm),M{mm}.numGparams,numobs(subcnt));
      %       [aicc(mm,jj,subcnt),bicc(mm,jj,subcnt)] = aicbic(pcmfit{jj}.Tcross.likelihood(subcnt,mm),M{mm}.numGparams,numobs(subcnt));
      LikGroup(mm,jj,subcnt)=pcmfit{jj}.Tgroup.likelihood(subcnt,mm);
    end
  end
  maicg=mean(aicg(:,jj,:),3);
  mmaicg=min(maicg);
  aicdif(:,jj)=maicg-mmaicg; % reject if >10
  aicrat(:,jj)=exp(.5*(mmaicg-maicg)); %reject if <.007
end
% addpath('D:/nutmeg_svn/'); % need nut_rownorm.m
% [srt,ind]=sort(nut_rownorm(aicdif(:,:)))
% figure;imagesc(aicdif(ind(1:end-1),:))


switch M{3}.type
  case 'component'
    Nceil=2;
    Nnull=17;
    close all;
    for jj=1:numtimesteps
      figure(1);subplot(5,5,jj);
      pcmfit{jj}.Tg = pcm_plotModelLikelihood(pcmfit{jj}.Tgroup,M,'upperceil',pcmfit{jj}.Tgroup.likelihood(:,Nceil),'normalize',0,'Nnull',Nnull,'Nceil',Nceil);title([num2str(jj*20) ' ms'])
      figure(4);subplot(5,5,jj);
      plot(mean(aicg(:,jj,:)-repmat(aicg(Nnull,jj,:),[17 1 1]),3),'o');ylim([-5 5]);title([num2str(jj*20) ' ms'])
    end
    for jj=1:numtimesteps
      try
        figure(2);subplot(5,5,jj);
        pcmfit{jj}.Tc = pcm_plotModelLikelihood(pcmfit{jj}.Tcross,M,'upperceil',pcmfit{jj}.Tgroup.likelihood(:,Nceil),'normalize',0,'Nnull',Nnull,'Nceil',Nceil);title([num2str(jj*20) ' ms'])
        figure(3);subplot(5,5,jj);
        plot(mean(pcmfit{jj}.Tc.likelihood_norm,1),'o');ylim([-.2 1]);title([num2str(jj*20) ' ms'])
      end
    end
  case 'feature'
end


likemat=nan(25,16);
for jj=1:numtimesteps
  for mm=1:length(M)
    lookmodel(jj,mm)=any(pcmfit{jj}.theta{mm}(1:M{mm}.numGparams)>-11);
    if lookmodel(jj,mm)
      likemat(jj,mm)=mean(pcmfit{jj}.Tgroup.likelihood(:,mm))/mean(pcmfit{jj}.Tgroup.likelihood(:,1));
    end
  end
end

%
