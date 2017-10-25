% Setup variables
clear
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
  idir='D:\audtac\sim_results\';
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
    idir='/home/zumerj/audtac/sim_results/';
  else % assume on VM linux of psychl-132432
    edir='/mnt/hgfs/D/audtac/eeg_data/';
    esdir='/mnt/hgfs/D/audtac/source_data/';
    ddir='/mnt/hgfs/D/audtac/legomagic/diaries/';
    bdir='/mnt/hgfs/D/audtac/behav_data/';
    sdir='/mnt/hgfs/D/audtac/spss_stuff/';
    fdir='/mnt/hgfs/D/audtac/figs/';
    mdir='/mnt/hgfs/D/audtac/structural_MRI/';
    pdir='/mnt/hgfs/D/audtac/polhemus/';
    idir='/mnt/hgfs/D/audtac/sim_results/';
  end
end
cd(edir)

% sub{100}='p01';
sub{1}='e01';
sub{2}='e02';
sub{3}='e03';
sub{4}='e04';
%%%  above is pilot, below is real
sub{5}='e05';
sub{6}='e06';
sub{7}='e07';
sub{8}='e08';
sub{9}='e09';
sub{10}='e10';
sub{11}='e11';
sub{12}='e12';
sub{13}='e13';
sub{14}='e14';
sub{15}='e15';
sub{16}='e16';
sub{17}='e17';
sub{18}='e18';
sub{19}='e19';
sub{20}='e20';
sub{21}='e21';
sub{22}='e22';
sub{23}='e23';
sub{24}='e24';
sub{25}='e25';
sub{26}='e26';
sub{27}='e27';
sub{28}='e28';
sub{29}='e29';
sub{30}='e30';
sub{31}='e31';
sub{32}='e32';

if ispc
  warning off
  rmpath(genpath('D:\matlab\spm8\external\fieldtrip\'))
  rmpath(genpath('D:\fieldtrip_svn\'))
  rmpath(genpath('D:\fieldtrip_git\'))
  warning on
  addpath('D:\fieldtrip_git\')
else
  warning off
  rmpath(genpath('/mnt/hgfs/D/matlab/spm8/external/fieldtrip/'))
  rmpath(genpath('/mnt/hgfs/D/fieldtrip_svn/'))
  rmpath(genpath('/mnt/hgfs/D/fieldtrip_git/'))
  warning on
  addpath('/mnt/hgfs/D/fieldtrip_git/')
end
which ft_defaults.m
ft_defaults;
addpath('D:\Matlab\from_mate');
% addpath('D:\Matlab\rsatoolbox\Engines');
addpath('D:\Matlab\rsatoolbox');
addpath('D:\Matlab\NearestSymmetricPositiveDefinite')
addpath('D:\Matlab\pcm_toolbox\'); % https://github.com/jdiedrichsen/pcm_toolbox

load([edir 'iikeep.mat'])
soades=[-.5 nan -.07 -.02 0 .02 .07 nan .5];

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

numtimesteps=25;
subcnt=0;
for ii=subuseall
  subcnt=subcnt+1;
  clear tlockdiff_slide data_condmean
  cd([edir sub{ii} ])
  
  % Load individual unaveraged data
  load(['tlock_pcm_' sub{ii} '_sleep' num2str(sleep) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat']);
  
  % Extract each difference across time windows of ERP, normalise.
  llcnt=0;
  trlcnt=0;
  numPart=5; % number of partitions for crossval
  for ll=soalist
    llcnt=llcnt+1;
    %     startind1=dsearchn(tlock_tacPaud{ll,tt,10}.time',timestart(ll));
    %     endind1=dsearchn(tlock_tacPaud{ll,tt,10}.time',timeend(ll));
    %     startind2=dsearchn(tlock_tacMSpN{ll,tt,10}.time',timestart(ll));
    %     endind2=dsearchn(tlock_tacMSpN{ll,tt,10}.time',timeend(ll));
    %     tlockdiff=tlock_tacPaud{ll,tt,10}.avg(:,startind1:endind1)-tlock_tacMSpN{ll,tt,10}.avg(:,startind2:endind2);
    for jj=1:numtimesteps; % every 20 ms
      tlockdiff_slide{llcnt}(:,:,jj)=mean(pcm_TPAmMSPN{ll,tt,ss}.trial(:,:,jj*20-19:jj*20+1),3);
    end
    designX{subcnt}(1+trlcnt:size(tlockdiff_slide{llcnt},1)+trlcnt,llcnt)=1;
    trlcnt=trlcnt+size(tlockdiff_slide{llcnt},1);
    data_condmean_squeeze{llcnt}=mean(tlockdiff_slide{llcnt},1);
    data_condmean{llcnt}=repmat(mean(tlockdiff_slide{llcnt},1),[size(tlockdiff_slide{llcnt},1) 1 1]);
  end
  datacat{subcnt}=cat(1,tlockdiff_slide{:}); % non-equal number of trials per condition; okay?
  partVec{subcnt}=[rem(1:size(datacat{subcnt},1),numPart)+1]';  % this is nonequal number of trials per partition & condition; okay?
%   partVec{subcnt}=ones(size(datacat{subcnt},1),1);
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
  
end % ii

save([edir 'pcm_datacat.mat'],'datacat*','designX','partVec');

%% create models

pcm_create_models;

%% call Model Group Fit

load([edir 'pcm_datacat.mat']);
load([edir 'pcm_models.mat']);
numtimesteps=size(datacat{1},3);

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

if 0
  for jj=1:numtimesteps
    % Y=datacat;
    % Y=datacat_meansubt;
    for ii=1:length(datacat_meansubt_noisenorm)
      Y{ii}=datacat_meansubt_noisenorm{ii}(:,:,jj);
    end
    [pcmfit{jj}.Tgroup,pcmfit{jj}.theta,pcmfit{jj}.Gpred] = pcm_fitModelGroup(Y,M,partVec,designX,'runEffect','fixed','fitScale',1);
    [pcmfit{jj}.Tcross,pcmfit{jj}.thetaCr,pcmfit{jj}.G_predcv] = pcm_fitModelGroupCrossval(Y,M,partVec,designX,'runEffect','fixed','groupFit',pcmfit{jj}.theta,'fitScale',1);
    save([edir 'pcm_datafit.mat'],'pcmfit');
  end
else % bluebear
  % audtac_pcm_run
  for jj=1:numtimesteps
    load([edir 'pcm_datafit_jj' num2str(jj) '.mat']); % using real partVec
    pcmfit{jj}=A;
  end
end


close all;figure
for jj=1:numtimesteps
  subplot(5,5,jj);
  pcmfit{jj}.T = pcm_plotModelLikelihood(pcmfit{jj}.Tcross,M,'upperceil',pcmfit{jj}.Tgroup.likelihood(:,2),'normalize',0,'Nnull',17,'Nceil',2);
end
for jj=1:numtimesteps
  subplot(5,5,jj);plot(mean(pcmfit{jj}.T.likelihood_norm,1),'o');ylim([-.2 1]);title([num2str(jj*20) ' ms'])
end

for jj=1:numtimesteps
  for mm=1:length(M)
    for subcnt=1:length(datacat)
      [~,bicg(mm,jj,subcnt)] = aicbic(pcmfit{jj}.Tgroup.likelihood(subcnt,mm),M{mm}.numGparams,numobs(subcnt));
      [~,bicc(mm,jj,subcnt)] = aicbic(pcmfit{jj}.Tcross.likelihood(subcnt,mm),M{mm}.numGparams,numobs(subcnt));
    end
  end
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
