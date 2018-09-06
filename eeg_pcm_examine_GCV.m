% eeg_pcm_examine_GCV

clear
% audtac_startup
sub{8}='e08';


sleep=0;
iter=27;
trialkc=-1;
soalist=[1 3 4 5 6 7 9];  % indicies of 7 different conditions
tt=3;
ss=10;

numtimesteps=numel(50:30:440);  % every 30 ms
subcnt=0;
% for ii=subuseall
ii=8; % just run for one subject
subcnt=subcnt+1;
load(['tlock_pcm_' sub{ii} '_sleep' num2str(sleep) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat']);
% Extract each difference across time windows of ERP, normalise.
llcnt=0;
trlcnt=0;
for ll=soalist
  llcnt=llcnt+1;
  for jj=1:numtimesteps; % every 30 ms over a window of 100ms duration
    tlockdiff_slide{llcnt}(:,:,jj)=mean(pcm_TPAmMSPN{ll,tt,ss}.trial(:,:,20+jj*30-49:20+jj*30+50),3);
  end
  designX{subcnt}(1+trlcnt:size(tlockdiff_slide{llcnt},1)+trlcnt,llcnt)=1;
  trlcnt=trlcnt+size(tlockdiff_slide{llcnt},1);
end
datacat{subcnt}=cat(1,tlockdiff_slide{:}); % non-equal number of trials per condition
partVec5{subcnt}=[rem(1:size(datacat{subcnt},1),5)+1]';  % this is nonequal number of trials per partition & condition; 
partVec2{subcnt}=[rem(1:size(datacat{subcnt},1),2)+1]';  % this is nonequal number of trials per partition & condition; 


% Multivariate noise normalisation
for iChannel = 1:size(datacat{subcnt},2)
    for jj=1:numtimesteps
        [~,~,r(:,iChannel,jj)] = regress(datacat{subcnt}(:,iChannel,jj),designX{subcnt});
    end
end
df = size(designX{subcnt},1) - size(designX{subcnt},2);
for jj=1:numtimesteps
    sigma(:,:,jj) = rsa.stat.covdiag(r(:,:,jj),df);  % RSA toolbox
end
for jj=1:numtimesteps
  datacat_noisenorm{subcnt}(:,:,jj)=datacat{subcnt}(:,:,jj)*(sigma(:,:,jj)^-0.5);
end

jj=14;
pcm_estGCrossval(datacat_noisenorm{subcnt}(:,:,jj),partVec5{subcnt},designX{subcnt})
% Note that there are several off-diagonal elements greater than either of
% the main-diagonal elements it is associate with (e.g. [1,2] or [5,6])

jj=2;
pcm_estGCrossval(datacat_noisenorm{subcnt}(:,:,jj),partVec5{subcnt},designX{subcnt})
% Note several negative values, including on main-diagonal!

