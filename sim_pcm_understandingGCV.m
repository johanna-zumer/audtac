% Simulated data for PCM based on EEG study with 64 channels and 7
% audio-tactile asynchronies
%
% Johanna Zumer and Remi Gau, 2017

clear
% Modify for your location
addpath('D:\Matlab\pcm_toolbox\')
addpath('D:\Matlab\rsatoolbox')  % used for noise normalisation

%% Main flags to change

plotflag=0;
loaddata=0; % =0 creates simulated data from scratch; =1 loads saved data
runflag='run'; % set to 'run' or 'load': PCM fit model to data


%%
% Model:      Y=Z*U + X*B + E;
D.numPart = 5;
D.numVox  = 64; % 64 channel EEG
ncond=7; % number of conditions
nsub=20; % number of subjects


%% Model creation based on features

% creates all models to test; model 17 has all features
saveflag=0;
create_sim3_featAllModels

for mm=1:numel(M)
%   M{mm}.fitAlgorithm='NR'; % NR or minimize
  M{mm}.fitAlgorithm='minimize'; % NR or minimize
end

%% Create/load simulated data with different weighting of features

switch loaddata
  case 1
    load('simulation1_data.mat');
  case 0
    
    % Activation U different for every subject, but same pattern.
    % 10 different models (of varying weights between 4 component G_ii)
    theta_feat=zeros(5,M{17}.numGparams);
    theta_feat(1,[1:3]) = 1;
    theta_feat(2,[2:3]) = 1;
    theta_feat(3,[1 4:6]) = 1;
    theta_feat(4,[2:10]) = 1;
    theta_feat(5,[1:10]) = 1;
    
    scale=abs(1 +10*randn(1,nsub)); % vary SNR 
    noise=abs(1 +100*randn(1,nsub));
    
    for tr=1:size(theta_feat,1)
      [Y{tr},part,conditions] = pcm_generateData(M{17},theta_feat(tr,:)',D,nsub,scale,noise);
    end
    
    Z=pcm_indicatorMatrix('identity',conditions);
    
    % add condition means and multivariate noise normalise
    X=Z;
    B=-ones(7,1);
    B=zeros(7,1);
    B=[1; 1.2; 1.4; 1.6; 1.4; 1.2; 1];
%     B=1000*ones(7,1);
    for tr=1:size(theta_feat,1)
      for ss=1:nsub
        % add condition mean
        Ym{tr}{ss}=Y{tr}{ss}+repmat(X*B,[1 64]);
        
        % multivariate noise normalise
        for iChannel = 1:D.numVox
          [~,~,r(:,iChannel)] = regress(Ym{tr}{ss}(:,iChannel),Z);
        end
        df = size(Z,1) - size(Z,2);
        sigma = rsa.stat.covdiag(r,df);
        Ym_nn{tr}{ss}=Ym{tr}{ss}*(sigma^-0.5);
      end
    end
    save('simulation1_data.mat','Ym_nn','part','conditions')
end

%%

for tr=1:size(theta_feat,1)
  for ss=1:nsub
    Gemp(:,:,tr,ss)=cov([Z'*Ym_nn{tr}{ss}]',1);
    GempCV(:,:,tr,ss)=pcm_estGCrossval(Ym_nn{tr}{ss},part,Z);
  end
end

nfig=2;
set(0,'DefaultFigurePosition',[15 450 1300 350+100*(nfig-3)])
figure;
for tr=1:size(theta_feat,1)
  subplot(nfig,size(theta_feat,1),tr);imagesc(mean(Gemp(:,:,tr,:),4));  %caxis([0.5*min(Gemp(:)) 0.5*max(Gemp(:))]);
  title([num2str(tr) ' thetaval'])
  if tr==1,ylabel(['Gemp ' ]);end
  subplot(nfig,size(theta_feat,1),size(theta_feat,1)+tr);imagesc(mean(GempCV(:,:,tr,:),4));  %caxis([0.5*min(GempCV(:)) 1*max(GempCV(:))]);
  if tr==1,ylabel(['GempCV ' ]);end
end


