% Simulated data for AT 7 asynchrony dataset for PCM
% Johanna Zumer, 2017

addpath('D:\Matlab\pcm_toolbox\'); % https://github.com/jdiedrichsen/pcm_toolbox
% addpath('D:\Matlab\NearestSymmetricPositiveDefinite')
% addpath('D:/nutmeg_svn/'); % just need nut_rownorm.m
% addpath('D:\Matlab\from_mate')
% addpath('D:\Matlab\rsatoolbox')
plotflag=0;
clear Y*

%%
% Model:      Y=Z*U + X*B + E;

Z=zeros(350,7); % assume 7 conditions with 50 trials per condition, same for all subjects
Z(1:50,1)=1;
Z(51:100,2)=1;
Z(101:150,3)=1;
Z(151:200,4)=1;
Z(201:250,5)=1;
Z(251:300,6)=1;
Z(301:350,7)=1;

sim1_RDM_G

% Activation U different for every subject, but same pattern.
% Take in to account modelling noise, whether AT vs A+T is difference than AT+N vs A+T
for ss=1:22  % pretend 22 subjects
  
  % %   if 0  % option 1: generate 4 separate conditions
  % %     % increased MS integration (e.g. AT+N greater than A+T as function of TIW, over all channels)
  % %     AAinit=demean(randn(7,64),2);
  % %     TTinit=demean(randn(7,64),2);
  % %     ATinit=demean(randn(7,64),2);
  % %     NNinit=demean(randn(7,64),2);
  % %
  % %     sig=3; % snr
  % %     % importnant that each condition has zero-mean over all channels, per
  % %     % condition (but can add mean later to simulate that aspect)
  % %     AA=[repmat(sig*[.0 .3  .3 .3  .3 .3  .3]',[1 32]) repmat(-sig*[.0 .3  .3 .3  .3 .3  .3]',[1 32])]+AAinit;
  % %     TT=[repmat(sig*[.1 .1  .1 .1  .1 .1  .0]',[1 32]) repmat(-sig*[.1 .1  .1 .1  .1 .1  .0]',[1 32])]+TTinit;
  % %     AT=[repmat(sig*[.1 .45 .5 .55 .5 .45 .3]',[1 32]) repmat(-sig*[.1 .45 .5 .55 .5 .45 .3]',[1 32])]+ATinit;
  % %     NN=demean(.1*randn(7,64),2)+NNinit;  % is this a valid 'null' model?
  % %
  % %     AAco=zeros(7);
  % %     AAco(2:7,2:7)=.3;
  % %     AAco(1,1)=.3;
  % %     AAco=AAco+.7*eye(7);
  % %     AAuse=AAco*AA;
  % %
  % %     TTco=zeros(7);
  % %     TTco(1:6,1:6)=.1;
  % %     TTco(7,7)=.1;
  % %     TTco=TTco+.9*eye(7);
  % %     TTuse=TTco*TT;
  % %
  % %     ATco=zeros(7);
  % %     ATco(2:6,2:6)=.3;
  % %     ATco(3:5,3:5)=ATco(3:5,3:5)+.1;
  % %     ATco(2,6)=.5;
  % %     ATco(6,2)=.5;
  % %     ATco(1:4,1:4)=ATco(1:4,1:4)+.3;
  % %     ATco(4:7,4:7)=ATco(4:7,4:7)+.3;
  % %     for cc=1:7,ATco(cc,cc)=1;end
  % %     ATuse=ATco*AT;
  % %
  % %     MScon=(AAuse+TTuse)-(ATuse+NN); % Correct Mulitsensory contrast
  % %     Ccon=(AAuse+TTuse)-(ATuse);  % Cecere et al. 2017 contrast
  
  % option 2: generate multisensory difference from the start
  
  
  % simulated model weights.  Gtotal = SUM_ii [exp(theta_ii)*G_ii]
  % order of G:  eye, TIW, Asymm, Sym_pair
  theta_real(1,:) = [0 -inf -inf -inf];
  theta_real(2,:) = [-inf 0 -inf -inf];
  theta_real(3,:) = [-inf -inf 0 -inf];
  theta_real(4,:) = [-inf -inf -inf 0];
  
  theta_real(5,:) = [0 0 -inf -inf];
  theta_real(6,:) = [-inf -inf 0 0 ];
  theta_real(7,:) = [-inf 0 0 -inf];
  
  theta_real(8,:) = [0 0 0 -inf];
  theta_real(9,:) = [0 2 0 -inf];
  theta_real(10,:) = [0 0 0 0];
  
  % X*B is condition mean
  X=Z;
  B=[1; 1.2; 1.4; 1.6; 1.4; 1.2; 1];
  
  % Average scaling factor (first column) and noise level (second column) for the G matrix
  Scale_noise=[100 1;   10 1;    1 1;    0.1 1];
  
  % Generate scaling factor for each subject (distrib around group mean)
  % take absolute values to make sure we get positive matrices in the end
  
  
  for tr=1:size(theta_real,1)
    G_sum(:,:,tr) = exp(theta_real(tr,1))*G_eye + exp(theta_real(tr,2))*G_1a + exp(theta_real(2))*G_1b + exp(theta_real(tr,3))*G_2a + exp(theta_real(tr,4))*G_3a;
    
    for sm=1:size(Scale_noise,1)
      theta_s(sm,ss)   = abs(Scale_noise(sm,1)+randn*.1);
      theta_sig(sm,ss) = abs(Scale_noise(sm,2)+randn*.1);
      
      V = Z*G_sum(:,:,tr)*theta_s(sm,ss)*Z';
      V = V + eye(size(V))*theta_sig(sm,ss);
      
      Y_ms_n{tr,sm}{ss}   =mvnrnd(zeros(1,350),V,64)';
      Y_ms_n_mc{tr,sm}{ss}=mvnrnd(X*B,         V,64)';
      %       MScon{tr}=mvnrnd(zeros(1,7),G_sum(:,:,tr),64)';
    end
  end
  
  %   end
  
  
  % % OLD way
  %   % with noise added
  %   sigmodel=[1 5 10];
  %   % with condition means added
  %   meanadd=[ones(50,64); 1.2*ones(50,64); 1.4*ones(50,64); 1.6*ones(50,64); 1.4*ones(50,64); 1.2*ones(50,64); ones(50,64)];
  %   for tr=1:length(MScon)
  %     Y_ms{tr}=Z*MScon{tr};
  %     %     Y_cc=Z*Ccon;
  %     for sm=1:length(sigmodel)
  %       Y_ms_n{tr,sm}{ss}=Y_ms{tr}+sigmodel(sm)*demean(randn(350,64),2);
  %       %       Y_cc_n{ss}=Y_cc+sigmodel*demean(randn(350,64),2);
  %       Y_ms_n_mc{tr,sm}{ss}=Y_ms_n{tr,sm}{ss}+meanadd;
  %       %         Y_cc_n_mc{ss}=Y_cc_n{ss}+meanadd;
  %     end
  %   end
  
  
  
end % ss

save([edir 'sim_Ydata.mat'],'Y*','Z','theta*');
return

%%

load([edir 'pcm_models.mat']);  % from pcm_create_models.m
load([edir 'sim_Ydata.mat']);

if 0
  
  % Fit the models on the group level
  % partVec=ones(350,1);
  % Still playing around with how best to partition for cross-validation
  % partVec=repmat([1 2],[1 175])';
  partVec=repmat([1 2 3 4 5],[1 70])';
  
  % delete(gcp('nocreate'))
  % mypool=parpool;
  % addAttachedFiles(mypool, [edir 'sim_pcm_output_cv.mat'])
  
  % values{2}=partVec;
  % values{3}=M;
  % snames{1}='ms_mr';
  % snames{2}='partVec';
  % snames{3}='M';
  
  for sm=1:size(Y_ms_n,2)
    for tr=1:size(Y_ms_n,1)
      % multisensory, mean-removed
      [ms_mr{tr,sm}.Tgroup,ms_mr{tr,sm}.theta,ms_mr{tr,sm}.G_pred,ms_mr{tr,sm}.G_hat] = pcm_fitModelGroup(Y_ms_n{tr,sm},M,partVec,Z,'runEffect','fixed','fitScale',1);
      [ms_mr{tr,sm}.Tcross,ms_mr{tr,sm}.thetaCr,ms_mr{tr,sm}.G_predcv] = pcm_fitModelGroupCrossval(Y_ms_n{tr,sm},M,partVec,Z,'runEffect','fixed','groupFit',ms_mr{tr,sm}.theta,'fitScale',1);
      %     save([edir 'sim_pcm_output.mat'],'ms_mr'); % using partVec=ones
      %     updateAttachedFiles(mypool)
      %   save([edir 'sim_pcm_output_cv.mat'],'ms_mr','partVec'); % using real partVec
      %     values{1}=ms_mr;
      %     parfor_save([edir 'sim_pcm_output_cv2_rownorm.mat'],values,snames)
      save([edir 'sim_pcm_output_cv2_rownorm_newY.mat'],'ms_mr','partVec','M'); % using real partVec
    end
  end
  
  for sm=1:size(Y_ms_n,2)
    for tr=1:size(Y_ms_n,1)
      % multisensory, mean condition included
      [ms_mc{tr,sm}.Tgroup,ms_mc{tr,sm}.theta,ms_mc{tr,sm}.G_pred] = pcm_fitModelGroup(Y_ms_n_mc{tr,sm},M,partVec,Z,'runEffect','fixed','fitScale',1);
      [ms_mc{tr,sm}.Tcross,ms_mc{tr,sm}.thetaCr,ms_mc{tr,sm}.G_predcv] = pcm_fitModelGroupCrossval(Y_ms_n_mc{tr,sm},M,partVec,Z,'runEffect','fixed','groupFit',ms_mc{tr,sm}.theta,'fitScale',1);
      save([edir 'sim_pcm_output_cv.mat'],'ms_mr','ms_mc','partVec');
    end
  end
  
  load([edir 'sim_pcm_output_cv.mat'],'ms_mr');
  
  % % Cecere subtraction, mean-removed
  % [cs_mr.Tgroup,cs_mr.theta,cs_mr.G_pred] = pcm_fitModelGroup(Y_cc_n,M,partVec,Z,'runEffect','fixed','fitScale',1);
  % [cs_mr.Tcross,cs_mr.thetaCr,cs_mr.G_predcv] = pcm_fitModelGroupCrossval(Y_cc_n,M,partVec,Z,'runEffect','fixed','groupFit',cs_mr.theta,'fitScale',1);
  % % Cecere subtraction, mean condition included
  % [cs_mc.Tgroup,cs_mc.theta,cs_mc.G_pred] = pcm_fitModelGroup(Y_cc_n_mc,M,partVec,Z,'runEffect','fixed','fitScale',1);
  % [cs_mc.Tcross,cs_mc.thetaCr,cs_mc.G_predcv] = pcm_fitModelGroupCrossval(Y_cc_n_mc,M,partVec,Z,'runEffect','fixed','groupFit',cs_mc.theta,'fitScale',1);
  
  % save([edir 'sim_pcm_output.mat'],'ms_mr','ms_mc','cs_mr','cs_mc');
  %
  % ms_mr.T = pcm_plotModelLikelihood(ms_mr.Tcross,M,'upperceil',ms_mr.Tgroup.likelihood(:,2));
  % ms_mr.T = pcm_plotModelLikelihood(ms_mr.Tcross,M,'upperceil',ms_mr.Tgroup.likelihood(:,2),'normalize',1);
  % mean(ms_mr.T.likelihood_norm,1)
  
else
  % OR, see results from running on BlueBear cluster
  for sm=1:size(Y_ms_n,2)
    for tr=1:size(Y_ms_n,1)
      load([idir 'sim_pcm_output_tr' num2str(tr) '_sm' num2str(sm) '.mat']); % using real partVec
      ms_mr{tr,sm}=A;
    end
  end
  
end



set_pcmsim_colors
for tr=1:size(Y_ms_n,1)
  if numel(M)==17
    Nnull=17;
    Nceil=2;
    colors{tr}(Nnull)=[];
    colors{tr}(Nceil)=[];
    for sm=1:size(Y_ms_n,2)
      figure(tr)
      try
        subplot(4,2,2*(sm-1)+1);
        ms_mr{tr,sm}.T = pcm_plotModelLikelihood(ms_mr{tr,sm}.Tcross,M,'upperceil',ms_mr{tr,sm}.Tcross.likelihood(:,2),'normalize',1,'Nnull',Nnull,'Nceil',Nceil,'colors',colors{tr});
        if sm==1,title('CV');end
%         ylim([-2 2]);
      end
      try
        subplot(4,2,2*(sm-1)+2);
        ms_mr{tr,sm}.T = pcm_plotModelLikelihood(ms_mr{tr,sm}.Tgroup,M,'upperceil',ms_mr{tr,sm}.Tgroup.likelihood(:,2),'normalize',1,'Nnull',Nnull,'Nceil',Nceil,'colors',colors{tr});
        if sm==1,title('non-CV');end
%         ylim([-2 2]);
      end
    end
    print(tr,[fdir 'simpcm_results_tr' num2str(tr) '.png'],'-dpng');
  elseif numel(M)==16
    Nnull=1;
    Nceil=2;
    colors{tr}(Nnull)=[];
    colors{tr}(Nceill)=[];
    for sm=1:size(Y_ms_n,2)
      figure(tr)
      subplot(3,2,2*(sm-1)+1);
      ms_mr{tr,sm}.T = pcm_plotModelLikelihood(ms_mr{tr,sm}.Tcross,M,'upperceil',ms_mr{tr,sm}.Tcross.likelihood(:,2),'normalize',1,'Nnull',Nnull,'Nceil',Nceil,'colors',colors{tr});
      ylim([-2 2]);
      subplot(3,2,2*(sm-1)+2);
      ms_mr{tr,sm}.T = pcm_plotModelLikelihood(ms_mr{tr,sm}.Tgroup,M,'upperceil',ms_mr{tr,sm}.Tgroup.likelihood(:,2),'normalize',1,'Nnull',Nnull,'Nceil',Nceil,'colors',colors{tr});
      ylim([-2 2]);
    end
  end
end



