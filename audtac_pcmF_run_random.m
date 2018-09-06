function audtac_pcmF_run(jj)
% function audtac_pcm_run(jj)
% for running on cluster
%
% jj = timestep

edir='/gpfs/bb/zumerj/nbu/audtac/eeg_data/';
sdir='/gpfs/bb/zumerj/nbu/audtac/sim_results/';
disp('path to remove:')
pathd=which('ft_defaults');
if ~isempty(pathd)
  rmpath(fileparts(pathd));
end
addpath('/gpfs/bb/zumerj/nbu/fieldtrip_git/')
addpath('/gpfs/bb/zumerj/nbu/pcm_toolbox/')
addpath('/gpfs/bb/zumerj/nbu/audtac/mfiles/')
addpath('/gpfs/bb/zumerj/nbu/rsatoolbox/')
disp('GROUP fit')


try
% load([sdir 'sim2_featmodels.mat']);  % from create_sim2_models.m
% load([sdir 'sim3_featmodels.mat']);  % from create_sim3_models.m
% load([sdir 'sim4_featmodels.mat']);  % from create_sim4_models.m
%  load([sdir 'sim4b_featmodels.mat']);  % from create_sim4_models.m
% load([sdir 'sim4bNR_featmodels.mat']);  % from create_sim4_models.m
%  load([sdir 'sim4dNR_featmodels.mat']);  % from create_sim4_models.m
%  load([sdir 'sim4eNR_featmodels.mat']);  % from create_sim4_models.m
  %load([sdir 'sim4fNR_featmodels.mat']);  % from create_sim4_models.m
  load([sdir 'sim4gNR_featmodels.mat']);  % from create_sim4_models.m
% load([sdir 'sim5_featmodels.mat']);  % from create_sim5_models.m
% load([sdir 'sim6_featmodels.mat']);  % 
% load([sdir 'sim7_featmodels.mat']);  % 
% load([sdir 'sim7NR_featmodels.mat']);  % 
% load([sdir 'sim7bNR_featmodels.mat']);  % 
catch
%	create_sim4b_featAllModels;
	%create_sim4bNR_featAllModels;
%	create_sim4dNR_featAllModels;
%	create_sim4eNR_featAllModels;
%	create_sim4fNR_featAllModels;
	create_sim4gNR_featAllModels;
%	create_sim5_featAllModels;
%	create_sim6_featAllModels;
%	create_sim7_featAllModels;
%	create_sim7NR_featAllModels;
%	create_sim7bNR_featAllModels;
end
% load([edir 'pcm_datacat.mat']);
% load([edir 'pcm_datacat_NumTime14.mat']);
% load([edir 'pcm_datacat_nt14_nc5_np4.mat']);
% load([edir 'pcm_datacat_nt17_nc15_np3.mat']);
 load([edir 'pcm_datacat_nt17_nc18_np2.mat']);
%  load([edir 'pcm_datacat_nt17_nc9_np4.mat']);

disp(['Iteration running: timestep ' num2str(jj) ])

for ii=1:length(datacat_chunk_nn)
%  Y{ii}=datacat_meansubt_noisenorm{ii}(:,:,jj);
  Y{ii}=datacat_chunk_nn{ii}(:,:,jj);
end

[A.Tgroup,A.theta,A.Gpred] = pcm_fitModelGroup(Y,M,partVecch,Xchunk,'runEffect','random','fitScale',1);
save([edir 'pcmF_datafit_TS17_Msim4gNRrand_nc18np2_jj' num2str(jj) '.mat'],'A','partVecch','M'); % using real partVec
[A.Tcross,A.thetaCr,A.G_predcv] = pcm_fitModelGroupCrossval(Y,M,partVecch,Xchunk,'runEffect','random','groupFit',A.theta,'fitScale',1);
save([edir 'pcmF_datafit_TS17_Msim4gNRrand_nc18np2_jj' num2str(jj) '.mat'],'A','partVecch','M'); % using real partVec


end

