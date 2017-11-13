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

load([sdir 'sim2_featmodels.mat']);  % from create_sim2_models.m
% load([edir 'pcm_datacat.mat']);
load([edir 'pcm_datacat_NumTime14.mat']);

disp(['Iteration running: timestep ' num2str(jj) ])

for ii=1:length(datacat_noisenorm)
%  Y{ii}=datacat_meansubt_noisenorm{ii}(:,:,jj);
  Y{ii}=datacat_noisenorm{ii}(:,:,jj);
end
[A.Tgroup,A.theta,A.Gpred] = pcm_fitModelGroup(Y,M,partVec2,designX,'runEffect','fixed','fitScale',1);
save([edir 'pcmF_datafit_p2_TS14_jj' num2str(jj) '.mat'],'A','partVec2','M'); % using real partVec
[A.Tcross,A.thetaCr,A.G_predcv] = pcm_fitModelGroupCrossval(Y,M,partVec2,designX,'runEffect','fixed','groupFit',A.theta,'fitScale',1);
save([edir 'pcmF_datafit_p2_TS14_jj' num2str(jj) '.mat'],'A','partVec2','M'); % using real partVec


end

