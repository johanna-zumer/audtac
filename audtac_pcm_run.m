function audtac_pcm_run(jj)
% function audtac_pcm_run(jj)
% for running on cluster
%
% jj = timestep

edir='/gpfs/bb/zumerj/nbu/audtac/eeg_data/';
disp('path to remove:')
pathd=which('ft_defaults');
if ~isempty(pathd)
  rmpath(fileparts(pathd));
end
addpath('/gpfs/bb/zumerj/nbu/fieldtrip_git/')
addpath('/gpfs/bb/zumerj/nbu/pcm_toolbox/')
addpath('/gpfs/bb/zumerj/nbu/audtac/mfiles/')

load([edir 'pcm_models.mat']);  % from pcm_create_models.m
load([edir 'pcm_datacat.mat']);

disp(['Iteration running: timestep ' num2str(jj) ])

for ii=1:length(datacat_meansubt_noisenorm)
  Y{ii}=datacat_meansubt_noisenorm{ii}(:,:,jj);
end
[A.Tgroup,A.theta,A.Gpred] = pcm_fitModelGroup(Y,M,partVec,designX,'runEffect','fixed','fitScale',1);
save([edir 'pcm_datafit_jj' num2str(jj) '.mat'],'A','partVec','M'); % using real partVec
[A.Tcross,A.thetaCr,A.G_predcv] = pcm_fitModelGroupCrossval(Y,M,partVec,designX,'runEffect','fixed','groupFit',A.theta,'fitScale',1);
save([edir 'pcm_datafit_jj' num2str(jj) '.mat'],'A','partVec','M'); % using real partVec


end

