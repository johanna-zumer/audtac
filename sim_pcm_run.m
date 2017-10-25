function sim_pcm_run(tr,sm)
% function sim_pcm_run(tr,sm)
% for running on cluster

sdir='/gpfs/bb/zumerj/nbu/audtac/sim_results/';
disp('path to remove:')
pathd=which('ft_defaults');
if ~isempty(pathd)
  rmpath(fileparts(pathd));
end
addpath('/gpfs/bb/zumerj/nbu/fieldtrip_git/')
addpath('/gpfs/bb/zumerj/nbu/pcm_toolbox/')
addpath('/gpfs/bb/zumerj/nbu/audtac/mfiles/')

load([sdir 'pcm_models.mat']);  % from pcm_create_models.m
load([sdir 'sim_Ydata.mat']);

partVec=repmat([1 2 3 4 5],[1 70])';

disp(['Iteration running: TR ' num2str(tr) ' SM ' num2str(sm)])

[A.Tgroup,A.theta,A.G_pred,A.G_hat] = pcm_fitModelGroup(Y_ms_n{tr,sm},M,partVec,Z,'runEffect','fixed','fitScale',1);
save([sdir 'sim_pcm_output_tr' num2str(tr) '_sm' num2str(sm) '.mat'],'A','partVec','M'); % using real partVec
[A.Tcross,A.thetaCr,A.G_predcv] = pcm_fitModelGroupCrossval(Y_ms_n{tr,sm},M,partVec,Z,'runEffect','fixed','groupFit',A.theta,'fitScale',1);
save([sdir 'sim_pcm_output_tr' num2str(tr) '_sm' num2str(sm) '.mat'],'A','partVec','M'); % using real partVec

end

