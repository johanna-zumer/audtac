function sim2_pcm_run(tr,sm)
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

load([sdir 'sim2_featmodels.mat']);  % from pcm_create_models.m
load([sdir 'sim2_Ydata.mat']);

partVec=repmat([1 2 3 4 5],[1 70])';

disp(['Iteration running: TR ' num2str(tr) ' SM ' num2str(sm)])

Y=Y_ms_n_mc_mvnn{tr,sm};

[A.Tgroup,A.theta,A.G_pred,A.G_hat] = pcm_fitModelGroup(Y,M,partVec,Z,'runEffect','fixed','fitScale',1);
save([sdir 'sim2_pcm_output_grpMC_tr' num2str(tr) '_sm' num2str(sm) '.mat'],'A','partVec','M'); % using real partVec
[A.Tcross,A.thetaCr,A.G_predcv] = pcm_fitModelGroupCrossval(Y,M,partVec,Z,'runEffect','fixed','groupFit',A.theta,'fitScale',1);
save([sdir 'sim2_pcm_output_grpMC_tr' num2str(tr) '_sm' num2str(sm) '.mat'],'A','partVec','M'); % using real partVec

clear A
[A.Tgroup,A.theta_hat,A.G_pred]=pcm_fitModelIndivid(Y,M,partVec,Z,'runEffect','fixed');
save([sdir 'sim2_pcm_output_indMC_tr' num2str(tr) '_sm' num2str(sm) '.mat'],'A','partVec','M'); % using real partVec
[A.D,A.Tcross,A.thetaCV]=pcm_fitModelIndividCrossval(Y,M,partVec,Z,'runEffect','fixed');
save([sdir 'sim2_pcm_output_indMC_tr' num2str(tr) '_sm' num2str(sm) '.mat'],'A','partVec','M'); % using real partVec


end

