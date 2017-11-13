function audtac_pcmF_ind_run(jj)
% function audtac_pcmF_ind_run(jj)
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

%Mind{1}=M{1}; % just test first 'freechol' model
Mind=M; % full 17 models
for ii=1:length(datacat_noisenorm)
%  Y{ii}=datacat_meansubt_noisenorm{ii}(:,:,jj);
  Y{ii}=datacat_noisenorm{ii}(:,:,jj);
end

% partVec5=partVec;clear partVec
% numPart=2;
% for subcnt=1:length(datacat)
%   partVec2{subcnt}=[rem(1:size(datacat{subcnt},1),numPart)+1]';  % this is nonequal number of trials per partition & condition; okay?
% end
numPart=3;
for subcnt=1:length(datacat)
  partVec3{subcnt}=[rem(1:size(datacat{subcnt},1),numPart)+1]';  % this is nonequal number of trials per partition & condition; okay?
end

[A{jj}.Tgroup,A{jj}.theta_hat,A{jj}.G_pred]=pcm_fitModelIndivid(Y,Mind,partVec5,designX,'runEffect','fixed');
save([edir 'pcmF_ind5_jj' num2str(jj) '.mat'],'A','partVec5','Mind'); % using real partVec
[A{jj}.D,A{jj}.Tcross,A{jj}.thetaCV]=pcm_fitModelIndividCrossval(Y,Mind,partVec5,designX,'runEffect','fixed');
save([edir 'pcmF_ind5_jj' num2str(jj) '.mat'],'A','partVec5','Mind'); % using real partVec

[A{jj}.Tgroup,A{jj}.theta_hat,A{jj}.G_pred]=pcm_fitModelIndivid(Y,Mind,partVec2,designX,'runEffect','fixed');
save([edir 'pcmF_ind2_jj' num2str(jj) '.mat'],'A','partVec2','Mind'); % using real partVec
[A{jj}.D,A{jj}.Tcross,A{jj}.thetaCV]=pcm_fitModelIndividCrossval(Y,Mind,partVec2,designX,'runEffect','fixed');
save([edir 'pcmF_ind2_jj' num2str(jj) '.mat'],'A','partVec2','Mind'); % using real partVec

[A{jj}.Tgroup,A{jj}.theta_hat,A{jj}.G_pred]=pcm_fitModelIndivid(Y,Mind,partVec3,designX,'runEffect','fixed');
save([edir 'pcmF_ind3_jj' num2str(jj) '.mat'],'A','partVec3','Mind'); % using real partVec
[A{jj}.D,A{jj}.Tcross,A{jj}.thetaCV]=pcm_fitModelIndividCrossval(Y,Mind,partVec3,designX,'runEffect','fixed');
save([edir 'pcmF_ind3_jj' num2str(jj) '.mat'],'A','partVec3','Mind'); % using real partVec

end

