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
addpath('/gpfs/bb/zumerj/nbu/rsatoolbox/')
disp('IND fit')

try
% load([sdir 'sim2_featmodels.mat']);  % from create_sim2_models.m
% load([sdir 'sim3_featmodels.mat']);  % from create_sim3_models.m
% load([sdir 'sim4_featmodels.mat']);  % from create_sim4_models.m
 load([sdir 'sim4b_featmodels.mat']);  % from create_sim4_models.m
%load([sdir 'sim4bNR_featmodels.mat']);  % from create_sim4_models.m
%load([sdir 'sim5_featmodels.mat']);  % from create_sim5_models.m
% load([sdir 'sim6_featmodels.mat']);  % from create_sim5_models.m
% load([sdir 'sim7_featmodels.mat']);  % from create_sim5_models.m
catch
	create_sim4b_featAllModels;
%	create_sim4bNR_featAllModels;
%	create_sim5_featAllModels;
%	create_sim6_featAllModels;
%	create_sim7_featAllModels;
end
% load([edir 'pcm_datacat.mat']);
% load([edir 'pcm_datacat_NumTime14.mat']);
% load([edir 'pcm_datacat_nt14_nc5_np4.mat']);
% load([edir 'pcm_datacat_nt17_nc15_np3.mat']);
load([edir 'pcm_datacat_nt17_nc18_np2.mat']);

disp(['Iteration running: timestep ' num2str(jj) ])

%Mind{1}=M{1}; % just test first 'freechol' model
Mind=M; % full all models
for ii=1:length(datacat_chunk_nn)
%  Y{ii}=datacat_meansubt_noisenorm{ii}(:,:,jj);
  Y{ii}=datacat_chunk_nn{ii}(:,:,jj);
end

if 1
%[A{jj}.Tgroup,A{jj}.theta_hat,A{jj}.G_pred]=pcm_fitModelIndivid(Y,Mind,partVecch,Xchunk,'runEffect','fixed');
%save([edir 'pcmF_ind17_Msim4b_nc18np2_jj' num2str(jj) '.mat'],'A','partVecch','Mind'); 
[A{jj}.D,A{jj}.Tcross,A{jj}.thetaCV]=pcm_fitModelIndividCrossval(Y,Mind,partVecch,Xchunk,'runEffect','fixed');
save([edir 'pcmF_ind17_Msim4bCV_nc18np2_jj' num2str(jj) '.mat'],'A','partVecch','Mind'); 

else

[A{jj}.Tgroup,A{jj}.theta_hat,A{jj}.G_pred]=pcm_fitModelIndivid(Y,Mind,partVecch,Xchunk,'runEffect','fixed');
save([edir 'pcmF_ind17_Msim4bNR_nc18np2_jj' num2str(jj) '.mat'],'A','partVecch','Mind'); 
[A{jj}.D,A{jj}.Tcross,A{jj}.thetaCV]=pcm_fitModelIndividCrossval(Y,Mind,partVecch,Xchunk,'runEffect','fixed');
save([edir 'pcmF_ind17_Msim4bNR_nc18np2_jj' num2str(jj) '.mat'],'A','partVecch','Mind'); 
for mm=1:length(M)
	Mind{mm}.theta0=A{jj}.thetaCV{mm};
end
[A{jj}.Tgroup,A{jj}.theta_hat,A{jj}.G_pred]=pcm_fitModelIndivid(Y,Mind,partVecch,Xchunk,'runEffect','fixed');
save([edir 'pcmF_ind17_Msim4b_nc18np2_jj' num2str(jj) '.mat'],'A','partVecch','Mind'); 
[A{jj}.D,A{jj}.Tcross,A{jj}.thetaCV]=pcm_fitModelIndividCrossval(Y,Mind,partVecch,Xchunk,'runEffect','fixed');
save([edir 'pcmF_ind17_Msim4b_nc18np2_jj' num2str(jj) '.mat'],'A','partVecch','Mind'); 

end


end

