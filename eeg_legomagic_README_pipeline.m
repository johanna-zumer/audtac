% EEG awake & sleep analysis pipeline

% code to call stimulation paradigm
call_lego_audtac % calls lego_audtac

% rt_legomagic_proc_part1 % behavioural RT results
% above made obsolete by update to 'concat_stim_files' below for including response-only_no-EEG session

eeg_legomagic_montage % converged on optimal montage for staging & feature viewing

% do sleep staging within these (or within ft_databrowser?)
spm_preproc_renamed 
spm_preproc2
spm_preproc_withTom % for final stage concensus results
% spm_prepareSleepStaging % used only for test data from Jo (not for main expt)

% analysis of sleep stage results
sleep_stage_results % generate iikeep.mat
% reran sleep_stage_results after getting CRC_*withTom.mat

% artifact & sleep feature detection
eeg_legomagic_eog
eeg_legomagic_muscle_artifact
eeg_legomagic_Kcomplex % calls crc_SWS_detect_jz;   DONE regarding using final sleep stages
eeg_legomagic_spindle % calls crc_sp_detect_jz;  DONE regarding using final sleep stages

% moderate both automatically and manually the features detected
eeg_legomagic_viewFeatures % view stages & features in ft_databrowser
% DONE regarding using final sleep stages
%%%% Do manual feature modification....(?) (Tom?)

% prepares stimulus information
concat_stim_files  % updated on 20 March, 2015 to process response-no-EEG session files as well

% epoching and assigning features above to each trial
% Don't need to call this independently (it will be called by erp_stats and freqanalsysis.
% eeg_legomagic_epoching2  % calls ft_trialfun_general_sleep
% call above as eeg_legomagic_epoching2(ii,sleep,1,saveflag) for erp & freq

% simulation_ERP_peak_height_jitter_SNR % just for test of effect of jitter on ERP

%% then independent paths for ERP and TFR

% % OLD for recreating awake-alone results in october 2014 (although
% unforutnately slightly modified to match newer way of doing it). top
% section is modified but bottom (stats and figures) is not modified as
% much.
% eeg_legomagic_erp_stats % calls eeg_legomagic_trialSelection1_wakeSleep;
% eeg_legomagic_erp_stats1 % calls eeg_legomagic_trialSelection2_wakeSleep
eeg_legomagic_erp_stats2_sepTacAud % calls eeg_legomagic_trialSelection2_wakeSleep_sepTacAud

eeg_legomagic_featurestats_sepTacAud % requires long trials but not save out to disk.
% above had called eeg_legomagic_epoching2_featurestats, but this is now
% deprecated and should call eeg_legomagic_epoching2(ii,sleep,0,saveflag).

% eeg_freqanalysis_sensor % calls eeg_legomagic_trialSelection_freqwide
eeg_freqanalysis_sensor1 % calls eeg_legomagic_trialSelection2_wakeSleep_sepTacAud
eeg_freqanalysis_nonlinearBootstrap_sensor % calls eeg_legomagic_trialSelection1

% RT analysis as well as RT-EEG correlations
eeg_legomagic_brainBehaviourCorrelations

% % Source analysis
% eeg_mriheadmodels
% eeg_erp_sourceloc
% eeg_tfr_sourceloc

%% obsolete
% eeg_legomagic_preproc2 % with ICA
% eeg_legomagic_preprocSleep2 % with ICA

