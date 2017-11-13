% Simulated data for PCM based on EEG study with 64 channels and 7
% audio-tactile asynchronies
%
% Johanna Zumer and Remi Gau, 2017

% Modify for your location
addpath('D:\Matlab\pcm_toolbox\')
addpath('D:\Matlab\NearestSymmetricPositiveDefinite')

ncond=7;
D.numPart = 5;
D.numVox  = 64; % 64 channel EEG
numSim = 22; % number of subjects

% Set .fitAlgorithm = 'minimize' or 'NR' inside this script
try
  load('pcm_models.mat');  % from create_sim1models
catch
  create_sim1models
end

theta_real(1,:) = [-inf -inf -inf -inf 0];
theta_real(2,:) = [0 0 -inf -inf -inf ];
theta_real(3,:) = [-inf -inf 0 -inf -inf ];
theta_real(4,:) = [-inf -inf -inf 0 -inf ];

theta_real(5,:) = [0 0 -inf -inf 0 ];
theta_real(6,:) = [-inf -inf 0 0 -inf ];
theta_real(7,:) = [0 0 0 -inf -inf ];

theta_real(8,:) = [0 0 0 -inf 0 ];
theta_real(9,:) = [2 2 0 -inf 0 ];
theta_real(10,:) = [0 0 0 0 0 ];

scale=abs(10+0.1*randn(1,numSim)); % vary this later
noise=abs(1 +0.1*randn(1,numSim));

% Model 16 is full model with all 5 components present; theta then weights them accordingly
for tr=1:size(theta_real,1)
  [Y{tr},part,conditions] = pcm_generateData(M{16},theta_real(tr,:)',D,numSim,scale,noise);
end

