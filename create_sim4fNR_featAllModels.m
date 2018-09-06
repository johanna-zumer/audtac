% set up Models for PCM
% Johanna Zumer, 2017

% plotflag=1;
clear M F
sim4f_featAll_G
plotflag=0;

%% 
% Model 1: Free model as Noise ceiling (upper limit)
M{1}.type       = 'freechol'; 
M{1}.numCond    = ncond;
M{1}.name       = 'noiseceiling'; 
M{1}           = pcm_prepFreeModel(M{1}); 

% Model 2: Complete null model assuming no activity or covariance anywhere (lower limit)
M{2}.type       = 'feature';
M{2}.numGparams = 1;
M{2}.Ac = zeros(ncond,1);
M{2}.name       = 'null_zero'; 

%%

% Model 3: Eye: All conditions independent variance
M{3}.type       = 'feature';
M{3}.numGparams = 7;
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{3}.name       = 'I'; 

% Model 4: TIW a and b; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 9;
M{end}.Ac = zeros(7,9,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,8)  = F{15};
M{end}.Ac(:,9,9)  = F{16};
M{end}.name       = 'TI'; 

% Model 5: Asymmetry; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 11;
M{end}.Ac = zeros(ncond,11,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,8)  = F{11};
M{end}.Ac(:,9,9)  = F{12};
M{end}.Ac(:,10,10)  = F{13};
M{end}.Ac(:,11,11)  = F{14};
M{end}.name       = 'AI'; 

% Model 6: Symmetric Pairs; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 10;
M{end}.Ac = zeros(ncond,10,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,8)  = F{8};
M{end}.Ac(:,9,9)  = F{9};
M{end}.Ac(:,10,10)  = F{10};
M{end}.name       = 'SI'; 

% Model 7: TIW a and b;  Asymmetry; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 13;
M{end}.Ac = zeros(ncond,13,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,8)  = F{11};
M{end}.Ac(:,9,9)  = F{12};
M{end}.Ac(:,10,10)  = F{13};
M{end}.Ac(:,11,11)  = F{14};
M{end}.Ac(:,12,12)  = F{15};
M{end}.Ac(:,13,13)  = F{16};
M{end}.name       = 'TAI'; 

% Model 8: TIW a and b;  Symmetric Pairs; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 12;
M{end}.Ac = zeros(ncond,12,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,8)  = F{8};
M{end}.Ac(:,9,9)  = F{9};
M{end}.Ac(:,10,10)  = F{10};
M{end}.Ac(:,11,11)  = F{15};
M{end}.Ac(:,12,12)  = F{16};
M{end}.name       = 'TSI'; 

% Model 9: Asymmetry;  Symmetric Pairs; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 14;
M{end}.Ac = zeros(ncond,14,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,8)  = F{8};
M{end}.Ac(:,9,9)  = F{9};
M{end}.Ac(:,10,10)  = F{10};
M{end}.Ac(:,11,11)  = F{11};
M{end}.Ac(:,12,12)  = F{12};
M{end}.Ac(:,13,13)  = F{13};
M{end}.Ac(:,14,14)  = F{14};
M{end}.name       = 'ASI'; 

% Model 10: TIW a and b;  Asymmetry;  Symmetric Pairs; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 16;
M{end}.Ac = zeros(ncond,16,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,8)  = F{8};
M{end}.Ac(:,9,9)  = F{9};
M{end}.Ac(:,10,10)  = F{10};
M{end}.Ac(:,11,11)  = F{11};
M{end}.Ac(:,12,12)  = F{12};
M{end}.Ac(:,13,13)  = F{13};
M{end}.Ac(:,14,14)  = F{14};
M{end}.Ac(:,15,15)  = F{15};
M{end}.Ac(:,16,16)  = F{16};
M{end}.name       = 'TASI'; 

%% With 17-20 instead of 11-14

% Model 11: Asymmetry; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 11;
M{end}.Ac = zeros(ncond,11,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,8)  = F{17};
M{end}.Ac(:,9,9)  = F{18};
M{end}.Ac(:,10,10)  = F{19};
M{end}.Ac(:,11,11)  = F{20};
M{end}.name       = 'AI'; 

% Model 12: TIW a and b;  Asymmetry; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 13;
M{end}.Ac = zeros(ncond,13,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,8)  = F{17};
M{end}.Ac(:,9,9)  = F{18};
M{end}.Ac(:,10,10)  = F{19};
M{end}.Ac(:,11,11)  = F{20};
M{end}.Ac(:,12,12)  = F{15};
M{end}.Ac(:,13,13)  = F{16};
M{end}.name       = 'TAI'; 

% Model 13: Asymmetry;  Symmetric Pairs; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 14;
M{end}.Ac = zeros(ncond,14,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,8)  = F{8};
M{end}.Ac(:,9,9)  = F{9};
M{end}.Ac(:,10,10)  = F{10};
M{end}.Ac(:,11,11)  = F{17};
M{end}.Ac(:,12,12)  = F{18};
M{end}.Ac(:,13,13)  = F{19};
M{end}.Ac(:,14,14)  = F{20};
M{end}.name       = 'ASI'; 

% Model 14: TIW a and b;  Asymmetry;  Symmetric Pairs; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 16;
M{end}.Ac = zeros(ncond,16,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,8)  = F{8};
M{end}.Ac(:,9,9)  = F{9};
M{end}.Ac(:,10,10)  = F{10};
M{end}.Ac(:,11,11)  = F{17};
M{end}.Ac(:,12,12)  = F{18};
M{end}.Ac(:,13,13)  = F{19};
M{end}.Ac(:,14,14)  = F{20};
M{end}.Ac(:,15,15)  = F{15};
M{end}.Ac(:,16,16)  = F{16};
M{end}.name       = 'TASI'; 

%% With 12,14,18,20 instead of 11-14

% Model 15: Asymmetry; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 11;
M{end}.Ac = zeros(ncond,11,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,8)  = F{12};
M{end}.Ac(:,9,9)  = F{14};
M{end}.Ac(:,10,10)  = F{18};
M{end}.Ac(:,11,11)  = F{20};
M{end}.name       = 'AI'; 

% Model 16: TIW a and b;  Asymmetry; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 13;
M{end}.Ac = zeros(ncond,13,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,8)  = F{12};
M{end}.Ac(:,9,9)  = F{14};
M{end}.Ac(:,10,10)  = F{18};
M{end}.Ac(:,11,11)  = F{20};
M{end}.Ac(:,12,12)  = F{15};
M{end}.Ac(:,13,13)  = F{16};
M{end}.name       = 'TAI'; 

% Model 17: Asymmetry;  Symmetric Pairs; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 14;
M{end}.Ac = zeros(ncond,14,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,8)  = F{8};
M{end}.Ac(:,9,9)  = F{9};
M{end}.Ac(:,10,10)  = F{10};
M{end}.Ac(:,11,11)  = F{12};
M{end}.Ac(:,12,12)  = F{14};
M{end}.Ac(:,13,13)  = F{18};
M{end}.Ac(:,14,14)  = F{20};
M{end}.name       = 'ASI'; 

% Model 18: TIW a and b;  Asymmetry;  Symmetric Pairs; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 16;
M{end}.Ac = zeros(ncond,16,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,8)  = F{8};
M{end}.Ac(:,9,9)  = F{9};
M{end}.Ac(:,10,10)  = F{10};
M{end}.Ac(:,11,11)  = F{12};
M{end}.Ac(:,12,12)  = F{14};
M{end}.Ac(:,13,13)  = F{18};
M{end}.Ac(:,14,14)  = F{20};
M{end}.Ac(:,15,15)  = F{15};
M{end}.Ac(:,16,16)  = F{16};
M{end}.name       = 'TASI'; 

%% With 11,13,17,19 instead of 11-14

% Model 19: Asymmetry; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 11;
M{end}.Ac = zeros(ncond,11,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,8)  = F{17};
M{end}.Ac(:,9,9)  = F{11};
M{end}.Ac(:,10,10)  = F{19};
M{end}.Ac(:,11,11)  = F{13};
M{end}.name       = 'AI'; 

% Model 20: TIW a and b;  Asymmetry; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 13;
M{end}.Ac = zeros(ncond,13,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,8)  = F{17};
M{end}.Ac(:,9,9)  = F{11};
M{end}.Ac(:,10,10)  = F{19};
M{end}.Ac(:,11,11)  = F{13};
M{end}.Ac(:,12,12)  = F{15};
M{end}.Ac(:,13,13)  = F{16};
M{end}.name       = 'TAI'; 

% Model 21: Asymmetry;  Symmetric Pairs; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 14;
M{end}.Ac = zeros(ncond,14,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,8)  = F{8};
M{end}.Ac(:,9,9)  = F{9};
M{end}.Ac(:,10,10)  = F{10};
M{end}.Ac(:,11,11)  = F{17};
M{end}.Ac(:,12,12)  = F{11};
M{end}.Ac(:,13,13)  = F{19};
M{end}.Ac(:,14,14)  = F{13};
M{end}.name       = 'ASI'; 

% Model 22: TIW a and b;  Asymmetry;  Symmetric Pairs; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 16;
M{end}.Ac = zeros(ncond,16,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,8)  = F{8};
M{end}.Ac(:,9,9)  = F{9};
M{end}.Ac(:,10,10)  = F{10};
M{end}.Ac(:,11,11)  = F{17};
M{end}.Ac(:,12,12)  = F{11};
M{end}.Ac(:,13,13)  = F{19};
M{end}.Ac(:,14,14)  = F{13};
M{end}.Ac(:,15,15)  = F{15};
M{end}.Ac(:,16,16)  = F{16};
M{end}.name       = 'TASI'; 

%%

% Model 23: TIW a and b;  Asymmetry;  Symmetric Pairs; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 19;
M{end}.Ac = zeros(ncond,19,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,8)  = F{8};
M{end}.Ac(:,9,9)  = F{9};
M{end}.Ac(:,10,10)  = F{10};
M{end}.Ac(:,11,11)  = F{11};
M{end}.Ac(:,12,12)  = F{12};
M{end}.Ac(:,13,13)  = F{13};
M{end}.Ac(:,14,14)  = F{14};
M{end}.Ac(:,15,15)  = F{15};
M{end}.Ac(:,16,16)  = F{16};
M{end}.Ac(:,17,17)  = F{21};
M{end}.Ac(:,18,18)  = F{22};
M{end}.Ac(:,19,19)  = F{23};
M{end}.name       = 'TASI'; 

% Model 24: TIW a and b;  Asymmetry;  Symmetric Pairs; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 19;
M{end}.Ac = zeros(ncond,19,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,8)  = F{8};
M{end}.Ac(:,9,9)  = F{9};
M{end}.Ac(:,10,10)  = F{10};
M{end}.Ac(:,11,11)  = F{17};
M{end}.Ac(:,12,12)  = F{18};
M{end}.Ac(:,13,13)  = F{19};
M{end}.Ac(:,14,14)  = F{20};
M{end}.Ac(:,15,15)  = F{15};
M{end}.Ac(:,16,16)  = F{16};
M{end}.Ac(:,17,17)  = F{21};
M{end}.Ac(:,18,18)  = F{22};
M{end}.Ac(:,19,19)  = F{23};
M{end}.name       = 'TASI'; 

% Model 25: TIW a and b;  Asymmetry;  Symmetric Pairs; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 19;
M{end}.Ac = zeros(ncond,19,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,8)  = F{8};
M{end}.Ac(:,9,9)  = F{9};
M{end}.Ac(:,10,10)  = F{10};
M{end}.Ac(:,11,11)  = F{12};
M{end}.Ac(:,12,12)  = F{14};
M{end}.Ac(:,13,13)  = F{18};
M{end}.Ac(:,14,14)  = F{20};
M{end}.Ac(:,15,15)  = F{15};
M{end}.Ac(:,16,16)  = F{16};
M{end}.Ac(:,17,17)  = F{21};
M{end}.Ac(:,18,18)  = F{22};
M{end}.Ac(:,19,19)  = F{23};
M{end}.name       = 'TASI'; 

% Model 26: TIW a and b;  Asymmetry;  Symmetric Pairs; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 17;
M{end}.Ac = zeros(ncond,17,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,8)  = F{8};
M{end}.Ac(:,9,9)  = F{9};
M{end}.Ac(:,10,10)  = F{10};
M{end}.Ac(:,11,11)  = F{12};
M{end}.Ac(:,12,12)  = F{14};
M{end}.Ac(:,13,13)  = F{18};
M{end}.Ac(:,14,14)  = F{20};
M{end}.Ac(:,15,15)  = F{15};
M{end}.Ac(:,16,16)  = F{16};
M{end}.Ac(:,17,17)  = F{21};
M{end}.name       = 'TASI'; 

% Model 27: TIW a and b;  Asymmetry;  Symmetric Pairs; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 19;
M{end}.Ac = zeros(ncond,19,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,8)  = F{8};
M{end}.Ac(:,9,9)  = F{9};
M{end}.Ac(:,10,10)  = F{10};
M{end}.Ac(:,11,11)  = F{17};
M{end}.Ac(:,12,12)  = F{11};
M{end}.Ac(:,13,13)  = F{19};
M{end}.Ac(:,14,14)  = F{13};
M{end}.Ac(:,15,15)  = F{15};
M{end}.Ac(:,16,16)  = F{16};
M{end}.Ac(:,17,17)  = F{21};
M{end}.Ac(:,18,18)  = F{22};
M{end}.Ac(:,19,19)  = F{23};
M{end}.name       = 'TASI'; 

% Model 28: TIW a and b;  Asymmetry;  Symmetric Pairs; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 18;
M{end}.Ac = zeros(ncond,18,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,8)  = F{8};
M{end}.Ac(:,9,9)  = F{9};
M{end}.Ac(:,10,10)  = F{10};
M{end}.Ac(:,11,11)  = F{17};
M{end}.Ac(:,12,12)  = F{11};
M{end}.Ac(:,13,13)  = F{19};
M{end}.Ac(:,14,14)  = F{13};
M{end}.Ac(:,15,15)  = F{15};
M{end}.Ac(:,16,16)  = F{16};
M{end}.Ac(:,17,17)  = F{22};
M{end}.Ac(:,18,18)  = F{23};
M{end}.name       = 'TASI'; 



%% 
for mm=1:numel(M)
  M{mm}.fitAlgorithm='NR'; % NR or minimize
%   M{mm}.fitAlgorithm='minimize'; % NR or minimize
end

if plotflag
  figure;
  for mm=3:18,subplot(2,16,mm-2);imagesc(sum(M{mm}.Ac,3));end
  for mm=3:18,
    theta0 = ones(size(M{mm}.Ac,3),1);
    Gmodel{mm}=pcm_calculateG(M{mm},theta0);
    subplot(2,16,mm-2+16);imagesc(Gmodel{mm});
  end
  
  
%   figure;
%   for mm=18:32,subplot(2,15,mm-17);imagesc(sum(M{mm}.Ac,3));end
%   for mm=18:32,
%     theta0 = ones(size(M{mm}.Ac,3),1);
%     Gmodel{mm}=pcm_calculateG(M{mm},theta0);
%     subplot(2,15,mm-17+15);imagesc(Gmodel{mm});
%   end
end

save([edir 'sim4fNR_featmodels.mat'],'M');


