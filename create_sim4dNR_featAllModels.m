% set up Models for PCM
% Johanna Zumer, 2017

% plotflag=1;
clear M F
sim4d_featAll_G
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

%%   With extra columm of ones for all conditions added

% Model 11: Eye: All conditions independent variance
M{end+1}.type       = 'feature';
M{end}.numGparams = 10;
M{end}.Ac = zeros(7,10,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
for ff=1:3
  M{end}.Ac(:,ff+7,ff+7)  = F{ff+16};
end
M{end}.name       = 'IO'; 

% Model 12: TIW a and b; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 12;
M{end}.Ac = zeros(7,12,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,8)  = F{15};
M{end}.Ac(:,9,9)  = F{16};
for ff=1:3
  M{end}.Ac(:,ff+9,ff+9)  = F{ff+16};
end
M{end}.name       = 'TIW1a1b_eye_ones'; 

% Model 13: Asymmetry; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 14;
M{end}.Ac = zeros(ncond,14,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,8)  = F{11};
M{end}.Ac(:,9,9)  = F{12};
M{end}.Ac(:,10,10)  = F{13};
M{end}.Ac(:,11,11)  = F{14};
M{end}.name       = 'Asym2_eye_ones'; 

% Model 14: Symmetric Pairs; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 13;
M{end}.Ac = zeros(ncond,13,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,8)  = F{8};
M{end}.Ac(:,9,9)  = F{9};
M{end}.Ac(:,10,10)  = F{10};
for ff=1:3
  M{end}.Ac(:,ff+10,ff+10)  = F{ff+16};
end
M{end}.name       = 'SymPair3_eye_ones'; 

% Model 15: TIW a and b;  Asymmetry; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 16;
M{end}.Ac = zeros(ncond,16,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,8)  = F{11};
M{end}.Ac(:,9,9)  = F{12};
M{end}.Ac(:,10,10)  = F{13};
M{end}.Ac(:,11,11)  = F{14};
M{end}.Ac(:,12,12)  = F{15};
M{end}.Ac(:,13,13)  = F{16};
for ff=1:3
  M{end}.Ac(:,ff+13,ff+13)  = F{ff+16};
end
M{end}.name       = 'TIW1a1b_Asym2_eye_ones'; 

% Model 16: TIW a and b;  Symmetric Pairs; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 15;
M{end}.Ac = zeros(ncond,15,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,8)  = F{8};
M{end}.Ac(:,9,9)  = F{9};
M{end}.Ac(:,10,10)  = F{10};
M{end}.Ac(:,11,11)  = F{15};
M{end}.Ac(:,12,12)  = F{16};
for ff=1:3
  M{end}.Ac(:,ff+12,ff+12)  = F{ff+16};
end
M{end}.name       = 'TIW1a1b_SymPair3_eye_ones'; 

% Model 17: Asymmetry;  Symmetric Pairs; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 17;
M{end}.Ac = zeros(ncond,17,M{end}.numGparams);
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
for ff=1:3
  M{end}.Ac(:,ff+14,ff+14)  = F{ff+16};
end
M{end}.name       = 'Asym2_SymPair3_eye_ones'; 

% Model 18: TIW a and b;  Asymmetry;  Symmetric Pairs; and eye
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
M{end}.Ac(:,15,15)  = F{15};  % TIW
M{end}.Ac(:,16,16)  = F{16};  % TIW
for ff=1:3
  M{end}.Ac(:,ff+16,ff+16)  = F{ff+16};
end
M{end}.name       = 'TIW1a1b_Asym2_SymPair3_eye_ones'; 





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
  
end

save([edir 'sim4dNR_featmodels.mat'],'M');


