% set up Models for PCM
% Johanna Zumer, 2017

% plotflag=1;
clear M F
sim5_featAll_G

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
M{3}.Ac = zeros(7,7,M{3}.numGparams);
for ff=1:7
  M{3}.Ac(:,:,ff)  = F{ff};
end
M{3}.name       = 'eye';

% Model 4: Asymmetry
M{4}.type       = 'feature';
M{4}.numGparams = 12;
M{4}.Ac = zeros(7,12,M{4}.numGparams);
for ff=1:12
  M{4}.Ac(:,ff,ff)  = F{ff+9}; % 10-21
end
M{4}.name       = 'Asym_OffTIW'; 

% Model 5: On TIW
M{5}.type       = 'feature';
M{5}.numGparams = 2;
M{5}.Ac = zeros(7,2,M{5}.numGparams);
M{5}.Ac(:,1,1)  = F{8};
M{5}.Ac(:,2,2)  = F{9};
M{5}.name       = 'OnTIW'; 

% Model 6: 
M{6}.type       = 'feature';
M{6}.numGparams = 14;
M{6}.Ac = zeros(7,14,M{6}.numGparams);
for ff=1:14
  M{6}.Ac(:,ff,ff)  = F{ff+7};
end
M{6}.name       = 'OnTIW_Asym_OffTIW'; 

M{7}.type       = 'feature';
M{7}.numGparams = 8;
M{7}.Ac = zeros(7,8,M{7}.numGparams);
M{7}.Ac(:,1,1)  = F{8};
M{7}.Ac(:,2,2)  = F{9};
M{7}.Ac(:,3,3)  = F{10};
M{7}.Ac(:,4,4)  = F{11};
M{7}.Ac(:,5,5)  = F{12};
M{7}.Ac(:,6,6)  = F{13};
M{7}.Ac(:,7,7)  = F{14};
M{7}.Ac(:,8,8)  = F{15};
M{7}.name       = 'Asym_OnTIW'; 

M{8}.type       = 'feature';
M{8}.numGparams = 8;
M{8}.Ac = zeros(7,8,M{8}.numGparams);
M{8}.Ac(:,1,1)  = F{8};
M{8}.Ac(:,2,2)  = F{9};
M{8}.Ac(:,3,3)  = F{16};
M{8}.Ac(:,4,4)  = F{17};
M{8}.Ac(:,5,5)  = F{18};
M{8}.Ac(:,6,6)  = F{19};
M{8}.Ac(:,7,7)  = F{20};
M{8}.Ac(:,8,8)  = F{21};
M{8}.name       = 'OnTIW_OffTIW'; 

M{9}.type       = 'feature';
M{9}.numGparams = 5;
M{9}.Ac = zeros(7,5,M{9}.numGparams);
M{9}.Ac(:,1,1)  = F{8};
M{9}.Ac(:,2,2)  = F{10};
M{9}.Ac(:,3,3)  = F{11};
M{9}.Ac(:,4,4)  = F{12};
M{9}.Ac(:,5,5)  = F{13};
M{9}.name       = 'All3'; 

M{10}.type       = 'feature';
M{10}.numGparams = 4;
M{10}.Ac = zeros(7,4,M{10}.numGparams);
M{10}.Ac(:,1,1)  = F{14};
M{10}.Ac(:,2,2)  = F{15};
M{10}.Ac(:,3,3)  = F{16};
M{10}.Ac(:,4,4)  = F{17};
M{10}.name       = 'All4'; 

% Model 11
M{end+1}.type       = 'feature';
M{end}.numGparams = 3;
M{end}.Ac = zeros(7,3,M{end}.numGparams);
M{end}.Ac(:,1,1)  = F{9};
M{end}.Ac(:,2,2)  = F{18};
M{end}.Ac(:,3,3)  = F{19};
M{end}.name       = 'All5'; 

% Model 12: 
M{end+1}.type       = 'feature';
M{end}.numGparams = 2;
M{end}.Ac = zeros(ncond,2,M{end}.numGparams);
M{end}.Ac(:,1,1)  = F{20};
M{end}.Ac(:,2,2)  = F{21};
M{end}.name       = 'All6'; 

% Model 13: 
M{end+1}.type       = 'feature';
M{end}.numGparams = 9;
M{end}.Ac = zeros(ncond,9,M{end}.numGparams);
M{end}.Ac(:,1,1)  = F{8};
M{end}.Ac(:,2,2)  = F{10};
M{end}.Ac(:,3,3)  = F{11};
M{end}.Ac(:,4,4)  = F{12};
M{end}.Ac(:,5,5)  = F{13};
M{end}.Ac(:,6,6)  = F{14};
M{end}.Ac(:,7,7)  = F{15};
M{end}.Ac(:,8,8)  = F{16};
M{end}.Ac(:,9,9)  = F{17};
M{end}.name       = 'all3_all4'; 

% Model 14: 
M{end+1}.type       = 'feature';
M{end}.numGparams = 7;
M{end}.Ac = zeros(ncond,7,M{end}.numGparams);
M{end}.Ac(:,1,1)  = F{9};
M{end}.Ac(:,2,2)  = F{18};
M{end}.Ac(:,3,3)  = F{19};
M{end}.Ac(:,4,4)  = F{14};
M{end}.Ac(:,5,5)  = F{15};
M{end}.Ac(:,6,6)  = F{16};
M{end}.Ac(:,7,7)  = F{17};
M{end}.name       = 'all4_all5'; 

% Model 15: 
M{end+1}.type       = 'feature';
M{end}.numGparams = 5;
M{end}.Ac = zeros(ncond,5,M{end}.numGparams);
M{end}.Ac(:,1,1)  = F{9};
M{end}.Ac(:,2,2)  = F{18};
M{end}.Ac(:,3,3)  = F{19};
M{end}.Ac(:,4,4)  = F{20};
M{end}.Ac(:,5,5)  = F{21};
M{end}.name       = 'all5_all6'; 

% Model 16: 
M{end+1}.type       = 'feature';
M{end}.numGparams = 12;
M{end}.Ac = zeros(ncond,12,M{end}.numGparams);
M{end}.Ac(:,1,1)  = F{8};
M{end}.Ac(:,2,2)  = F{10};
M{end}.Ac(:,3,3)  = F{11};
M{end}.Ac(:,4,4)  = F{12};
M{end}.Ac(:,5,5)  = F{13};
M{end}.Ac(:,6,6)  = F{14};
M{end}.Ac(:,7,7)  = F{15};
M{end}.Ac(:,8,8)  = F{16};
M{end}.Ac(:,9,9)  = F{17};
M{end}.Ac(:,10,10)  = F{9};
M{end}.Ac(:,11,11)  = F{18};
M{end}.Ac(:,12,12)  = F{19};
M{end}.name       = 'all3_all4_all5'; 

% Model 17: 
M{end+1}.type       = 'feature';
M{end}.numGparams = 9;
M{end}.Ac = zeros(ncond,9,M{end}.numGparams);
M{end}.Ac(:,1,1)  = F{9};
M{end}.Ac(:,2,2)  = F{18};
M{end}.Ac(:,3,3)  = F{19};
M{end}.Ac(:,4,4)  = F{14};
M{end}.Ac(:,5,5)  = F{15};
M{end}.Ac(:,6,6)  = F{16};
M{end}.Ac(:,7,7)  = F{17};
M{end}.Ac(:,8,8)  = F{20};
M{end}.Ac(:,9,9)  = F{21};
M{end}.name       = 'all4_all5_all6'; 


%%   With Identity added to all models above (4-17) are now 

% Model 18: 
M{end+1}.type       = 'feature';
M{end}.numGparams = 19;
M{end}.Ac = zeros(7,19,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
for ff=1:12
  M{end}.Ac(:,ff+7,ff+7)  = F{ff+9}; % 10-21
end
M{end}.name       = 'Asym_OffTIW_eye';


% Model 19: On TIW
M{end+1}.type       = 'feature';
M{end}.numGparams = 9;
M{end}.Ac = zeros(7,9,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,9)  = F{8};
M{end}.Ac(:,8,9)  = F{9};
M{end}.name       = 'OnTIW_eye'; 

% Model 20
M{end+1}.type       = 'feature';
M{end}.numGparams = 21;
M{end}.Ac = zeros(7,21,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
for ff=1:14
  M{end}.Ac(:,ff+7,ff+7)  = F{ff+7};
end
M{end}.name       = 'Asym_OnTIW_OffTIW_eye'; 

% Model 21
M{end+1}.type       = 'feature';
M{end}.numGparams = 15;
M{end}.Ac = zeros(7,15,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
for ff=1:8
  M{end}.Ac(:,ff+7,ff+7)  = F{ff+7};
end
M{end}.name       = 'Asym_OnTIW_eye'; 

% Model 22
M{end+1}.type       = 'feature';
M{end}.numGparams = 15;
M{end}.Ac = zeros(7,15,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,8)  = F{8};
M{end}.Ac(:,9,9)  = F{9};
M{end}.Ac(:,10,10)  = F{16};
M{end}.Ac(:,11,11)  = F{17};
M{end}.Ac(:,12,12)  = F{18};
M{end}.Ac(:,13,13)  = F{19};
M{end}.Ac(:,14,14)  = F{20};
M{end}.Ac(:,15,15)  = F{21};
M{end}.name       = 'OnTIW_OffTIW_eye'; 

% Model 23
M{end+1}.type       = 'feature';
M{end}.numGparams = 12;
M{end}.Ac = zeros(7,12,M{end}.numGparams);
M{end}.Ac(:,1,1)  = F{8};
M{end}.Ac(:,2,2)  = F{10};
M{end}.Ac(:,3,3)  = F{11};
M{end}.Ac(:,4,4)  = F{12};
M{end}.Ac(:,5,5)  = F{13};
for ff=1:7
  M{end}.Ac(:,6:end,ff+5)  = F{ff};
end
M{end}.name       = 'All3_eye'; 

% Model 24
M{end+1}.type       = 'feature';
M{end}.numGparams = 11;
M{end}.Ac = zeros(7,11,M{end}.numGparams);
M{end}.Ac(:,1,1)  = F{14};
M{end}.Ac(:,2,2)  = F{15};
M{end}.Ac(:,3,3)  = F{16};
M{end}.Ac(:,4,4)  = F{17};
for ff=1:7
  M{end}.Ac(:,5:end,ff+4)  = F{ff};
end
M{end}.name       = 'all4_eye'; 

% Model 25
M{end+1}.type       = 'feature';
M{end}.numGparams = 10;
M{end}.Ac = zeros(7,10,M{end}.numGparams);
M{end}.Ac(:,1,1)  = F{9};
M{end}.Ac(:,2,2)  = F{18};
M{end}.Ac(:,3,3)  = F{19};
for ff=1:7
  M{end}.Ac(:,4:end,ff+3)  = F{ff};
end
M{end}.name       = 'all5_eye'; 

% Model 26
M{end+1}.type       = 'feature';
M{end}.numGparams = 9;
M{end}.Ac = zeros(7,9,M{end}.numGparams);
M{end}.Ac(:,1,1)  = F{20};
M{end}.Ac(:,2,2)  = F{21};
for ff=1:7
  M{end}.Ac(:,3:end,ff+2)  = F{ff};
end
M{end}.name       = 'all6_eye'; 

% Model 27
M{end+1}.type       = 'feature';
M{end}.numGparams = 16;
M{end}.Ac = zeros(7,16,M{end}.numGparams);
M{end}.Ac(:,1,1)  = F{8};
M{end}.Ac(:,2,2)  = F{10};
M{end}.Ac(:,3,3)  = F{11};
M{end}.Ac(:,4,4)  = F{12};
M{end}.Ac(:,5,5)  = F{13};
M{end}.Ac(:,6,6)  = F{14};
M{end}.Ac(:,7,7)  = F{15};
M{end}.Ac(:,8,8)  = F{16};
M{end}.Ac(:,9,9)  = F{17};
for ff=1:7
  M{end}.Ac(:,10:end,ff+9)  = F{ff};
end
M{end}.name       = 'all3_all4_eye'; 

% Model 28
M{end+1}.type       = 'feature';
M{end}.numGparams = 14;
M{end}.Ac = zeros(7,14,M{end}.numGparams);
M{end}.Ac(:,1,1)  = F{9};
M{end}.Ac(:,2,2)  = F{18};
M{end}.Ac(:,3,3)  = F{19};
M{end}.Ac(:,4,4)  = F{14};
M{end}.Ac(:,5,5)  = F{15};
M{end}.Ac(:,6,6)  = F{16};
M{end}.Ac(:,7,7)  = F{17};
for ff=1:7
  M{end}.Ac(:,8:end,ff+7)  = F{ff};
end
M{end}.name       = 'all4_all5_eye'; 

% Model 29
M{end+1}.type       = 'feature';
M{end}.numGparams = 12;
M{end}.Ac = zeros(7,12,M{end}.numGparams);
M{end}.Ac(:,1,1)  = F{9};
M{end}.Ac(:,2,2)  = F{18};
M{end}.Ac(:,3,3)  = F{19};
M{end}.Ac(:,4,4)  = F{20};
M{end}.Ac(:,5,5)  = F{21};
for ff=1:7
  M{end}.Ac(:,6:end,ff++5)  = F{ff};
end
M{end}.name       = 'all5_all6_eye'; 

% Model 30
M{end+1}.type       = 'feature';
M{end}.numGparams = 19;
M{end}.Ac = zeros(7,19,M{end}.numGparams);
M{end}.Ac(:,1,1)  = F{8};
M{end}.Ac(:,2,2)  = F{10};
M{end}.Ac(:,3,3)  = F{11};
M{end}.Ac(:,4,4)  = F{12};
M{end}.Ac(:,5,5)  = F{13};
M{end}.Ac(:,6,6)  = F{14};
M{end}.Ac(:,7,7)  = F{15};
M{end}.Ac(:,8,8)  = F{16};
M{end}.Ac(:,9,9)  = F{17};
M{end}.Ac(:,10,10)  = F{9};
M{end}.Ac(:,11,11)  = F{18};
M{end}.Ac(:,12,12)  = F{19};
for ff=1:7
  M{end}.Ac(:,13:end,ff+12)  = F{ff};
end
M{end}.name       = 'all3_all4_all5_eye'; 

% Model 31
M{end+1}.type       = 'feature';
M{end}.numGparams = 16;
M{end}.Ac = zeros(7,16,M{end}.numGparams);
M{end}.Ac(:,1,1)  = F{9};
M{end}.Ac(:,2,2)  = F{18};
M{end}.Ac(:,3,3)  = F{19};
M{end}.Ac(:,4,4)  = F{14};
M{end}.Ac(:,5,5)  = F{15};
M{end}.Ac(:,6,6)  = F{16};
M{end}.Ac(:,7,7)  = F{17};
M{end}.Ac(:,8,8)  = F{20};
M{end}.Ac(:,9,9)  = F{21};
for ff=1:7
  M{end}.Ac(:,10:end,ff+9)  = F{ff};
end
M{end}.name       = 'all4_all5_all6_eye'; 





%% 
for mm=1:numel(M)
%   M{mm}.fitAlgorithm='NR'; % NR or minimize
  M{mm}.fitAlgorithm='minimize'; % NR or minimize
end

if 0 % this needs updating
if plotflag
  figure;
  for mm=3:17,subplot(2,15,mm-2);imagesc(sum(M{mm}.Ac,3));end
  for mm=3:17,
    theta0 = ones(size(M{mm}.Ac,3),1);
    Gmodel{mm}=pcm_calculateG(M{mm},theta0);
    subplot(2,15,mm-2+15);imagesc(Gmodel{mm});
  end
  
  
  figure;
  for mm=18:32,subplot(2,15,mm-17);imagesc(sum(M{mm}.Ac,3));end
  for mm=18:32,
    theta0 = ones(size(M{mm}.Ac,3),1);
    Gmodel{mm}=pcm_calculateG(M{mm},theta0);
    subplot(2,15,mm-17+15);imagesc(Gmodel{mm});
  end
end
end

save([edir 'sim5_featmodels.mat'],'M');


