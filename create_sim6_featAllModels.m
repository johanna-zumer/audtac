% set up Models for PCM
% Johanna Zumer, 2017

% plotflag=1;
clear M F
sim5_featAll_G  % Yes, using same features as sim5, but in different arrangement here

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

% Model 3: 
M{end+1}.type       = 'feature';
M{end}.numGparams = 7;
M{end}.Ac = zeros(7,7,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.name       = 'eye';

% Model 4: On TIW
M{end+1}.type       = 'feature';
M{end}.numGparams = 9;
M{end}.Ac = zeros(7,9,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
M{end}.Ac(:,8,9)  = F{8};
M{end}.Ac(:,8,9)  = F{9};
M{end}.name       = 'OnTIW_eye'; 

% Model 5: 
M{end+1}.type       = 'feature';
M{end}.numGparams = 15;
M{end}.Ac = zeros(7,15,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
for ff=1:8
  M{end}.Ac(:,ff+7,ff+7)  = F{ff+7}; 
end
M{end}.name       = 'OnTIW_Asym_eye';

% Model 6: 
M{end+1}.type       = 'feature';
M{end}.numGparams = 17;
M{end}.Ac = zeros(7,17,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
for ff=1:10
  M{end}.Ac(:,ff+7,ff+7)  = F{ff+7}; 
end
M{end}.name       = 'OnTIW_Asym_Green_eye';

% Model 7: 
M{end+1}.type       = 'feature';
M{end}.numGparams = 19;
M{end}.Ac = zeros(7,19,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
for ff=1:12
  M{end}.Ac(:,ff+7,ff+7)  = F{ff+7}; 
end
M{end}.name       = 'OnTIW_Asym_Green_Cyan_eye';

% Model 8: 
M{end+1}.type       = 'feature';
M{end}.numGparams = 13;
M{end}.Ac = zeros(7,13,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
for ff=1:6
  M{end}.Ac(:,ff+7,ff+7)  = F{ff+9}; 
end
M{end}.name       = 'Asym_eye';

% Model 9: 
M{end+1}.type       = 'feature';
M{end}.numGparams = 15;
M{end}.Ac = zeros(7,15,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
for ff=1:8
  M{end}.Ac(:,ff+7,ff+7)  = F{ff+9}; 
end
M{end}.name       = 'Asym_Green_eye';

% Model 10: 
M{end+1}.type       = 'feature';
M{end}.numGparams = 17;
M{end}.Ac = zeros(7,17,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
for ff=1:10
  M{end}.Ac(:,ff+7,ff+7)  = F{ff+9}; 
end
M{end}.name       = 'Asym_Green_Cyan_eye';

% Model 11: 
M{end+1}.type       = 'feature';
M{end}.numGparams = 9;
M{end}.Ac = zeros(7,9,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
for ff=1:2
  M{end}.Ac(:,ff+7,ff+7)  = F{ff+15}; 
end
M{end}.name       = 'Green_eye';


% Model 12: 
M{end+1}.type       = 'feature';
M{end}.numGparams = 11;
M{end}.Ac = zeros(7,11,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
for ff=1:4
  M{end}.Ac(:,ff+7,ff+7)  = F{ff+15}; 
end
M{end}.name       = 'Green_Cyan_eye';

% Model 13: 
M{end+1}.type       = 'feature';
M{end}.numGparams = 9;
M{end}.Ac = zeros(7,9,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
for ff=1:2
  M{end}.Ac(:,ff+7,ff+7)  = F{ff+17}; 
end
M{end}.name       = 'Cyan_eye';

% Model 14: 
M{end+1}.type       = 'feature';
M{end}.numGparams = 11;
M{end}.Ac = zeros(7,11,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
for ff=1:2
  M{end}.Ac(:,ff+7,ff+7)  = F{ff+7}; 
end
for ff=1:2
  M{end}.Ac(:,ff+9,ff+9)  = F{ff+15}; 
end
M{end}.name       = 'OnTIW_Green_eye';


% Model 15: 
M{end+1}.type       = 'feature';
M{end}.numGparams = 11;
M{end}.Ac = zeros(7,11,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
for ff=1:2
  M{end}.Ac(:,ff+7,ff+7)  = F{ff+7}; 
end
for ff=1:2
  M{end}.Ac(:,ff+9,ff+9)  = F{ff+17}; 
end
M{end}.name       = 'OnTIW_Cyan_eye';


% Model 16: 
M{end+1}.type       = 'feature';
M{end}.numGparams = 13;
M{end}.Ac = zeros(7,13,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
for ff=1:2
  M{end}.Ac(:,ff+7,ff+7)  = F{ff+7}; 
end
for ff=1:2
  M{end}.Ac(:,ff+9,ff+9)  = F{ff+15}; 
end
for ff=1:2
  M{end}.Ac(:,ff+11,ff+11)  = F{ff+17}; 
end
M{end}.name       = 'OnTIW_Green_Cyan_eye';

% Model 17: 
M{end+1}.type       = 'feature';
M{end}.numGparams = 15;
M{end}.Ac = zeros(7,15,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
for ff=1:6
  M{end}.Ac(:,ff+7,ff+7)  = F{ff+9}; 
end
for ff=1:2
  M{end}.Ac(:,ff+13,ff+13)  = F{ff+17}; 
end
M{end}.name       = 'Asym_Cyan_eye';

% Model 18: 
M{end+1}.type       = 'feature';
M{end}.numGparams = 17;
M{end}.Ac = zeros(7,17,M{end}.numGparams);
for ff=1:7
  M{end}.Ac(:,1:ncond,ff)  = F{ff};
end
for ff=1:6
  M{end}.Ac(:,ff+7,ff+7)  = F{ff+9}; 
end
for ff=1:2
  M{end}.Ac(:,ff+13,ff+13)  = F{ff+17}; 
end
for ff=1:2
  M{end}.Ac(:,ff+15,ff+15)  = F{ff+7}; 
end
M{end}.name       = 'Asym_Cyan_OnTIW_eye';


%% 
% Model 19
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

% Model 20
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

% Model 21
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

% Model 22
M{end+1}.type       = 'feature';
M{end}.numGparams = 9;
M{end}.Ac = zeros(7,9,M{end}.numGparams);
M{end}.Ac(:,1,1)  = F{20};
M{end}.Ac(:,2,2)  = F{21};
for ff=1:7
  M{end}.Ac(:,3:end,ff+2)  = F{ff};
end
M{end}.name       = 'all6_eye'; 

% Model 23
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

% Model 24
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

% Model 25
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

% Model 26
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

% Model 27
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

save([edir 'sim6_featmodels.mat'],'M');


