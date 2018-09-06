% set up Models for PCM
% Johanna Zumer, 2017

% plotflag=1;
clear M F

%% Features

% sim2_featAll_G
% 'Component 1': Diagonal
for ff=1:ncond
  F{ff}=zeros(ncond);
  F{ff}(ff,ff)=1;
end
% 'Component 2': Symmetric Pairs
F{8}=[0 0 1 0 1 0 0]';
F{9}=[0 1 0 0 0 1 0]';
F{10}=[1 0 0 0 0 0 1]';
% 'Component 3': Asymmetry
F{11}=[1 1 1 1 0 0 0]';
F{12}=[1 1 1 0 0 0 0]';
F{13}=[0 0 0 1 1 1 1]';
F{14}=[0 0 0 0 1 1 1]';
% 'Component 4': Integration Window
F{15}=[0 0 1 1 1 0 0]';
F{16}=[0 1 1 1 1 1 0]';


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
  M{3}.Ac(:,:,ff) = F{ff};
end
M{3}.name       = 'eye'; 

% Model 4: TIW a and b
M{4}.type       = 'feature';
M{4}.numGparams = 2;
M{4}.Ac = zeros(7,2,M{4}.numGparams);
M{4}.Ac(:,1,1)  = F{15};
M{4}.Ac(:,2,2)  = F{16};
M{4}.name       = 'TIW1a1b'; 

% Model 5: Asymmetry
M{5}.type       = 'feature';
M{5}.numGparams = 4;
M{5}.Ac = zeros(7,4,M{5}.numGparams);
M{5}.Ac(:,1,1)  = F{11};
M{5}.Ac(:,2,2)  = F{12};
M{5}.Ac(:,3,3)  = F{13};
M{5}.Ac(:,4,4)  = F{14};
M{5}.name       = 'Asym2'; 

% Model 6: Symmetric Pairs
M{6}.type       = 'feature';
M{6}.numGparams = 3;
M{6}.Ac = zeros(7,3,M{6}.numGparams);
M{6}.Ac(:,1,1)  = F{8};
M{6}.Ac(:,2,2)  = F{9};
M{6}.Ac(:,3,3)  = F{10};
M{6}.name       = 'SymPair3'; 

% Model 7: TIW a and b;  Asymmetry
M{7}.type       = 'feature';
M{7}.numGparams = 6;
M{7}.Ac = zeros(7,6,M{7}.numGparams);
M{7}.Ac(:,1,1)  = F{15};
M{7}.Ac(:,2,2)  = F{16};
M{7}.Ac(:,3,3)  = F{11};
M{7}.Ac(:,4,4)  = F{12};
M{7}.Ac(:,5,5)  = F{13};
M{7}.Ac(:,6,6)  = F{14};
M{7}.name       = 'TIW1a1b_Asym2'; 

% Model 8: TIW a and b;  Symmetric Pairs
M{8}.type       = 'feature';
M{8}.numGparams = 5;
M{8}.Ac = zeros(7,5,M{8}.numGparams);
M{8}.Ac(:,1,1)  = F{8};
M{8}.Ac(:,2,2)  = F{9};
M{8}.Ac(:,3,3)  = F{10};
M{8}.Ac(:,4,4)  = F{15};
M{8}.Ac(:,5,5)  = F{16};
M{8}.name       = 'TIW1a1b_SymPair3'; 

% Model 9: Asymmetry;  Symmetric Pairs
M{9}.type       = 'feature';
M{9}.numGparams = 7;
M{9}.Ac = zeros(7,7,M{9}.numGparams);
M{9}.Ac(:,1,1)  = F{8};
M{9}.Ac(:,2,2)  = F{9};
M{9}.Ac(:,3,3)  = F{10};
M{9}.Ac(:,4,4)  = F{11};
M{9}.Ac(:,5,5)  = F{12};
M{9}.Ac(:,6,6)  = F{13};
M{9}.Ac(:,7,7)  = F{14};
M{9}.name       = 'Asym2_SymPair3'; 

% Model 10: TIW a and b;  Asymmetry;  Symmetric Pairs
M{10}.type       = 'feature';
M{10}.numGparams = 9;
M{10}.Ac = zeros(7,9,M{9}.numGparams);
M{10}.Ac(:,1,1)  = F{8};
M{10}.Ac(:,2,2)  = F{9};
M{10}.Ac(:,3,3)  = F{10};
M{10}.Ac(:,4,4)  = F{11};
M{10}.Ac(:,5,5)  = F{12};
M{10}.Ac(:,6,6)  = F{13};
M{10}.Ac(:,7,7)  = F{14};
M{10}.Ac(:,8,8)  = F{15};
M{10}.Ac(:,9,9)  = F{16};
M{10}.name       = 'TIW1a1b_Asym2_SymPair3'; 

% Model 11: TIW a and b; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 9;
M{end}.Ac = zeros(7,9,M{end}.numGparams);
for ff=1:ncond
  M{end}.Ac(:,1:ncond,ff) = F{ff};
end
M{end}.Ac(:,8,8)  = F{15};
M{end}.Ac(:,9,9)  = F{16};
M{end}.name       = 'TIW1a1b_eye'; 

% Model 12: Asymmetry; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 11;
M{end}.Ac = zeros(ncond,11,M{end}.numGparams);
for ff=1:ncond
  M{end}.Ac(:,1:ncond,ff) = F{ff};
end
M{end}.Ac(:,8,8)  = F{11};
M{end}.Ac(:,9,9)  = F{12};
M{end}.Ac(:,10,10)  = F{13};
M{end}.Ac(:,11,11)  = F{14};
M{end}.name       = 'Asym2_eye'; 

% Model 13: Symmetric Pairs; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 10;
M{end}.Ac = zeros(ncond,10,M{end}.numGparams);
for ff=1:ncond
  M{end}.Ac(:,1:ncond,ff) = F{ff};
end
M{end}.Ac(:,8,8)  = F{8};
M{end}.Ac(:,9,9)  = F{9};
M{end}.Ac(:,10,10)  = F{10};
M{end}.name       = 'SymPair3_eye'; 

% Model 14: TIW a and b;  Asymmetry; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 13;
M{end}.Ac = zeros(ncond,13,M{end}.numGparams);
for ff=1:ncond
  M{end}.Ac(:,1:ncond,ff) = F{ff};
end
M{end}.Ac(:,8,8)  = F{11};
M{end}.Ac(:,9,9)  = F{12};
M{end}.Ac(:,10,10)  = F{13};
M{end}.Ac(:,11,11)  = F{14};
M{end}.Ac(:,12,12)  = F{15};
M{end}.Ac(:,13,13)  = F{16};
M{end}.name       = 'TIW1a1b_Asym2_eye'; 

% Model 15: TIW a and b;  Symmetric Pairs; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 12;
M{end}.Ac = zeros(ncond,12,M{end}.numGparams);
for ff=1:ncond
  M{end}.Ac(:,1:ncond,ff) = F{ff};
end
M{end}.Ac(:,8,8)  = F{8};
M{end}.Ac(:,9,9)  = F{9};
M{end}.Ac(:,10,10)  = F{10};
M{end}.Ac(:,11,11)  = F{15};
M{end}.Ac(:,12,12)  = F{16};
M{end}.name       = 'TIW1a1b_SymPair3_eye'; 

% Model 16: Asymmetry;  Symmetric Pairs; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 14;
M{end}.Ac = zeros(ncond,14,M{end}.numGparams);
for ff=1:ncond
  M{end}.Ac(:,1:ncond,ff) = F{ff};
end
M{end}.Ac(:,8,8)  = F{8};
M{end}.Ac(:,9,9)  = F{9};
M{end}.Ac(:,10,10)  = F{10};
M{end}.Ac(:,11,11)  = F{11};
M{end}.Ac(:,12,12)  = F{12};
M{end}.Ac(:,13,13)  = F{13};
M{end}.Ac(:,14,14)  = F{14};
M{end}.name       = 'Asym2_SymPair3_eye'; 

% Model 17: TIW a and b;  Asymmetry;  Symmetric Pairs; and eye
M{end+1}.type       = 'feature';
M{end}.numGparams = 16;
M{end}.Ac = zeros(ncond,16,M{end}.numGparams);
for ff=1:ncond
  M{end}.Ac(:,1:ncond,ff) = F{ff};
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
M{end}.name       = 'TIW1a1b_Asym2_SymPair3_eye'; 

for mm=1:numel(M)
%   M{mm}.fitAlgorithm='NR'; % NR or minimize
  M{mm}.fitAlgorithm='minimize'; % NR or minimize
end

if plotflag
  figure;
  for mm=3:17,subplot(2,15,mm-2);imagesc(sum(M{mm}.Ac,3));end
  for mm=3:17,
    if mm>1
      theta0 = ones(size(M{mm}.Ac,3),1);
      Gmodel{mm}=pcm_calculateG(M{mm},theta0);
    end
    subplot(2,15,mm-2+15);imagesc(Gmodel{mm});
  end
end

save([edir 'sim2_featmodels.mat'],'M');


