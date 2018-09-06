% set up Models for PCM
% Johanna Zumer, 2017

% plotflag=1;
clear M F
sim7_featAll_G  % Yes, using same features as sim5, but in different arrangement here

% addpath('D:/Matlab/permn/')

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

fint=8:15;  % indicies of features of interest

% 8choose1
featall=nchoosek(1:8,1);
for mm=1:size(featall,1)
  featuse=featall(mm,:);
  M{end+1}.type       = 'feature';
  numpadd=0;for mp=1:length(featuse),numpadd=numpadd+size(F{fint(featuse(mp))},2);end
  M{end}.numGparams = 7+numpadd;
  M{end}.Ac = zeros(7,M{end}.numGparams,M{end}.numGparams);
  for ff=1:7
    M{end}.Ac(:,1:ncond,ff)  = F{ff};
  end
  ffadd=0;
  for ff=1:length(featuse)
    for cc=1:size(F{fint(featuse(ff))},2)
      ffadd=ffadd+1;
      M{end}.Ac(:,7+ffadd,7+ffadd)  = F{fint(featuse(ff))}(:,cc);
    end
  end
  M{end}.name       = ['1Feat_' num2str(featuse(1)) ];
  M{end}.featincl=featuse;
end


% 8choose2
featall=nchoosek(1:8,2);
for mm=1:size(featall,1)
  featuse=featall(mm,:);
  M{end+1}.type       = 'feature';
  numpadd=0;for mp=1:length(featuse),numpadd=numpadd+size(F{fint(featuse(mp))},2);end
  M{end}.numGparams = 7+numpadd;
  M{end}.Ac = zeros(7,M{end}.numGparams,M{end}.numGparams);
  for ff=1:7
    M{end}.Ac(:,1:ncond,ff)  = F{ff};
  end
  ffadd=0;
  for ff=1:length(featuse)
    for cc=1:size(F{fint(featuse(ff))},2)
      ffadd=ffadd+1;
      M{end}.Ac(:,7+ffadd,7+ffadd)  = F{fint(featuse(ff))}(:,cc);
    end
  end
  M{end}.name       = ['2Feat_' num2str(featuse(1)) '_' num2str(featuse(2))];
  M{end}.featincl=featuse;
end


% 8choose3
featall=nchoosek(1:8,3);
for mm=1:size(featall,1)
  featuse=featall(mm,:);
  M{end+1}.type       = 'feature';
  numpadd=0;for mp=1:length(featuse),numpadd=numpadd+size(F{fint(featuse(mp))},2);end
  M{end}.numGparams = 7+numpadd;
  M{end}.Ac = zeros(7,M{end}.numGparams,M{end}.numGparams);
  for ff=1:7
    M{end}.Ac(:,1:ncond,ff)  = F{ff};
  end
  ffadd=0;
  for ff=1:length(featuse)
    for cc=1:size(F{fint(featuse(ff))},2)
      ffadd=ffadd+1;
      M{end}.Ac(:,7+ffadd,7+ffadd)  = F{fint(featuse(ff))}(:,cc);
    end
  end
  M{end}.name       = ['3Feat_' num2str(featuse(1)) '_' num2str(featuse(2)) '_' num2str(featuse(3))];
  M{end}.featincl=featuse;
end

% 8choose4
featall=nchoosek(1:8,4);
for mm=1:size(featall,1)
  featuse=featall(mm,:);
  M{end+1}.type       = 'feature';
  numpadd=0;for mp=1:length(featuse),numpadd=numpadd+size(F{fint(featuse(mp))},2);end
  M{end}.numGparams = 7+numpadd;
  M{end}.Ac = zeros(7,M{end}.numGparams,M{end}.numGparams);
  for ff=1:7
    M{end}.Ac(:,1:ncond,ff)  = F{ff};
  end
  ffadd=0;
  for ff=1:length(featuse)
    for cc=1:size(F{fint(featuse(ff))},2)
      ffadd=ffadd+1;
      M{end}.Ac(:,7+ffadd,7+ffadd)  = F{fint(featuse(ff))}(:,cc);
    end
  end
  M{end}.name       = ['4Feat_' num2str(featuse(1)) '_' num2str(featuse(2)) '_' num2str(featuse(3)) '_' num2str(featuse(4))];
  M{end}.featincl=featuse;
end


% 8choose5
featall=nchoosek(1:8,5);
for mm=1:size(featall,1)
  featuse=featall(mm,:);
  M{end+1}.type       = 'feature';
  numpadd=0;for mp=1:length(featuse),numpadd=numpadd+size(F{fint(featuse(mp))},2);end
  M{end}.numGparams = 7+numpadd;
  M{end}.Ac = zeros(7,M{end}.numGparams,M{end}.numGparams);
  for ff=1:7
    M{end}.Ac(:,1:ncond,ff)  = F{ff};
  end
  ffadd=0;
  for ff=1:length(featuse)
    for cc=1:size(F{fint(featuse(ff))},2)
      ffadd=ffadd+1;
      M{end}.Ac(:,7+ffadd,7+ffadd)  = F{fint(featuse(ff))}(:,cc);
    end
  end
  M{end}.name       = ['5Feat_' num2str(featuse(1)) '_' num2str(featuse(2)) '_' num2str(featuse(3)) '_' num2str(featuse(4)) '_' num2str(featuse(5))];
  M{end}.featincl=featuse;
end

% 8choose6
featall=nchoosek(1:8,6);
for mm=1:size(featall,1)
  featuse=featall(mm,:);
  M{end+1}.type       = 'feature';
  numpadd=0;for mp=1:length(featuse),numpadd=numpadd+size(F{fint(featuse(mp))},2);end
  M{end}.numGparams = 7+numpadd;
  M{end}.Ac = zeros(7,M{end}.numGparams,M{end}.numGparams);
  for ff=1:7
    M{end}.Ac(:,1:ncond,ff)  = F{ff};
  end
  ffadd=0;
  for ff=1:length(featuse)
    for cc=1:size(F{fint(featuse(ff))},2)
      ffadd=ffadd+1;
      M{end}.Ac(:,7+ffadd,7+ffadd)  = F{fint(featuse(ff))}(:,cc);
    end
  end
  M{end}.name       = ['6Feat_' num2str(featuse(1)) '_' num2str(featuse(2)) '_' num2str(featuse(3)) '_' num2str(featuse(4)) '_' num2str(featuse(5)) '_' num2str(featuse(6))];
  M{end}.featincl=featuse;
end

% 8choose7
featall=nchoosek(1:8,7);
for mm=1:size(featall,1)
  featuse=featall(mm,:);
  M{end+1}.type       = 'feature';
  numpadd=0;for mp=1:length(featuse),numpadd=numpadd+size(F{fint(featuse(mp))},2);end
  M{end}.numGparams = 7+numpadd;
  M{end}.Ac = zeros(7,M{end}.numGparams,M{end}.numGparams);
  for ff=1:7
    M{end}.Ac(:,1:ncond,ff)  = F{ff};
  end
  ffadd=0;
  for ff=1:length(featuse)
    for cc=1:size(F{fint(featuse(ff))},2)
      ffadd=ffadd+1;
      M{end}.Ac(:,7+ffadd,7+ffadd)  = F{fint(featuse(ff))}(:,cc);
    end
  end
  M{end}.name       = ['7Feat_' num2str(featuse(1)) '_' num2str(featuse(2)) '_' num2str(featuse(3)) '_' num2str(featuse(4)) '_' num2str(featuse(5)) '_' num2str(featuse(6)) '_' num2str(featuse(7))];
  M{end}.featincl=featuse;
end

% 8choose8
featall=nchoosek(1:8,8);
for mm=1:size(featall,1)
  featuse=featall(mm,:);
  M{end+1}.type       = 'feature';
  numpadd=0;for mp=1:length(featuse),numpadd=numpadd+size(F{fint(featuse(mp))},2);end
  M{end}.numGparams = 7+numpadd;
  M{end}.Ac = zeros(7,M{end}.numGparams,M{end}.numGparams);
  for ff=1:7
    M{end}.Ac(:,1:ncond,ff)  = F{ff};
  end
  ffadd=0;
  for ff=1:length(featuse)
    for cc=1:size(F{fint(featuse(ff))},2)
      ffadd=ffadd+1;
      M{end}.Ac(:,7+ffadd,7+ffadd)  = F{fint(featuse(ff))}(:,cc);
    end
  end
  M{end}.name       = ['8Feat_' num2str(featuse(1)) '_' num2str(featuse(2)) '_' num2str(featuse(3)) '_' num2str(featuse(4)) '_' num2str(featuse(5)) '_' num2str(featuse(6)) '_' num2str(featuse(7)) '_' num2str(featuse(8))];
  M{end}.featincl=featuse;
end

%% 
for mm=1:numel(M)
  M{mm}.fitAlgorithm='NR'; % NR or minimize
%   M{mm}.fitAlgorithm='minimize'; % NR or minimize
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

save([edir 'sim7NR_featmodels.mat'],'M');


