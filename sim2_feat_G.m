% create RDM and G for aud-tac simulations and real data analysis
% Johanna Zumer, 2017


sim2_featAll_G

% Model 1: Eye: All conditions independent variance
Mg{1}.type       = 'feature';
Mg{1}.numGparams = 7;
for ff=1:7
  Mg{1}.Ac(:,1:ncond,ff) = F{ff};
end
Mg{1}.name       = 'eye'; 

% Model 17: TIW a and b;  Asymmetry;  Symmetric Pairs; and eye
Mg{end+1}.type       = 'feature';
Mg{end}.numGparams = 16;
Mg{end}.Ac = zeros(ncond,16,Mg{end}.numGparams);
for ff=1:ncond
  Mg{end}.Ac(:,1:ncond,ff) = F{ff};
end
Mg{end}.Ac(:,8,8)  = F{8};
Mg{end}.Ac(:,9,9)  = F{9};
Mg{end}.Ac(:,10,10)  = F{10};
Mg{end}.Ac(:,11,11)  = F{11};
Mg{end}.Ac(:,12,12)  = F{12};
Mg{end}.Ac(:,13,13)  = F{13};
Mg{end}.Ac(:,14,14)  = F{14};
Mg{end}.Ac(:,15,15)  = F{15};
Mg{end}.Ac(:,16,16)  = F{16};
Mg{end}.name       = 'TIW1a1b_Asym2_SymPair3_eye'; 


theta=zeros(Mg{end}.numGparams,1);
theta(1:7)=1;
G_eye = pcm_calculateG(Mg{2},theta);  % haha, what a complicated way to simply say G_eye=eye(ncond)!


% Temporal integration window
% 1) TIW, varying width and centre in 2 different sub-models
theta=zeros(Mg{end}.numGparams,1);
theta(15:16)=1;
G_1 = pcm_calculateG(Mg{2},theta);  % haha, what a complicated way to simply say G_eye=eye(ncond)!

% 2) Asymmetry: AT will be similar to each other but different from TA (vice versa)
theta=zeros(Mg{end}.numGparams,1);
theta(11:14)=1;
G_2 = pcm_calculateG(Mg{2},theta);  % haha, what a complicated way to simply say G_eye=eye(ncond)!

% 3) Symmetric pairs: TA70 will be more like itself and AT70 than others
theta=zeros(Mg{end}.numGparams,1);
theta(8:10)=1;
G_3 = pcm_calculateG(Mg{2},theta);  % haha, what a complicated way to simply say G_eye=eye(ncond)!

