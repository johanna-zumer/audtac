% sim1_featAll_G
%
% Each F is a feature
% One or more features come together to make a component

clear F
ncond=7;

% 'Component 1' Identity
% 1-7 constitute 'eye' (each condition has some individual variance
for ff=1:ncond
  F{ff}=zeros(ncond);
  F{ff}(ff,ff)=1;
end
% F{1}=zeros(ncond,ncond);
% F{1}(1,1)=1;
% F{1}(7,7)=1;
% 
% F{2}=zeros(ncond,ncond);
% F{2}(2,2)=1;
% F{2}(6,6)=1;
% 
% F{3}=zeros(ncond,ncond);
% F{3}(3,3)=1;
% F{3}(5,5)=1;
% 
% F{4}=zeros(ncond,ncond);
% F{4}(4,4)=1;

% Coded for Asymmetric or Symmetric
% TIW 3, 4, 5, or 6

% 'Component 3': On-centre Temporal Integration Window
F{8}=[0 0 1 1 1 0 0]';  % S 3

F{9}=[0 1 1 1 1 1 0]';  % S 5

% 'Component 2': Asymmetry of different sizes 
F{10}=[1 1 1 0 0 0 0; 0 0 0 0 1 1 1 ]';  % A 3

F{11}=[0 1 1 1 0 0 0; 0 0 0 1 1 1 0]';  % A 3

F{12}=[1 1 1 1 0 0 0; 0 0 0 1 1 1 1]';  % A 4

% 'Component 4': Off-centre TIW
F{13}=[0 1 1 1 1 0 0; 0 0 1 1 1 1 0]';  % A 4

F{14}=[1 1 1 1 1 0 0; 0 0 1 1 1 1 1]';  % A 5

F{15}=[1 1 1 1 1 1 0; 0 1 1 1 1 1 1]';  % A 6

