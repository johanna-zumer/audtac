% sim1_featAll_G
%
% Each F is a feature
% One or more features come together to make a component

ncond=7;

% 'Component 1' Identity
% 1-7 constitute 'eye' (each condition has some individual variance
% for ff=1:ncond
%   F{ff}=zeros(ncond);
%   F{ff}(ff,ff)=1;
% end
F{1}=zeros(ncond,ncond);
F{1}(1,1)=1;
F{1}(7,7)=1;

F{2}=zeros(ncond,ncond);
F{2}(2,2)=1;
F{2}(6,6)=1;

F{3}=zeros(ncond,ncond);
F{3}(3,3)=1;
F{3}(5,5)=1;

F{4}=zeros(ncond,ncond);
F{4}(4,4)=1;

% Component 0: Background everywhere
F{5}=ones(ncond,1);

% 'Component 2': Symmetric Pairs
% F{8}(3,5)=1;
% F{8}(5,3)=1;
% F{9}(2,6)=1;
% F{9}(6,2)=1;
% F{10}(1,7)=1;
% F{10}(7,1)=1;
F{8}=[0 0 1 0 1 0 0]';
F{9}=[0 1 0 0 0 1 0]';
F{10}=[1 0 0 0 0 0 1]';

% 'Component 3': Asymmetry
F{11}=[1 1 1 1 0 0 0]';
F{12}=[1 1 1 0 0 0 0]';
F{13}=[0 0 0 1 1 1 1]';
F{14}=[0 0 0 0 1 1 1]';

% 'Component 4': Temporal Integration Window
F{15}=[0 0 1 1 1 0 0]';
F{16}=[0 1 1 1 1 1 0]';

%%


