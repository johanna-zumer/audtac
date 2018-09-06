% Each F is a feature
% One or more features come together to make a component

clear F
ncond=7;

% 'Component 1' Identity
% for ff=1:ncond
%   F{ff}=zeros(ncond);
%   F{ff}(ff,ff)=1;
% end
F{1}=eye(7);


% 'Component 2': Symmetric Pairs
F{8}=[0 0 1 0 1 0 0]';
F{9}=[0 1 0 0 0 1 0]';
F{10}=[1 0 0 0 0 0 1]';

% 'Component 3': Asymmetry
F{11}=[1 1 1 1 0 0 0]';
F{13}=[0 0 0 1 1 1 1]';


% 'Component 4': Temporal Integration Window
F{15}=[0 0 1 1 1 0 0]';
F{16}=[0 1 1 1 1 1 0]';

% Component 0: Background everywhere
F{21}=[1 1 1 1 1 1 1]';

%%


