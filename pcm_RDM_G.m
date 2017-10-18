% create RDM and G for aud-tac simulations and real data analysis


% start with generating G matrices: start with dissimilarity matrix and convert to G (second moment)
G_eye=eye(7);


% Temporal integration window
% 1) TIW, varying width and centre
rdmodel_1a=1-eye(7);
rdmodel_1a(3:5,3:5)=0;
%     G_1a=nearestSPD(rdm2G(rdmodel_1a,plotflag,'model 1a'));
% normbyrownorm recommended based on Diedcrichsen paper 2017 just after equation 30
G_1a=nearestSPD(rdm2G(jz_normbyrownorm(rdmodel_1a),plotflag,'model 1a'));

rdmodel_1b=1-eye(7);
rdmodel_1b(2:6,2:6)=0;
G_1b=nearestSPD(rdm2G(jz_normbyrownorm(rdmodel_1b),plotflag,'model 1b'));


rdmodel_1c=1-eye(7);
rdmodel_1c(2:5,2:5)=0;
G_1c=nearestSPD(rdm2G(jz_normbyrownorm(rdmodel_1c),plotflag,'model 1c'));

rdmodel_1d=1-eye(7);
rdmodel_1d(3:6,3:6)=0;
G_1d=nearestSPD(rdm2G(jz_normbyrownorm(rdmodel_1d),plotflag,'model 1d'));

rdmodel_1e=1-eye(7);
rdmodel_1e(2:4,2:4)=0;
G_1e=nearestSPD(rdm2G(jz_normbyrownorm(rdmodel_1e),plotflag,'model 1e'));

rdmodel_1f=1-eye(7);
rdmodel_1f(4:6,4:6)=0;
G_1f=nearestSPD(rdm2G(jz_normbyrownorm(rdmodel_1f),plotflag,'model 1f'));

% 2) Asymmetry: AT will be similar to each other but different from TA (vice versa)
rdmodel_2a=1-eye(7);
rdmodel_2a(1:4,1:4)=0;
rdmodel_2a(4:7,4:7)=0;
rdmodel_2a_norm = jz_normbyrownorm(rdmodel_2a);
rdmodel_2a_norm(4,:)=eps;
G_2a=nearestSPD(rdm2G(rdmodel_2a_norm,plotflag,'model 2'));



% 3) Symmetric pairs: TA70 will be more like itself and AT70 than others
rdmodel_3a=1-eye(7);
rdmodel_3a(1,7)=0;
rdmodel_3a(2,6)=0;
rdmodel_3a(3,5)=0;
rdmodel_3a(4,4)=0;
rdmodel_3a(5,3)=0;
rdmodel_3a(6,2)=0;
rdmodel_3a(7,1)=0;
G_3a=nearestSPD(rdm2G(jz_normbyrownorm(rdmodel_3a),plotflag,'model 3'));




% Interaction terms:
interact_1a2a=rdmodel_1a.*rdmodel_2a;
interact_1b2a=rdmodel_1b.*rdmodel_2a;
interact_1a3a=rdmodel_1a.*rdmodel_3a;
interact_1b3a=rdmodel_1b.*rdmodel_3a;
interact_2a3a=rdmodel_2a.*rdmodel_3a;
G_1a2a=nearestSPD(rdm2G(interact_1a2a,plotflag,'model 1a X 2'));
G_1b2a=nearestSPD(rdm2G(interact_1b2a,plotflag,'model 1b X 2'));
G_1a3a=nearestSPD(rdm2G(interact_1a3a,plotflag,'model 1a X 3'));
G_1b3a=nearestSPD(rdm2G(interact_1b3a,plotflag,'model 1b X 3'));
G_2a3a=nearestSPD(rdm2G(interact_2a3a,plotflag,'model 2 X 3'));


