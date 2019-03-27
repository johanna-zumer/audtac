function design=set_cfg_design_depT(nsub)

design=zeros(2,2*nsub);
design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
design(2,:)=[1:nsub 1:nsub];
