function design=set_cfg_design_depF(nsub)

design=zeros(2,7*nsub);
design(1,:)=[ones(1,nsub) 2*ones(1,nsub) 3*ones(1,nsub) 4*ones(1,nsub) 5*ones(1,nsub) 6*ones(1,nsub) 7*ones(1,nsub)];
design(2,:)=[1:nsub 1:nsub 1:nsub 1:nsub 1:nsub 1:nsub 1:nsub];
