% MRI head models

clear all
close all
if ispc
  edir='D:\audtac\eeg_data\';
  esdir='D:\audtac\source_data\';
  ddir='D:\audtac\legomagic\diaries\';
  bdir='D:\audtac\behav_data\';
  sdir='D:\audtac\spss_stuff\';
  mdir='D:\audtac\structural_MRI\';
  pdir='D:\audtac\polhemus\';
else
  [~,hostname]=system('hostname');
  if ~isempty(strfind(hostname,'les')) | ~isempty(strfind(hostname,'LES')) % either COLLES-151401 or LES-LINUX_FS3
    edir='/home/zumerj/audtac/eeg_data/';
    esdir='/home/zumerj/audtac/source_data/';
    ddir='/home/zumerj/audtac/legomagic/diaries/';
    bdir='/home/zumerj/audtac/behav_data/';
    sdir='/home/zumerj/audtac/spss_stuff/';
    mdir='/home/zumerj/audtac/structural_MRI/';
    pdir='/home/zumerj/audtac/polhemus/';
  else % assume on VM linux of psychl-132432
    edir='/mnt/hgfs/D/audtac/eeg_data/';
    esdir='/mnt/hgfs/D/audtac/source_data/';
    ddir='/mnt/hgfs/D/audtac/legomagic/diaries/';
    bdir='/mnt/hgfs/D/audtac/behav_data/';
    sdir='/mnt/hgfs/D/audtac/spss_stuff/';
    mdir='/mnt/hgfs/D/audtac/structural_MRI/';
    pdir='/mnt/hgfs/D/audtac/polhemus/';
  end
end

sub{5}='e05'; % ma.a. 25/06/14
sub{6}='e06'; % e.u.  01/07/14
sub{7}='e07'; % a.s.  07/07/14
sub{8}='e08'; % k.t.  09/07/14
sub{9}='e09';% d.a.  14/07/14
sub{10}='e10';% k.l.  15/07/14
sub{11}='e11';% ab.m.  16/07/14  % from here on, had EOG/EMG
sub{12}='e12';% b.s.  17/07/14
sub{13}='e13';% d.t.  21/07/14
sub{14}='e14';% f.g.  22/07/14
sub{15}='e15';% r.m.  23/07/14
sub{16}='e16';% t.p.  24/07/14 % from here on, had attempt Polhemus
sub{17}='e17';% t.t.  28/07/14
sub{18}='e18';% k.n.v.  29/07/14
sub{19}='e19';% j.b.  30/07/14
sub{20}='e20';% n.m.  31/07/14
sub{21}='e21';% l.c.  04/08/14
sub{22}='e22';% a.b.  05/08/14
sub{23}='e23';% r.c.
sub{24}='e24';% a.d.
sub{25}='e25';% j.c.
sub{26}='e26';% r.s.
sub{27}='e27';% a.p.
sub{28}='e28';% w.p.
sub{29}='e29';% i.r.
sub{30}='e30';% o.y.l.
sub{31}='e31';% r.b.
sub{32}='e32';% i.f.

if ispc
  rmpath(genpath('D:\matlab\spm8\external\fieldtrip\'))
  rmpath(genpath('D:\fieldtrip_svn\'))
  addpath('D:\fieldtrip_svn\')
else
  rmpath(genpath('/mnt/hgfs/D/matlab/spm8/external/fieldtrip/'))
  rmpath(genpath('/mnt/hgfs/D/fieldtrip_svn/'))
  addpath('/mnt/hgfs/D/fieldtrip_svn/')
end
which ft_defaults.m
ft_defaults;

fwd=fileparts(which('ft_defaults.m'))
if ispc
  fwd=[fwd '\'];
else
  fwd=[fwd '/'];
end

load([edir 'iikeep.mat'])
iiuse=setdiff(union(iiBuse,iiSuse),[3 11]);

%% Overview

% Precursor: get individual MRI realigned & warped to MNI
%
% Need 3 components for forward model:
% 1) Volume conductor model: Use precomputed standard BEM or computed from individual MRI
% 2) Electrode positions (polhemus or standard) aligned to MRI
% 3) Grid points: use MNI grid, warped to subject-specific MRI)
%
% Finally, compute lead field


%% Precursor: Realign/normalise individual MRI to standard template MNI
for ii=iiuse
  cd([mdir sub{ii} ])
  try
    load mrinorm.mat
  catch
    if ispc
      mri_orig_name=[mdir sub{ii} '\' sub{ii} 'struct.nii'   ];
    else
      mri_orig_name=[mdir sub{ii} '/' sub{ii} 'struct.nii'   ];
    end
    mri = ft_read_mri(mri_orig_name);
    
    cfg=[];
    cfg.coordsys='ras';
    % % MUST this occur linear, not nonlinear so that fiducials match up
    % later?
    cfg.nonlinear='no';
    mrinormlin=ft_volumenormalise(cfg,mri); % I've modified spm_defaults to expand .bb
    mrinormlin=ft_convert_units(mrinormlin,'mm');
    
    % The above 'mrinormlin' should be considered subject-specific MRI (but
    % in Affine-coregistered to MNI space
    
    cfg=[];
    cfg.coordsys='ras';
    % % MUST this occur linear, not nonlinear so that fiducials match up
    % later?
    cfg.nonlinear='yes';
    %   mrinormwarp=ft_volumenormalise(cfg,mri); % I've modified spm_defaults to expand .bb
    mrinormwarp=ft_volumenormalise(cfg,mrinormlin); % I've modified spm_defaults to expand .bb
    mrinormwarp=ft_convert_units(mrinormwarp,'mm');
    close all
    save('mrinorm','mrinormlin','mrinormwarp');
    
    params=mrinormwarp.params;
    save('params_sn.mat','-struct','params'); % to be used later with John's Gems
  end
end

%% Compute BEM for each subject

% 3 options:  (in order of preference, do try/catch sequentially until last
% one succeeds)
%
% 1) create segmentation and BEM from individual MRI.  Try first, test if
% fails (only works for 5 participants)
% 2) use existing meshes from existing BEM but warp shape to individual
%    head, then recompute vol  (works for everyone)
% 3) just use existing BEM vol but with warped pnts. (no need as option 2
% works)

bemoption=0;
for ii=setdiff(iiuse,[5])
  cd([mdir sub{ii} ])
  
  % Option 1
  try % create segmentation and BEM from individual MRI
    % not using b/c segmentation often fails, but worth a try
    try
      load segmentedmri.mat
    catch
      
      load('mrinorm');
      
      cfg           = [];
      cfg.output    = {'brain','skull','scalp'};
      segmentedmri  = ft_volumesegment(cfg, mrinormlin);
      
      % % %
      keyboard; % see test_bug1954 for an improvement on this.
      if 0
        scalp = imfill(segmentedmri.scalp, [128 128 128]); % these numbers may need to change?
        skull = imfill(segmentedmri.skull, [128 128 128]);
        brain = imfill(segmentedmri.brain, [128 128 128]);
        
        % the brain and skull go too far down at the bottom
        skull(:,:,1:5) = 0;
        brain(:,:,1:10) = 0;
        
        % ensure that they don't overlap
        skull = skull & imerode(scalp, strel_bol(2));
        brain = brain & imerode(skull, strel_bol(2));
        
        segmentedmri.scalp = scalp;
        segmentedmri.skull = skull;
        segmentedmri.brain = brain;
      end
      % % %
      
      save segmentedmri segmentedmri
    end
    
    cfg=[];
    cfg.tissue={'brain','skull','scalp'};
    cfg.numvertices = [3000 2000 1000];
    cfg.method='projectmesh';
    bnd=ft_prepare_mesh(cfg,segmentedmri);
    
    cfg        = [];
    cfg.method ='dipoli';
    voldipolimm     = ft_prepare_headmodel(cfg, ft_convert_units(bnd,'mm'));
    try
      voldipolimm.mat;
      disp('got 1')
    catch
      
      % fix bnd
      ft_hastoolbox('iso2mesh',1)
      cd([fwd 'private'])
      bnd = decouplesurf(bnd);    % decouplesurf is an unimplemented subfunction temporarily stashed in prepare_mesh_segmentation
      
      % Mesh repairs - Not yet implemented in FT
      % % Decimate
      % [bnd(1).pnt, bnd(1).tri] = meshresample(bnd(1).pnt, bnd(1).tri, 1000/size(bnd(1).pnt,1));
      % [bnd(2).pnt, bnd(2).tri] = meshresample(bnd(2).pnt, bnd(2).tri, 2000/size(bnd(2).pnt,1));
      % [bnd(3).pnt, bnd(3).tri] = meshresample(bnd(3).pnt, bnd(3).tri, 3000/size(bnd(3).pnt,1));
      % [bnd(4).pnt, bnd(4).tri] = meshresample(bnd(4).pnt, bnd(4).tri, 3000/size(bnd(4).pnt,1));
      
      % Check and repair individual meshes using iso2mesh
      for bb = 1:length(bnd)
        [bnd(bb).pnt, bnd(bb).tri] = meshcheckrepair(bnd(bb).pnt, bnd(bb).tri, 'dup');
        [bnd(bb).pnt, bnd(bb).tri] = meshcheckrepair(bnd(bb).pnt, bnd(bb).tri, 'isolated');
        [bnd(bb).pnt, bnd(bb).tri] = meshcheckrepair(bnd(bb).pnt, bnd(bb).tri, 'deep');
        [bnd(bb).pnt, bnd(bb).tri] = meshcheckrepair(bnd(bb).pnt, bnd(bb).tri, 'meshfix');
      end
      
      % Ensure no overlaps
      cd([fwd 'private'])
      bnd = decouplesurf(bnd);    % decouplesurf is an unimplemented subfunction temporarily stashed in prepare_mesh_segmentation
      
      cfg        = [];
      cfg.method ='dipoli';
      voldipolimm     = ft_prepare_headmodel(cfg, ft_convert_units(bnd,'mm'));
    end
    cd([mdir sub{ii} ])
    voldipolimm.mat;
    disp('got 2')
    bemoption=1;
    save vol1.mat voldipolimm bemoption
    save bnd1.mat bnd bemoption
  end
  
  % no point in this, as it's just based on 'subject01' from FT.
  %   try
  %     % try existing segmented MRI from FT
  %     load([fwd 'data/segmentedmri.mat'])
  %     mri01 = ft_read_mri('Subject01.mri');
  %
  %     segmentendmri.anatomy=mri01.anatomy;
  %
  %
  %     cfg=[];
  %     cfg.tissue={'brain','skull','scalp'};
  %     cfg.numvertices = [3000 2000 1000];
  %     cfg.method='projectmesh';
  %     bnd=ft_prepare_mesh(cfg,segmentedmri);
  %
  %     cfg        = [];
  %     cfg.method ='dipoli';
  %     voldipolimm     = ft_prepare_headmodel(cfg, ft_convert_units(bnd,'mm'));
  %
  %   end
  
  % Option 2 (do this for everyone in case we decide to use it for
  % everyone)
  cd([mdir sub{ii} ])
  % use existing meshes from existing BEM but warp shape to individual
  % head, then recompute vol
  tmp=load('standard_bem');
  cd([mdir sub{ii} ])
  params=load('params_sn');
  bnd=tmp.vol.bnd;
  bnd(1).pnt=ft_warp_apply(params,tmp.vol.bnd(1).pnt,'sn2individual');
  bnd(2).pnt=ft_warp_apply(params,tmp.vol.bnd(2).pnt,'sn2individual');
  bnd(3).pnt=ft_warp_apply(params,tmp.vol.bnd(3).pnt,'sn2individual');
  
  cfg        = [];
  cfg.method ='dipoli';
  voldipolimm     = ft_prepare_headmodel(cfg, ft_convert_units(bnd,'mm'));
  
  voldipolimm.mat;
  bemoption=2;
  
  save vol.mat voldipolimm bemoption
  save bnd.mat bnd bemoption
  
  %   ft_plot_mesh(bnd(3), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
  %   hold on;
  %   ft_plot_mesh(bnd_ind(1),'edgecolor','none','facecolor',[0.4 0.6 0.4]);
  
  % Option 3
  %       % realistically, bemoption 2 works for everyone, so this is never reached.
  %       disp(ME.message);
  %       % just use existing BEM vol but with warped pnts.
  %       tmp=load('standard_bem');
  %       cd([mdir sub{ii} ])
  %       params=load('params_sn');
  %       bnd=tmp.vol.bnd;
  %       bnd(1).pnt=ft_warp_apply(params,tmp.vol.bnd(1).pnt,'sn2individual');
  %       bnd(2).pnt=ft_warp_apply(params,tmp.vol.bnd(2).pnt,'sn2individual');
  %       bnd(3).pnt=ft_warp_apply(params,tmp.vol.bnd(3).pnt,'sn2individual');
  %
  %       voldipolimm=tmp.vol;
  %       voldipolimm.bnd=bnd;
  %
  %       bemoption=3;
  %       save vol3.mat voldipolimm bemoption
  %       save bnd3.mat bnd bemoption
  
  
  
  
  %   ft_plot_mesh(bnd(3), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
  %   hold on;
  %   ft_plot_mesh(bnd(2),'edgecolor','none','facealpha',0.4);
  %   hold on;
  %   ft_plot_mesh(bnd(1),'edgecolor','none','facecolor',[0.4 0.6 0.4]);
  %
  % ft_plot_mesh(bnd(3), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
  % hold on;
  % ft_plot_mesh(bnd_ind(3),'edgecolor','none','facecolor',[0.4 0.6 0.4]);
  
  %   for ii=iiuse
  %   cfg        = [];
  %   cfg.method ='concentricspheres';
  %   volcs        = ft_prepare_headmodel(cfg, ft_convert_units(bnd,'mm'));
  %
  %   % dipoli requires Linux or Mac to run a mex file (not Windows)
  %   cfg        = [];
  %   cfg.method ='bemcp';
  %   volbemcpmm        = ft_prepare_headmodel(cfg, ft_convert_units(bnd,'mm'));
  %
  %   cfg        = [];
  %   cfg.method ='bemcp';
  %   volbemcpindmm        = ft_prepare_headmodel(cfg, ft_convert_units(bnd_ind,'mm'));
  %
  %   cfg        = [];
  %   cfg.method ='bemcp';
  %   volbemcpstdmm        = ft_prepare_headmodel(cfg, ft_convert_units(tmp.vol.bnd,'mm'));
  %
  %   cfg        = [];
  %   cfg.method ='bemcp';
  %   volbemcpcm        = ft_prepare_headmodel(cfg, ft_convert_units(bnd,'cm'));
  %
  %   cfg        = [];
  %   cfg.method ='bemcp';
  %   volbemcpm     = ft_prepare_headmodel(cfg, ft_convert_units(bnd,'m'));
  %
  %     cfg        = [];
  %     cfg.method ='dipoli';
  %     voldipolimm     = ft_prepare_headmodel(cfg, ft_convert_units(bnd_ind,'mm'));
  %
  %     ft_plot_mesh(vol.bnd(3), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
  %     hold on;
  %     ft_plot_mesh(vol.bnd(2),'edgecolor','none','facealpha',0.4);
  %     hold on;
  %     ft_plot_mesh(vol.bnd(1),'edgecolor','none','facecolor',[0.4 0.6 0.4]);
  %
  %     save vol vol*
  %
  % cd(esdir)
end

% Note:  which bemoption?
% bemoption 1: [8 9 17 25 27 ]
% bemoption 2: [5 6 7 10 12 13 14 15 16 18 20 21 22 23 24 26 28 29 30 31 32]

% Only needed as temporary hack, since code above didn't look like it does
% now originally.
%
% for ii=[8 9 17 25 27 ]
%   cd([mdir sub{ii} ])
%   !mv vol.mat vol1.mat
%   !mv bnd.mat bnd1.mat
%       tmp=load('standard_bem');
%       params=load('params_sn');
%       bnd=tmp.vol.bnd;
%       bnd(1).pnt=ft_warp_apply(params,tmp.vol.bnd(1).pnt,'sn2individual');
%       bnd(2).pnt=ft_warp_apply(params,tmp.vol.bnd(2).pnt,'sn2individual');
%       bnd(3).pnt=ft_warp_apply(params,tmp.vol.bnd(3).pnt,'sn2individual');
%
%       cfg        = [];
%       cfg.method ='dipoli';
%       voldipolimm     = ft_prepare_headmodel(cfg, ft_convert_units(bnd,'mm'));
%
%       voldipolimm.mat;
%       bemoption=2;
%
%       save vol.mat voldipolimm bemoption
%       save bnd.mat bnd bemoption
% end


%% Read polhemus and convert to elec sens file

if 0 % needed to be done once, but not called again after having done once
  cd(pdir);
  files=dir('e*');
  for ff=1:length(files)
    lff=str2num(files(ff).name(2:3));
    elecfile{lff}=files(ff).name;
  end
  
  for ii=iiuse
    if any(ii==setdiff(16:32,[16 17 20 21 22 23 27])) % even though files exist for 16 onwards, they were still messed up/distorted for some
      elec=ft_read_sens(elecfile{ii});
      poluse=1;
    else
      elec=ft_read_sens('standard_1005.elc');
      poluse=0;
    end
    elec=ft_convert_units(elec,'mm');
    save([esdir 'elec_' sub{ii} '.mat'],'elec','poluse');
    %   ft_plot_sens(elec);
  end % ii
end


%% Following http://fieldtrip.fcdonders.nl/tutorial/headmodel_eeg?s[]=standard&s[]=bem#align_the_electrodes

% for ii=setdiff(iiuse,[18:32])
% for ii=setdiff(iiuse,[5:17 18:23])
for ii=[20 21 22 23 27]
  
  load([esdir 'elec_' sub{ii} '.mat']); % this is either polhemus or standard file
  cd([mdir sub{ii} ])
  try
    load vol1.mat
  catch
    load vol.mat
  end
  
  
  
  
  if poluse==0
    params=load('params_sn');
    
    elec_warp2ind=elec;
    %   elec_warp2ind.chanpos=get_orig_coord2(elec.chanpos, 'params_sn.mat');
    %   elec_warp2ind.elecpos=get_orig_coord2(elec.elecpos, 'params_sn.mat');
    elec_warp2ind.chanpos=ft_warp_apply(params,elec.chanpos, 'sn2individual');
    elec_warp2ind.elecpos=ft_warp_apply(params,elec.elecpos, 'sn2individual');
    elec_use=elec_warp2ind;
    elecuse=0;
    
    %     cfg=[];
    %     cfg.locationcoordinates='head';
    %     cfg.location=elec.chanpos(1,:);
    %     ft_sourceplot(cfg,mrinormwarp);
    %     cfg.location=elec.chanpos(2,:);
    %     ft_sourceplot(cfg,mrinormwarp);
    %     cfg.location=elec.chanpos(3,:);
    %     ft_sourceplot(cfg,mrinormwarp);
    %     % -> check if elec_warp2ind.chanpos(1:3,:)  matches onto fiducial spots on MRI
    %     cfg=[];
    %     cfg.locationcoordinates='head';
    %     cfg.location=elec_warp2ind.chanpos(1,:);
    %     ft_sourceplot(cfg,mrinormlin);
    %     cfg.location=elec_warp2ind.chanpos(2,:);
    %     ft_sourceplot(cfg,mrinormlin);
    %     cfg.location=elec_warp2ind.chanpos(3,:);
    %     ft_sourceplot(cfg,mrinormlin);
  elseif poluse==1
    load mrinorm.mat
    
    elec.chanpos=elec.chanpos([1:64 68:70],:);
    elec.elecpos=elec.elecpos([1:64 68:70],:);
    elec.label=elec.label([1:64 68:70]);
    
    %     % presumably this alignment is very far off.
    %     cfg=[];
    %     cfg.locationcoordinates='head';
    %     cfg.location=elec.chanpos(65,:);
    %     ft_sourceplot(cfg,mrinormlin);
    %     cfg.location=elec.chanpos(66,:);
    %     ft_sourceplot(cfg,mrinormlin);
    %     cfg.location=elec.chanpos(67,:);
    %     ft_sourceplot(cfg,mrinormlin);
    
    % Mark fiducials interactively on mrinormlin
    try % maybe it was already done for this participant
      load('fiducial_linvox');
    catch
      cfg=[];
      cfg.method='interactive';
      cfg.coordsys='ctf'; %This doesn't matter if 'ctf' or 'spm' if just extracting voxel-based coord of fiducials
      mrinormlinfid = ft_volumerealign(cfg, mrinormlin);
      fiducial_linvox=mrinormlinfid.cfg.fiducial;
      save('fiducial_linvox.mat','fiducial_linvox'); % to be used later with ft_warp_apply
    end
    nas=fiducial_linvox.nas;
    lpa=fiducial_linvox.lpa;
    rpa=fiducial_linvox.rpa;
    
    % use original mri transform
    nas=ft_warp_apply(mrinormlin.transform,nas, 'homogenous');
    lpa=ft_warp_apply(mrinormlin.transform,lpa, 'homogenous');
    rpa=ft_warp_apply(mrinormlin.transform,rpa, 'homogenous');
    
    % create a structure similar to a template set of electrodes
    fid.chanpos       = [nas; lpa; rpa];       % ctf-coordinates of fiducials
    fid.label         = elec.label(65:67);     % same labels as in elec
    fid.unit          = mrinormlin.unit;       % same units as mri
    
    % MNI space is consistently 'too big' compared to native space.  Scale
    % up the electrode cloud first before alignment
    elec_scale=elec;
    elec_scale.chanpos=1.09*elec.chanpos;
    elec_scale.elecpos=1.09*elec.elecpos;
    
    % alignment
    cfg               = [];
    cfg.method        = 'fiducial';
    cfg.template      = fid;                   % see above
    cfg.fiducial      = elec.label(65:67);  % labels of fiducials in fid and in elec
    cfg.warp          = 'rigidbody';
    cfg.elec          = elec;
    elec_aligned_notsc= ft_electroderealign(cfg);
    cfg.elec          = elec_scale;
    elec_aligned      = ft_electroderealign(cfg);
    
    
    % -> check if elec_warp2ind.chanpos(1:3,:)  matches onto fiducial spots on MRI
    cfg=[];
    cfg.locationcoordinates='head';
    cfg.location=elec_aligned.chanpos(65,:);
    ft_sourceplot(cfg,mrinormlin);
    cfg.location=elec_aligned.chanpos(66,:);
    ft_sourceplot(cfg,mrinormlin);
    cfg.location=elec_aligned.chanpos(67,:);
    ft_sourceplot(cfg,mrinormlin);
    keyboard
    close all
    
    figure;
    % head surface (scalp)
    ft_plot_mesh(voldipolimm.bnd(1), 'edgecolor','none','facealpha',0.8,'facecolor',[0.6 0.6 0.8]);
    hold on;
    % electrodes
    ft_plot_sens(elec_aligned,'style', 'sk');
    
    
    %     a=1.1;
    %     elec_manrescale=elec_aligned;
    %     elec_warp2ind.chanpos=ft_warp_apply(params,elec.chanpos, 'homogenous');
    %     elec_warp2ind.elecpos=ft_warp_apply(params,elec.elecpos, 'homogenous');
    %     elec_use=elec_warp2ind;
    %
    
    % do a bit of final registration (ideally only shift/rotation i.e. rigidbody)
    cfg               = [];
    cfg.method        = 'interactive';
    cfg.elec          = elec_aligned;
    cfg.headshape     = voldipolimm.bnd(1);
    elec_interactive  = ft_electroderealign(cfg);
    
    figure;
    % head surface (scalp)
    ft_plot_mesh(voldipolimm.bnd(1), 'edgecolor','none','facealpha',0.8,'facecolor',[0.6 0.6 0.8]);
    hold on;
    % electrodes
    ft_plot_sens(elec_interactive,'style', 'sk');
    
    cfg = [];
    cfg.method='template';
    cfg.elec=elec_interactive;
    cfg.headshape=voldipolimm.bnd(1);
    cfg.warp='rigidbody';
    elec_headshape1=ft_electroderealign(cfg);
    
    figure;
    % head surface (scalp)
    ft_plot_mesh(voldipolimm.bnd(1), 'edgecolor','none','facealpha',0.8,'facecolor',[0.6 0.6 0.8]);
    hold on;
    % electrodes
    ft_plot_sens(elec_headshape1,'style', 'sk');
    
    
    
    
    
    %     cfg = [];
    %     cfg.method='template';
    %     cfg.warp='globalrescale';
    %     cfg.elec=elec_aligned;
    %     cfg.template='standard_1005.elc';
    %     elec_templatelin=ft_electroderealign(cfg);
    %
    %     figure;
    %     % head surface (scalp)
    %     ft_plot_mesh(voldipolimm.bnd(1), 'edgecolor','none','facealpha',0.8,'facecolor',[0.6 0.6 0.8]);
    %     hold on;
    %     % electrodes
    %     ft_plot_sens(elec_templatelin,'style', 'sk');
    %
    %     cfg               = [];
    %     cfg.method        = 'interactive';
    %     cfg.elec          = elec_templatelin;
    %     cfg.headshape     = voldipolimm.bnd(1);
    %     elec_interactive  = ft_electroderealign(cfg);
    %
    %     cfg = [];
    %     cfg.method='template';
    %     cfg.elec=elec_interactive;
    %     cfg.headshape=voldipolimm.bnd(1);
    %     cfg.warp='rigidbody';
    %     elec_headshape3=ft_electroderealign(cfg);
    %
    %     cfg = [];
    %     cfg.method='template';
    %     cfg.elec=elec_aligned;
    %     cfg.headshape=voldipolimm.bnd(1);
    %     cfg.warp='rigidbody';
    %     elec_headshape1=ft_electroderealign(cfg);
    %     cfg.warp='globalrescale';
    %     elec_headshape2=ft_electroderealign(cfg);
    %     cfg.warp='traditional';
    %     elec_headshape3=ft_electroderealign(cfg);
    %
    %     figure;
    %     % head surface (scalp)
    %     ft_plot_mesh(voldipolimm.bnd(1), 'edgecolor','none','facealpha',0.8,'facecolor',[0.6 0.6 0.8]);
    %     hold on;
    %     % electrodes
    %     ft_plot_sens(elec_headshape3,'style', 'sk');
    %     hold on;
    %     ft_plot_sens(elec_headshape1,'style', 'sk');
    
    
    elecuse=inputdlg('Use elec_aligned (1) or elec_interactive (2) or elec_headshape1 (3)?');
    if strcmp(elecuse,'1')
      elec_use=elec_aligned;
    elseif strcmp(elecuse,'2')
      elec_use=elec_interactive;
    elseif strcmp(elecuse,'3')
      elec_use=elec_headshape1;
    else
      error('wrong answer to question')
    end
    
  end
  
  %   figure;
  %   % head surface (scalp)
  %   ft_plot_mesh(voldipolimm.bnd(1), 'edgecolor','none','facealpha',0.8,'facecolor',[0.6 0.6 0.8]);
  %   hold on;
  %   % electrodes
  %   ft_plot_sens(elec,'style', 'sk');
  %
  %   figure;
  %   % head surface (scalp)
  %   ft_plot_mesh(voldipolimm.bnd(1), 'edgecolor','none','facealpha',0.8,'facecolor',[0.6 0.6 0.8]);
  %   hold on;
  %   % electrodes
  %   ft_plot_sens(elec_use,'style', 'sk');
  
  keyboard
  close all
  
  % any still need work? ii=[8 9 17  ]
  if exist('elec.mat','file')
    !mv elec.mat elec_bak.mat
  end
%   save elec.mat elec_use elec_aligned elec_interactive elec_headshape1 elecuse
  save elec.mat elec*
  
end

% %% Following 'Align EEG electrode positions to MRI and BEM headmodel' on wiki
% % http://fieldtrip.fcdonders.nl/example/align_eeg_electrode_positions_to_bem_headmodel
%
% for ii=iiuse
%
%   load([esdir 'elec_' sub{ii} '.mat']);
%
%
%   cd([mdir sub{ii} ])
%   mri_orig_name=[mdir sub{ii} '\' sub{ii} 'struct.nii'   ];
%   mri = ft_read_mri(mri_orig_name);
%
%
%   cfg=[];
%   cfg.method='interactive';
%   cfg.coordsys='ctf'; %This doesn't matter if 'ctf' or 'spm' if just extracting voxel-based coord of fiducials
%   mrire=ft_volumerealign(cfg,mri);
%   %       cfg=[];
%   %       cfg.method='fiducial';
%   %       cfg.fiducial=mrire.cfg.fiducial;
%   %       cfg.coordsys='spm'; %But why, when mrinorml is in 'spm'?
%   %       mrirespm=ft_volumerealign(cfg,mri);
%   fiducial_vox=mrire.cfg.fiducial;
%
%   % transform the voxel indices to headcoordinates in mm; use original mri.transform
%   head_nas= ft_warp_apply(mri.transform, fiducial_vox.nas); %, 'homogenous'); % nasion
%   head_lpa= ft_warp_apply(mri.transform, fiducial_vox.lpa); %, 'homogenous'); % Left preauricular
%   head_rpa= ft_warp_apply(mri.transform, fiducial_vox.rpa); %, 'homogenous'); % Right preauricular
%
%
%
%   cfg = [];
%   cfg.downsample = 2;
%   seg = ft_volumesegment(cfg, mrinormlre);
%
%   cfg = [];
%   cfg.method = 'singleshell';
%   vol = ft_prepare_headmodel(cfg, seg);
%
%
% end % ii
%

%% Fixing accidental channel label typo

cd(esdir)
for ii=intersect(setdiff(16:32,[16 17 20 21 22 23 27]),iiSuse)
  cd([mdir sub{ii} ]);
  load elec.mat
  elec_use.label{match_str(elec_use.label,'FPz')}='Fpz';
  save elec.mat elec*  
end


%% Get and use template grid, warped to individual MRI
% http://fieldtrip.fcdonders.nl/example/create_single-subject_grids_in_individual_head_space_that_are_all_aligned_in_mni_space

load standard_sourcemodel3d8mm;
% ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:));

sourcemodel=ft_convert_units(sourcemodel,'mm');


for ii=iiuse
  
  cd([mdir sub{ii} ])
  load('mrinorm');
  cfg = [];
  cfg.grid.warpmni   = 'yes';
  cfg.grid.template  = sourcemodel;
  cfg.grid.nonlinear = 'yes'; % use non-linear normalization
  cfg.grid.unit      = 'mm';
  cfg.mri            = mrinormlin;
  grid               = ft_prepare_sourcemodel(cfg);
  
  load vol
  figure; % plot skull boundary
  ft_plot_mesh(voldipolimm.bnd(2), 'edgecolor', 'none'); alpha 0.4;
  ft_plot_mesh(grid.pos(grid.inside,:));
  
  keyboard
  close all
  
  save grid.mat grid
  
  %   segmentedmri.anatomy   = mrinormlin.anatomy;
  %
  %   cfg = [];
  %   cfg.funparameter = 'brain';
  %   ft_sourceplot(cfg,segmentedmri); %segmented gray matter on top
  %   cfg.funparameter = 'skull';
  %   ft_sourceplot(cfg,segmentedmri); %segmented white matter on top
  %   cfg.funparameter = 'scalp';
  %   ft_sourceplot(cfg,segmentedmri); %segmented csf matter on top
  
end

%% Finally, compute leadfield

% load('standard_bem');
% voldipolimm_standard=vol;
ii=5;
tfile=dir([edir sub{ii} '/tlock_diff*']);
load([edir sub{ii} '/' tfile(2).name],'tlock_tacAll')

% for ii=setdiff(iiuse,[5])
% for ii=[20 21 22 23 27]
for ii=intersect(setdiff(16:32,[16 17 20 21 22 23 27]),iiSuse)
  
  cd([mdir sub{ii} ]);
  load vol.mat
  load grid.mat
  load elec.mat
  %   cd([edir sub{ii} ]);
  
  %   cfg=[];
  %   cfg.vol=vol;
  %   cfg.grid=grid;
  %   cfg.elec=elec_use;
  %   cfg.channel=tlock_tacAll{2,10}.label;
  %   lf=ft_prepare_leadfield(cfg);
  %
  %   gridcm=ft_convert_units(grid,'cm');
  %   eleccm=ft_convert_units(elec_use,'cm');
  %   cfg=[];
  %   cfg.vol=volbemcpcm;
  %   cfg.grid=gridcm;
  %   cfg.elec=eleccm;
  %   cfg.channel=tlock_tacAll{2,10}.label;
  %   lfcm=ft_prepare_leadfield(cfg);
  %
  %   gridm=ft_convert_units(grid,'m');
  %   elecm=ft_convert_units(elec_use,'m');
  %   cfg=[];
  %   cfg.vol=volbemcpm;
  %   cfg.grid=gridm;
  %   cfg.elec=elecm;
  %   cfg.channel=tlock_tacAll{2,10}.label;
  %   lfm=ft_prepare_leadfield(cfg);
  %
  %   gridmm=ft_convert_units(grid,'mm');
  %   elecmm=ft_convert_units(elec_use,'mm');
  %   cfg=[];
  %   cfg.vol=volbemcpmm;
  %   cfg.grid=gridmm;
  %   cfg.elec=elecmm;
  %   cfg.channel=tlock_tacAll{2,10}.label;
  %   lfmm=ft_prepare_leadfield(cfg);
  %
  %   gridmm=ft_convert_units(grid,'mm');
  %   elecmm=ft_convert_units(elec_use,'mm');
  %   cfg=[];
  %   cfg.vol=volbemcpstdmm;
  %   cfg.grid=gridmm;
  %   cfg.elec=elecmm;
  %   cfg.channel=tlock_tacAll{2,10}.label;
  %   lfsmm=ft_prepare_leadfield(cfg);
  %
  %   gridmm=ft_convert_units(grid,'mm');
  %   elecmm=ft_convert_units(elec_use,'mm');
  %   cfg=[];
  %   cfg.vol=volbemcpindmm;
  %   cfg.grid=gridmm;
  %   cfg.elec=elecmm;
  %   cfg.channel=tlock_tacAll{2,10}.label;
  %   lfimm=ft_prepare_leadfield(cfg);
  %
  %   gridmm=ft_convert_units(grid,'mm');
  %   elecmm=ft_convert_units(elec_use,'mm');
  %   cfg=[];
  %   cfg.vol=voldipolimm_standard;
  %   cfg.grid=gridmm;
  %   cfg.elec=elecmm;
  %   cfg.channel=tlock_tacAll{2,10}.label;
  %   lfdmm=ft_prepare_leadfield(cfg);
  
  cfg=[];
  cfg.vol=voldipolimm;
  cfg.grid=ft_convert_units(grid,'mm');
  cfg.elec=ft_convert_units(elec_use,'mm');
  cfg.channel=tlock_tacAll{2,10}.label;
  lf=ft_prepare_leadfield(cfg);
  
  %   gridmm=ft_convert_units(grid,'mm');
  %   elecmm=ft_convert_units(elec_use,'mm');
  %   cfg=[];
  %   cfg.vol=volcs;
  %   cfg.grid=gridmm;
  %   cfg.elec=elecmm;
  %   cfg.channel=tlock_tacAll{2,10}.label;
  %   lfcmm=ft_prepare_leadfield(cfg);
  
  tlock=[];
  tlock.time=[1 2 3];
  tlock.avg=lf.leadfield{dsearchn(grid.pos,[-20 0 50])};
  tlock.label=lf.cfg.channel;
  tlock.dimord='chan_time';
  
  cfg=[];
  cfg.layout='elec1010.lay';
  cfg.xlim=[0.9 1.1];
  figure;
  ft_topoplotER(cfg,tlock);
  cfg.xlim=[1.9 2.1];
  figure;
  ft_topoplotER(cfg,tlock);
  cfg.xlim=[2.9 3.1];
  figure;
  ft_topoplotER(cfg,tlock);
  
  keyboard
  close all
  save([esdir 'lf_' sub{ii} '.mat'],'lf');

  % redo/inspect: 20 21 22 23 ~27 ~28 ~29
  
end

% these the lf was 'wobbly' so maybe electrode positions not right
for ii=[20 21 22 23 27 28 29]
  cd([mdir sub{ii} ]);
  load vol.mat
  load elec.mat

    figure;
    % head surface (scalp)
    ft_plot_mesh(voldipolimm.bnd(1), 'edgecolor','none','facealpha',0.8,'facecolor',[0.6 0.6 0.8]);
    hold on;
    % electrodes
    ft_plot_sens(elec_interactive,'style', 'sk');
    
    % bail on polhemus, use std_elec:  20 21 22 23 27
    % polhemus ok: 28 29
  
end

% scroll through to check all
cd(esdir)
for ii=iiuse
  load(['lf_' sub{ii} '.mat']);
  tlock=[];
  tlock.time=[1 2 3];
  tlock.avg=lf.leadfield{dsearchn(lf.pos,[-20 0 50])};
  tlock.label=lf.cfg.channel;
  tlock.dimord='chan_time';
  
  cfg=[];
  cfg.layout='elec1010.lay';
  cfg.xlim=[0.9 1.1];
  figure;
  ft_topoplotER(cfg,tlock);
  cfg.xlim=[1.9 2.1];
  figure;
  ft_topoplotER(cfg,tlock);
  cfg.xlim=[2.9 3.1];
  figure;
  ft_topoplotER(cfg,tlock);
  keyboard
  close all
  
end

% tlockcm=[];
% tlockcm.time=[1 2 3];
% tlockcm.avg=lfcm.leadfield{dsearchn(grid.pos,[-20 0 50])};
% tlockcm.label=lfcm.cfg.channel;
% tlockcm.dimord='chan_time';
%
% tlockm=[];
% tlockm.time=[1 2 3];
% tlockm.avg=lfm.leadfield{dsearchn(grid.pos,[-20 0 50])};
% tlockm.label=lfm.cfg.channel;
% tlockm.dimord='chan_time';
%
% tlockmm=[];
% tlockmm.time=[1 2 3];
% tlockmm.avg=lfmm.leadfield{dsearchn(grid.pos,[-20 0 50])};
% tlockmm.label=lfmm.cfg.channel;
% tlockmm.dimord='chan_time';
%
% tlocksmm=[];
% tlocksmm.time=[1 2 3];
% tlocksmm.avg=lfsmm.leadfield{dsearchn(grid.pos,[-20 0 50])};
% tlocksmm.label=lfsmm.cfg.channel;
% tlocksmm.dimord='chan_time';
%
% tlockimm=[];
% tlockimm.time=[1 2 3];
% tlockimm.avg=lfimm.leadfield{dsearchn(grid.pos,[-20 0 50])};
% tlockimm.label=lfimm.cfg.channel;
% tlockimm.dimord='chan_time';
%
% tlockdmm=[];
% tlockdmm.time=[1 2 3];
% tlockdmm.avg=lfdmm.leadfield{dsearchn(grid.pos,[-20 0 50])};
% tlockdmm.label=lfdmm.cfg.channel;
% tlockdmm.dimord='chan_time';
%
% tlockdmmn=[];
% tlockdmmn.time=[1 2 3];
% tlockdmmn.avg=lfdmm_new.leadfield{dsearchn(grid.pos,[-20 0 50])};
% tlockdmmn.label=lfdmm_new.cfg.channel;
% tlockdmmn.dimord='chan_time';
%
% tlockcmm=[];
% tlockcmm.time=[1 2 3];
% tlockcmm.avg=lfcmm.leadfield{dsearchn(grid.pos,[-20 0 50])};
% tlockcmm.label=lfcmm.cfg.channel;
% tlockcmm.dimord='chan_time';

% cfg=[];
% cfg.layout='elec1010.lay';
% cfg.xlim=[0.9 1.1];
% figure;
% ft_topoplotER(cfg,tlockcm);
% cfg.xlim=[1.9 2.1];
% figure;
% ft_topoplotER(cfg,tlockcm);
% cfg.xlim=[2.9 3.1];
% figure;
% ft_topoplotER(cfg,tlockcm);
%
% cfg=[];
% cfg.layout='elec1010.lay';
% cfg.xlim=[0.9 1.1];
% figure;
% ft_topoplotER(cfg,tlockm);
% cfg.xlim=[1.9 2.1];
% figure;
% ft_topoplotER(cfg,tlockm);
% cfg.xlim=[2.9 3.1];
% figure;
% ft_topoplotER(cfg,tlockm);
%
% cfg=[];
% cfg.layout='elec1010.lay';
% cfg.xlim=[0.9 1.1];
% figure;
% ft_topoplotER(cfg,tlockmm);
% cfg.xlim=[1.9 2.1];
% figure;
% ft_topoplotER(cfg,tlockmm);
% cfg.xlim=[2.9 3.1];
% figure;
% ft_topoplotER(cfg,tlockmm);
%
% cfg=[];
% cfg.layout='elec1010.lay';
% cfg.xlim=[0.9 1.1];
% figure;
% ft_topoplotER(cfg,tlocksmm);
% cfg.xlim=[1.9 2.1];
% figure;
% ft_topoplotER(cfg,tlocksmm);
% cfg.xlim=[2.9 3.1];
% figure;
% ft_topoplotER(cfg,tlocksmm);
%
% cfg=[];
% cfg.layout='elec1010.lay';
% cfg.xlim=[0.9 1.1];
% figure;
% ft_topoplotER(cfg,tlockimm);
% cfg.xlim=[1.9 2.1];
% figure;
% ft_topoplotER(cfg,tlockimm);
% cfg.xlim=[2.9 3.1];
% figure;
% ft_topoplotER(cfg,tlockimm);
%
% cfg=[];
% cfg.layout='elec1010.lay';
% cfg.xlim=[0.9 1.1];
% figure;
% ft_topoplotER(cfg,tlockdmm); % works
% cfg.xlim=[1.9 2.1];
% figure;
% ft_topoplotER(cfg,tlockdmm);
% cfg.xlim=[2.9 3.1];
% figure;
% ft_topoplotER(cfg,tlockdmm);
%
% cfg=[];
% cfg.layout='elec1010.lay';
% cfg.xlim=[0.9 1.1];
% figure;
% ft_topoplotER(cfg,tlockdmmn); %
% cfg.xlim=[1.9 2.1];
% figure;
% ft_topoplotER(cfg,tlockdmmn);
% cfg.xlim=[2.9 3.1];
% figure;
% ft_topoplotER(cfg,tlockdmmn);
%
% cfg=[];
% cfg.layout='elec1010.lay';
% cfg.xlim=[0.9 1.1];
% figure;
% ft_topoplotER(cfg,tlockcmm); % workds
% cfg.xlim=[1.9 2.1];
% figure;
% ft_topoplotER(cfg,tlockcmm);
% cfg.xlim=[2.9 3.1];
% figure;
% ft_topoplotER(cfg,tlockcmm);









