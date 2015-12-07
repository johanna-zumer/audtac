%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                          header                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;
close all;
% addpath(genpath('C:\r2agui'))
% r2agui;
% NOTE read from DICOM is saver in terms of correctness NIFTI is sometimes
% wrong
addpath('C:\fieldtrip-20140601');
ft_defaults;
addpath(genpath('C:\headmodel_tools'))
addpath(genpath('C:\fieldtrip-20140601\external\freesurfer'));
addpath(genpath('C:\fieldtrip-20140601\external\iso2mesh'));
addpath(genpath('C:\spm8'));
rmpath(genpath('C:\\spm8\external\fieldtrip'));

%% Preface

path = 'C:\data_reinst\Subjects';
sj = 'Sub006';
folder = [path '\' sj '\MRI\headmodeling'];

mri_orig_name = [folder '\' sj '.nii'];
mri_spm_name = [folder '\' sj '_spm.nii'];
elec_orig = [path '\' sj '\Epos\' sj '.elc'];

mri_source = [path '\' sj '\MRI\' sj '.nii'];

if ~exist(folder)
    mkdir(folder)
end
copyfile(mri_source, mri_orig_name); 

cd(folder);

%% lets start by changing our nifti file to spm

mri = ft_read_mri(mri_orig_name);

%realign the original mri now to spm space
cfg = [];
cfg.coordsys = 'spm';
cfg.method = 'interactive';
mri_spm = ft_volumerealign(cfg, mri);

%% reslice the MRI so it can actually be stored in SPM coordsys
cfg              = [];

% ... in case we want to cut because we have to much black space in the mri
%cfg.xrange     = [-78 81]; %, in physical units
cfg.yrange     = [-140 110]; %, in physical units
cfg.zrange     = [-141 108]; %in physical units

% reslice
mri_reslice    = ft_volumereslice(cfg,mri_spm);

% write the new mri into a NIFTI file 
ft_write_mri(mri_spm_name, mri_reslice.anatomy, 'dataformat', 'nifti', 'transform', mri_reslice.transform);

%% Align the electrodes to our new MRI
mri_file = mri_spm_name;
mri = ft_read_mri(mri_file);
%realign the original mri now to ctf space (all we want are the fiducials in voxelspace)

%read the electrode positions
elec = ft_read_sens(elec_orig);   
% bring into ctf system (not really we just want the original fid_coordinates)
cfg = [];
cfg.coordsys = 'ctf';
cfg.method = 'interactive';
mri_ctf = ft_volumerealign(cfg, mri);

% get the fiducial coordinates in voxelspace of the original MRI
vox_Nas = mri_ctf.cfg.fiducial.nas;  % fiducials saved in mri structure
vox_Lpa = mri_ctf.cfg.fiducial.lpa;     
vox_Rpa = mri_ctf.cfg.fiducial.rpa;

% use the original transformation Matrix 
Transform = mri.transform;

% transform the voxel indices of SPM - MRI to SPM - headcoordinates in mm
head_Nas          = ft_warp_apply(Transform, vox_Nas, 'homogenous'); % nasion 
head_Lpa          = ft_warp_apply(Transform, vox_Lpa, 'homogenous'); % Left preauricular
head_Rpa          = ft_warp_apply(Transform, vox_Rpa, 'homogenous'); % Right preauricular

% align electrodes to the new headPositions 
cfg = [];
cfg.method = 'template';
cfg.elec = elec;
cfg.fiducial = elec.label(1:3);
cfg.template.chanpos(1,:) = head_Nas;
cfg.template.chanpos(2,:) = head_Lpa;
cfg.template.chanpos(3,:) = head_Rpa;
cfg.template.label = elec.label(1:3);
cfg.template.unit = 'mm';
elec_aligned = ft_electroderealign(cfg);

save(fullfile(folder,['elec_aligned_' sj '.mat']), 'elec_aligned', '-v7.3');

%% now lets do the segmentation
mri_file = mri_spm_name;
mri = ft_read_mri(mri_file);
% segment here
seg = create_bem_segmentation('t1filename', mri_file,'writeseg','yes');

%fake knowledge
seg.transform = mri.transform;
seg.dim = mri.dim;
seg.unit = mri.unit;
seg.coordsys = 'spm';

%save
save(fullfile(folder,['seg_' sj '.mat']), 'seg', '-v7.3');

%% Plot masks to check
figure, ft_plot_slice(seg.scalp,'location',[88 128 128],'orientation',[0 0 1]), hold on
ft_plot_slice(seg.skull,'location',[88 128 128],'orientation',[0 1 0]);
ft_plot_slice(seg.csf+seg.brain,'location',[88 128 128],'orientation',[1 0 0]);

%% Generate meshes
clear cfg
cfg.tissue = {'scalp', 'skull', 'csf', 'brain'};
cfg.method = 'iso2mesh';    % 'projectmesh';
cfg.numvertices = 10000;    % We'll decimate later
bnd = ft_prepare_mesh(cfg, seg);

% Mesh repairs - Not yet implemented in FT
% Decimate
[bnd(1).pnt, bnd(1).tri] = meshresample(bnd(1).pnt, bnd(1).tri, 1000/size(bnd(1).pnt,1));
[bnd(2).pnt, bnd(2).tri] = meshresample(bnd(2).pnt, bnd(2).tri, 2000/size(bnd(2).pnt,1));
[bnd(3).pnt, bnd(3).tri] = meshresample(bnd(3).pnt, bnd(3).tri, 3000/size(bnd(3).pnt,1));
[bnd(4).pnt, bnd(4).tri] = meshresample(bnd(4).pnt, bnd(4).tri, 3000/size(bnd(4).pnt,1));

% Check and repair individual meshes using iso2mesh
for ii = 1:length(bnd)
    [bnd(ii).pnt, bnd(ii).tri] = meshcheckrepair(bnd(ii).pnt, bnd(ii).tri, 'dup');
    [bnd(ii).pnt, bnd(ii).tri] = meshcheckrepair(bnd(ii).pnt, bnd(ii).tri, 'isolated');
    [bnd(ii).pnt, bnd(ii).tri] = meshcheckrepair(bnd(ii).pnt, bnd(ii).tri, 'deep');
    [bnd(ii).pnt, bnd(ii).tri] = meshcheckrepair(bnd(ii).pnt, bnd(ii).tri, 'meshfix');
end

% Ensure no overlaps
bnd = decouplesurf(bnd);    % decouplesurf is an unimplemented subfunction temporarily stashed in prepare_mesh_segmentation

% Save
save(fullfile(folder,['bnd_' sj '.mat']), 'bnd', '-v7.3');

%% Plot headmodel together with electrodes
figure;
ft_plot_mesh(bnd(1),'facecolor','none'); %scalp
hold;
ft_plot_mesh(bnd(2), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
hold ;
ft_plot_mesh(bnd(3),'edgecolor','none','facealpha',0.4);
hold;
ft_plot_mesh(bnd(4),'edgecolor','none','facecolor',[0.4 0.6 0.4], 'facealpha', 0.7);
hold;

hold;
ft_plot_sens(elec_aligned ,'style', 'b*');

%% Fine tuning
cfg           = [];
cfg.method    = 'interactive';
cfg.elec      = elec_aligned;
cfg.headshape = bnd(1);
elec_fin  = ft_electroderealign(cfg);

%save
save(fullfile(folder,['elec_fin_' sj '.mat']), 'elec_fin', '-v7.3');

%% just one last time...

%plot headmodel together with electrodes
figure;
ft_plot_mesh(bnd(1),'facecolor','none'); %scalp
hold;
ft_plot_mesh(bnd(2), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
hold ;
ft_plot_mesh(bnd(3),'edgecolor','none','facealpha',0.4);
hold;
ft_plot_mesh(bnd(4),'edgecolor','none','facecolor',[0.4 0.6 0.4], 'facealpha', 0.7);
hold;
ft_plot_sens(elec_fin ,'style', 'b*');

%%
copyfile([folder '\elec_fin_' sj '.mat'], [path '\' sj '\Epos']);
if ~exist([path '\' sj '\Source']) mkdir([path '\' sj '\Source']);end
copyfile([folder '\bnd_' sj '.mat'], [path '\' sj '\Source']);
copyfile([folder '\elec_fin_' sj '.mat'], [path '\' sj '\Source']);
copyfile([folder '\' sj '_spm.nii'], [path '\' sj '\Source']);
%% cool views..  Plot meshes on top of MRI
figure;
ft_plot_mesh(bnd(1), 'facealpha', 0.25), hold on
ft_plot_mesh(bnd(2), 'facealpha', 0.25, 'facecolor', 'blue')
ft_plot_mesh(bnd(3), 'facealpha', 0.25, 'facecolor', 'yellow')
ft_plot_mesh(bnd(4), 'facecolor', 'red')
ft_plot_ortho(mri.anatomy,'transform',mri.transform,'style','intersect')


hold;
ft_plot_sens(elec_fin ,'style', 'g*');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% THIS IS THE END....                                           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                       FINALLY                                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                   ... BUT FOR HOW LONG???     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


