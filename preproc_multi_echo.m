%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   MULTI - ECHO    %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

addpath('MB');

main_dir = '/network/lustre/iss02/cenir/analyse/irm/users/benoit.beranger/SOUNDfMRI/nifti';

e = exam(main_dir, 'SOUND'); % all subjects with multi-echo


%% Common

e.addSerie('T1w$', 'anat_T1' );
e.getSerie('anat').addVolume('^v_.*nii','v',1);


%% Pilotes

e_tmp = e.getExam('Pilote01');

% Func
run_list = {
    'Staircase'
    'Block_0_A'
    'Block_1_A'
    'Block_2_A'
    'Block_3_A'
    'Block_4_A'
    'Block_5_A'
    'Block_6_A'
    };
for r = 1 : length(run_list)
    run_name = run_list{r};
    e_tmp.addSerie([run_name           '$'], ['run_' run_name], 1);
    e_tmp.addSerie([run_name '_PhysioLog$'], ['phy_' run_name], 1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e_tmp = e.getExam('Pilote02');

% Func
run_list = {
    'Staircase'
    'Block_0_A'
    'Block_1_A'
    'Block_2_A'
    'Block_3_A'
    'Block_4_A'
    'Block_5_A'
    };
for r = 1 : length(run_list)
    run_name = run_list{r};
    e_tmp.addSerie([run_name           '$'], ['run_' run_name], 1);
    e_tmp.addSerie([run_name '_PhysioLog$'], ['phy_' run_name], 1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e_tmp = e.getExam('Pilote03');

% Func
run_list = {
    'Staircase'
    'Block_0_A'
    'Block_1_A'
    'Block_2_A'
    'Block_3_A'
    };
for r = 1 : length(run_list)
    run_name = run_list{r};
    e_tmp.addSerie([run_name           '$'], ['run_' run_name], 1);
    e_tmp.addSerie([run_name '_PhysioLog$'], ['phy_' run_name], 1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e_tmp = e.getExam('Pilote04');

% Func
run_list = {
    'Staircase'
    'Block_4_A'
    'Block_5_A'
    'Block_6_A'
    'Block_7_A'
    'Block_8_A'
    'Block_0_P'
    'Block_1_P'
    };
for r = 1 : length(run_list)
    run_name = run_list{r};
    e_tmp.addSerie([run_name           '$'], ['run_' run_name], 1);
    e_tmp.addSerie([run_name '_PhysioLog$'], ['phy_' run_name], 1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e_tmp = e.getExam('DEV2_039_01_SOUNDFMRI_Pilote05');

% Func
run_list = {
    'Staircase'
    'Block_0_A'
    'Block_1_A'
    'Block_2_A'
    'Block_3_A'
    'Block_4_A'
    'Block_5_A'
    'Block_6_A'
    };
for r = 1 : length(run_list)
    run_name = run_list{r};
    e_tmp.addSerie([run_name           '$'], ['run_' run_name], 1);
    e_tmp.addSerie([run_name '_PhysioLog$'], ['phy_' run_name], 1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e_tmp = e.getExam('DEV2_039_02_SOUNDFMRI_Pilote05');

% Func
run_list = {
    'Staircase'
    'Block_0_P'
    'Block_1_P'
    'Block_2_P'
    'Block_3_P'
    'Block_4_P'
    'Block_5_P'
    'Block_6_P'
    'Block_7_P'
    'Block_8_P'
    'Block_9_P'
    };
for r = 1 : length(run_list)
    run_name = run_list{r};
    e_tmp.addSerie([run_name           '$'], ['run_' run_name], 1);
    e_tmp.addSerie([run_name '_PhysioLog$'], ['phy_' run_name], 1);
end


%% SUJET05

e_SUJET05_ACTIVE = e.getExam('2022_07_28_SOUNDFMRI_SUJET05_ACTIVE');

% Func
run_list = {
    'Staircase'
    'Block_0_A'
    'Block_1_A'
    'Block_2_A'
    'Block_3_A'
    'Block_4_A'
    'Block_5_A'
    'Block_6_A'
    'Block_7_A'
    'Block_8_A'
    'Block_9_P' % +++
    'Block_10_P'% +++
    };
for r = 1 : length(run_list)
    run_name = run_list{r};
    e_SUJET05_ACTIVE.addSerie([run_name           '$'], ['run_' run_name], 1);
    e_SUJET05_ACTIVE.addSerie([run_name '_PhysioLog$'], ['phy_' run_name], 1);
end


e_SUJET05_PASSIVE = e.getExam('2022_07_26_SOUNDFMRI_SUJET05_PASSIVE');

% Func
run_list = {
    'Staircase'
    'Block_0_P'
    'Block_1_P'
    'Block_2_P'
    'Block_3_P'
    'Block_4_P'
    'Block_5_P'
    'Block_6_P'
    'Block_7_P'
    'Block_8_P'
    % 'Block_9_P'
    % 'Block_10_P'
    };
for r = 1 : length(run_list)
    run_name = run_list{r};
    e_SUJET05_PASSIVE.addSerie([run_name           '$'], ['run_' run_name], 1);
    e_SUJET05_PASSIVE.addSerie([run_name '_PhysioLog$'], ['phy_' run_name], 1);
end


%% exception

exception_list = e_SUJET05_ACTIVE + e_SUJET05_PASSIVE;


%% specific

e_passive = e.getExam('PASSIVE');
e_passive = e_passive - exception_list;

% Func
run_list = {
    'Staircase'
    'Block_0_P'
    'Block_1_P'
    'Block_2_P'
    'Block_3_P'
    'Block_4_P'
    'Block_5_P'
    'Block_6_P'
    'Block_7_P'
    'Block_8_P'
    'Block_9_P'
    'Block_10_P'
    };
for r = 1 : length(run_list)
    run_name = run_list{r};
    e_passive.addSerie([run_name           '$'], ['run_' run_name], 1);
    e_passive.addSerie([run_name '_PhysioLog$'], ['phy_' run_name], 1);
end


e_active = e.getExam('ACTIVE');
e_active = e_active - exception_list;

% Func
run_list = {
    'Staircase'
    'Block_0_A'
    'Block_1_A'
    'Block_2_A'
    'Block_3_A'
    'Block_4_A'
    'Block_5_A'
    'Block_6_A'
    'Block_7_A'
    'Block_8_A'
    };
for r = 1 : length(run_list)
    run_name = run_list{r};
    e_active.addSerie([run_name           '$'], ['run_' run_name], 1);
    e_active.addSerie([run_name '_PhysioLog$'], ['phy_' run_name], 1);
end

%% common

e.getSerie('run').addVolume('^v_.*nii$',   'v', 3);
e.getSerie('phy').addPhysio(     'dcm$', 'dcm', 1);

% ep2d_se : distortion correction
e.addSerie('ep2d_se_PA_forward_SBRef$', 'se_forward', 1);
e.addSerie('ep2d_se_PA_reverse_SBRef$', 'se_reverse', 1);
e.getSerie('se_').addVolume('^v_.*nii$', 'v', 1);

e.reorderSeries('name');

% e.explore


%% Cluster ?

CLUSTER = 0;


%% segment cat12 #CAT12->SPM12

anat = e.gser('anat_T1').gvol('^v');

clear par
if CLUSTER
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
else
    par.run = 1;
    par.sge = 0;
end
par.redo    = 0;
par.display = 0;

par.subfolder = 0;

par.GM        = [1 0 1 0]; % warped_space_Unmodulated (wp1*)     / warped_space_modulated (mwp1*)     / native_space (p1*)     / native_space_dartel_import (rp1*)
par.WM        = [1 0 1 0]; %                          (wp2*)     /                        (mwp2*)     /              (p2*)     /                            (rp2*)
par.CSF       = [1 0 1 0]; %                          (wp3*)     /                        (mwp3*)     /              (p3*)     /                            (rp3*)
par.TPMC      = [0 0 0 0]; %                          (wp[456]*) /                        (mwp[456]*) /              (p[456]*) /                            (rp[456]*)   This will create other probalities map (p4 p5 p6)

par.label     = [1 1 0] ;  % native (p0*)  / normalize (wp0*)  / dartel (rp0*)       This will create a label map : p0 = (1 x p1) + (3 x p2) + (1 x p3)
par.bias      = [1 1 0] ;  % native (ms*)  / normalize (wms*)  / dartel (rms*)       This will save the bias field corrected  + SANLM (global) T1
par.las       = [0 0 0] ;  % native (mis*) / normalize (wmis*) / dartel (rmis*)       This will save the bias field corrected  + SANLM (local) T1

par.warp      = [1 1];     % warp fields  : native->template (y_*) / native<-template (iy_*)

par.doSurface = 0;
par.jacobian  = 0;         % write jacobian determinant in normalize space
par.doROI     = 0;         % will compute the volume in each atlas region

job_do_segmentCAT12(anat,par);


%% Sort echos #MATLAB/matvol

clear par
par.run  = 1;
par.fake = 0;
par.sge  = 0;

par.redo = 0;
meinfo = job_sort_echos( e.getSerie('run') , par );


%% job_afni_proc_multi_echo #ANFI/afni_proc.py

clear par
if CLUSTER
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
else
    par.run = 1;
    par.sge = 0;
end
par.fake = 0;
par.redo = 0;

par.seperate = 1;
par.write_nifti = 1;

par.blocks  = {'tshift', 'volreg', 'blip'};
par.blip.forward = e.getSerie('se_forward').getVolume();
par.blip.reverse = e.getSerie('se_reverse').getVolume();

afni_prefix = char(par.blocks); % {'despike', 'tshift', 'volreg'}
afni_prefix = afni_prefix(:,1)';
afni_prefix = fliplr(afni_prefix);   % 'vtd'

afni_subdir = ['afni_' afni_prefix];
par.subdir = afni_subdir;
job_afni_proc_multi_echo( meinfo, par );


%% do_fsl_robust_mask_epi #FSL

fin  = e.getSerie('run').getVolume(['^' afni_prefix 'e1']);

clear par
if CLUSTER
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
else
    par.run = 1;
    par.sge = 0;
end
par.fake  = 0;
par.redo  = 0;
par.fsl_output_format = 'NIFTI_GZ';
do_fsl_robust_mask_epi( fin, par );

% Checkpoint & unzip
par.jobname = 'unzip_and_keep__bet';
e.getSerie('run').getVolume(['^bet_Tmean_' afni_prefix 'e1$']).removeEmpty().unzip_and_keep(par);



%% TEDANA #Python

clear par
if CLUSTER
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
else
    par.run = 1;
    par.sge = 0;
end
par.fake  = 0;
par.redo  = 0;
par.pct = 0;

% cluster
par.walltime = '12:00:00';      % HH:MM:SS
par.mem      = '16G';           % ICA is very memory consuming
par.sge_nb_coeur = 2;           % I dont't know why, but 2 CPU increase the "stability" of the job on the cluster

tedana_subdir = ['tedana0011_' afni_prefix];
job_tedana_0011( meinfo, afni_prefix, tedana_subdir, ['bet_Tmean_' afni_prefix 'e1_mask.nii.gz'], par );


% Checkpoint & unzip
par.jobname = 'unzip_and_keep__tedana';
e.getSerie('run').getVolume('ts_OC').removeEmpty().unzip_and_keep(par);


%% Coregister TEDANA outputs to Anat #SPM12

clear par
if CLUSTER
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
else
    par.run = 1;
    par.sge = 0;
end
par.redo  = 0;
par.type  = 'estimate';

src = e.getSerie('run').removeEmpty().getVolume(['^bet_Tmean_' afni_prefix 'e1$']);
oth = e.getSerie('run').removeEmpty().getVolume('^ts_OC');
ref = e.getSerie('run').removeEmpty().getExam.getSerie('anat_T1').getVolume('^p0');

par.jobname = 'spm_coreg_epi2anat';
job_coregister(src,ref,oth,par);


%% Normalize TEDANA outputs #SPM12

clear par
if CLUSTER
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
else
    par.run = 1;
    par.sge = 0;
end
par.redo = 0;
par.vox = [2.5 2.5 2.5]; % IMPORTANT keep original EPI voxel size
img = e.getSerie('run').getVolume('^ts_OC').removeEmpty();
y   = img.getExam.getSerie('anat_T1').getVolume('^y');
par.jobname = 'spm_normalize_epi';
job_apply_normalize(y,img,par);

% Nomalize Tmean, used later for PhysIO Noise ROI
img = e.getSerie('run').removeEmpty().getVolume(['^bet_Tmean_' afni_prefix 'e1$']);
y   = img.getExam.getSerie('anat_T1').getVolume('^y');
par.jobname = 'spm_normalize_meanepi';
job_apply_normalize(y,img,par);


%%  Smooth TEDANA outputs #SPM12

clear par
if CLUSTER
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
else
    par.run = 1;
    par.sge = 0;
end
par.redo = 0;

img = e.getSerie('run').getVolume('^wts_OC').removeEmpty();

par.smooth   = [5 5 5];
par.prefix   = 's5';
job_smooth(img,par);

par.smooth   = [8 8 8];
par.prefix   = 's8';
job_smooth(img,par);


%% coregister WM & CSF on functionnal (using the warped mean) #SPM12
% This will be used for TAPAS:PhysIO

clear par
if CLUSTER
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
else
    par.run = 1;
    par.sge = 0;
end
ref = e.getSerie('run');
ref = ref(:,1).getVolume(['wbet_Tmean_' afni_prefix 'e1']);
src = e.getSerie('anat_T1').getVolume('^wp2');
oth = e.getSerie('anat_T1').getVolume('^wp3');
par.type = 'estimate_and_write';
par.jobname = 'spm_coreg_WMCSF2wEPI';
job_coregister(src,ref,oth,par);


%% rp afni2spm #matlab/matvol

% input
dfile = e.getSerie('run').getRP('rp_afni').removeEmpty();

% output
output_dir = fullfile( dfile.getSerie().getPath(), tedana_subdir );

% go
job_rp_afni2spm(dfile, output_dir);


%% extract physio from special dicom

% https://github.com/CMRR-C2P/MB

e.getSerie('phy').getPhysio('dcm').extract()

% e.getSerie('phy').getPhysio('phy').check() % takes a bit of time, use it once to verify your data


%% PhysIO nuisance regressor generation #matlab/TAPAS-PhysIO
%% Prepare files

info = e.getSerie('phy').removeEmpty().getPhysio('info');
puls = e.getSerie('phy').removeEmpty().getPhysio('puls');
resp = e.getSerie('phy').removeEmpty().getPhysio('resp');
run  = e.getSerie('run').removeEmpty();

%--------------------------------------------------------------------------
% in AFNI, outside the mask is NaN
% in SPM, outise the mask is 0
% so convert NaN to 0

clear par
if CLUSTER
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
else
    par.run = 1;
    par.sge = 0;
end
job_afni_remove_nan( run.getVolume('^wts_OC'), par );
%--------------------------------------------------------------------------

if CLUSTER
    run.addVolume('tedana','^nwts_OC','nwts_OC',1);
end
volume = run.getVolume('^nwts_OC');

outdir = volume.getDir();

rp = run.getRP('rp_spm');

mask = run.getExam.getSerie('anat').getVolume('^rwp[23]');


%% Prepare job

clear par
%----------------------------------------------------------------------------------------------------------------------------------------------------
% ALWAYS MANDATORY
%----------------------------------------------------------------------------------------------------------------------------------------------------

par.physio   = 1;
par.noiseROI = 1;
par.rp       = 1;

par.TR     = 1.660;
par.nSlice = 60;

par.volume = volume;
par.outdir = outdir;

%----------------------------------------------------------------------------------------------------------------------------------------------------
% Physio
%----------------------------------------------------------------------------------------------------------------------------------------------------

par.physio_Info = info;
par.physio_PULS = puls;
par.physio_RESP = resp;

par.physio_RETROICOR        = 1;
par.physio_HRV              = 1;
par.physio_RVT              = 1;
par.physio_logfiles_vendor  = 'Siemens_Tics'; % Siemens CMRR multiband sequence, only this one is coded yet
par.physio_logfiles_align_scan = 'last';         % 'last' / 'first'
% Determines which scan shall be aligned to which part of the logfile.
% Typically, aligning the last scan to the end of the logfile is beneficial, since start of logfile and scans might be shifted due to pre-scans;
par.physio_slice_to_realign    = 'middle';       % 'first' / 'middle' / 'last' / sliceNumber (integer)
% Slice to which regressors are temporally aligned. Typically the slice where your most important activation is expected.


%----------------------------------------------------------------------------------------------------------------------------------------------------
% noiseROI
%----------------------------------------------------------------------------------------------------------------------------------------------------

par.noiseROI_mask   = mask;
par.noiseROI_volume = volume;

par.noiseROI_thresholds   = [0.95 0.70];     % keep voxels with tissu probabilty >= 95%
par.noiseROI_n_voxel_crop = [2 1];           % crop n voxels in each direction, to avoid partial volume
par.noiseROI_n_components = 10;              % keep n PCA componenets


%----------------------------------------------------------------------------------------------------------------------------------------------------
% Realignment Parameters
%----------------------------------------------------------------------------------------------------------------------------------------------------

par.rp_file = rp;

par.rp_order     = 24;   % can be 6, 12, 24
% 6 = just add rp, 12 = also adds first order derivatives, 24 = also adds first + second order derivatives
par.rp_method    = 'FD'; % 'MAXVAL' / 'FD' / 'DVARS'
par.rp_threshold = +Inf;  % Threshold above which a stick regressor is created for corresponding volume of exceeding value


%----------------------------------------------------------------------------------------------------------------------------------------------------
% Other
%----------------------------------------------------------------------------------------------------------------------------------------------------
par.print_figures = 0; % 0 , 1 , 2 , 3

% classic matvol
if CLUSTER
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
else
    par.run = 1;
    par.sge = 0;
end
par.display  = 0;
par.redo     = 0;

% cluster
par.jobname  = 'spm_physio';
par.walltime = '04:00:00';
par.mem      = '4G';

job_physio_tapas( par );


%% save

save e e

