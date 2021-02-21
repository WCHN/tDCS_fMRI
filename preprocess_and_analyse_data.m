% Davide Nardo @ MRC-CBU 2020
%
% This script preprocesses and analyses the shared data in the same way as
% described in the Methods and Results sections.
%
% Make sure you set the folder containing the data as your working
% directory in MATLAB before launching this script.
%
% The following pipeline assumes co-registration between functional and
% structural images.

% Requirements:
% * SPM: https://www.fil.ion.ucl.ac.uk/spm/
% * MarsBar: http://marsbar.sourceforge.net/

% SET VARIABLES
clear
root_path = pwd;

spm_path = spm('Dir');
spm('defaults','fmri');
spm_jobman('initcfg');

marsbar_path = fullfile(spm_path,'toolbox','marsbar-0.44'); % <-- change accordingly
addpath(marsbar_path);

folders = {'No_tDCS', 'Sham_tDCS', 'Anodal_tDCS'};
FMdef = fullfile(spm_path,'toolbox','FieldMap','FIL','pm_defaults_Prisma.m');

weight      = 0; % weight vector for std (0 or 1)
TR          = 3.36; % <-- change accordingly
coordinates = [-52 18 16; 21 59 11; 2 -56 51]; % x y z cortical projections for: FC5; FP2; PZ <-- change accordingly
radius      = 40; % sphere radius in mm for ROIs
roi_names   = {'FC5','FP2','PZ'};

outDir      = fullfile(root_path,'ROIs'); % directory for ROIs
spm_mkdir(outDir);
repository  = fullfile(outDir,'resized'); % directory for resized ROIs
spm_mkdir(repository);

% file prefixes and extensions
EPI     = '^f.*\.nii$';
fMAPS   = '^s.*\.nii$';
VDM     = '^v.*\.nii$';
ufEPI   = '^uf.*\.nii$';
rp_file = '^rp.*\.txt$';
anat    = '^s.*\.nii$';

% SEGMENT MPRAGE
anatomical = cellstr(spm_select('FPList',fullfile(root_path,'data','MPRAGE'), anat));
matlabbatch{1}.spm.spatial.preproc.channel.vols = anatomical;
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {fullfile(spm_path, 'tpm', 'TPM.nii,1')};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {fullfile(spm_path, 'tpm', 'TPM.nii,2')};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {fullfile(spm_path, 'tpm', 'TPM.nii,3')};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {fullfile(spm_path, 'tpm', 'TPM.nii,4')};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {fullfile(spm_path, 'tpm', 'TPM.nii,5')};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {fullfile(spm_path, 'tpm', 'TPM.nii,6')};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.vox = NaN;
matlabbatch{1}.spm.spatial.preproc.warp.bb = [NaN NaN NaN;NaN NaN NaN];
spm_jobman('run',matlabbatch);
clear matlabbatch

% CREATE ROIs
% reslice GM to EPI space
GM = spm_select('FPList', fullfile(root_path, 'data', 'MPRAGE'), 'c1');
EPI_img = spm_select('FPList', fullfile(root_path, 'data', folders{1}, 'fMRI'), 'f');
EPI_img = EPI_img(1,:);
copyfile(GM, repository)
GM = spm_select('FPList', repository, 'c1');
Vi = {EPI_img; GM}; % model space; space to change
spm_reslice(Vi, struct('mean', false, 'which', 1, 'prefix', ''));
spm_imcalc(GM,GM,'i1>0.05'); % binarize

% make spheres
for r = 1 : size(coordinates,1)
    current_coord = coordinates(r,:);
    roi_label = roi_names{r};
    sphere = maroi_sphere(struct('centre', current_coord, 'radius', radius));
    output_name = fullfile(outDir, [roi_label '_roi']);
    saveroi(sphere, [output_name '.mat']); % save MarsBaR _roi (.mat) file
    save_as_image(sphere, [output_name '.nii']); % save Nifti (.nii) file
    
    % reslice sphere to EPI space
    ROI_img = spm_select('FPList', outDir, [roi_label '_roi.nii']);
    copyfile(ROI_img, repository)
    ROI_img = spm_select('FPList', repository, [roi_label '_roi.nii']);
    Vi = {EPI_img; ROI_img}; % model space; space to change
    spm_reslice(Vi, struct('mean', false, 'which', 1, 'prefix', ''));
    spm_imcalc(ROI_img,ROI_img,'i1>0.05'); % binarize
    
    % merge ROI & GM
    matlabbatch{1}.spm.util.imcalc.input = {GM;ROI_img};
    matlabbatch{1}.spm.util.imcalc.output = ['GM_' roi_label];
    matlabbatch{1}.spm.util.imcalc.outdir = {repository};
    matlabbatch{1}.spm.util.imcalc.expression = 'i1&i2';
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    spm_jobman('run',matlabbatch);
    clear matlabbatch
end

% PREPROCESS DATA & ANALYSE T-SCORE
for fold = 1 : numel(folders)
    current_folder = fullfile(root_path, 'data', folders{fold});
    
    % fieldmap
    epiDir = cellstr(fullfile(current_folder,'fMRI'));
    datafiles{1} = cellstr(spm_select('FPList', epiDir, EPI));
    fMaps = cellstr(spm_select('FPList',fullfile(current_folder,'FieldMap'), fMAPS));
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase = {fMaps{3}};
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = {fMaps{1}};
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsfile = {FMdef};
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session(1).epi = datafiles{1,1}(1,:);
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 1;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'session';
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 1;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.anat = '';
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 0;
    spm_jobman('run', matlabbatch);
    clear matlabbatch epiDir datafiles fMaps
    
    % realign & unwarp
    epiDir = cellstr(fullfile(current_folder,'fMRI'));
    datafiles{1} = cellstr(spm_select('FPList', epiDir, EPI));
    fMaps = cellstr(spm_select('FPList',fullfile(current_folder,'FieldMap'), VDM));
    matlabbatch{1}.spm.spatial.realignunwarp.data(1).scans = datafiles{1,1};
    matlabbatch{1}.spm.spatial.realignunwarp.data(1).pmscan = {fMaps{1}};
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.quality = 1;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.sep = 4;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.rtm = 0;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.einterp = 2;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.weight = '';
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.jm = 0;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.sot = [];
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.rem = 1;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.noi = 5;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.mask = 1;
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';
    spm_jobman('run', matlabbatch);
    clear matlabbatch epiDir datafiles fMaps
    
    % specify GLM
    output_folder = fullfile(current_folder,'GLM');
    epiDir = fullfile(current_folder,'fMRI');
    datafiles = cellstr(spm_select('FPList', epiDir, ufEPI));
    motion_params = spm_select('FPList', epiDir, rp_file);
    matlabbatch{1}.spm.stats.fmri_spec.dir = {output_folder};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = datafiles;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {motion_params};
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 1;
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    spm_jobman('run',matlabbatch);
    clear matlabbatch epiDir datafiles motion_params
    
    % estimate GLM
    estimate_file = spm_select('FPlist',output_folder,'SPM.mat');
    matlabbatch{1}.spm.stats.fmri_est.spmmat = {estimate_file};
    spm_jobman('run',matlabbatch);
    clear matlabbatch
    
    % compute contrast (T-score of the mean)
    contrast = [zeros(1,6) 1];
    matlabbatch{1}.spm.stats.con.spmmat = {estimate_file};
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 't_score_mean';
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = contrast;
    spm_jobman('run',matlabbatch);
    clear matlabbatch output_folder
end

% COMPUTE & PLOT METRICS
% combine GM & mask image for GM ROI
mask = spm_select('FPList', fullfile(root_path, 'data', folders{1}, 'GLM'), 'mask');
mask_hdr = spm_vol(mask);
mask_img = spm_read_vols(mask_hdr);
GM_hdr = spm_vol(GM);
GM_img = spm_read_vols(GM_hdr);
GM_img = GM_img.*mask_img;
output_hdr = GM_hdr;
output_hdr.fname = fullfile(root_path, 'ROIs', 'resized', 'GM_brain.nii');
spm_write_vol(output_hdr,GM_img);

% set variables
ROIname = {'GM_brain','GM_FC5','GM_FP2','GM_PZ'};
metrics = {'spmT_0001','fpm'};
conditions = {'No','Sham','Anodal'};
ROIlabels = {'GM','FC5','FP2','PZ'};

% extract values
for r = 1 : numel(ROIname)
    ROI = spm_vol(fullfile(root_path, 'ROIs', 'resized', [ROIname{r} '.nii']));
    ROIdata = spm_read_vols(ROI);
    for m = 1 : numel(metrics)
        for c = 1 : numel(folders)
            % Loop over metrics and conditions
            switch m
                case 1
                    data = spm_vol(fullfile(root_path, 'data', folders{c}, 'GLM', [metrics{m} '.nii']));
                case 2
                    data = spm_vol(spm_select('FPList', fullfile(root_path, 'data', folders{c}, 'FieldMap'), metrics{m}));
            end
            summary{r,c,m}.metricVals = spm_summarise(data, ROI);
        end
    end
end

% display
figure
nPoints = 15;
bins(1,:) = linspace(0,1500,nPoints);   % t-score
bins(2,:) = linspace(-50,100,nPoints);  % Frequency, Hz
metricLabel = {'T-score of mean','Frequency (Hz)'};
for r = 1 : numel(ROIname)
    for m = 1 : numel(metrics)
        subplot(numel(metrics),numel(ROIname),(m-1)*numel(ROIname)+r)
        for cond = 1 : numel(folders)
            h = hist(summary{r,cond,m}.metricVals(:), bins(m,:));
            plot(bins(m,:), h, 'LineWidth',2)
            hold on
        end
        hold off; grid minor; set(gca, 'XLim', [bins(m,1) bins(m,end)])
        xlabel(metricLabel{m}); ylabel('Number of Voxels')
    end
end
for r = 1 : numel(ROIname)
    subplot(numel(metrics),numel(ROIname),r)
    title(ROIlabels{r})
end
subplot(numel(metrics),numel(ROIname),1);
legend(conditions, 'Location','Best')
