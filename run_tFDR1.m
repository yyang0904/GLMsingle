function run_tFDR1(subject)

%% Varibles

alphaList = 10 .^ (-(0:10)); % list of alpha values

%% Directories

% in LAB laptop
% paths.main = '/Users/yyang3077/Desktop/AMPB';
% paths.software = '/Users/yyang3077/Documents/Toolbox';
% addpath(genpath(fullfile('/Applications/freesurfer/8.0.0-beta/matlab')));
% addpath(genpath(fullfile(paths.software,'spm12-main')));
% paths.data = fullfile(paths.main, 'glm_results'); 
% paths.glm  = fullfile(paths.data, subject);

% in pc
paths.main = '/Volumes/T7/AMPB';
paths.software = '/Users/yang/Downloads';
addpath(genpath(fullfile(paths.software, 'freesurfer', 'matlab')));
addpath(genpath(fullfile(paths.software,'spm12-main')));
    rmpath(genpath('/Users/yang/Downloads/spm12-main/external/fieldtrip'));
paths.data = fullfile(paths.main, 'glm_results'); 
paths.glm  = fullfile(paths.data, subject);

% addpath(genpath(fullfile(paths.main, 'scripts', 'helpers'))); 

%% Calculate FDR Thresholded T-Statistics Maps

% tFiles = dir(fullfile(paths.glm, '*_tstat.mgz')); 
% tFiles = fullfile(paths.glm, {tFiles.name}); 

tFiles = dir(fullfile(paths.glm, '*_tstat.mgz')); 
tFiles = tFiles(~startsWith({tFiles.name}, '._')); % Exclude invalid files
tFiles = fullfile(paths.glm, {tFiles.name});

for i = 1:length(tFiles) % for each t-stat file
    [~,fname,ext] = fileparts(tFiles{i}); 
    fprintf('Processing: %s\n', [fname, ext]); 
    tCurr = squeeze(MRIread(tFiles{i}).vol); % current t-stat file
    pCurr = squeeze(MRIread(regexprep(tFiles{i}, '_tstat', '_pval')).vol); 
    [pSort,indx] = sort(pCurr); % sort p-values
    nv = length(pSort); % number of vertices

    tFDR = NaN(length(tCurr), length(alphaList)); % initialize
    for i2 = 1:length(alphaList) % for each alpha level
        %%% compute FDR threshold for each vertex
        fdr = alphaList(i2) ./ (nv - (0:(nv - 1)));
        h = pSort(:) < fdr(:); % vertices that are significant
        thr = find(h == 0, 1); % find first non-signficant vertex
        sigIndx = indx(1:(thr-1)); % significant vertex indices

        %%% create FDR thresholded t-statistics map
        tFDR(sigIndx,i2) = tCurr(sigIndx);
    end

    clear mri; mri.vol = shiftdim(tFDR, -2); % assign to surface
    mri.nframes = length(alphaList); % number of alpha / frames
    saveName = regexprep(tFiles{i}, '_tstat', '_tfdr');
    MRIwrite(mri, saveName); [~,fname,ext] = fileparts(saveName); 
    fprintf('Saved: %s\n', [fname, ext]); 
end

end