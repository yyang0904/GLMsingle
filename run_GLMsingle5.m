%% input run_GLMsingle5('sub-NSxLxIUx1994'); to run the script, only for Localizer

function run_GLMsingle5(subject)

%% Variables
TR = 2; % repetition time in seconds

%%% define contrasts by task
contrasts.ptlocal = 'motion = 1 : stationary = -1';
contrasts.mtlocal = 'motion = 1 : stationary = -1';

%% Directories

    paths.base = '/Volumes/T7/AMPB';

    if contains(subject, 'EB')
        session = 'ses-01b';
    elseif contains(subject,'NS')
        session = 'ses-01';
    else
        error('subject vision unclear')
    end

    paths.main = fullfile(paths.base, 'derivatives', 'fmriprep', subject, session, 'func');
    paths.software = '/Users/yang/Downloads';
    addpath(genpath(fullfile(paths.software, 'freesurfer', 'matlab'))); % FreeSurfer MATLAB functions
    addpath(genpath(fullfile(paths.software, 'GLMsingle-main')));       % GLMsingle functions
    addpath(genpath(fullfile(paths.software, 'gifti-main')));
    addpath(genpath(fullfile(paths.software, 'matStats-master')));
    addpath(genpath(fullfile(paths.software, 'SML-main')));

    % Data directories
    paths.data = paths.main; 
    paths.func = paths.main;
    paths.save = fullfile(paths.base, 'glm_results', subject); % Folder to save GLM output
    

    % % Find BOLD files (GIFTI format)
    % boldFiles  = dir(fullfile(paths.func, '*_bold.func.gii')); 
    % boldFiles  = fullfile(paths.func, {boldFiles.name}); 
    
    % List all GIFTI files ending with '_bold.func.gii'
    files = dir(fullfile(paths.func, '*_bold.func.gii'));
    validFiles = files(~startsWith({files.name}, '._'));     % Filter out any files starting with '._'
    boldFiles = fullfile(paths.func, {validFiles.name});     % Construct full paths to the valid files

    
    % Identify task types from the BOLD file names
    taskList = regexprep(boldFiles, '.+_task-(\w+)_.+', '$1'); 
    
    clear list; % clear previous task list
    list.task = unique(taskList); % get unique task list
    list.hemi = {'hemi-L', 'hemi-R'}; 
    
    % Generate all combinations of tasks and hemispheres
    tasks = list.task;
    hemis = list.hemi;
    
    % Initialize an array of conditions
    conds = struct('task', {}, 'hemi', {});

    % Create combinations of tasks and hemispheres
    for t = 1:length(tasks)
        for h = 1:length(hemis)
            conds(end+1).task = tasks{t};  % Task
            conds(end).hemi = hemis{h};    % Hemisphere
        end
    end
    
    % Define base path for event files
    % basePath = '/Volumes/cos-lab-wpark78/AMPB'; 
    basePath = '/Volumes/Extreme SSD/AMPB/data';
 
    %% for each task and hemisphere
    for i = 1:length(conds) 
        
        %%% Subset BOLD files by task types
        indx = contains(boldFiles, conds(i).task) & ...
            contains(boldFiles, conds(i).hemi); 
        currBold = boldFiles(indx); 
    
        % Generate corresponding event files dynamically based on the number of BOLD runs
        eventFiles = cell(1, length(currBold));  % initialize as a cell array of the same length as BOLD files
        for j = 1:length(currBold)
            eventFiles{j} = fullfile(basePath, subject, session, 'func', ...
                sprintf('%s_%s_task-%s_run-%d_events.tsv', subject, session, conds(i).task, j));
        end
    
        clear opt; % clear previous options
        opt.sessionindicator = str2double(regexprep(currBold, '.+_ses-(\d+)b?_.+', '$1'));
    
        % Debugging print for the number of runs in data and design
        disp(['Number of runs in data: ', num2str(length(currBold))]);
    
        [data, t] = create_data_matrix(currBold, TR); 
        disp(['Number of runs in data: ', num2str(length(data))]);
        [design, stimdur, condnames] = create_design_matrix(eventFiles, t); 
    
        % Debugging print for event files and design matrix length
        disp(['Number of event files: ', num2str(length(eventFiles))]);
        disp(['Number of runs in design: ', num2str(length(design))]);
    
        %%% Check if design and data lengths match
        if length(design) ~= length(data)
            error('Mismatch between the number of runs in the design matrix and the data matrix.');
        end
    
        %%% Call GLMsingle and save results in task-specific directories
        paths.glm = fullfile(paths.save, conds(i).task, conds(i).hemi); 
        if ~exist(paths.glm, 'dir')
            mkdir(paths.glm);
        end

        GLMestimatesingletrial(design, data, stimdur, TR, paths.glm, opt); 
        save(fullfile(paths.glm, 'CONDITION_NAMES.mat'), 'condnames');
    
        %%% Load GLMsingle results (beta weights) and single-trial design matrix
        clear modelmd designSINGLE; % clear previous results
        results = load(fullfile(paths.glm, 'TYPED_FITHRF_GLMDENOISE_RR.mat')); 
        load(fullfile(paths.glm, 'DESIGNINFO.mat'), 'designSINGLE'); 
    
        %%% Unscale beta weights that were saved as percent signal change
        %%% -> beta_scaled = <beta> / mean * 100 
        modelmd = squeeze(results.modelmd) ./ 100; % remove singleton dimensions
        modelmd = bsxfun(@times, modelmd, results.meanvol); % apply amplitude scaling
    
        %%% Write surface maps (.mgz), task-dependent
        switch conds(i).task

           % case 'ampb' % Separate by runs and save
           %      for i2 = 1:length(currBold) % for each run
           %          % Extract beta weights for the current run
           %          beta_map = modelmd(:, :, :, i2); % Adjust dimensions if necessary
           % 
           %          % Prepare MRI structure for saving
           %          clear mri; mri.vol = shiftdim(beta_map, -2); % Adjust dimensions as needed
           %          mri.nframes = size(mri.vol, 4); % number of frames
           % 
           %          % Generate the filename and save
           %          [~, fname] = fileparts(currBold{i2}); % current run filename
           %          saveName = regexprep(fname, '_bold.func', '_beta.mgz'); 
           %          MRIwrite(mri, fullfile(paths.save, saveName));
           %          fprintf('Saved beta map for run: %s\n', saveName);
           %      end

            case {'mtlocal', 'ptlocal'}
                % Identify each trial in designSINGLE for indexing the beta weights
                [~, indx] = cellfun(@(x) find(x), designSINGLE, 'UniformOutput', false);

            
          %% Save beta maps for each run
                for i2 = 1:length(currBold) % for each run
                    % Extract beta weights for the current run using the indices
                    beta_map = modelmd(:, indx{i2}); % Adjust dimensions if necessary
            
                    % Prepare MRI structure for saving
                    clear mri;
                    mri.vol = shiftdim(beta_map, -2); % Adjust dimensions as needed
                    mri.nframes = size(mri.vol, 4); % number of frames
            
                    % Generate the filename and save
                    [~, fname] = fileparts(currBold{i2}); % current run filename
                    saveName = regexprep(fname, '_bold.func', '_beta.mgz');
                    MRIwrite(mri, fullfile(paths.save, saveName));
                    fprintf('Saved beta map for run: %s\n', saveName);
                end
            
          %% Collapse across runs and perform contrast

                %%% Assign condition names to single-trial events
                indx = cellfun(@find_condition_index, design, 'UniformOutput', false);
                betaConditions = condnames(cat(1, indx{:}));
                nt = size(data{1}, 2); % number of time points
    
                %%% Create contrast vector
                c = zeros(length(betaConditions), 1); % initialize
                contrastStr = contrasts.(conds(i).task); % task-based contrast
                contrastStr = strsplit(regexprep(contrastStr, '\s', ''), ':');
                contrastStr = regexp(contrastStr, '(?<name>\w+)=(?<value>[\-\w]+)', 'names');
                contrastStr = cat(1, contrastStr{:}); % uncell structures
                for i2 = 1:length(contrastStr) % for each contrast
                    indx = strcmp(betaConditions, contrastStr(i2).name);
                    c(indx) = str2double(contrastStr(i2).value);
                end
    
                %%% Compute design matrix for each hrf in the HRF library
                hrfs = num2cell(getcanonicalhrflibrary(stimdur, TR), 2); % HRF library
                d = sum(cat(3, designSINGLE{:}), 3); % collapse single-trial design matrix
                X = cellfun(@(x) convn(d, x(:)), hrfs, 'UniformOutput', false);
                X = cellfun(@(x) TR .* x(1:nt,:), X, 'UniformOutput', false); 
    
                %%% Compute variance in the residuals for each vertex
                q = 1 - (results.R2 ./ 100); % variance NOT explained, proportion
                e = q .* var(cat(3, data{:}), [], 2:3, 'omitnan'); % variance of residuals    
    
                %%% Compute numerator (contrast-weighted beta weights) and
                %%% denominator (contrast-weighted design matrix) for each hrf
                w = c' * modelmd'; % weighted contrast
                terror = cellfun(@(x) c' * pinv(x' * x) * c, X); 
    
                %%% Calculate t-statistic and p-value for each vertex
                t = w(:) ./ sqrt(e .* terror(results.HRFindex));
                p = 2 .* tcdf(-abs(t), abs(length(c) - nt));
    
                %%% Save contrast t-statistic map
                clear mri; mri.vol = shiftdim(t, -2); % assign to surface
                contrastDesc = strjoin(lower({contrastStr.name}), 'X');
                saveName = sprintf('%s_task-%s_%s_space-fsnative_desc-%s_tstat.mgz', ...
                    subject, conds(i).task, conds(i).hemi, contrastDesc);
                MRIwrite(mri, fullfile(paths.save, saveName)); 
                fprintf('Saved: %s\n', saveName);
    
                %%% Save contrast p-value map
                clear mri; mri.vol = shiftdim(p, -2); % assign to surface 
                saveName = sprintf('%s_task-%s_%s_space-fsnative_desc-%s_pval.mgz', ...
                    subject, conds(i).task, conds(i).hemi, contrastDesc);
                MRIwrite(mri, fullfile(paths.save, saveName)); 
                fprintf('Saved: %s\n', saveName); 



        end

    end
    
end
    

    %% Helper Functions
    % create_data_matrix
    function [data, t] = create_data_matrix(boldFiles, TR)
        data = cell(1, length(boldFiles)); % initialize data cell array
        for i = 1:length(boldFiles) % for each BOLD file
            [~, ~, ext] = fileparts(boldFiles{i});
            if contains(ext, 'gii')
                bold = gifti(boldFiles{i}).cdata; 
            else
                bold = MRIread(boldFiles{i}).vol; % load MRI data
            end
            data{i} = single(bold); % store data in the cell array
            if i == 1
                t = ((1:size(bold, ndims(bold))) - 1) .* TR; % create time vector
            end
        end
    end

    % create_design_matrix
    function [design, stimdur, allConditions] = create_design_matrix(eventFiles, t)
        %%% Determine all conditions across event files
        allConditions = cellfun(@(x) unique(cellstr(tdfread(x).trial_type)), ...
            eventFiles, 'UniformOutput', false);
        allConditions = setdiff(unique(cat(1, allConditions{:})), {'blank'});
        
        design = cell(1, length(eventFiles)); % initialize design matrix cell array
        for i = 1:length(eventFiles) % for each event file
            %%% Read events files and extract trial conditions
            events = tdfread(eventFiles{i}); % read events file
            conditions = cellstr(events.trial_type);
    
            %%% Create design matrix per run
            X = zeros(length(t), length(allConditions)); stimdur = NaN; % initialize
            for i2 = 1:length(conditions) % for each condition
                k = strcmp(allConditions, conditions{i2}); 
                if any(k) % if task condition exists
                    if isnan(stimdur); stimdur = events.duration(i2); end % get stimdur
                    tIndx = t >= events.onset(i2) & ...
                        t < (events.onset(i2) + events.duration(i2));
                    indx = find(tIndx, 1); % only select the onset TR
                    X(indx, k) = 1; % assign condition to design matrix
                end
            end
    
            %%% Assign current run's design matrix into cell
            design{i} = X; 
        end
    end


    function [cindx] = find_condition_index(design)
        [nt, cindx] = find(design); [~, indx] = sort(nt); cindx = cindx(indx); 
    end
