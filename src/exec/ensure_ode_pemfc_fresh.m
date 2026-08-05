function ensure_ode_pemfc_fresh()
% ENSURE_ODE_PEMFC_FRESH  Rebuild the ode_PEMFC MEX file if stale.


    thisFile = mfilename('fullpath');
    root = thisFile;
    while ~isfolder(fullfile(root, 'bin')) && ~strcmp(root, fileparts(root))
        root = fileparts(root);
    end

    if ~isfolder(fullfile(root, 'bin'))
        error('[ensure_ode_pemfc_fresh] Could not locate a "bin" folder walking up from %s. Check project layout.', thisFile);
    end

    % Search recursively so this works regardless of exact MEX output location
    mexFiles = dir(fullfile(root, '**', ['ode_PEMFC.' mexext]));

    if isempty(mexFiles)
        fprintf('[ensure_ode_pemfc_fresh] No compiled MEX found -- building it now.\n');
        rebuild();
        return;
    end

    mexPath = fullfile(mexFiles(1).folder, mexFiles(1).name);
    mexTime = mexFiles(1).datenum;

    % Directories whose .m files can affect the compiled physics:
    % - src/model   (ode_PEMFC.m and everything it calls: channels,
    %                catalyst layers, GDLs, membrane, solid part, etc.)
    % - src/param   (mod_param_PEMFC.m -- default parameter values,
    %                including p.counter_flow)
    % - src/state   (state2struct/struct2state helpers used inside p)
    % - src/util    (vec2struct/struct2vec helpers used inside p)
    watchDirs = {'src/model', 'src/param', 'src/state', 'src/util'};
    newestSrcTime = 0;
    newestSrcFile = '';

    for d = 1:numel(watchDirs)
        fullDir = fullfile(root, watchDirs{d});
        if ~isfolder(fullDir)
            warning('[ensure_ode_pemfc_fresh] Watch directory not found: %s', fullDir);
            continue;
        end
        files = dir(fullfile(fullDir, '*.m'));
        for f = 1:numel(files)
            if files(f).datenum > newestSrcTime
                newestSrcTime = files(f).datenum;
                newestSrcFile = fullfile(files(f).folder, files(f).name);
            end
        end
    end

    if newestSrcTime > mexTime
        fprintf(['[ensure_ode_pemfc_fresh] Stale MEX detected.\n', ...
                 '  Newest source change: %s (%s)\n', ...
                 '  Compiled MEX:         %s (%s)\n', ...
                 '  Rebuilding...\n'], ...
                 newestSrcFile, datetime(newestSrcTime,'ConvertFrom','datenum'), ...
                 mexPath, datetime(mexTime,'ConvertFrom','datenum'));
        rebuild();
    else
        fprintf('[ensure_ode_pemfc_fresh] MEX is up to date (%s). Skipping rebuild.\n', ...
                datetime(mexTime,'ConvertFrom','datenum'));
    end
end

function rebuild()
    ode_PEMFC_CoderScript();
    clear ode_PEMFC;   % drop any cached handle to the old MEX from memory
end