function ensure_ode_pemfc_fresh()
% ENSURE_ODE_PEMFC_FRESH  Rebuild the ode_PEMFC MEX file only if stale.
%
% Call this once at the top of any exec script that relies on ode_PEMFC
% being executed via its compiled MEX version (i.e. any script calling
% ode_PEMFC(t,x,u) with 3 args, which MATLAB will route to the MEX file
% instead of the .m source if one exists on the path).
%
% Compares the MEX file's modification time against every .m file under
% src/ that could feed into the compiled call tree (model equations +
% parameters). If any source file is newer than the MEX file -- or the
% MEX file doesn't exist yet -- it reruns ode_PEMFC_CoderScript() to
% rebuild it. Otherwise it does nothing (fast path, no recompile).
%
% This gives you the same safety as recompiling on every run, without
% paying the codegen cost when nothing relevant has changed.
    thisFile = mfilename('fullpath');
    root = thisFile;
    while ~isfolder(fullfile(root, 'bin')) && ~strcmp(root, fileparts(root))
    root = fileparts(root);
    end
    mexFiles = dir(fullfile(root, 'bin', 'ode_PEMFC.mex*'));

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
        files = dir(fullfile(root, watchDirs{d}, '*.m'));
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
                 newestSrcFile, datetime(newestSrcTime), ...
                 mexPath, datetime(mexTime));
        rebuild();
    else
        fprintf('[ensure_ode_pemfc_fresh] MEX is up to date (%s). Skipping rebuild.\n', ...
                datetime(mexTime));
    end
end

function rebuild()
    ode_PEMFC_CoderScript();
    clear ode_PEMFC;   % drop any cached handle to the old MEX from memory
end
