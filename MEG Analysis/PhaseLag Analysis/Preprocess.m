%% Preprocessing.m
% Converts raw MaxFiltered MEG recordings into SPM format, epochs them,
% detects bad channels/trials, and merges the three main-task runs per
% subject. Produces the merged, artefact-flagged files (cesub-*.mat) that
% Analysis.m operates on.
%
% Run top to bottom in one session (the four stages share variables such
% as no_loc, which is computed during Conversion and used in every stage
% after it). Each %% section can also be run/evaluated on its own (e.g.
% via "Run Section" or selecting and evaluating a block of code), since
% nothing here depends on functions defined elsewhere in the file.

%% Setup
% Paths, channel selection, participants, conditions, time windows and
% frequency bands used throughout this script.

% `clear` alone only removes workspace variables, not persistent variables
% inside functions (e.g. spm.m caches which SPM installation is
% "active"), globals from a previous SPM/OSL initialisation, or leftover
% figure windows. Stale state like this can cause SPM/FieldTrip errors
% that have nothing to do with this script's own logic (e.g. spm_version
% failing to read its own revision info, or FieldTrip functions resolving
% inconsistently) if SPM or OSL was already initialised earlier in the
% same MATLAB session. `clear all` and `close all hidden` reset that
% properly; if problems persist, a full MATLAB restart is the most
% reliable fix.
close all hidden
clear all

% Paths
addpath '/imaging/henson/TIME/LatestMEGscripts'
addpath /imaging/local/software/spm_toolbox/osl/osl-core
addpath('/imaging/henson/TIME/MEG_scripts')

% Initialise OSL once, at the start, before any SPM/FieldTrip calls. OSL
% bundles its own complete, internally-consistent SPM12 + FieldTrip.
% Functions used throughout this pipeline (not just osl_detect_artefacts)
% depend on it -- e.g. spm_eeg_load itself calls FieldTrip's
% ft_datatype_sens internally, which is missing from at least one
% separately-installed standalone SPM12 on this cluster. Do not also
% addpath a different SPM12 install alongside OSL's: having two
% different SPM/FieldTrip installations on the path at once causes
% FieldTrip's internal functions to resolve inconsistently between them.
osl_startup('/imaging/local/software/spm_toolbox/osl', 'shared')

% Initialise SPM's M/EEG defaults without opening its GUI. `spm eeg`
% launches the Menu, Interactive and Graphics windows, and publish()
% captures every open figure, so those windows end up embedded in the HTML
% output. spm('defaults', ...) configures exactly the same defaults with no
% windows, and is safe to call repeatedly.
spm('defaults', 'EEG');

% Close any SPM windows that initialisation opened regardless, for the same
% reason. SPM tags its windows, so they can be removed without touching the
% figures this script creates itself.
delete(findall(0, 'Type', 'figure', 'Tag', 'Menu'));
delete(findall(0, 'Type', 'figure', 'Tag', 'Interactive'));
delete(findall(0, 'Type', 'figure', 'Tag', 'Graphics'));

BIDSdir    = '/imaging/henson/TIME/timeBIDS/';
megBIDSdir = '/imaging/henson/TIME/timeBIDS/derivatives/meg-derivatives/';
outdir     = '/imaging/henson/TIME/meg-derivatives/';

% Channels
load('meg_only_chan_names');   % loads "chan_names"
ChanType = 'MEGPLANAR';        % use gradiometers throughout

% Participants and conditions
Subs  = 1:32;
Conds = {'SyncTheta', 'AsyncTheta', 'NoFlicker', 'SyncDelta'};

transdef = '';   % set to 'td' to use the head-motion-transformed files instead

% Time windows (all in seconds unless stated)
TWin = [-1 4];        % epoch window
EWin = [0.75 2.75];   % entrainment window used for condition comparisons and the phase-lag analysis
BWin = [-0.5 -0.3];   % baseline-correction window
SWin = [0.75 3];      % window used for the trial-level theta summary measures

% Frequency bands
frequencies    = [1.2 1.7 2:.5:8];   % 15 frequencies, ~evenly spaced, including 4 Hz and 1.7 Hz exactly
freqres        = .25;                % Hilbert filter half-bandwidth (Hz)
ThetaHz        = 4;
DeltaHz        = 1.7;
ThetaBand      = [3 8];               % broad theta range
RelThetaBound  = [3.5 4.5];           % flanking frequencies for the "relative 4 Hz" trial-level power normalisation (distinct from the weight-estimation flanking band used in Analysis.m)

Ind4Hz    = find(frequencies == ThetaHz);
Ind3Hz    = find(frequencies == ThetaBand(1));
Ind8Hz    = find(frequencies == ThetaBand(2));
Ind3_5Hz  = find(frequencies == RelThetaBound(1));
Ind4_5Hz  = find(frequencies == RelThetaBound(2));

%% Convert raw .fif files to SPM format
% Both the three main-task runs and the auditory localiser run are
% converted the same way (build an spm_eeg_convert input struct, run it).

no_loc = [];   % subjects missing a localiser recording

parfor iSub = 1:length(Subs)

    sub_nam = sprintf('sub-%02d', Subs(iSub));
    sub_out = fullfile(outdir, sub_nam, 'meg');
    if ~exist(sub_out, 'dir')
        mkdir(sub_out);
    end
    cd(sub_out)

    sub_in = fullfile(megBIDSdir, sub_nam, 'meg');

    for iRun = 1:3
        run_nam = sprintf('%s_task-TIME_run-%d_meg_proc-sss%s', sub_nam, iRun, transdef);
        S = [];
        S.dataset  = fullfile(sub_in, [run_nam '.fif']);
        S.mode     = 'continuous';
        S.outfile  = fullfile(sub_out, run_nam);
        S.channels = chan_names;
        spm_eeg_convert(S);
    end

    run_nam    = sprintf('%s_task-loc_meg_proc-sss%s', sub_nam, transdef);
    loc_infile = fullfile(sub_in, [run_nam '.fif']);
    if exist(loc_infile, 'file')
        S = [];
        S.dataset  = loc_infile;
        S.mode     = 'continuous';
        S.outfile  = fullfile(sub_out, run_nam);
        S.channels = chan_names;
        spm_eeg_convert(S);
    else
        warning('Subject %s has no localiser FIF file: %s', sub_nam, loc_infile);
        no_loc = [no_loc iSub];
    end
end

%% Epoching
% Trials are cut from -1000 to 4000 ms relative to stimulus onset, with a
% 17 ms correction for display latency (one refresh at 60 Hz).

epoch_times = [-1000 4000];   % ms

% Stimulus delivery delays, applied to the recorded trigger times.
% The visual delay is the projector's display latency, about one refresh at
% 60 Hz. It applies to the main task, where the trigger marks the frame
% request rather than the frame appearing.
%
% It does NOT apply to the auditory localiser: sound does not pass through
% the projector, so shifting those epochs by a display latency introduces an
% error rather than correcting one. The original script applied the visual
% delay to both. The effect on the auditory weights is small, since they come
% from a power spectrum over a 2.5 s window and a 17 ms shift barely changes
% spectral magnitude, but it is wrong in principle and is separated here.
% If the auditory delivery latency through the tube phones has been measured,
% put it in aud_delay; otherwise 0 is the honest default.
vis_delay   = 17;              % ms, ~1 refresh at 60 Hz
aud_delay   = 0;               % ms, auditory delivery latency (0 unless measured)

ntrls = cell(1, length(Subs));

parfor iSub = 1:length(Subs)

    sub_nam = sprintf('sub-%02d', Subs(iSub));
    sub_out = fullfile(outdir, sub_nam, 'meg');
    cd(sub_out)
    ntrls{iSub} = nan(1, 4);   % columns: run 1, run 2, run 3, localiser

    for iRun = 1:3
        infile = sprintf('sub-%02d_task-TIME_run-%d_meg_proc-sss%s.mat', Subs(iSub), iRun, transdef);
        D = spm_eeg_load(infile);

        trl = spm_load(fullfile(BIDSdir, sub_nam, 'meg', ...
            sprintf('sub-%02d_task-TIME_run-%d_events.tsv', Subs(iSub), iRun)));
        trl.sample = round(trl.sample / (1000 / D.fsample));   % rescale if downsampled

        ntrls{iSub}(iRun) = size(trl.sample, 1);
        swin = round(epoch_times * D.fsample / 1000);

        S = [];
        S.D  = D;
        S.bc = 1;
        onsets  = trl.sample + round(vis_delay * D.fsample / 1000);
        S.trl   = [onsets + swin(1), onsets + swin(2), repmat(swin(1), ntrls{iSub}(iRun), 1)];  % FieldTrip trl format
        S.conditionlabels = trl.trial_type;

        D = spm_eeg_epochs(S);
        ntrls{iSub}(iRun) = D.ntrials;
    end

    if ~ismember(iSub, no_loc)
        infile = sprintf('sub-%02d_task-loc_meg_proc-sss%s.mat', Subs(iSub), transdef);
        D = spm_eeg_load(infile);

        trl = spm_load(fullfile(BIDSdir, sub_nam, 'meg', ...
            sprintf('sub-%02d_task-loc_events.tsv', Subs(iSub))));
        trl.sample = round(trl.sample / (1000 / D.fsample));

        ntrls{iSub}(4) = size(trl.sample, 1);
        swin = round(epoch_times * D.fsample / 1000);

        S = [];
        S.D  = D;
        S.bc = 1;
        onsets = trl.sample + round(aud_delay * D.fsample / 1000);   % auditory, not visual
        S.trl  = [onsets + swin(1), onsets + swin(2), repmat(swin(1), ntrls{iSub}(4), 1)];
        S.conditionlabels = repmat({'AudLocSyncTheta'}, size(S.trl, 1), 1);

        D = spm_eeg_epochs(S);
        ntrls{iSub}(4) = D.ntrials;
    end
end

ntrls = cat(1, ntrls{:});

%% Bad-channel and bad-trial detection
% How osl_detect_artefacts decides, from the OSL source
% (github.com/OHBA-analysis/osl-core/blob/master/osl_detect_artefacts.m):
%
%   metric      standard deviation, taken per channel across all samples of
%               the good trials for bad channels, and per trial across all
%               samples of the good channels for bad trials
%   test        generalised extreme studentised deviate (GESD), alpha 0.05
%   iteration   alternates channel and trial passes, up to 5 rounds, stopping
%               when a round finds nothing
%   caps        at most 10 new bad channels per modality in total, and at most
%               20% of remaining trials per trial pass
%   scope       for epoched data, trials are marked bad across all modalities
%               together, whereas channels are handled within modality
%   carry-over  anything already marked bad stays bad
%
% Nothing is deleted; trials and channels are flagged, and every later step
% selects on 'GOOD'.
%
% Two consequences worth recording. The criterion is variance alone, so a
% trial is rejected for being unusually large or small overall, not for
% containing any particular artefact. And the 20% cap means rejection is
% relative rather than absolute: with five iterations the ceiling is high, but
% a participant whose whole recording is noisy loses no more trials than one
% whose recording is clean, because the test is against that participant's own
% distribution.
% Uses OSL's automatic artefact detection on each run separately (before
% merging), so that a channel or trial marked bad in one run does not
% discard good data from the others. OSL was already initialised once in
% Setup above.
%
% *No component-based artefact removal is performed.* osl_detect_artefacts
% flags whole channels and whole trials whose amplitude or variance is
% atypical; it does not decompose the data and subtract ocular, cardiac or
% muscle components. Wang et al. (2026), whose analysis this pipeline is
% partly replicating, ran ICA and removed components identified as eye
% movement, heartbeat and muscle before epoching.
%
% This matters most for anterior sensors, where ocular activity projects
% strongly, and therefore for any effect whose topography is frontally
% weighted: such an effect cannot be distinguished from ocular contamination
% with the present preprocessing. Note also that chan_names retains MEG
% channels only, so EOG and ECG are discarded at conversion and are not
% available either as ICA classifiers or as a covariate for checking
% contamination after the fact.
%
% Adding ICA would require re-converting with the EOG and ECG channels kept,
% then running OSL's AFRICA (osl_africa) or an equivalent between epoching
% and artefact detection. Until then, a frontally weighted result should be
% reported with that limitation stated explicitly.

bad_chans  = cell(1, length(Subs));
bad_trials = cell(1, length(Subs));
nTrials    = cell(1, length(Subs));
RunPresent = false(length(Subs), 3);
LocPresent = false(1, length(Subs));
DataRank   = nan(1, length(Subs));
nGoodChan  = nan(1, length(Subs));

parfor iSub = 1:length(Subs)

    sub_nam = sprintf('sub-%02d', Subs(iSub));
    sub_in  = fullfile(outdir, sub_nam, 'meg');
    cd(sub_in)

    bad_chans{iSub}  = cell(4, 1);
    bad_trials{iSub} = cell(4, 1);
    nTrials{iSub}    = nan(1, 4);

    % Accumulated locally and written to RunPresent once below. parfor
    % requires every access to a sliced variable to use the same subscripts,
    % so it cannot be indexed as (iSub, iRun) here and (iSub, :) later.
    runOK = false(1, 3);

    for iRun = 1:3
        runFile = sprintf('e%s_task-TIME_run-%d_meg_proc-sss%s.mat', sub_nam, iRun, transdef);
        if exist(fullfile(sub_in, runFile), 'file') ~= 2
            % Absent runs are logged rather than skipped silently, so that a
            % missing file cannot be mistaken for a run with no bad trials.
            continue
        end
        runOK(iRun) = true;
        D = spm_eeg_load(runFile);
        D = osl_detect_artefacts(D);
        bad_chans{iSub}{iRun}  = D.badchannels;
        bad_trials{iSub}{iRun} = D.badtrials;
        nTrials{iSub}(iRun)    = D.ntrials;
        D.save;
    end

    locFile = sprintf('e%s_task-loc_meg_proc-sss%s.mat', sub_nam, transdef);
    LocPresent(iSub) = exist(fullfile(sub_in, locFile), 'file') == 2;
    if LocPresent(iSub)
        D = spm_eeg_load(fullfile(sub_in, locFile));
        D = osl_detect_artefacts(D);
        bad_chans{iSub}{4}  = D.badchannels;
        bad_trials{iSub}{4} = D.badtrials;
        nTrials{iSub}(4)    = D.ntrials;
        D.save;
    end

    % Rank after tSSS, from the good gradiometers of the first available run.
    % SSS projects onto a basis of far fewer components than there are
    % channels, so this is well below the channel count and is the number that
    % matters for how much independent information the data carry.
    RunPresent(iSub, :) = runOK;
    firstRun = find(runOK, 1);
    if ~isempty(firstRun)
        Dr = spm_eeg_load(sprintf('e%s_task-TIME_run-%d_meg_proc-sss%s.mat', ...
            sub_nam, firstRun, transdef));
        gi = indchantype(Dr, 'MEGPLANAR', 'GOOD');
        gt = indtrial(Dr, Conds, 'GOOD');
        if ~isempty(gi) && ~isempty(gt)
            ev = sort(eig(cov(reshape(Dr(gi, :, gt), numel(gi), [])')), 'descend');
            DataRank(iSub)  = sum(ev > ev(1) * 1e-6);
            nGoodChan(iSub) = numel(gi);
        end
    end
end

save(fullfile(outdir, 'bad_things'), 'bad_chans', 'bad_trials', 'nTrials')

% Summary for the methods section.
save(fullfile(outdir, 'preproc_report'), 'bad_chans', 'bad_trials', 'nTrials', ...
     'RunPresent', 'LocPresent', 'DataRank', 'nGoodChan')

fprintf('\n%-9s %8s %8s %8s %10s %9s %9s %7s\n', 'Subject', ...
    'run-1', 'run-2', 'run-3', 'localiser', 'good trls', 'bad chans', 'rank');
for iSub = 1:length(Subs)
    pres = cell(1, 3);
    for iRun = 1:3
        if ~RunPresent(iSub, iRun)
            pres{iRun} = 'MISSING';
        else
            pres{iRun} = sprintf('%d/%d', nTrials{iSub}(iRun) - numel(bad_trials{iSub}{iRun}), ...
                nTrials{iSub}(iRun));
        end
    end
    if LocPresent(iSub)
        locStr = sprintf('%d/%d', nTrials{iSub}(4) - numel(bad_trials{iSub}{4}), nTrials{iSub}(4));
    else
        locStr = 'MISSING';
    end
    goodMain = 0; presMain = 0;
    for iRun = 1:3
        if RunPresent(iSub, iRun)
            goodMain = goodMain + nTrials{iSub}(iRun) - numel(bad_trials{iSub}{iRun});
            presMain = presMain + nTrials{iSub}(iRun);
        end
    end
    medChan = median(cellfun(@numel, bad_chans{iSub}(RunPresent(iSub, :))));
    fprintf('%-9s %8s %8s %8s %10s %9d %9.0f %7d\n', sprintf('sub-%02d', Subs(iSub)), ...
        pres{1}, pres{2}, pres{3}, locStr, goodMain, medChan, DataRank(iSub));
end

nMissRun = sum(~RunPresent(:));
fprintf(['\nMain-task runs absent: %d, affecting %d participants (%s).\n'], ...
    nMissRun, sum(any(~RunPresent, 2)), ...
    strjoin(arrayfun(@(k) sprintf('sub-%02d', Subs(k)), find(any(~RunPresent, 2)), ...
    'UniformOutput', false), ', '));
fprintf('Localiser absent for %d participants (%s).\n', sum(~LocPresent), ...
    strjoin(arrayfun(@(k) sprintf('sub-%02d', Subs(k)), find(~LocPresent), ...
    'UniformOutput', false), ', '));

goodAll = nan(1, length(Subs)); presAll = nan(1, length(Subs));
for iSub = 1:length(Subs)
    g = 0; p = 0;
    for iRun = 1:3
        if RunPresent(iSub, iRun)
            g = g + nTrials{iSub}(iRun) - numel(bad_trials{iSub}{iRun});
            p = p + nTrials{iSub}(iRun);
        end
    end
    goodAll(iSub) = g; presAll(iSub) = p;
end
fprintf(['Good main-task trials: median %d, range %d to %d. ' ...
         'Rejected %.1f%% on average (range %.1f to %.1f).\n'], ...
    median(goodAll), min(goodAll), max(goodAll), ...
    100 * mean(1 - goodAll ./ presAll), 100 * min(1 - goodAll ./ presAll), ...
    100 * max(1 - goodAll ./ presAll));

allBadChan = [];
for iSub = 1:length(Subs)
    allBadChan = [allBadChan, cellfun(@numel, bad_chans{iSub}(RunPresent(iSub, :)))']; %#ok<AGROW>
end
fprintf('Bad channels per run: median %d, range %d to %d, of %d gradiometers.\n', ...
    median(allBadChan), min(allBadChan), max(allBadChan), median(nGoodChan, 'omitnan'));
fprintf('Rank after tSSS: median %d, range %d to %d.\n', ...
    median(DataRank, 'omitnan'), min(DataRank), max(DataRank));


%% Merging runs
% Note on sampling rate. These data are at 500 Hz, having been downsampled
% by two from the acquired 1 kHz upstream of this pipeline, which is why the
% epoching section rescales event samples relative to 1000 Hz. That happens
% to match the rate Wang et al. resample to, so no further downsampling is
% needed. It also conditions the narrow phase-extraction filters well: a
% 1 Hz band at 4 Hz is stable at the requested order 4 here, where at 1 kHz
% it would have to be reduced.

% The three main-task runs are concatenated into a single file per
% subject. The per-run files are deleted afterward, in a separate step
% below, once every subject has been merged.

ntrls = cell(1, length(Subs));

for iSub = 1:length(Subs)

    sub_nam = sprintf('sub-%02d', Subs(iSub));
    sub_in  = fullfile(outdir, sub_nam, 'meg');
    cd(sub_in)

    fn = cell(1, 3);
    for iRun = 1:3
        fn{iRun} = fullfile(sub_in, sprintf('e%s_task-TIME_run-%d_meg_proc-sss%s.mat', sub_nam, iRun, transdef));
    end

    % Check every expected file exists before calling spm_eeg_merge, so a
    % missing/incomplete file fails with a clear, specific message rather
    % than SPM's generic "Trouble reading files".
    missingFiles = fn(cellfun(@(f) exist(f, 'file') ~= 2, fn));
    if ~isempty(missingFiles)
        error('Subject %s: missing expected epoched run file(s), cannot merge:\n%s', ...
            sub_nam, strjoin(missingFiles, newline));
    end

    S = [];
    S.D = char(fn);   % full paths, rather than bare filenames relying on the current directory
    D = spm_eeg_merge(S);

    for c = 1:length(D.condlist)
        ntrls{iSub}(c) = length(indtrial(D, D.condlist{c}, 'GOOD'));
    end

    if ~ismember(iSub, no_loc)
        D_loc = spm_eeg_load(fullfile(sub_in, sprintf('e%s_task-loc_meg_proc-sss%s.mat', sub_nam, transdef)));
        ntrls{iSub}(end + 1) = length(indtrial(D_loc, D_loc.condlist, 'GOOD'));
    else
        ntrls{iSub}(end + 1) = 0;
    end
end

ntrls = cat(1, ntrls{:});
save(fullfile(outdir, 'ntrials_minus_bad'), 'ntrls')

% Delete the per-run epoched files now that every subject has been merged
% successfully, since only the merged (cesub-*) file is used from this
% point on. Kept separate from the merge loop above so a failed/missing
% merge for one subject never risks deleting another subject's per-run
% files, and so nothing is deleted until the whole merge step has
% completed without error. Additionally, the merged file's existence is
% checked per subject before deleting anything for that subject, so a
% merge that silently failed to save doesn't leave you with neither the
% per-run files nor a usable merged file.
parfor iSub = 1:length(Subs)
    sub_nam = sprintf('sub-%02d', Subs(iSub));
    sub_in  = fullfile(outdir, sub_nam, 'meg');

    mergedFile = fullfile(sub_in, sprintf('cesub-%02d_task-TIME_run-1_meg_proc-sss.mat', Subs(iSub)));
    if exist(mergedFile, 'file') ~= 2
        error('Subject %s: merged file not found, refusing to delete per-run files:\n%s', ...
            sub_nam, mergedFile);
    end

    for iRun = 1:3
        fn = fullfile(sub_in, sprintf('e%s_task-TIME_run-%d_meg_proc-sss%s.mat', sub_nam, iRun, transdef));
        delete(spm_eeg_load(fn));
    end
end