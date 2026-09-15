% TIME_TF.m
%
% Setup and time-frequency decomposition. Run this once, then Analysis.m.
%
% The two are separate because this script needs a parallel pool and reads
% every subject's raw data, while Analysis.m works entirely from what this
% one writes. Re-run this only when something upstream changes: the
% frequencies vector, the epoch or baseline windows, ChanType, or the
% preprocessing itself.
%
% Writes to resultdir: TF_output.mat, holding ThetaBand_pSub, bThetaBand_pSub,
% DeltaBand_pSub, TF_pSub and bTF_pSub.
%
% Analysis.m repeats the setup below, so the two must be kept in step. The
% parameters are identical in both by design; if one changes, change both.

%% Setup
% Paths, channel selection, participants, conditions, time windows and
% frequency bands used throughout this script. Kept in sync with the
% Setup section at the top of Preprocessing.m -- if you change one,
% change the other to match.

close all hidden
clear all

% Small helper used when reporting binary outcomes in text
ternary_str = @(cond, a, b) subsref({b, a}, struct('type', '{}', 'subs', {{cond + 1}}));

% Paths
addpath '/imaging/henson/TIME/LatestMEGscripts'
addpath /imaging/local/software/spm_toolbox/osl/osl-core
addpath('/imaging/henson/TIME/MEG_scripts')

% Initialise OSL once, at the start, before any SPM/FieldTrip calls. OSL
% bundles its own complete, internally-consistent SPM12 + FieldTrip. Do
% not also addpath a different SPM12 install alongside OSL's: having two
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
resultdir = '/imaging/henson/TIME/LatestMEGscripts/MEGxBehavioural';   % where all outputs are written

% Channels
load('meg_only_chan_names');   % loads "chan_names"
ChanType = 'MEGPLANAR';        % use gradiometers throughout

% Participants and conditions
Subs  = 1:32;
Conds = {'SyncTheta', 'AsyncTheta', 'NoFlicker', 'SyncDelta'};

transdef = '';   % set to 'td' to use the head-motion-transformed files instead

% Time windows (all in seconds unless stated)
TWin = [-1 4];        % epoch window
EWin = [0.75 2.75];   % entrainment window: used for the trial-level power measures and
                      % for estimating the sensor weights. The phase analysis uses its
                      % own window, PhaseWin, matching Wang et al.
BWin = [-0.5 -0.3];   % baseline-correction window

% Frequency bands
frequencies    = [1.2 1.7 2:.5:8];   % 15 frequencies, ~evenly spaced, including 4 Hz and 1.7 Hz exactly
freqres        = .25;                % Hilbert filter half-bandwidth (Hz)
ThetaHz        = 4;
DeltaHz        = 1.7;
ThetaBand      = [3 8];               % frequency span of ThetaBand_pSub, not a measure

Ind4Hz    = find(frequencies == ThetaHz);
Ind3Hz    = find(frequencies == ThetaBand(1));
Ind8Hz    = find(frequencies == ThetaBand(2));
if any(cellfun(@isempty, {Ind4Hz, Ind3Hz, Ind8Hz}))
    error(['Theta frequencies %g, %g, %g are not all in the frequencies ' ...
           'vector; every trial-level power measure indexes these bins directly.'], ...
        ThetaHz, ThetaBand(1), ThetaBand(2));
end

% Position of 4 Hz within ThetaBand_pSub's frequency dimension, which spans ThetaBand rather than all of frequencies.
Ind4HzInBand = Ind4Hz - Ind3Hz + 1;

DeltaBand = [1.2 2];   % frequency span of DeltaBand_pSub, not a measure
Ind1_7Hz  = find(frequencies == DeltaHz);
Ind1_2Hz  = find(frequencies == DeltaBand(1));
Ind2Hz    = find(frequencies == DeltaBand(2));
if isempty(Ind1_7Hz) || isempty(Ind1_2Hz) || isempty(Ind2Hz)
    error('Delta frequencies %g, %g, %g are not all in the frequencies vector.', ...
        DeltaHz, DeltaBand(1), DeltaBand(2));
end

Ind1_7HzInBand = Ind1_7Hz - Ind1_2Hz + 1;

% Sensor layout (2D coordinates for topographies), used later in this
% script's Diagnostics section. Loaded from a non-transformed subject
% (sub-30) because the transformed ("td") .fif files have updated xy
% coordinates that produce distorted topographies.
D  = spm_eeg_load(fullfile(outdir, 'sub-30/meg/esub-30_task-loc_meg_proc-sss.mat'));
XY = D.coor2D;
[~, XYinds] = intersect(indchantype(D, {'MEG', 'MEGPLANAR'}), indchantype(D, ChanType));
XY = XY(:, XYinds);


%% Time-frequency analysis
% Power is estimated per trial with a Hilbert transform at the 15
% frequencies defined in the Setup section, then baseline-corrected two
% ways: subtractive (|Pow|) and percentage-change (|bPow|). Two summaries
% are kept per trial: the channel-averaged time-frequency map
% (|TF|/|bTF|) and the band-limited, source-window-averaged power for theta
% (|ThetaBandPow|/|bThetaBandPow|) and delta (|DeltaBandPow|). The band
% arrays keep every channel, so Analysis.m can weight channels itself; TF
% is already averaged across channels.
%
% *Fix relative to the original script:* |TF_pSub|/|bTF_pSub|
% (channel-averaged, frequency x time x trial) are saved here and used
% directly in the next section, instead of a separate |Pow_pSub|
% variable that was referenced later in the original script but never
% actually computed or saved.

ThetaBand_pSub  = cell(1, length(Subs));
bThetaBand_pSub = cell(1, length(Subs));
DeltaBand_pSub  = cell(1, length(Subs));
TF_pSub         = cell(1, length(Subs));
bTF_pSub        = cell(1, length(Subs));

nThetaFreqs = length(Ind3Hz:Ind8Hz);
nDeltaFreqs = length(Ind1_2Hz:Ind2Hz);

parfor iSub = 1:length(Subs)

    sub_nam = sprintf('sub-%02d', Subs(iSub));
    sub_in  = fullfile(outdir, sub_nam, 'meg');
    cd(sub_in)

    D = spm_eeg_load(fullfile(sub_in, sprintf('cesub-%02d_task-TIME_run-1_meg_proc-sss.mat', Subs(iSub))));

    % Restricted to the analysed sensor type. The original script used
    % indchantype(D, 'ALL', 'GOOD'), which averages magnetometers and
    % gradiometers together. Those carry different units (T against T/m) and
    % differ in numerical magnitude by roughly two orders, so the average was
    % dominated by one type in a way that has nothing to do with signal, and
    % mixed two sensor geometries in a single number.
    ChanInds = indchantype(D, ChanType, 'GOOD');

    nTimeSamples  = length(D.indsample(TWin(1)):D.indsample(TWin(2)));
    ThetaBandPow  = nan(D.size(1), nThetaFreqs, D.size(3));
    bThetaBandPow = nan(D.size(1), nThetaFreqs, D.size(3));
    DeltaBandPow  = nan(D.size(1), nDeltaFreqs, D.size(3));
    TF            = nan(length(frequencies), nTimeSamples, D.size(3));
    bTF           = nan(length(frequencies), nTimeSamples, D.size(3));

    for iTrial = 1:size(D, 3)

        Y = ft_specest_hilbert(D(:, :, iTrial), D.time, ...
            'freqoi', frequencies, 'width', freqres, 'filttype', 'but', ...
            'filtorder', 2, 'filtdir', 'twopass', 'polyremoval', 1, 'verbose', 0);

        Pow = abs(Y);

        BLine = mean(Pow(:, :, D.indsample(BWin(1)):D.indsample(BWin(2))), 3);
        bPow  = 100 * (Pow ./ repmat(BLine, [1 1 D.nsamples]) - 1);
        Pow   = Pow - repmat(BLine, [1 1 D.nsamples]);

        ThetaBandPow(:, :, iTrial)  = squeeze(mean(Pow(:, Ind3Hz:Ind8Hz, D.indsample(EWin(1)):D.indsample(EWin(2))), 3));
        bThetaBandPow(:, :, iTrial) = squeeze(mean(bPow(:, Ind3Hz:Ind8Hz, D.indsample(EWin(1)):D.indsample(EWin(2))), 3));

        DeltaBandPow(:, :, iTrial)  = squeeze(mean(Pow(:, Ind1_2Hz:Ind2Hz, D.indsample(EWin(1)):D.indsample(EWin(2))), 3));

        TF(:, :, iTrial)  = squeeze(mean(Pow(ChanInds, :, D.indsample(TWin(1)):D.indsample(TWin(2))), 1));
        bTF(:, :, iTrial) = squeeze(mean(bPow(ChanInds, :, D.indsample(TWin(1)):D.indsample(TWin(2))), 1));
    end

    ThetaBand_pSub{iSub}  = ThetaBandPow;
    bThetaBand_pSub{iSub} = bThetaBandPow;
    DeltaBand_pSub{iSub}  = DeltaBandPow;
    TF_pSub{iSub}         = TF;
    bTF_pSub{iSub}        = bTF;
end

cd(resultdir)
% -v7.3 because the combined arrays exceed the 2 GB limit of the default format.
save('TF_output.mat', 'ThetaBand_pSub', 'bThetaBand_pSub', 'DeltaBand_pSub', ...
     'TF_pSub', 'bTF_pSub', '-v7.3')