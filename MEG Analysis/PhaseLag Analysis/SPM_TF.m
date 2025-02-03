%% Time-Frequency Analysis with SPM
% This script extracts frequency power and saves them as .mat format
%
% The script is written to be compatible for MRC CBU computing cluster.
% Change paths as necessary to to link data and spm.

clear

% set paths and initialize spm
BIDSdir = '/imaging/henson/TIME/timeBIDS/';
megBIDSdir = '/imaging/henson/TIME/timeBIDS/derivatives/meg-derivatives/';
outdir  = '/imaging/henson/TIME/meg-derivatives/';
addpath('/imaging/henson/TIME/MEG_scripts')

addpath '/imaging/henson/TIME/LatestMEGscripts'
addpath /imaging/local/software/spm_toolbox/osl/osl-core
addpath /imaging/local/software/spm_cbu_svn/releases/spm12_latest/
spm eeg

% set channel type
load('meg_only_chan_names') % contains "chan_names"
%chantype = 'MEGMAG';
chantype = 'MEGPLANAR'; % Use gradiometers
in.chantype = chantype;

% Set time windows: 0 is the onset of stimuli
TWin = [-1 4]; % whole epoch time window in secs
EWin = [0 3]; % entraiment window in secs
BWin = [-500 -300]; % baseline correction time window. Needs to be in ms

% Note that if using transdef, xy coords seem updated, so topos better if
% use xy coords from non-transed data file, so take, eg sub32
D = spm_eeg_load(fullfile(outdir,sprintf('sub-30/meg/esub-30_task-loc_meg_proc-sss.mat')));
XY = D.coor2D;
[~,XYinds] = intersect(indchantype(D,{'MEG','MEGPLANAR'}),indchantype(D,chantype));
XY = XY(:,XYinds);

% Set participants
Subs = 1:32;
% Subs = setdiff(Subs,11); % s11 cannot be maxfiltered unless downsampled
% Subs = setdiff(Subs,21);

% Set frequencies
ThetaHz = 4;
RelThetaBound = [3.5 4.5];
ThetaBroad = [3 8];
DeltaHz = 1.7;

% Set conditions
Conds = {'SyncTheta','AsyncTheta','NoFlicker','SyncDelta'};

%% Convert from raw .fif

transdef = ''
% transdef = 'td'

no_loc = [];
parfor iSub = 1:length(Subs)
    
    sub_nam = sprintf('sub-%02d', Subs(iSub));
    sub_out = fullfile(outdir, sub_nam);
    try mkdir(sub_out); end
    
    sub_out = fullfile(sub_out,'meg');
    try mkdir(sub_out); end
    cd(sub_out)
    
    sub_in = fullfile(megBIDSdir,sub_nam,'meg');
    
    for iRun=1:3
        run_nam = sprintf('%s_task-TIME_run-%d_meg_proc-sss%s',sub_nam,iRun,transdef)
        S = [];
        S.dataset = fullfile(sub_in,[run_nam '.fif']);
        S.mode = 'continuous';
        S.outfile = fullfile(sub_out,run_nam); % don't need "spmeeg_*"
        S.channels = chan_names;
        D = spm_eeg_convert(S);
    end
    
    %% Localiser
    run_nam = sprintf('%s_task-loc_meg_proc-sss%s',sub_nam,transdef)
    
    if exist(fullfile(sub_in,[run_nam '.fif']),'file')
        S = [];
        S.dataset = fullfile(sub_in,[run_nam '.fif']);
        S.mode = 'continuous';
        S.outfile = fullfile(sub_out,run_nam); % don't need "spmeeg_*"
        S.channels = chan_names;
        D = spm_eeg_convert(S);
    else
        warning('Subject %s has no localiser FIF file %s?',sub_nam,fullfile(sub_in,[run_nam '.fif']));
        no_loc = [no_loc iSub];
    end
end


%% Epoch

epoch_times = [-1000 4000];
vis_delay = 17;  % In ms. Assume 1 refresh at 60Hz?

ntrls = {};

parfor iSub = 1:length(Subs)
    
    sub_nam = sprintf('sub-%02d', Subs(iSub));
    sub_out = fullfile(outdir, sub_nam,'meg');
    cd(sub_out)
    ntrls{iSub} = nan(1,4);
    
    for iRun=1:3
        S = []; 
        S.D = fullfile(sprintf('sub-%02d_task-TIME_run-%d_meg_proc-sss%s.mat',Subs(iSub),iRun,transdef));
    
        D = spm_eeg_load(S.D);
        
        trl = spm_load(fullfile(BIDSdir,sprintf('sub-%02d',Subs(iSub)),'meg',sprintf('sub-%02d_task-TIME_run-%d_events.tsv',Subs(iSub),iRun)));
        
        trl.sample = round(trl.sample/(1000/D.fsample)); % If data downsampled relative to original 1000Hz
        
        ntrls{iSub}(iRun) = size(trl.sample,1);
        swin = round(epoch_times * D.fsample/1000);
        
        S = []; S.D = D;
        S.bc = 1;
        S.trl = trl.sample + round(vis_delay * D.fsample/1000);
        S.trl = [S.trl+swin(1) S.trl+swin(2) repmat(swin(1),ntrls{iSub}(iRun),1)]; % FieldTrip's "trl" format = [trlbeg trlend trloff]
        S.conditionlabels = trl.trial_type;
        D = spm_eeg_epochs(S);
        ntrls{iSub}(iRun) = D.ntrials;
        
    end
    
    %% Do Localiser
    if ~ismember(iSub,no_loc)
        S = [];
        S.D = fullfile(sprintf('sub-%02d_task-loc_meg_proc-sss%s.mat',Subs(iSub),transdef));
        
        D = spm_eeg_load(S.D);
        
        trl = spm_load(fullfile(BIDSdir,sprintf('sub-%02d',Subs(iSub)),'meg',sprintf('sub-%02d_task-loc_events.tsv',Subs(iSub))));
        
        trl.sample = round(trl.sample/(1000/D.fsample)); % If data downsampled relative to original 1000Hz
        
        ntrls{iSub}(4) = size(trl.sample,1);
        swin = round(epoch_times * D.fsample/1000);
        
        S = []; S.D = D;
        S.bc = 1;
        S.trl = trl.sample + round(vis_delay * D.fsample/1000);
        S.trl = [S.trl+swin(1) S.trl+swin(2) repmat(swin(1),ntrls{iSub}(4),1)]; % FieldTrip's "trl" format = [trlbeg trlend trloff]
        %         S.conditionlabels = trl.trial_type;
        S.conditionlabels = repmat({'AudLocSyncTheta'},size(S.trl,1),1);
        D = spm_eeg_epochs(S);
        ntrls{iSub}(4) = D.ntrials;
        
    end
    
end

ntrls = cat(1,ntrls{:})


%% Detect bad trials and channels
osl_startup('/imaging/local/software/spm_toolbox/osl','shared')


bad_chans = cell(1,length(Subs));
bad_trials = cell(1,length(Subs));
parfor iSub = 1:length(Subs)
    
    sub_nam = sprintf('sub-%02d', Subs(iSub));
    sub_in = fullfile(outdir,sub_nam,'meg');
    cd(sub_in)
    
    bad_chans{iSub} = cell(4,1);
    bad_trials{iSub} = cell(4,1);
    
    for iRun=1:3
        D = spm_eeg_load(sprintf('e%s_task-TIME_run-%d_meg_proc-sss%s.mat',sub_nam,iRun,transdef));
        D = osl_detect_artefacts(D);
        bad_chans{iSub}{iRun} = D.badchannels;
        bad_trials{iSub}{iRun} = D.badtrials;
        D.save;
    end
    
    if ~ismember(iSub,no_loc)
        D = spm_eeg_load(fullfile(sub_in,sprintf('e%s_task-loc_meg_proc-sss%s.mat',sub_nam,transdef)));
        D = osl_detect_artefacts(D);
        bad_chans{iSub}{4} = D.badchannels;
        bad_trials{iSub}{4} = D.badtrials;
        D.save;
    end
end
bad_chans{:}
bad_trials{:}
save(fullfile(outdir,'bad_things'),'bad_chans','bad_trials')

%% Merge (concatenate) runs

ntrls = {};

parfor iSub = 1:length(Subs)
    sub_nam = sprintf('sub-%02d', Subs(iSub));
    sub_in = fullfile(outdir,sub_nam,'meg');
    cd(sub_in)
    
    S = []; fn = {};
    for iRun=1:3
        fn{iRun} = fullfile(sprintf('e%s_task-TIME_run-%d_meg_proc-sss%s.mat',sub_nam,iRun,transdef));
    end
    S.D = strvcat(fn);
    D = spm_eeg_merge(S);
    
    for c = 1:length(D.condlist)
        ntrls{iSub}(c) = length(indtrial(D, D.condlist{c}, 'GOOD'));
    end
    
    if ~ismember(iSub,no_loc)
        D = spm_eeg_load(fullfile(sub_in,sprintf('e%s_task-loc_meg_proc-sss%s.mat',sub_nam,transdef)));
        ntrls{iSub}(end+1) = length(indtrial(D, D.condlist, 'GOOD'));
    else
        ntrls{iSub}(end+1) = 0;
    end
    
    for iRun=1:3
        delete(spm_eeg_load(fn{iRun}));
    end
end

ntrls = cat(1,ntrls{:})
save(fullfile(outdir,'ntrials_minus_bad'),'ntrls')



%% TF
% TF with hilbert and baseline correction with spm_rescale 
Theta4Hz_All = [];
RelTheta4Hz_All = [];
BroadTheta_All = [];

parfor iSub = 1:length(Subs)
    
    sub_nam = sprintf('sub-%02d', Subs(iSub));
    sub_in = fullfile(outdir,sub_nam,'meg');
    cd(sub_in)
    
    S = [];
    % Currently uses down-sampled, epoched, and runs-combined MEG files
    S.D = spm_eeg_load(fullfile(sub_in,sprintf('cesub-%02d_task-TIME_run-1_meg_proc-sss.mat',Subs(iSub))));
    
    S.frequencies = [1.2 1.7 2:.5:8]; 
    S.freqres = .25;
    
    S.settings.frequencies = S.frequencies;
    S.settings.freqres = S.freqres;
    
    S.method = [];
    S.method = 'hilbert';
    S.phase = 1;
    S.settings.filter.type = 'but';
    S.settings.filter.order = 2;
    S.settings.filter.dir = 'twopass';
    S.settings.polyorder = 1;
    S.channels = 'MEGPLANAR';
    S.settings.taper = 'hanning'
    S.taper = 'hanning'
    
    D = spm_eeg_tf(S);
    
    if ~ismember(iSub,no_loc)
        S = [];
        % Currently uses down-sampled, epoched, and runs-combined MEG files
        S.D = spm_eeg_load(fullfile(sub_in,sprintf('esub-%02d_task-loc_meg_proc-sss.mat',Subs(iSub))));
        
        S.frequencies = [1.2 1.7 2:.5:8];
        S.freqres = .25;
        
        S.settings.frequencies = S.frequencies;
        S.settings.freqres = S.freqres;
        
        S.method = [];
        S.method = 'hilbert';
        S.phase = 1;
        S.settings.filter.type = 'but';
        S.settings.filter.order = 2;
        S.settings.filter.dir = 'twopass';
        S.settings.polyorder = 1;
        S.channels = 'MEGPLANAR';
        
        S.settings.taper = 'hanning'
        S.taper = 'hanning'
        
        D = spm_eeg_tf(S);
    end
end

%% Baseline correction
parfor iSub = 1:length(Subs)
    
    sub_nam = sprintf('sub-%02d', Subs(iSub));
    sub_in = fullfile(outdir,sub_nam,'meg');
    cd(sub_in)

    S = [];
    S.D = spm_eeg_load(fullfile(sub_in,sprintf('tf_cesub-%02d_task-TIME_run-1_meg_proc-sss.mat',Subs(iSub))));
    S.method = 'Rel';
    S.timewin = BWin;
    S.pooledbaseline = 0; % Do it trialwise
    
    D = spm_eeg_tf_rescale(S);
    
end

%% Plot power spectra
AllA = nan(15,2001,4,32);
% each subject
lc = {'r','b','g','k'}; % 
for iSub = 1:length(Subs)
    iSub
    f1 = figure; hold on
    for iCond = 1:length(Conds)
        
        sub_nam = sprintf('sub-%02d', Subs(iSub));
        sub_in = fullfile(outdir,sub_nam,'meg');
        cd(sub_in)
        
        D = spm_eeg_load(fullfile(sub_in,sprintf('rtf_cesub-%02d_task-TIME_run-1_meg_proc-sss.mat',Subs(iSub))));
        
        A = squeeze(mean(D(:,:,1:indsample(D,3),indtrial(D, Conds{iCond},'GOOD')),4));
        A = squeeze(mean(A,1));
        AllA(:,:,iCond,iSub) = A;
        
        % Time-Freq
%         figure
%         contourf(D.time(1:indsample(D,3)),D.frequencies, A)
%         title([Conds{iCond} sub_nam chantype])
%         colorbar
        
        % Freq
        figure(f1)
        plot(D.frequencies, mean(A(:,indsample(D,0):end),2),lc{iCond})
        title([sub_nam chantype])
        
    end
    legend(Conds);
%     pause
end

% average per condition
f1 = figure; hold on
for iCond = 1:length(Conds)
TimeFreqs = squeeze(mean(AllA(:,indsample(D,0):indsample(D,3),iCond,:),2,'omitnan'))';
TimeFreqs = TimeFreqs([1:28, 30],:);
% Time-Freq
figure
contourf(D.time(1:indsample(D,3)),D.frequencies, mean(AllA(:,:,iCond,:),4,'omitnan'))
title([Conds{iCond} ' N=30 ' chantype])
xlabel('Time')
ylabel('Freq')
colorbar

% Freq
figure(f1)
% plot(D.frequencies, mean(AllA(:,indsample(D,0):end,iCond,:),[2 4],'omitnan'),lc{iCond})
shadedErrorBar(D.frequencies, TimeFreqs,{@mean,@std},'lineprops', lc(iCond),'patchSaturation',0.1)
end
legend(Conds);

%% Plot topos

% individual
ChanInds = indchantype(D,chantype);
Topo4Hz = [];
Topo1_7Hz = [];

for iSub = 1:length(Subs)
    iSub
    for iCond = 1:length(Conds)
        
        sub_nam = sprintf('sub-%02d', Subs(iSub));
        sub_in = fullfile(outdir,sub_nam,'meg');
        cd(sub_in)
        D = spm_eeg_load(fullfile(sub_in,sprintf('tf_cesub-%02d_task-TIME_run-1_meg_proc-sss.mat',Subs(iSub))));
        
        Topo4Hz(:,iCond,iSub) = mean(D(:,find(D.frequencies == ThetaHz),indsample(D,EWin(1)):indsample(D,EWin(2)),indtrial(D, Conds{iCond}, 'GOOD')),[2 3 4]);
        
        spm_eeg_plotScalpData(Topo4Hz(:,iCond,iSub),XY,D.chanlabels(ChanInds),in);
        title([sub_nam ' ' Conds{iCond} ' ' chantype]) 
        pause
        
        if strcmp(Conds{iCond}, 'SyncDelta')
           Topo1_7Hz(:,iSub) = mean(D(:,find(D.frequencies == DeltaHz),indsample(D,EWin(1)):indsample(D,EWin(2)),indtrial(D, Conds{iCond}, 'GOOD')),[2 3 4]);
        
%             spm_eeg_plotScalpData(Topo1_7Hz(:,iSub),XY,D.chanlabels(ChanInds),in);
%             title([sub_nam ' ' Conds{iCond} ' ' chantype]) 
%             pause
        end
    end
end

% !!! Sub-21 channel weights seem to have an error
Subs = 1:32;
Subs = setdiff(Subs,21);

% plot averages
for iCond = 1:length(Conds)
    
    if strcmp(Conds{iCond}, 'SyncDelta')
        
        spm_eeg_plotScalpData(mean(Topo1_7Hz(:,Subs),2),XY,D.chanlabels(ChanInds),in);
    else
        spm_eeg_plotScalpData(mean(Topo4Hz(:,iCond,Subs),[2 3]),XY,D.chanlabels(ChanInds),in);        
    end
    title([Conds{iCond} ' ' chantype])
end

%% Extract and save theta power
% Theta 4 Hz: Only 4 Hz power
% Relative 4 Hz: 4 Hz power relative to 3.5 and 4.5 Hz power
% Broad Theta power: average power over 3 to 8 Hz freqs

% !!! Sub-21 channel weights seem to have an error
Subs = 1:32;
Subs = setdiff(Subs,21);

for iSub = 1:length(Subs)
    
    % Load baseline corrected TF data
    sub_nam = sprintf('sub-%02d', Subs(iSub));
    sub_in = fullfile(outdir,sub_nam,'meg');
    cd(sub_in)
    
    D = spm_eeg_load(fullfile(sub_in,sprintf('rtf_cesub-%02d_task-TIME_run-1_meg_proc-sss.mat',Subs(iSub))));
    
    % Omit bad channels and trials
    ChanInds = indchantype(D,chantype,'GOOD');
    [~,badinds] = intersect(indchantype(D,chantype),D.badchannels);
    TrialInds = indtrial(D, Conds, 'GOOD');
    
    %% Get Theta power 
    Table_info = {'ParticipantID', 'TrialNum', 'Theta4Hz'};
    
    Theta4Hz_pTrial = cell2table(cell(length(TrialInds)...
        ,size(Table_info,2)),'VariableNames',Table_info);
    
    Theta4HzP = squeeze(mean(D(ChanInds,find(D.frequencies == ThetaHz),indsample(D,EWin(1)):indsample(D,EWin(2)),TrialInds), [1 2 3]));
    
    Theta4Hz_pTrial.Theta4Hz = Theta4HzP;
        
    Theta4Hz_pTrial.TrialNum = TrialInds';   
    
    % save to table
    
    Theta4Hz_pTrial.ParticipantID = repmat(Subs(iSub),[length(Theta4Hz_pTrial.TrialNum) 1]);

    Theta4Hz_All{iSub} = Theta4Hz_pTrial;
    
    %% Get Relative to 3.5 and 4.5 Theta
    Table_info = {'ParticipantID', 'TrialNum', 'RelTheta4Hz'};
    
    RelTheta4Hz_pTrial = cell2table(cell(length(TrialInds)...
        ,size(Table_info,2)),'VariableNames',Table_info);
    
    p3_5Hz = squeeze(mean(D(ChanInds,find(D.frequencies == RelThetaBound(1)),indsample(D,EWin(1)):indsample(D,EWin(2)),TrialInds), [1 2 3]));
    p4_5Hz = squeeze(mean(D(ChanInds,find(D.frequencies == RelThetaBound(2)),indsample(D,EWin(1)):indsample(D,EWin(2)),TrialInds), [1 2 3]));
    
    RelTheta4HzP = Theta4HzP ./ (mean([p3_5Hz; p4_5Hz],1));
    
    RelTheta4Hz_pTrial.RelTheta4Hz = RelTheta4HzP;
        
    RelTheta4Hz_pTrial.TrialNum = TrialInds';
    
    % save to table
    
    RelTheta4Hz_pTrial.ParticipantID = repmat(Subs(iSub),[length(RelTheta4Hz_pTrial.TrialNum) 1]);

    RelTheta4Hz_All{iSub} = RelTheta4Hz_pTrial;
    
    %% Get broad theta [3 8]
    Table_info = {'ParticipantID', 'TrialNum', 'BroadTheta'};
    
    BroadTheta_pTrial = cell2table(cell(length(TrialInds)...
        ,size(Table_info,2)),'VariableNames',Table_info);
    
    Ind3Hz = find(D.frequencies == ThetaBroad(1));
    Ind8Hz = find(D.frequencies == ThetaBroad(2));
    
    BroadTheta = squeeze(mean(D(ChanInds,Ind3Hz:Ind8Hz,indsample(D,EWin(1)):indsample(D,EWin(2)),TrialInds), [1 2 3]));
    
    BroadTheta_pTrial.BroadTheta = BroadTheta;
        
    BroadTheta_pTrial.TrialNum = TrialInds';
    
    % save to table
    
    BroadTheta_pTrial.ParticipantID = repmat(Subs(iSub),[length(BroadTheta_pTrial.TrialNum) 1]);

    BroadTheta_All{iSub} = BroadTheta_pTrial;
    
end

Theta4Hz_pSub = [];
RelTheta4Hz_pSub = [];
BroadTheta_pSub = [];

for iSub = 1:length(Subs)
    
    Theta4Hz_pSub = [Theta4Hz_pSub; Theta4Hz_All{iSub}];
    RelTheta4Hz_pSub = [RelTheta4Hz_pSub; RelTheta4Hz_All{iSub}];
    BroadTheta_pSub = [BroadTheta_pSub; BroadTheta_All{iSub}];
    
end

cd /imaging/henson/TIME/LatestMEGscripts/MEGxBehavioural
save('Theta4Hz.mat','Theta4Hz_pSub')
save('RelTheta4Hz.mat','RelTheta4Hz_pSub')
save('BroadTheta.mat','BroadTheta_pSub')


%% File deleter
% !Careful! removes data that are not needed anymore

% for iSub = 1:length(Subs)
%     
%     sub_nam = sprintf('sub-%02d', Subs(iSub));
%     sub_in = fullfile(outdir,sub_nam,'meg');
%     cd(sub_in)
%     % remove unnecessary files
%     delete *.mat *.dat osl.conf
%     
% end