%% Sensor-space Source Localisation
% The part of the script calculates auditory and visual activity weights
% from the sensor-space and saves them as .mat format.
% 
% Most plotting codes are commented. Uncomment as necessary.
%
% See the second part for the application of the weights to calculate the
% phase lag between visual and auditory activity.
%
% The script is written to be compatible for MRC CBU computing cluster.
% Change paths as necessary to to link data and spm.

clear

% set paths and initialize spm
outdir  = '/imaging/henson/TIME/meg-derivatives/';
addpath('/imaging/henson/TIME/MEG_scripts')

addpath '/imaging/henson/TIME/LatestMEGscripts'
addpath /imaging/local/software/spm_cbu_svn/releases/spm12_latest/
spm eeg

%chantype = 'MEGMAG';
chantype = 'MEGPLANAR';
in.chantype = chantype;

% Set time windows: 0 is the onset of stimuli
TWin = [-1 4]; % whole epoch time window in secs
EWin = [0.5 3]; % entraiment window in secs

% Set frequencies
ThetaHz = 4;
RelThetaBound = [3 5];
DeltaHz = 1.7;

% change to use trans default or only trans maxfiltered data
%transdef = 'td'
transdef = ''

% Note that if using transdef, xy coords seem updated, so topos better if
% use xy coords from non-transed data file, so take, eg sub32
if ~isempty(transdef)
    D = spm_eeg_load(fullfile(outdir,sprintf('sub-32/meg/esub-32_task-loc_meg_proc-sss.mat')));
    XY = D.coor2D;
    [~,XYinds] = intersect(indchantype(D,{'MEG','MEGPLANAR'}),indchantype(D,chantype));
    XY = XY(:,XYinds);
end

Subs = 1:32;
Subs = setdiff(Subs,11); % s11 cannot be maxfiltered unless downsampled

%% Get channel weights for visual and auditory activity
% auditory weights from localizer
% visual from encoding phase of main runs

figa=figure;
figv=figure;
figap=figure;
figvp=figure;

maxFreq = nan(length(Subs),2);
AudioWeights = {};
VisuoAudioWeights = {};
VisuoOrthAudioWeights = {};

for iSub = 1:length(Subs)
    
    % change to subject directory
    sub_nam = sprintf('sub-%02d', Subs(iSub));
    sub_in = fullfile(outdir,sub_nam,'meg');
    cd(sub_in)
    
    % Extract Auditory weights if there is localizer
    % subs 4:11 do not have the localizer
    
    % load audio localizer
    loc_file = sprintf('esub-%02d_task-loc_meg_proc-sss%s.mat',Subs(iSub),transdef);
    
    if exist(loc_file)
        
        % load data
        D = spm_eeg_load(loc_file);
        
        % get trial and channel indices
        % two different names for the trials due to change during collection
        TrialInds = indtrial(D, {'sound','AudLocSyncTheta'}, 'GOOD');
        ChanInds = indchantype(D,chantype);
        
        % index data
        Y = D(ChanInds,(indsample(D,EWin(1))+1):indsample(D,EWin(2)),TrialInds);
        % concetanete trials for freq resolution
        Y = reshape(Y,size(Y,1),size(Y,2)*size(Y,3))';
        
        %figure,plot(y);
        
        % use the custom fft function
        [f, p, a] = pow_spec(Y,1/D.fsample,0);
        %        figure(fap),clf,plot(f,p),xlim([3 5]);
        
        % average over channels to check if 4Hz has peak power in theta band
        wp = mean(p,2);
        %         figure(figap),clf,plot(f,wp),xlim([3 5]);
        %         title(['Sub-' num2str(Subs(iSub)) ' Audio PowSpec'])
        
        indFreq = find(f>RelThetaBound(1) & f<RelThetaBound(2)); [~,indMaxFreq] = max(wp(indFreq));
        maxFreq(Subs(iSub),1) = f(indFreq(indMaxFreq));
        
        df = round(f,3);
        
        indFreq = find(df == ThetaHz);
        
        % Get the 4Hz power relative to the frequencies between 3 and 4
        AudioWeights{Subs(iSub)} = mean(p(indFreq,:),1) ./ mean(p(find(df>RelThetaBound(1) & df<RelThetaBound(2) & df~=ThetaHz),:),1);
        
        wp = p*AudioWeights{Subs(iSub)}';
        %         figure(fap),clf,plot(f,wp),xlim([3 5]); % biased to show 4 Hz!
        
        %         if isempty(transdef)
        %             XY = D.coor2D;
        %             [~,XYinds] = intersect(indchantype(D,{'MEG','MEGPLANAR'}),ChanInds);
        %             XY = XY(:,XYinds);
        %         end
        %         in.f = figa.Number; figure(in.f); clf
        %         [ZI,f] = spm_eeg_plotScalpData(AudioWeights{Subs(iSub)}',XY,D.chanlabels(ChanInds),in);
        %         title(['Sub-' num2str(Subs(iSub)) ' Audio ' chantype])
    else
        figure(figa); clf; figure(figap); clf;
    end
    
    % Load main run for localising visual activity
    D = spm_eeg_load(sprintf('cedsub-%02d_task-TIME_run-1_meg_proc-sss%s.mat',Subs(iSub),transdef));
    TrialInds = indtrial(D, {'SyncTheta'; 'AsyncTheta'}, 'GOOD');
    ChanInds = indchantype(D,chantype);
    
    Y = D(ChanInds,(indsample(D,EWin(1))+1):indsample(D,EWin(2)),TrialInds);
    Y = reshape(Y,size(Y,1),size(Y,2)*size(Y,3))';
    %figure,plot(y);
    
    % use the custom fft function
    [f,p,a] = pow_spec(Y,1/D.fsample,0);
    %     %figure(fvp),clf,plot(f,p),xlim([3 5]);
    
    % average over channels to check if 4Hz has peak power in theta band
    wp = mean(p,2);
    %     figure(figvp),clf,plot(f,wp),xlim([3 5]);
    %     title(['Sub-' num2str(Subs(iSub)) ' Video PowSpec'])
    
    indFreq = find(f>3 & f<5); [~,indMaxFreq] = max(wp(indFreq));
    maxFreq(Subs(iSub),2) = f(indFreq(indMaxFreq));
    
    df = round(f,3);
    
    indFreq = find(df == ThetaHz);
    
    % Get the 4Hz power relative to the frequencies between 3 and 4
    VisuoAudioWeights{Subs(iSub)} = mean(p(indFreq,:),1) ./ mean(p(find(df>RelThetaBound(1) & df<RelThetaBound(2) & df~=ThetaHz),:),1);
    
    wp = p*VisuoAudioWeights{Subs(iSub)}';
    %     figure(fvp),clf,plot(f,wp),xlim([3 5]); % biased to show 4 Hz!
    
    %     if isempty(transdef)
    %         XY = D.coor2D;
    %         [~,XYinds] = intersect(indchantype(D,{'MEG','MEGPLANAR'}),ChanInds);
    %         XY = XY(:,XYinds);
    %     end
    %     in.f = figv.Number; figure(in.f); clf
    %     [ZI,f] = spm_eeg_plotScalpData(VisuoAudioWeights{Subs(iSub)}',XY,D.chanlabels(ChanInds),in);
    %     title(['Sub-' num2str(Subs(iSub)) ' Visuo ' chantype])
    
    %     drawnow
    %     pause
end

maxFreq

% Subjects with no Loc 4Hz (or no Loc at all)
find(round(maxFreq(:,1)*100)/100 ~= ThetaHz)
% Subjects with main 4Hz
find(round(maxFreq(:,2)*100)/100 ~= ThetaHz)

% Normalise for total power
for iSub = 1:length(Subs)
    AudioWeights{Subs(iSub)} = AudioWeights{Subs(iSub)}/sum(AudioWeights{Subs(iSub)});
    VisuoAudioWeights{Subs(iSub)} = VisuoAudioWeights{Subs(iSub)}/sum(VisuoAudioWeights{Subs(iSub)});
end

mAudioWeights = mean(cat(1,AudioWeights{:}),1);

%% Plot subject-average weights
figa=figure;
figv=figure;
% Use coords from subject 32 non-trans
D = spm_eeg_load(fullfile(outdir,sprintf('sub-30/meg/esub-30_task-loc_meg_proc-sss.mat')));
XY = D.coor2D;
ChanInds = indchantype(D,chantype);
[~,XYinds] = intersect(indchantype(D,{'MEG','MEGPLANAR'}),ChanInds);
XY = XY(:,XYinds);
in.f = figa.Number; figure(in.f); clf
[ZI,f] = spm_eeg_plotScalpData(mAudioWeights',XY,D.chanlabels(ChanInds),in);
title(['Mean Auditory Weights ' chantype])

mVisuoAudioWeights = mean(cat(1,VisuoAudioWeights{:}),1);
in.f = figv.Number; figure(in.f); clf
[ZI,f] = spm_eeg_plotScalpData(mVisuoAudioWeights',XY,D.chanlabels(ChanInds),in);
title(['Mean Visual Weights ' chantype])

% orthgonal seems to make it worse
mVisuoOrthAudioWeights = orthog(mVisuoAudioWeights', mAudioWeights')';
mVisuoOrthAudioWeights = mVisuoOrthAudioWeights - min(mVisuoOrthAudioWeights); % remove negative
in.f = figv.Number; figure(in.f); clf
[ZI,f] = spm_eeg_plotScalpData(mVisuoOrthAudioWeights',XY,D.chanlabels(ChanInds),in);
title(['Mean audioOrth ' chantype])

% Care - Orth weights can now be negative...

for iSub = 1:length(Subs)
    if ~isempty(transdef)
        
        AudioWeights{Subs(iSub)} = mAudioWeights;
        VisuoAudioWeights{Subs(iSub)} = mVisuoAudioWeights;
        VisuoOrthAudioWeights{Subs(iSub)} = mVisuoOrthAudioWeights;
        
    else
        
        if isempty(AudioWeights{Subs(iSub)})
            AudioWeights{Subs(iSub)} = mAudioWeights;
        end
        
        %     xy = D.coor2D;
        %     [~,xyinds] = intersect(indchantype(D,{'MEG','MEGPLANAR'}),chinds);
        %     xy = xy(:,xyinds);
        %     in.f = fa.Number; figure(in.f); clf
        %     [ZI,f] = spm_eeg_plotScalpData(AudioWeights{iSub}',xy,D.chanlabels(chinds),in);
        %     title(['Sub-' num2str(sub_nums(iSub)) ' Audio ' chantype])
        
        VisuoOrthAudioWeights{Subs(iSub)} = orthog(VisuoAudioWeights{Subs(iSub)}', AudioWeights{Subs(iSub)}')';
        
        % %     xy = D.coor2D;
        % %     [~,xyinds] = intersect(indchantype(D,{'MEG','MEGPLANAR'}),chinds);
        % %     xy = xy(:,xyinds);
        % %     in.f = fa.Number; figure(in.f); clf
        % %     [ZI,f] = spm_eeg_plotScalpData(VisuoAudioWeights{iSub}',xy,D.chanlabels(chinds),in);
        % %     title(['Sub-' num2str(sub_nums(iSub)) ' Visuo ' chantype])
        
        %     in.f = fv.Number; figure(in.f); clf
        %     [ZI,f] = spm_eeg_plotScalpData(VisuoOrthAudioWeights{iSub}',xy,D.chanlabels(chinds),in);
        %     title(['Sub-' num2str(sub_nums(iSub)) ' Visuo Orth Audio ' chantype])
        
        %     drawnow
        %     pause(1)
    end
end

cd('/imaging/henson/TIME/meg-derivatives')
save(sprintf('AudioWeights%s',chantype),'AudioWeights')
save(sprintf('VisuoAudioWeights%s',chantype),'VisuoAudioWeights')
save(sprintf('VisuoOrthAudioWeights%s',chantype),'VisuoOrthAudioWeights')

%%  Apply Channel Weights to data and calculate phase diff
clear

outdir  = '/imaging/henson/TIME/meg-derivatives/';

Subs = 1:32;
Subs = setdiff(Subs,11); % s11 cannot be maxfiltered unless downsampled
% Subs = [1:3 13:32];
transdef = '';

chantype = 'MEGPLANAR';
in.chantype = chantype;

% load weights
load('/imaging/henson/TIME/meg-derivatives/AudioWeightsMEGPLANAR.mat')
load('/imaging/henson/TIME/meg-derivatives/VisuoAudioWeightsMEGPLANAR.mat')
load('/imaging/henson/TIME/meg-derivatives/VisuoOrthAudioWeightsMEGPLANAR.mat')

Conds = {'SyncTheta','AsyncTheta'};

% time windows
Twin = [-1 4];
EWin = [0.5 3];

HilbFreqBand = [3.5 4.5];

% figah = figure;
% figvh = figure;
% figvah = figure;
% figap = figure;
% figvp = figure;

FFTphs = [];
VAphsdif = nan(length(Subs),2);
SApd = [];
WithinModality = [];

Table_info = {'ParticipantID', 'TrialNum', 'Synchronicity'};

for iSub = 1:length(Subs)
    iSub
    % change to subject directory
    sub_nam = sprintf('sub-%02d', Subs(iSub));
    sub_in = fullfile(outdir,sub_nam,'meg');
    cd(sub_in)
    
    % load data
    D = spm_eeg_load(sprintf('cedsub-%02d_task-TIME_run-1_meg_proc-sss%s.mat',Subs(iSub),transdef));
    
    ChanInds = indchantype(D,chantype,'GOOD');
    [~,badinds] = intersect(indchantype(D,chantype),D.badchannels);
    
    TrialInds = indtrial(D, Conds,'GOOD');
    
    if Subs(iSub) == 31
        ay = D(ChanInds,(indsample(D,Twin(1))+1):indsample(D,Twin(2)),TrialInds);
        TrialInds = find(sum(ay,[1 2]))';
    end
    
    Sync_pTrial = cell2table(cell(length(TrialInds)...
        ,size(Table_info,2)),'VariableNames',Table_info);
    
    hfAud = {}; hfVis = {}; VAphsdif_pSub = {};
    
    signflip1 = []; % IMPORTANT, ie any flip applied to all conditions
    signflip2 = []; % IMPORTANT, ie any flip applied to all conditions
    
    for iCond = 1:length(Conds)
        
        TrialInds = indtrial(D, Conds{iCond}, 'GOOD');
        ay = D(ChanInds,(indsample(D,Twin(1))+1):indsample(D,Twin(2)),TrialInds);
        
        if Subs(iSub) == 31
            TrialInds = find(sum(ay,[1 2]))';
        end
        
        AudS = []; VisS = []; AudPw = []; AudPA = []; VisPw = []; VisPA = [];
        for iTrial = 1:length(TrialInds)
            
            Y = squeeze(ay(:,:,iTrial))';
            
            % Currently selecting the average sensor weights
            if any(Subs(iSub) == 4:12)
                
                ChanWeights = mean(cat(1,AudioWeights{:}),1)';
            else
                %                 ChanWeights = AudioWeights{Subs(iSub)}';
                ChanWeights = mean(cat(1,AudioWeights{:}),1)';
            end
            ChanWeights(badinds) = [];
            ChanWeights = ChanWeights/sum(abs(ChanWeights));
            [~,iMaxchan] = max(ChanWeights);
            
            % uncomment this to use the heighest weight channel only rather
            % than a distribution of weights across channels
            ChanWeights = zeros(length(ChanWeights),1);
            ChanWeights(iMaxchan) = 1;
            
            % no need to flip channels when using max weight chan
            %             for iChan = 1:length(ChanWeights)
            %                 % correlate each channel with the maximum weight channel
            %                 ChanSign = sign(corr(Y(:,iChan),Y(:,iMaxchan)));
            %                 % reverse the weight of the channel if correlation is negative
            %                 ChanWeights(iChan) = ChanSign*ChanWeights(iChan);
            %             end
            
            
            AudS(:,iTrial) = Y * ChanWeights;
            
            %             [f,AudPw(:,iTrial),AudPA(:,iTrial)] = pow_spec(AudS(:,iTrial),1/D.fsample,0);
            
            % Currently selecting the average sensor weights
%             ChanWeights = VisuoAudioWeights{Subs(iSub)}';
            % uncomment to use average weight for each participant
            %             ChanWeights = mean(cat(1,VisuoAudioWeights{:}),1)';
            
            ChanWeights = VisuoOrthAudioWeights{iSub}';
            ChanWeights(badinds) = [];
            ChanWeights = ChanWeights/sum(abs(ChanWeights));
            [~,iMaxchan] = max(ChanWeights);
            
            % uncomment this to use the heighest weight channel only rather
            % than a distribution of weights across channels
            ChanWeights = zeros(length(ChanWeights),1);
            ChanWeights(iMaxchan) = 1;
            
            % no need to flip channels when using max weight chan
            %             for iChan = 1:length(ChanWeights)
            %                 ChanSign = sign(corr(Y(:,iChan),Y(:,iMaxchan)));
            %                 ChanWeights(iChan) = ChanSign*ChanWeights(iChan);
            %             end
            
            VisS(:,iTrial) = Y * ChanWeights;
            
            %             if isempty(signflip1) % CARE: if not trial-averaged, then will be taken from first trial only...
            %                 AudViscor = corr(AudS(:,iTrial),VisS(:,iTrial))
            %                 signflip1 = sign(AudViscor);
            %             end
            %             VisS(:,iTrial) = signflip1*VisS(:,iTrial);
            
            %             [f,VisPw(:,iTrial),VisPA(:,iTrial)] = pow_spec(VisS(:,iTrial),1/D.fsample,0);
            
        end
        
        %         figure(figap),subplot(2,1,iCond), plot(f,mean(AudPw,2)),xlim([1.5 9]);
        %         title(sprintf('Aud: %s',Conds{iCond}))
        %         df = round(100*f)/100; indFreq = find(df==4); f(indFreq);
        %        FFTphs(iSub,1,c) = angle(mean(exp(1i*apa(indFreq,:))));
        
        %         figure(figvp),subplot(2,1,iCond), plot(f,mean(VisPw,2)),xlim([1.5 9]);
        %         title(sprintf('Vis: %s',Conds{iCond}))
        %         df = round(100*f)/100; indFreq = find(df==4); f(indFreq);
        %        FFTphs(iSub,2,c) = angle(mean(exp(1i*vpa(indFreq,:))));
        %         PhaseDiff = VisPA(indFreq,:) - AudPA(indFreq,:);
        
        %         a = mean(a,2);
        %         v = mean(v,2);
        
        fAud = ft_preproc_bandpassfilter(AudS',D.fsample,HilbFreqBand,5,'but','twopass','reduce')';
        fVis = ft_preproc_bandpassfilter(VisS',D.fsample,HilbFreqBand,5,'but','twopass','reduce')';
        
        % sign flip 2
        if isempty(signflip2) % CARE: if not trial-averaged, then will be taken from first trial only...
            
            AudViscor = corr(mean(fAud(:,:),2),mean(fVis(:,:),2));
            signflip2 = sign(AudViscor);
            SignFlips(Subs(iSub),iCond) = signflip2;
        end
        fVis = signflip2*fVis;
        
        %         if iCond==1
        %             figure(figah),plot([mean(fAud,2) mean(fVis,2)]); title(sprintf('Sub-%d, Sync',sub_nums(iSub))); drawnow; %xlim([1000 4000]);
        %         elseif iCond==2
        %             figure(figvh),plot([mean(fAud,2) mean(fVis,2)]); title(sprintf('Sub-%d, Async',sub_nums(iSub))); drawnow; %xlim([1000 4000]);
        %         end
        %         drawnow
        
        hfAud{iCond} = hilbert(fAud);
        hfVis{iCond} = hilbert(fVis);
        
        hfAud{iCond} = hfAud{iCond}((indsample(D,EWin(1))+1):indsample(D,EWin(2)),:);
        hfVis{iCond} = hfVis{iCond}((indsample(D,EWin(1))+1):indsample(D,EWin(2)),:);
        
        PhaseDiff = angle(hfVis{iCond}) - angle(hfAud{iCond});
        PhaseDiff = angle(exp(1i*PhaseDiff)); % unwrap
        
        %         figure(figvah),subplot(2,1,iCond), plot(PhaseDiff); ylim([-pi pi]); title(Conds{iCond}) %compass(exp(1i*pd))
        
        if length(TrialInds) > 1
            VAphsdif(iSub,iCond) = angle(mean(exp(1i*PhaseDiff),[1 2]));
            VAphsdif_pSub{iCond} = angle(mean(exp(1i*PhaseDiff), 1));
            AllTrialInds{iCond} = TrialInds;
            
        else
            VAphsdif(iSub,iCond) = angle(mean(exp(1i*PhaseDiff)));
        end
        
    end
    
    % save to table
    
    if length(TrialInds) > 1
        Sync_pTrial.TrialNum = [AllTrialInds{1}'; AllTrialInds{2}'];
        Sync_pTrial.ParticipantID = repmat(Subs(iSub),[length(Sync_pTrial.TrialNum) 1]);
        Sync_pTrial.Synchronicity = [VAphsdif_pSub{1}'; VAphsdif_pSub{2}'];
        
        if iSub == 1
            Sync_pSub = Sync_pTrial;
        else
            Sync_pSub = [Sync_pSub;Sync_pTrial];
        end
    end
    %     pause
end

% save synchronicity values
cd /imaging/henson/TIME/LatestMEGscripts/MEGxBehavioural
save('Synchronicity.mat','Sync_pSub')

SignFlips

VAphsdif
mmpd = angle(mean(exp(1i*VAphsdif)))*180/pi
mtd = (1000/4)*(mmpd/360)

%% compass plots
labels = {'Sync:Vis-Aud','Asyc:Vis-Aud'};
for iCond = 1:length(labels)
    z = VAphsdif(:,iCond);
    %     z(end+1,:) = exp(angle(mean(z))*1i);
    figure
    %     CFig = compass(z);
    %     MeanLine = CFig(end);
    %     MeanLine.LineWidth = 2;
    %     MeanLine.Color = 'r';
    polarhistogram(z,20);
    title([labels{iCond} " Mean= " num2str(round(mmpd(iCond),2))])
end