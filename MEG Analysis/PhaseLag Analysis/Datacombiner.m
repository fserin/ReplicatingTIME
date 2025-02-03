%% Combine behavioural data with MEG data for R analysis
clear

cd /imaging/henson/TIME/LatestMEGscripts/MEGxBehavioural

load Theta4Hz.mat
load RelTheta4Hz.mat
load Synchronicity.mat
load BroadTheta.mat

BehDir = '/imaging/henson/TIME/BehavioralAnalysisScripts/CombinedMemoryData';

d = dir([BehDir '/*.csv']);

for iSub = 1:length(d)
    
    % read sub data
    Sub_Data = readtable([BehDir '/' d(iSub).name]);
    
    % Add relative theta power
    
    RelThetaRows = RelTheta4Hz_pSub.TrialNum(RelTheta4Hz_pSub.ParticipantID == iSub);
    RelThetas = RelTheta4Hz_pSub.RelTheta4Hz(RelTheta4Hz_pSub.ParticipantID == iSub);
    
    Sub_Data.('RelTheta4Hz') = nan(size(Sub_Data,1),1);
    
    Sub_Data.RelTheta4Hz(RelThetaRows) = RelThetas;
    
    % Add regular theta power
    
    ThetaPowRows = Theta4Hz_pSub.TrialNum(Theta4Hz_pSub.ParticipantID == iSub);
    ThetaPows = Theta4Hz_pSub.Theta4Hz(Theta4Hz_pSub.ParticipantID == iSub);
    
    Sub_Data.('Theta4Hz') = nan(size(Sub_Data,1),1);
    
    Sub_Data.Theta4Hz(ThetaPowRows) = ThetaPows;
    
    % Add broad theta power
    
    BroadThetaPowRows = BroadTheta_pSub.TrialNum(BroadTheta_pSub.ParticipantID == iSub);
    BroadThetaPows = BroadTheta_pSub.BroadTheta(BroadTheta_pSub.ParticipantID == iSub);
    
    Sub_Data.('BroadTheta') = nan(size(Sub_Data,1),1);
    
    Sub_Data.BroadTheta(BroadThetaPowRows) = BroadThetaPows;
    
    % Add Synchronicity
    
    SyncRows = Sync_pSub.TrialNum(Sync_pSub.ParticipantID == iSub);
    Synchronicity = Sync_pSub.Synchronicity(Sync_pSub.ParticipantID == iSub);
    
    Sub_Data.('Synchronicity') = nan(size(Sub_Data,1),1);
    
    Sub_Data.Synchronicity(SyncRows) = Synchronicity;
    
    % combine all subjects into one .csv
    if iSub == 1
        Subs_Data = Sub_Data;
    else
        Subs_Data = [Subs_Data; Sub_Data];
    end
end

cd /imaging/henson/TIME/LatestMEGscripts/MEGxBehavioural
writetable(Subs_Data,'DatawTheta.csv')