% Analysis of MISC channels
clear
close all

%% Select Subjects
% latest recovered #8 did not have proper MEG settings loaded.
% #9 has very noisy diode, possibly misplaced diode.
% #22 has same issue as #9
Subs = [1:7 10:21 23:32];

% Set relevant time window and freq window to filter the data
Twin = [0 3];
FreqBand = [3.9 4.1];

% load MISC data
load('MISCData_Trialavg.mat')
load('MISCData_pTrial.mat')

% reset variables
meanPhaseDiff = []; pTrial_PhaseDiffs = [];
TrialAvg_PhaseDiff = []; TrialAvg_meanPhaseDiff = [];
Phasediff_pTrial = [];

% loop through subs
for iSub = 1:length(Subs)
    
    MISC_pTrial = MISCData_pTrial{Subs(iSub)};
    MISC_Trialavg = MISCData_Trialavg{Subs(iSub)};
    Conds = MISC_pTrial.condlist;
    
    for iCond = 1:2 % only sync and async
        
        % Get trial indices and index relevant epochs
        TrialInds = indtrial(MISC_pTrial,Conds{iCond});
        pTrial_VisSignal = squeeze(MISC_pTrial(indchannel(MISC_pTrial,'MISC008'),:,TrialInds));
        pTrial_AudSignal = squeeze(MISC_pTrial(indchannel(MISC_pTrial,'MISC006'),:,TrialInds));
        
        TrialInds = indtrial(MISC_Trialavg,Conds{iCond});
        meanRawSignal = squeeze(MISC_Trialavg(:,:,TrialInds));
        
        if length(TrialInds) > 1
            meanRawSignal = reshape(meanRawSignal,[size(meanRawSignal,1) size(meanRawSignal,2)*size(meanRawSignal,3)]);
        end
        meanRawSignal = meanRawSignal';
        
        % Take the envelope of only the audio as visual signal is already the envelope
        pTrial_AudSignal = abs(hilbert(pTrial_AudSignal));
        
        SignalEnv(:,1) = abs(hilbert(meanRawSignal(:,1)));
        SignalEnv(:,2) = meanRawSignal(:,2);
        
        % Add tube delay (10 ms)!!!!!!!!!
        TubeDelay = round(MISC_pTrial.fsample*10/1000);
        SoundDelayVec = circshift(pTrial_AudSignal(:,:), TubeDelay, 1);
        SoundDelayVec(1:TubeDelay,:) = 0;
        pTrial_AudSignal = SoundDelayVec;
        
        SoundDelayVec = circshift(SignalEnv(:,1), TubeDelay, 1);
        SoundDelayVec(1:TubeDelay,:) = 0;
        SignalEnv(:,1) = SoundDelayVec;
        
        % Band-pass filter data
        fSignalEnv = lowpass(SignalEnv,FreqBand(2),MISC_pTrial.fsample,'ImpulseResponse','iir');
        fSignalEnv = highpass(fSignalEnv,FreqBand(1),MISC_pTrial.fsample,'ImpulseResponse','iir');
        
        fpTrial_VisSignal = lowpass(pTrial_VisSignal,FreqBand(2),MISC_pTrial.fsample,'ImpulseResponse','iir');
        fpTrial_VisSignal = highpass(fpTrial_VisSignal,FreqBand(1),MISC_pTrial.fsample,'ImpulseResponse','iir');
        
        fpTrial_AudSignal = lowpass(pTrial_AudSignal,FreqBand(2),MISC_pTrial.fsample,'ImpulseResponse','iir');
        fpTrial_AudSignal = highpass(fpTrial_AudSignal,FreqBand(1),MISC_pTrial.fsample,'ImpulseResponse','iir');
        
        % get the relevant timepoints
        TimeInds = [indsample(MISC_pTrial,Twin(1)):indsample(MISC_pTrial,Twin(2))]';
        ati = TimeInds;
        if length(TrialInds) > 1
            for t = 2:length(TrialInds)
                ati = [ati; TimeInds + MISC_Trialavg.nsamples];
            end
        end
        
        % hilbert transform
        hfpTrial_VisSignal = hilbert(fpTrial_VisSignal);
        hfpTrial_AudSignal = hilbert(fpTrial_AudSignal);
        hfSignalEnv = hilbert(fSignalEnv);
        
        % index relevant timepoints
        meanRawSignal = meanRawSignal(ati,:);
        hfpTrial_VisSignal = hfpTrial_VisSignal(TimeInds,:);
        hfpTrial_AudSignal = hfpTrial_AudSignal(TimeInds,:);
        hfSignalEnv = hfSignalEnv(ati,:);
        
        % Take phase diff
        TrialAvg_PhaseDiff(iSub,:,iCond) = angle(hfSignalEnv(:,2)) - angle(hfSignalEnv(:,1));
        TrialAvg_meanPhaseDiff(iSub,iCond) = angle(mean(exp(sqrt(-1)*TrialAvg_PhaseDiff(iSub,:,iCond))))*180/pi;
        
        pTrial_PhaseDiffs = angle(hfpTrial_VisSignal) - angle(hfpTrial_AudSignal);
        
        for iTrial = 1:size(pTrial_PhaseDiffs,2)
            Phasediff_pTrial(iTrial,iCond,iSub) = angle(mean(exp(sqrt(-1)*pTrial_PhaseDiffs(:,iTrial))))*180/pi; % mean phase difference
        end
        
        meanPhaseDiff(iSub,iCond) = angle(mean(exp(sqrt(-1)*Phasediff_pTrial(:,iCond,iSub)*pi/180)))*180/pi;
        
        if length(TrialInds) > 1
            pst = 1:length(ati);
        else
            pst = MISC_Trialavg.time(TimeInds);
        end
        
        % plot raw signal
%         figure(fnum(iCond)),clf
%         plot(pst,zscore(meanRawSignal)),legend('Aud','Vis'),title([sub_nam ' ' Conds{iCond}])
%         % plot phase info
%         figure(fnum(iCond+2)),clf
%         plot(pst,zscore(hfd)),legend('Aud','Vis'),title([sub_nam ' ' Conds{iCond}])
    end
    % pause to inspect plots
%     drawnow
%     pause
end

% get grand means
TrialAvg_gmean_PhaseDiff = angle(mean(exp(sqrt(-1)*TrialAvg_meanPhaseDiff*pi/180)))*180/pi
TrialAvg_meanTimeDiff = (1000/4)*(TrialAvg_gmean_PhaseDiff/360)

gmean_PhaseDiff = angle(mean(exp(sqrt(-1)*meanPhaseDiff*pi/180)))*180/pi
meanTimeDiff = (1000/4)*(gmean_PhaseDiff/360)

% Compass plots
for iCond = 1:2
    z = exp(sqrt(-1)*TrialAvg_meanPhaseDiff(:,iCond)*pi/180);
    figure,compass(z)
    for s = 1:length(Subs)
        line([0 real(z(s))]', [0 imag(z(s))]')
    end
    line([0 real(mean(z))]', [0 imag(mean(z))]','Color',[1 0 0],'LineWidth',2)
    title(Conds{iCond})
end

% histogram of each participants' trials
% for iSub = 1:length(Subs)
%     for iCond = 1:2
%         figure, histogram(Phasediff_pTrial(:,iCond,iSub))
%         title([num2str(iSub) ' ' Conds{iCond}])
%     end
% end

% plot phase against time
% pdss = squeeze(mean(TrialAvg_PhaseDiff(:,:,1),1));
% figure, bar(ati,angle(exp(sqrt(-1)*pdss))*180/pi), title('sync')
% 
% pdss = squeeze(mean(TrialAvg_PhaseDiff(:,:,2),1));
% figure, bar(ati,angle(exp(sqrt(-1)*pdss))*180/pi), title('async')