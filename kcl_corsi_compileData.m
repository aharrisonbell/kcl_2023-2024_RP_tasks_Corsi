%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kcl_corsi_compileData.m %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% started by AHB, Feb 2024
% v1.0 - first draft


%% Load data (Daria - you will need to update these folders OR make sure you run this program from within the same directory as the dat files are located
if ispc % if this program is run on a windows PC
    rootdir='C:\Users\K...\OneDrive - King''s College London\MATLAB';
else % if this program is run on a macbook
    rootdir='~/OneDrive - King''s College London/MATLAB/currentProjects/proj_Corsi';
end

%% find datafiles
datafiles = dir('kcl_corsi_ppt_*.mat'); % find all datafiles that are NOT intermediate

disp(['Found ',num2str(numel(datafiles)), ' datafiles'])
summaryData = []; % initialise empty matrix
% SummaryData Structure
% 1)  Participant Number
% 2)  BlockNumber
% 3)  TrialNumber
% 4)  SequenceLength
% 5)  NumberOfCorrectBlocks
% 6)  Trial Performance (% of total sequence correct)
% 

% 11)  TimetoSelectBlock1
% 12)  TimetoSelectBlock2
% 13)  TimetoSelectBlock3
% 14)  TimetoSelectBlock4
% 15) TimetoSelectBlock5
% 16) TimetoSelectBlock6
% 17) TimetoSelectBlock7
% 18) TimetoSelectBlock8



%% block number and sequence number, number of blocks per sequence, number of blocks correctly reproduced per sequence, time taken to reproduce each sequence. Is this all you need from us and does this look okay?
for dd = 1:numel(datafiles)
    clear respMat temp_*
    load(datafiles(dd).name)
    temp_numTrials = numel(respMat); % total number of trials in session
    temp_summaryData = nan(temp_numTrials, 13);
    
    temp_summaryData(:,1) = respMat.participantNumber;

    for tt = 1:temp_numTrials
        temp_summaryData(tt,2) = respMat(tt).currBlock;
        temp_summaryData(tt,3) = respMat(tt).totalTrialExperiment;
        temp_summaryData(tt,4) = respMat(tt).trialSequenceLength;
        if respMat(tt).errorIndex == 0 % perfect
            temp_summaryData(tt,5) = temp_summaryData(tt,4);
        else
            temp_summaryData(tt,5) = respMat(tt).errorIndex - 1;
        end

        temp_summaryData(tt,6) = temp_summaryData(tt,5) / respMat(tt).trialSequenceLength * 100;

        % Click Times
        temp_numClicks = length(respMat(tt).clickTimes);
        temp_summaryData(tt,11:11+temp_numClicks-1) = respMat(tt).clickTimes/1000;
    end

    summaryData = [summaryData; temp_summaryData];

end