%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kcl_corsi_compileData.m %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% started by AHB, Feb 2024
% v1.0 - first draft
% v1.1 - added some descriptives


%% Load data (Esha/Nikita - you will need to update these folders OR make sure you run this program from within the same directory as the dat files are located
if ispc % if this program is run on a windows PC
    rootdir='C:\Users\K...\OneDrive - King''s College London\MATLAB';
else % if this program is run on a macbook
    rootdir='~/OneDrive - King''s College London/MATLAB/currentProjects/proj_Corsi/corsiData';
end

%% find datafiles
datafiles = dir([rootdir, filesep, 'kcl_corsi_ppt_*.mat']); % find all datafiles that are NOT intermediate
disp(['Found ',num2str(numel(datafiles)), ' datafiles'])

%% Create large summaryData matrix with all GOOD data

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
rawData = []; % initialise empty matrix
for dd = 1:numel(datafiles)
    clear respMat temp_*
    load([rootdir, filesep, datafiles(dd).name])
    temp_numTrials = numel(respMat); % total number of trials in session
    temp_summaryData = nan(temp_numTrials, 13);

    % Correct Participant Number
    %temp_participantNumber = str2num(datafiles(dd).name(14:22));
    %if respMat(1,1) ~= temp_participantNumber
    %    disp([datafiles(dd).name,' - Participant numbers don''t match. Switching to number in FILENAME'])
    %    respMat(:,1) = temp_participantNumber;
    %end
    % Participant Number
    disp([datafiles(dd).name,' - Participant Number: ', num2str(respMat(1).participantNumber)])

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

    rawData = [rawData; temp_summaryData];
end

%% Save Raw Data
writematrix(rawData, [rootdir, filesep, 'kcl_corsi_RawData_', date, '_', num2str(size(datafiles, 1)), 'files.csv']) %#ok<*DATE>
writecell({datafiles.name}', [rootdir, filesep, 'kcl_corsi_RawData_', date, '_list_of_files.csv'])

%% Scroll through each participant to generate summary statistics
individualParticipants = unique(rawData(:,1));


% Structure of summaryData
% 1)  Participant Number
% 2)  GROUP NUMBER (Daria - you will need to sort this out)
% 3)  Block Number
% 4)  average per trial performance
% 5)  std per trial performance
% 6,7) mean/std % Correct/Block : Sequence Length 5
% 8,9) mean/std % Correct/Block : Sequence Length 5
% 10,11)mean/std % Correct/Block : Sequence Length 5
% 12,13)mean/std % Correct/Block : Sequence Length 5


%%
summaryData = []; % initialise matrix
for pp = 1:length(individualParticipants)
    temp_summaryData =nan(4, 13);
    disp(['Analysing participant ',num2str(individualParticipants(pp)),'...'])
    tempData = rawData(rawData(:,1) == individualParticipants(pp),:);
    
    temp_summaryData(:,1) = tempData(1,1); % paste participant number
    for bb = 1:length(unique(tempData(:,2))) % scroll through each block
        temp_summaryData(bb, 3) = bb; % block number
        temp_summaryData(bb, 4) = mean(tempData(tempData(:,2) == bb, 6));
        temp_summaryData(bb, 4) = std (tempData(tempData(:,2) == bb, 6));
        
        % per sequence length (difficulty)

        seqlengths = unique(tempData(:,4));
        for cc = 1:length(seqlengths) % scroll through each difficulty (coherence) level
            temp_summaryData(bb, 4+(cc*2)) = mean(tempData(tempData(:,2) == bb & tempData(:,4) == seqlengths(cc), 6));
            temp_summaryData(bb, 5+(cc*2)) = std (tempData(tempData(:,2) == bb & tempData(:,4) == seqlengths(cc), 6));
        end
    end

    summaryData = [summaryData; temp_summaryData];
    clear temp_summaryData

end
summaryData(1:4,3:end)

%% Plot Data
figure
hold on
errorbar(5:8, mean(summaryData(summaryData(:,3)==1, [6, 8, 10, 12])), ...
    sem(summaryData(summaryData(:,3)==1, [6, 8, 10, 12])), 'rs-', 'LineWidth', 2); % block 1
errorbar(5:8, mean(summaryData(summaryData(:,3)==2, [6, 8, 10, 12])), ...
    sem(summaryData(summaryData(:,3)==2, [6, 8, 10, 12])), 'bs-', 'LineWidth', 2); % block 2
errorbar(5:8, mean(summaryData(summaryData(:,3)==3, [6, 8, 10, 12])), ...
    sem(summaryData(summaryData(:,3)==3, [6, 8, 10, 12])), 'gs-', 'LineWidth', 2); % block 3
errorbar(5:8, mean(summaryData(summaryData(:,3)==4, [6, 8, 10, 12])), ...
    sem(summaryData(summaryData(:,3)==4, [6, 8, 10, 12])), 'ms-', 'LineWidth', 2); % block 4
legend('Block 1','Block 2',' Block 3','Block 4')

title({'Performance as a function of BLOCK and Sequence Length','(mean +/- sem)'}, 'FontSize', 14)
xlabel('Difficulty (Sequence Length)', 'FontSize', 12)
ylabel('Performance (% Correct)', 'FontSize', 12)

writematrix(summaryData, [rootdir, filesep, 'kcl_corsi_SummaryData_', date, '_', num2str(size(datafiles, 1)), 'files.csv']) %#ok<*DATE>

