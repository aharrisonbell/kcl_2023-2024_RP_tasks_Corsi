%% kcl_CorsiTask
clearvars;
% by AHB, started Oct 2023
% Version 1.3 - Feb 6, 2024 (post dress rehearsal)
%task_version = 1.4; % Version 1.4 - Feb 16, 2024 (post dress rehearsal. post tweaking)
task_version = 1.5; % Version 1.5 - Feb 20, 2024 (changed debrief)
% Developed for BSc Psychology/Neuroscience and Psychology Research Project
% Uses Psychtoolbox V3
% Corsi trials are presented in 4 blocks, 20-30 trials per block

%% Clear workspace and screen
sca;
close all;


%% Setup PTB with some default values
PsychDefaultSetup(2);
KbName('UnifyKeyNames');
KbCheck;
% Set the screen number to the external secondary monitor if connected
screenNumber = max(Screen('Screens'));

%% Set display parameters (THESE MAY NEED TO CHANGE DEPENDING ON SETUP)
display.dist = 50;  % cm
display.width = 30; % cm
display.skipChecks = 1; % avoid Screen's timing checks and verbosity

% Generate a brief display to acquire parameters
display = OpenWindow(display);
screenResX = display.resolution(1);
screenResY = display.resolution(2);
Screen('CloseAll');

%% Set up stimulus and task parameters
% colours
red    = [255 0   0  ];
blue   = [100 175 255]; % light blue
green  = [0   255 0  ];
magenta= [255 0   255];
cyan   = [0   255 255];
yellow = [255 255 0  ];
orange = [255 128 0  ];
white =  [255 255 255];
neutralColour = yellow;
flipColour = blue;
selectColour = orange;

% Colour codings
zone1and2 = blue;
zone3and4 = red;
zone5and6 = magenta; % m = magenta
zone7and8 = green; % g = green
allzones = [blue; red; magenta; green];

% Establish Task Parameters
seqlength = 5:8; % possible length of sequences
total_number_of_blocks = 4;
total_trials_block = 28;
seqlength_by_blockType = [repmat(5,1,7) repmat(6,1,7) repmat(7,1,7) repmat(8,1,7);... % 7 repeats of seqlengths between 5-8
    repmat(5,1,7) repmat(6,1,7) repmat(7,1,7) repmat(8,1,7);...
    repmat(5,1,7) repmat(6,1,7) repmat(7,1,7) repmat(8,1,7);...
    repmat(5,1,7) repmat(6,1,7) repmat(7,1,7) repmat(8,1,7)];

% Timings
displayArrayDuration = 2; % display initial array for x s
displaySequenceDuration = .750; % display each element in the sequence array for x s
delayAfterSequence = .750; % delay (in s) to wait after displaying sequence % should this be randomised?


% Possible Stimulus Locations
% Setup stimulus location (based on screen size)
% Oct 23, 2023 - will need to modify this in final version of the task
possibleLocations = [...
    1 1 .10*screenResX .10*screenResY;... % should change these to a random range between .10 and .40
    1 2 .20*screenResX .10*screenResY;...
    1 3 .30*screenResX .10*screenResY;...
    1 4 .40*screenResX .10*screenResY;...
    2 1 .10*screenResX .30*screenResY;...
    2 2 .20*screenResX .30*screenResY;...
    2 3 .30*screenResX .30*screenResY;...
    2 4 .40*screenResX .30*screenResY;...
    3 1 .10*screenResX .50*screenResY;...
    3 2 .20*screenResX .50*screenResY;...
    3 3 .30*screenResX .50*screenResY;...
    3 4 .40*screenResX .50*screenResY;...
    4 1 .10*screenResX .70*screenResY;...
    4 2 .20*screenResX .70*screenResY;...
    4 3 .30*screenResX .70*screenResY;...
    4 4 .40*screenResX .70*screenResY;...
    5 1 .60*screenResX .10*screenResY;...
    5 2 .70*screenResX .10*screenResY;...
    5 3 .80*screenResX .10*screenResY;...
    5 4 .90*screenResX .10*screenResY;...
    6 1 .60*screenResX .30*screenResY;...
    6 2 .70*screenResX .30*screenResY;...
    6 3 .80*screenResX .30*screenResY;...
    6 4 .90*screenResX .30*screenResY;...
    7 1 .60*screenResX .50*screenResY;...
    7 2 .70*screenResX .50*screenResY;...
    7 3 .80*screenResX .50*screenResY;...
    7 4 .90*screenResX .50*screenResY;...
    8 1 .60*screenResX .70*screenResY;...
    8 2 .70*screenResX .70*screenResY;...
    8 3 .80*screenResX .70*screenResY;...
    8 4 .90*screenResX .70*screenResY];

%% Trial Parameters
breaktime = 3; % countdown variable at the start of the block

%% Prompt screen to enter ppt info to be written in logfile name
question = {'Participant Number:'; 'Block Number'};
title = 'Experiment Setup';
NumOfLines = [1 75; 1 75];
prompt= inputdlg(question,title,NumOfLines);
participantNumber  = str2double(prompt(1));
startingBlock = str2double(prompt(2));

if isempty(startingBlock)||isnan(startingBlock)
    startingBlock = 1;
end

% Do dummy calls to GetSecs, WaitSecs, KbCheck
% KbCheck;
% WaitSecs(0.1);
% GetSecs;

% Initialise a response matrix
respMat = struct('participantNumber', []);

%% Initialise Task
KbReleaseWait; % make sure person is not pushing down any buttons
mouseID = GetMouseIndices; % get mouse address
mouseID = mouseID(1);
escapeKey = KbName('Q');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXPERIMENTAL (TRIAL) LOOP (possibly)
% Single Trial:
% 1) Display Array for X ms
% 2) Display Sequence by flipping colour (each block for X ms)
% 3) Pause for X ms
% 4) Collect click location and response time until a mistake is made
% 5) If mistake is made - display feedback
% 6) If no mistake is made - display feedback "CORRECT!"
% 7) If no mistake is made, increase sequence length by 1 (to maximum of 8)
% 8) If mistake is made, decrease sequence length by 1 (to minimum of 3)
%%% Need to confirm staircase method

display = OpenWindow(display);
center  = display.resolution/2;
boxModifier = display.resolution(1)/50;
totalTrialExperiment = 1; % initialise total experiment trial counter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for currBlock = startingBlock:total_number_of_blocks % allows for manually starting at any block
    if currBlock == 1 % display instructions for block 1
        instructionimg = imread('kcl_corsi_promptScreen_block1.jpg');
        texI = Screen('MakeTexture', display.windowPtr, instructionimg);
        % rect = [50 -250 1600 1250];
        Screen('DrawTexture', display.windowPtr, texI) %, rect);
        Screen('Flip',display.windowPtr);
        KbWait();

    elseif currBlock == 2 % display instructions for block 2
        instructionimg = imread('kcl_corsi_promptScreen_block2.jpg');
        texI = Screen('MakeTexture', display.windowPtr, instructionimg);
        % rect = [50 -250 1600 1250];
        Screen('DrawTexture', display.windowPtr, texI) %, rect);
        Screen('Flip',display.windowPtr);
        KbWait();

    elseif currBlock == 3 % display instructions for block 3
        instructionimg = imread('kcl_corsi_promptScreen_block3.jpg');
        texI = Screen('MakeTexture', display.windowPtr, instructionimg);
        % rect = [50 -250 1600 1250];
        Screen('DrawTexture', display.windowPtr, texI) %, rect);
        Screen('Flip',display.windowPtr);
        KbWait();

    else % display instructions for block 4
        instructionimg = imread('kcl_corsi_promptScreen_block4_part1.jpg');
        texI = Screen('MakeTexture', display.windowPtr, instructionimg);
        % rect = [50 -250 1600 1250];
        Screen('DrawTexture', display.windowPtr, texI) %, rect);
        Screen('Flip',display.windowPtr);
        KbStrokeWait;

        instructionimg = imread('kcl_corsi_promptScreen_block4_part2.jpg'); 
        texI = Screen('MakeTexture', display.windowPtr, instructionimg);
        % rect = [50 -250 1600 1250];
        Screen('DrawTexture', display.windowPtr, texI) %, rect);
        Screen('Flip',display.windowPtr);
        KbWait();
    end

    %% generate ITIs
    %ITIs = [ones(total_trials_block/4,1)*.75;ones(total_trials_block/4,1)*1;ones(total_trials_block/4,1)*1.25];% create vector for ITIs
    %ITI = Shuffle(ITIs);% shuffle ITIs

    % Generate list of seq_lengths
    blockseq_lengths = Shuffle(seqlength_by_blockType(currBlock,:));

    %% Starting Trial Block
    for t = 1:total_trials_block
        HideCursor;

        %% 00) Establish Trial Parameters and eequence
        tr_seqlength = blockseq_lengths(t); % TEMPORARY - for testing purposes sequence_length for trial t
        tr_corsiPattern = nan(tr_seqlength,2);

        % Generate trial sequence
        clear tr_zoneOrder rect cp
        if tr_seqlength <= 8
            tr_zoneOrder = randperm(8, tr_seqlength);
        elseif tr_seqlength > 8
            tr_zoneOrder = [randperm(8, 8) randperm(8, tr_seqlength - 8)];
        end
        rect = nan(4, tr_seqlength);
        for cp = 1:tr_seqlength % set through each stim in the corsi pattern (cp)
            temp = randperm(4,1); % randomiser for stimulus within zone location (1-4)

            tr_corsiPattern(cp,[1,2]) = possibleLocations(possibleLocations(:,1)==tr_zoneOrder(cp) & ... % select zone
                possibleLocations(:,2)==temp,[3,4]); % select random configuration of stimuli ~~ stim location within zone

            % Prepare matrix containing block size and location
            rect(1, cp) =  tr_corsiPattern(cp,1) - boxModifier;
            rect(2, cp) =  tr_corsiPattern(cp,2) - boxModifier;
            rect(3, cp) =  tr_corsiPattern(cp,1) + boxModifier;
            rect(4, cp) =  tr_corsiPattern(cp,2) + boxModifier;
        end


        %% 1) Display Array
        Screen('FillRect', display.windowPtr, neutralColour, rect) % prepare complete array display
        Screen('Flip',display.windowPtr, [], 1); % display array
        pause(displayArrayDuration) % display array for DISPLAYARRAYDURATION seconds

        %% 2) Display Sequence
        for cp = 1:tr_seqlength
            % Flip to NEW colour (flipColour)
            Screen('FillRect', display.windowPtr,flipColour, rect(:,cp))
            Screen('Flip',display.windowPtr, [], 1);
            pause(displaySequenceDuration)

            % Flip BACK to old colour (neutralColour)
            Screen('FillRect', display.windowPtr,neutralColour, rect(:,cp))
            Screen('Flip',display.windowPtr, [], 1);
        end

        %% 3) Delay after presenting sequence
        pause(delayAfterSequence)

        %% 4) Collect response
        % display mouse cursor
        ShowCursor('Hand', display.windowPtr, mouseID)
        SetMouse(center(1),center(2));

        % Generate correctSequence_Locations matrix
        % (correctSequence_Locations = x_start, x_end, y_start, y_end)
        temp_correctSequence_Locations = nan(tr_seqlength, 4);
        for scp = 1:tr_seqlength
            temp_correctSequence_Locations(scp, 1) = tr_corsiPattern(scp, 1) - boxModifier;
            temp_correctSequence_Locations(scp, 2) = tr_corsiPattern(scp, 1) + boxModifier;
            temp_correctSequence_Locations(scp, 3) = tr_corsiPattern(scp, 2) - boxModifier;
            temp_correctSequence_Locations(scp, 4) = tr_corsiPattern(scp, 2) + boxModifier;
        end

        errorIndex = 0; % initialise error index

        % Initialise Mouse Trackers
        x=nan(tr_seqlength,1); y=nan(tr_seqlength,1); clicks=nan(tr_seqlength, 1); clickSecs=nan(tr_seqlength, 1);

        % Allow user input and check for errors
        startTime = GetSecs;
        for userClicks = 1:tr_seqlength % maximum number of clicks
            [clicks(userClicks), x(userClicks), y(userClicks), clickSecs(userClicks)] = GetClicks(display.windowPtr, 0);

            % Figure out which block was selected and change its colour
            for j_cp = 1:tr_seqlength
                if x(userClicks) > rect(1, j_cp) && x(userClicks) < rect(3, j_cp) && ...
                        y(userClicks) > rect(2, j_cp) && y(userClicks) < rect(4, j_cp)
                    % This is where the CLUE is inserted


                    if currBlock > 1 && j_cp < tr_seqlength
                        switch tr_zoneOrder(j_cp+1) % depends on the NEXT stimulus in the sequence
                            case 1
                                selectColour = zone1and2;
                                % disp("zone1and2")
                            case 2
                                selectColour = zone1and2;
                                % disp("zone1and2")
                            case 3
                                selectColour = zone3and4;
                                % disp("zone3and4")
                            case 4
                                selectColour = zone3and4;
                                % disp("zone3and4")
                            case 5
                                selectColour = zone5and6;
                                % disp("zone5and6")
                            case 6
                                selectColour = zone5and6;
                                % disp("zone5and6")
                            case 7
                                selectColour = zone7and8;
                                % disp("zone7and8")
                            case 8
                                selectColour = zone7and8;
                                % disp("zone7and8")
                        end
                    else % block = 1
                        selectColour = allzones(randperm(4,1), :); % randomly select one of the four colours
                    end
                    Screen('FillRect', display.windowPtr, selectColour, rect(:, j_cp))
                    Screen('Flip',display.windowPtr, [], 1);
                    break
                end % which block
            end % scroll through each block (figure out which one was clicked)

            % Check if location is correct
            if x(userClicks) < temp_correctSequence_Locations(userClicks, 1) || x(userClicks) > temp_correctSequence_Locations(userClicks, 2) || ...
                    y(userClicks) < temp_correctSequence_Locations(userClicks, 3) || y(userClicks) > temp_correctSequence_Locations(userClicks, 4)
                errorIndex = userClicks; % the click number they got wrong
                break % leave loop and go to error feedback
            end
        end % user clicks
        clickTimes = clickSecs - startTime;
        Screen('Flip',display.windowPtr);
        if errorIndex > 0
            Screen(display.windowPtr,'DrawText','Incorrect :(' ,center(1),center(2),[255,255,255]);
        else
            Screen(display.windowPtr,'DrawText','Correct! Well done!' ,center(1),center(2),[255,255,255]);
        end
        Screen('Flip',display.windowPtr);
        pause(1);

        %% Write trial data to response matrix
        respMat(totalTrialExperiment).participantNumber = participantNumber; % participant number
        respMat(totalTrialExperiment).totalTrialExperiment = totalTrialExperiment; % trial number within experiment
        respMat(totalTrialExperiment).currBlock = currBlock; % current block
        respMat(totalTrialExperiment).trialWithinBlock = t; % trial within block
        respMat(totalTrialExperiment).trialSequenceLength = tr_seqlength; % sequence length 
        respMat(totalTrialExperiment).tr_zoneOrder = tr_zoneOrder;
        respMat(totalTrialExperiment).tr_corsiPattern = tr_corsiPattern;
        respMat(totalTrialExperiment).correctSequence_Locations = temp_correctSequence_Locations;
        respMat(totalTrialExperiment).clickTimes = clickTimes;
        respMat(totalTrialExperiment).errorIndex = errorIndex;
        respMat(totalTrialExperiment).task_version = task_version;
        % what else?
        totalTrialExperiment = totalTrialExperiment + 1;

    end % individual block loop
    saveAfterCrash_Corsi; % save data after each block ... just in case
    %% DISPLAY END OF BLOCK SCREEN
    if currBlock == 1 % display instructions for block 1
        instructionimg = imread('kcl_corsi_promptScreen_block1end.jpg');
        texI = Screen('MakeTexture', display.windowPtr, instructionimg);
        % rect = [50 -250 1600 1250];
        Screen('DrawTexture', display.windowPtr, texI) %, rect);
        Screen('Flip',display.windowPtr);
        KbStrokeWait;
    elseif currBlock == 2 % display instructions for block 2
        instructionimg = imread('kcl_corsi_promptScreen_block2end.jpg');
        texI = Screen('MakeTexture', display.windowPtr, instructionimg);
        % rect = [50 -250 1600 1250];
        Screen('DrawTexture', display.windowPtr, texI) %, rect);
        Screen('Flip',display.windowPtr);
        KbStrokeWait;
    elseif currBlock == 3 % display instructions for block 3
        instructionimg = imread('kcl_corsi_promptScreen_block3end.jpg');
        texI = Screen('MakeTexture', display.windowPtr, instructionimg);
        % rect = [50 -250 1600 1250];
        Screen('DrawTexture', display.windowPtr, texI) %, rect);
        Screen('Flip',display.windowPtr);
        KbStrokeWait;
    % else % display instructions for block 4
    %     instructionimg = imread('kcl_corsi_promptScreen_block4_part.jpg');
    %     texI = Screen('MakeTexture', display.windowPtr, instructionimg);
    %     % rect = [50 -250 1600 1250];
    %     Screen('DrawTexture', display.windowPtr, texI) %, rect);
    %     Screen('Flip',display.windowPtr);
    %     KbWait();
    end


end % all block loop end

%% save data
save(['kcl_corsi_ppt_' , num2str(participantNumber), '_', datestr(now,'mmmm-dd-yyyy_HH-MM-SS AM'), '.mat'],'respMat') %#ok<*TNOW1,*DATST,*DLMWT>

debriefimg = imread('kcl_corsi_promptScreen_debrief.jpg');
texI = Screen('MakeTexture', display.windowPtr, debriefimg);
rect = [50 -250 1600 1250];
Screen('DrawTexture', display.windowPtr, texI) %, rect);
Screen('Flip',display.windowPtr);
KbWait();

Screen('CloseAll');
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IF THE KEYBOARD DOESN'T WORK AFTER PTB CRASHES TRY RUNNING THE FOLLOWING LINE OF CODE
% BY HIGHLIGHTING THIS TEXT AND RIGHT-CLICKING -> EVALUATE SELECTION IN COMMAND WINDOW
%ListenChar(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IF THE PROGRAM CRASHES BETWEEN BLOCKS - DO THE FOLLOWING TO SAVE THE EXISTING DATA
% BY HIGHLIGHTING THIS TEXT AND RIGHT-CLICKING -> EVALUATE SELECTION IN COMMAND WINDOW
%saveAfterCrash_Corsi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

