function display = OpenWindow(display)

if ~exist('display','var')
    display.screenNum = 0;
end

if ~isfield(display,'screenNum')
    display.screenNum = 0;
end

if ~isfield(display,'bkColor')
    display.bkColor = [0,0,0]; %black
end

if ~isfield(display,'skipChecks')
    display.skipChecks = 0;
end

if display.skipChecks
    Screen('Preference', 'Verbosity', 0);
    Screen('Preference', 'SkipSyncTests',1);
    Screen('Preference', 'VisualDebugLevel',0);
end

%Open the window
[display.windowPtr,res]=Screen('OpenWindow',display.screenNum,display.bkColor);

%Set the display parameters 'frameRate' and 'resolution'
display.frameRate = 1/Screen('GetFlipInterval',display.windowPtr); %Hz

if ~isfield(display,'resolution')
    display.resolution = res([3,4]);
end