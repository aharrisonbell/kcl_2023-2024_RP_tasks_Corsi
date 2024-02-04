function   display = drawFixation(display)
global feedback

%Deal with default values
if ~isfield(display,'fixation')
    display.fixation = [];
end

%Size
if ~isfield(display.fixation,'size')
    display.fixation.size = .5; %degrees
end

%Mask
if ~isfield(display.fixation,'mask')
    display.fixation.mask = 2;  %degrees
end

%Color
if ~isfield(display.fixation,'color')
    display.fixation.color = {[255,255,255],[0,0,0]};
end

%Flip
if ~isfield(display.fixation,'flip')
    display.fixation.flip = 1;  %flip by default
end

global center
center = display.resolution/2;

%Calculate size of boxes in screen-coordinates
sz(1) = angle2pix(display,display.fixation.size/2);
sz(2) = angle2pix(display,display.fixation.size/4);
sz(3) = angle2pix(display,display.fixation.mask/2);

%Calculate the rectangles in screen-coordinates [l,t,r,b]
for i=1:3
    rect{i}= [-sz(i)+center(1),-sz(i)+center(2),sz(i)+center(1),sz(i)+center(2)];
end

%Mask (background color)
Screen('FillOval', display.windowPtr, display.bkColor,rect{3});
%Outer rectangle (default is white)
Screen('FillRect', display.windowPtr, display.fixation.color{1},rect{1});
%Inner rectangle (default is black)
Screen('FillRect', display.windowPtr, display.fixation.color{2},rect{2});

if feedback == 0
    feedbackimg = imread('Fix.png');
elseif feedback == 1
    feedbackimg = imread('Smiley.png');
elseif feedback == 2
    feedbackimg = imread('Frowney.png');
end
texA = Screen('MakeTexture', display.windowPtr, feedbackimg); 
Screen('DrawTexture', display.windowPtr, texA, [], [], 0);

if display.fixation.flip
     Screen('Flip',display.windowPtr);
end