function movingDots(display,dots,duration)
 
leftKey   = KbName('Z');
rightKey  = KbName('M');
escapeKey = KbName('Q');

rt = [];
global RT 
global runs
RT = [];
global response
global tarresp
global feedback
%Calculate total number of dots across fields
nDots = sum([dots.nDots]);

%Zero out the color and size vectors
colors = zeros(3,nDots);
sizes  = zeros(1,nDots);

%Generate a random order to draw the dots so that one field won't occlude
%another field.
order=  randperm(nDots);

%Intitialize the dot positions and define some other initial parameters
tStart = GetSecs;
count = 1;
for i=1:length(dots) %Loop through the fields

    %Calculate the left, right top and bottom of each aperture (in degrees)
    l(i) = dots(i).center(1)-dots(i).apertureSize(1)/2;
    r(i) = dots(i).center(1)+dots(i).apertureSize(1)/2;
    b(i) = dots(i).center(2)-dots(i).apertureSize(2)/2;
    t(i) = dots(i).center(2)+dots(i).apertureSize(2)/2;

    %Generate random starting positions
    dots(i).x = (rand(1,dots(i).nDots)-.5)*dots(i).apertureSize(1) + dots(i).center(1);
    dots(i).y = (rand(1,dots(i).nDots)-.5)*dots(i).apertureSize(2) + dots(i).center(2);

    %Create a direction vector for a given coherence level
    direction = rand(1,dots(i).nDots)*360;
    nCoherent = ceil(dots(i).coherence*dots(i).nDots);  %Start w/ all random directions
    direction(1:nCoherent) = dots(i).direction;  %Set the 'coherent' directions
    direction = Shuffle(direction);
    
    %Calculate dx and dy vectors in real-world coordinates
    dots(i).dx = dots(i).speed*sin(direction*pi/180)/display.frameRate;
    dots(i).dy = -dots(i).speed*cos(direction*pi/180)/display.frameRate;
    dots(i).life =    ceil(rand(1,dots(i).nDots)*dots(i).lifetime);
    
    %Fill in the 'colors' and 'sizes' vectors for this field
    id = count:(count+dots(i).nDots-1);  %index into the nDots length vector for this field
    nColCoherent   = int64(dots.nDots*dots.col.coherence);
    nColIncoherent = int64(dots.nDots-(dots.nDots*dots.col.coherence));
    Distr     = [randperm(255);randperm(255);randperm(255)];
    SelDistr  = Distr(:,1:nColIncoherent);
    Ndistr    = int64(nColIncoherent/2);
    dots(i).distr1.color = [255 0 0];
    dots(i).distr2.color = [0 255 0];
    
    t1 = repmat(dots(i).color(:),1,nColCoherent);
    t2 = repmat(dots(i).distr1.color(:),1,nColIncoherent/2);
    t3 = repmat(dots(i).distr2.color(:),1,nColIncoherent/2);
    colors = [t1,t2,t3];
    
    sizes(order(id)) = repmat(dots(i).size,1,dots(i).nDots);
    count = count+dots(i).nDots;
end

%Zero out the screen position vectors and the 'goodDots' vector
pixpos.x = zeros(1,nDots);
pixpos.y = zeros(1,nDots);
goodDots = false(zeros(1,nDots));

%Calculate total number of temporal frames
nFrames = secs2frames(display,duration);

%Loop through the frames

for frameNum=1:nFrames
    count = 1;
    for i=1:length(dots)  %Loop through the fields

        %Update the dot position's real-world coordinates
        dots(i).x = dots(i).x + dots(i).dx;
        dots(i).y = dots(i).y + dots(i).dy;

        %Move the dots that are outside the aperture back one aperture width.
        dots(i).x(dots(i).x<l(i)) = dots(i).x(dots(i).x<l(i)) + dots(i).apertureSize(1);
        dots(i).x(dots(i).x>r(i)) = dots(i).x(dots(i).x>r(i)) - dots(i).apertureSize(1);
        dots(i).y(dots(i).y<b(i)) = dots(i).y(dots(i).y<b(i)) + dots(i).apertureSize(2);
        dots(i).y(dots(i).y>t(i)) = dots(i).y(dots(i).y>t(i)) - dots(i).apertureSize(2);

        %Increment the 'life' of each dot
        dots(i).life = dots(i).life+1;

        %Find the 'dead' dots
        deadDots = mod(dots(i).life,dots(i).lifetime)==0;

        %Replace the positions of the dead dots to random locations
        dots(i).x(deadDots) = (rand(1,sum(deadDots))-.5)*dots(i).apertureSize(1) + dots(i).center(1);
        dots(i).y(deadDots) = (rand(1,sum(deadDots))-.5)*dots(i).apertureSize(2) + dots(i).center(2);

        %Calculate the index for this field's dots into the whole list of
        %dots.  Using the vector 'order' means that, for example, the first
        %field is represented not in the first n values, but rather is
        %distributed throughout the whole list.
        id = order(count:(count+dots(i).nDots-1));
        
        %Calculate the screen positions for this field from the real-world coordinates
        pixpos.x = angle2pix(display,dots(i).x)+ display.resolution(1)/2;
        pixpos.y = angle2pix(display,dots(i).y)+ display.resolution(2)/2;

        %Determine which of the dots in this field are outside this field's
        %elliptical apertured
        goodDots = (dots(i).x-dots(i).center(1)).^2/(dots(i).apertureSize(1)/2)^2 + ...
            (dots(i).y-dots(i).center(2)).^2/(dots(i).apertureSize(2)/2)^2 < 1;
        
        count = count+dots(i).nDots;
    end
    
    %Draw all fields at once
    Screen('DrawDots',display.windowPtr,[pixpos.x(goodDots);pixpos.y(goodDots)], sizes(goodDots), colors(:,goodDots),[0,0],1);
   
    %Draw the fixation point (and call Screen's Flip')
    drawFixation(display);
    
    [keyIsDown,secs, keyCode] = KbCheck;
    
    %Collect Resp    
        if     keyCode(leftKey)
               rt = [rt;GetSecs - tStart];
               response = 1;
               frameNum = nFrames
               if tarresp == 1 
                   feedback = 1;
               elseif tarresp == 2
                   feedback = 2;
               end
        elseif keyCode(rightKey)
               rt = [rt;GetSecs - tStart];
               response = 2;
               frameNum = nFrames
               if tarresp == 2 
                   feedback = 1;
               elseif tarresp == 1
                   feedback = 2;
               end               
        elseif keyCode(escapeKey)
               ShowCursor;
               sca;
               ListenChar(1);
               return
        end                 
end
if     length(rt) > 1
       rt = rt(1);
elseif length(rt) <= 1
       rt = NaN;
       response = 0;
end
RT = [rt];
feedback = 0;
%clear the screen and leave the fixation point
drawFixation(display);
