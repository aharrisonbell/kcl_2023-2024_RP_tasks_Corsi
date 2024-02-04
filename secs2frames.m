function frames=secs2frames(display,secs)
%converts time in seconds to frames 
frames = round(secs*display.frameRate);