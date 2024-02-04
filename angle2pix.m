function pix = angle2pix(display,ang)
%Calculate pixel size
pixSize = display.width/display.resolution(1);
sz = 2*display.dist*tan(pi*ang/(2*180));
pix = round(sz/pixSize);
return

angle2pix(display,ang)