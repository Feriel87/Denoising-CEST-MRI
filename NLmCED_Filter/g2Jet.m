function [fs, fsx, fsy, fsxx, fsxy, fsyy] = g2Jet( f, scale )
% Second order Gaussian jet

fs = gD( f, scale, 0, 0 );
fsx = gD( f, scale, 1, 0 );
fsy = gD( f, scale, 0, 1 );
fsxx = gD( f, scale, 2, 0 );
fsxy = gD( f, scale, 1, 1 );
fsyy = gD( f, scale, 0, 2 );