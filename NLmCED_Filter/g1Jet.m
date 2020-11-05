function [fs, fsx, fsy] = g1Jet( f, scale )
% First order Gaussian jet
fs = gD( f, scale, 0, 0 );
fsx = gD( f, scale, 1, 0 );
fsy = gD( f, scale, 0, 1 );