function out = My3DConv( I, gx, gy, gz )
% Computes a separable 3D convolution
gx  = gx(:);
gx  = permute(gx,[2,1,3]);
gy  = gy(:);
gz  = gz(:);
gz  = permute(gz,[3,2,1]);
I   = convn( I, gx, 'same' );
I   = convn( I, gy, 'same' );
out = convn( I, gz, 'same' );
return;