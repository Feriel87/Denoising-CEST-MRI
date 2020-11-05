function g=gD(f,scale,ox,oy)

k = ceil(3*scale);
x = -k:k;
Gs = exp(-x.^2/2*scale^2);
Gs = Gs/sum(Gs);
Gsx = gDerivative( ox, x, Gs, scale );
Gsy = gDerivative( oy, x, Gs, scale );

g = convSepBrd( f, Gsx, Gsy );