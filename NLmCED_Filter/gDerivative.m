function r=gDerivative( order, x, Gs, scale )
switch order
case 0
r = Gs;
case 1
r = -x/(scale^2) .* Gs;
case 2
r = (x.^2-scale^2)/(scale^4) .* Gs;
    otherwise
    error('only derivatives up to second order are supported');
end