function hi = hypsometricIntegral(z)
%HYPSOMETRICINTEGRAL Equal-area raster HI = (mean(z)-min(z))/(max(z)-min(z)).
arguments
    z {mustBeNumeric}
end
z=double(z(isfinite(z)));
if isempty(z) || max(z)==min(z)
    hi=NaN;
else
    hi=(mean(z)-min(z))/(max(z)-min(z));
end
end
