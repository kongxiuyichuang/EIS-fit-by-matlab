function Z= Zw(w,options)
%Return the impedance of Zw(warburg),
% complex size(1:length(w))
%semi-infinite linear diffusion
%Inputs
%---------------
%w = Angular frequency [1/s]
%options is a structure with fields:
%sigma= Constant phase element [ohm/s^(-0.5)];
% sigma=2RT/(n^2F^2*sqrt(2D)*c),
%  C is the concentration of the redox species assuming that this concentration...
%  is the same as the one in the bulk solution and that Cox = CRed = C.
% Ref:----------------
%  https://doi.org/10.1021/acsmeasuresciau.2c00070

arguments
    w (1,:) double
    options.sigma (1,1) double=[];   
end
sigma=options.sigma;
Z = sigma*(w.^(-0.5))-1j*sigma*(w.^(-0.5));
end