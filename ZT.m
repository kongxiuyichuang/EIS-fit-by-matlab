function Z= ZT(w,options)
%Return the impedance of ZT(diffusion impendance),
% complex size(1:length(w))
%ZT Generalized Finite Space Warburg element (GFSW).
% also published as Finite Length Warburg with Open Circuit terminated
%Inputs
%---------------
%w = Angular frequency [1/s]
%options is a structure with fields:
%Rw=RT/nFc_0*D [ohm cm^2], polarization resistance;
%tau=L^2/D [s], time constant;
% L[cm] is the thickness of the diffusion layer and 
% D [cm2*s–1] is the diffusion coefficient.
%n_w=0~0.5[-]
% Ref:--------------
% https://doi.org/10.1021/acsmeasuresciau.2c00070
% Illig J. Physically based impedance modelling of lithium-ion cells[M]. KIT Scientific Publishing, 2014.

arguments
    w (1,:) double
    options.Rw (1,1) double=[]; 
    options.tau (1,1) double=[];
    options.n_w (1,1) double=[0.5];

end
Rw=options.Rw;
tau=options.tau;
n_w=options.n_w;

Z = Rw.*coth((1j*w*tau).^n_w)./((1j*w*tau).^n_w);
end