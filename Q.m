function Z= Q(w,options)
%Return the impedance of Q(constant phase element),complex size(1:length(w))
%Inputs
%---------------
%w = Angular frequency [1/s]
%options is a structure with fields:
%Q = Constant phase element [s^n/ohm]
%n = Constant phase elelment exponent [-]
arguments
    w (1,:) double
    options.Q (1,1) double=[];
    options.n (1,1) double=[];
end
Q=options.Q;
n=options.n;
Z = 1./(Q*(w*1j).^n);
end