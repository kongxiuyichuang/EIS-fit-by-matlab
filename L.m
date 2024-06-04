function Z= L(w,options)
%Return the impedance of an ideal inductance  
%complex size(1:length(w))
%Inputs
%---------------
%w = Angular frequency [1/s]
%options is a structure with fields:
%L = inductance [H]  
arguments
    w (1,:) double
    options.L (1,1) double=[];
end
L=options.L;
Z = 1j*w*L;
end