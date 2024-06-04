function Z= C(w,options)
%Return the impedance of an ideal capacitor 
%Inputs
%---------------
%w = Angular frequency [1/s]
%options is a structure with fields:
%C = Capacitance [F]  
arguments
    w (1,:) double
    options.C (1,1) double=[];
end
C=options.C;
Z = 1./(C*(w*1i));
end