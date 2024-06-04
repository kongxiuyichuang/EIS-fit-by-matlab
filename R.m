function z = R(w,options)
%Return the impedance of R
%complex size(1:length(w))
%Inputs
%----------
%w = Angular frequency [1/s]
%options is a structure with fields:
%R = resistance [Ohm]

arguments
    w (1,:) double
    options.R (1,1) double=[];
end
R=options.R;
z = R *ones(1, length(w));
end