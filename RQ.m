function z = RQ(w, options)
% Return the impedance of an RQ circuit.
% complex size(1:length(w))
% Inputs
% ----------------
% w = Angular frequency [1/s]
% options is a structure with fields:
% R = Resistance [Ohm]
% Q = Constant phase element [s^n/ohm]
% n = Constant phase elelment exponent [-]
% fs = Summit frequency of RQ circuit [Hz]

arguments
    w (1,:) double
    options.R  double = []
    options.Q  double = []
    options.n  double = []
    options.fs  double = []
end

R=options.R;
Q=options.Q;
n=options.n;
fs=options.fs;

if isempty(R)
    R = (1/(Q.*(2*pi*fs).^n));
elseif isempty(Q)
    Q = (1/(R.*(2*pi*fs).^n));
elseif isempty(n)
    n = log(Q.*R)./log(1./(2*pi*fs));
end
% z = (R./(1+R.*Q.*(w*1i).^n));
z=1./(1./R+Q.*(w*1j).^n);
end


% function z = RQ(w, varargin)
%     Return the impedance of an Rs-RQ circuit.
%     Inputs
%     w = Angular frequency [1/s]
%     varargin is a series of name-value pairs:
%     'R' = Resistance [Ohm]
%     'Q' = Constant phase element [s^n/ohm]
%     'n' = Constant phase elelment exponent [-]
%     'fs' = Summit frequency of RQ circuit [Hz]
%
%     Create an inputParser instance
%     p = inputParser;
%
%     Define default values for the options
%     defaultR = [];
%     defaultQ = [];
%     defaultn = [];
%     defaultfs = [];
%
%     Add optional name-value pair arguments
%     addParameter(p, 'R', defaultR, @isnumeric);
%     addParameter(p, 'Q', defaultQ, @isnumeric);
%     addParameter(p, 'n', defaultn, @isnumeric);
%     addParameter(p, 'fs', defaultfs, @isnumeric);
%
%     Parse the input options
%     parse(p, varargin{:});
%
%     Extract the values from the inputParser
%     R = p.Results.R;
%     Q = p.Results.Q;
%     n = p.Results.n;
%     fs = p.Results.fs;
%
%     if isempty(R)
%         R = (1./(Q.*(2*pi*fs).^n));
%     elseif isempty(Q)
%         Q = (1./(R.*(2*pi*fs).^n));
%     elseif isempty(n)
%         n = log(Q.*R)./log(1./(2*pi.*fs));
%     end
%     z = (R./(1+R.*Q.*(w*1i).^n));
% end
