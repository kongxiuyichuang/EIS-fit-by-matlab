function result=KK_link(f,Z_exp,options)
%LINKK Test impedance spectra for causality, time-invariance, and linearity.
% Validate an impedance spectrum against the Kramers–Kronig relations by
% fitting an LTI model comprising a number of series-connected Voigt 
% (parallel resistor-capacitor) to the spectrum.
% Implements an algorithm that automatically selectes the number of
% elements (N) to avoid under- and over-fitting.
%
% In circuit form, the model is represented as follows:
% 
%                                   C1                  CM
%                               |---||----|         |---||----|
% --\/\/\-----)()()(-----||-----|         |-- ... --|         |--
%    R0         L0       C0     |--\/\/\--|         |--\/\/\--|
%                    (optional)     R1                  RN
%                                          (N times)
% 
% where the elements may take on either positive or negative values
% (negative values do not have physical meaning, but still satisfy the
% the KK relations).
% Adapted to MATLAB from https://github.com/whileman133/LMB_Toolkit.git
%   which is based on the method presented in:
%   Schönleber, M. et al. A Method for Improving the Robustness of
%   linear Kramers-Kronig Validity Tests. 
%   Electrochimica Acta 131, 20–27 (2014)
%   doi: 10.1016/j.electacta.2014.01.034
%   https://doi.org/10.1016/j.electacta.2014.01.034
%----------------------
% input:
% f：frequency [Hz]
% Z_exp：the experiment impedance data, a complex matrix
% options is structure with fields:
% maxN: the max number of RC elements
% C_switch:includes a series integrator (capacitor) in 
% the model for fitting the low-frequency response of a
% battery cell. The default is 'on','off' means no capacitor is included.
% c:sets the threshold for the under/over-fitting factor (mu) to c instead
% of the default 0.8, c must be between 0 and 1.
%----------------------
% output:
% result is structure with fields:
% f:frequency [Hz]
% Z_exp:the experiment impedance data
% Z_model:the model impedance
% elements: The values ​​of the model elements,elements is struct
% residuals: Complex residuals, normalized to the magnitude of the true
% measured impedance ||Z||
% mu:Under/over-fitting factor for the final fit model.
% N: The number of series RC element

arguments
    f double;
    Z_exp double;
    options.maxN (1,1) double =100;
    options.C_switch string='on';
    options.c (1,1) {mustBeNumeric} = 0.8;

end
z=Z_exp;
C_switch=options.C_switch;
c=options.c;
maxN=options.maxN;
tauL = 1/2/pi/max(f);
tauH = 1/2/pi/min(f);

for N = 1:maxN
    % Logarithmically space time-constants of the RC circuits.
    tau = logspace(log10(tauL),log10(tauH),N);
    result = fitKK(f,z,tau,C_switch);
    if result.mu <= c
        break;
    end
end

result.N = N;
plot_result(z,result.Z_model,f,result.residuals)
end



function result = fitKK(f,z,tau,C_switch)
tau = tau(:)';       % Force row vector.
z = z(:); f = f(:);  % Force column vector.
magZ = abs(z);
omega = 2*pi*f;
s = 1j*omega;
M = length(tau);
coeffRn = 1./(1+s*tau);

% The assumed parameter vector is: x = [R0, L0, 1/C0, R1, R2, ... RM]'.

% Build measurement vector (y).
y = [real(z)./magZ; imag(z)./magZ];

% Build measurement matrix (H). Columns correspond to elements
% in the circuit model.
switch C_switch
    case 'on'
        Hrealz = zeros(length(z),3+M);   Himagz = zeros(length(z),3+M);
        Hrealz(:,1) = 1;                 Himagz(:,1) = 0;           % R0
        Hrealz(:,2) = 0;                 Himagz(:,2) = omega;       % L0
        Hrealz(:,3) = 0;                 Himagz(:,3) = -1./omega;   % 1/C0
        Hrealz(:,4:end) = real(coeffRn); Himagz(:,4:end) = imag(coeffRn); % R1...RM
        H = [Hrealz./magZ; Himagz./magZ];
    case 'off'
        Hrealz = zeros(length(z),2+M);   Himagz = zeros(length(z),2+M);
        Hrealz(:,1) = 1;                 Himagz(:,1) = 0;           % R0
        Hrealz(:,2) = 0;                 Himagz(:,2) = omega;       % L0
        Hrealz(:,3:end) = real(coeffRn); Himagz(:,3:end) = imag(coeffRn); % R1...RM
        H = [Hrealz./magZ; Himagz./magZ];
    otherwise
        error('C_switch is not a valid string, legal strings are "on" or "off"')
end

% Find parameter vector by least-squares regression.
% H*x=y
x = H\y;

% Compute impedance predicted by the fit model.
switch C_switch
    case 'on'
    elem.R0 = x(1);
    elem.L0 = x(2);
    elem.C0inv = x(3);
    elem.Rn = x(4:end);
    Z_model = elem.R0 + s*elem.L0 + elem.C0inv./s + coeffRn*elem.Rn;
    case 'off'
    elem.R0 = x(1);
    elem.L0 = x(2);
    elem.Rn = x(3:end);
    Z_model = elem.R0 + s*elem.L0 + coeffRn*elem.Rn;
    otherwise

end

% Compute under/over-fitting factor and complex residuals.
Rn = elem.Rn;
mu = 1 - sum(abs(Rn(Rn<0)))/sum(abs(Rn(Rn>=0)));
if isnan(mu) || mu < 0
    mu = 1;
end
residuals = (Z_model-z)./abs(z);

% Store result in output structure.
result.f=f;
result.Z_exp=z;
result.Z_model = Z_model;
result.elements = elem;
result.residuals = residuals;
result.mu = mu;


end

function plot_result(Z_exp,Z_modle,f,residuals)
z=Z_exp;
Zmodel=Z_modle;
screenSize = get(0, 'ScreenSize');
pos = [screenSize(3)/2 - 400, screenSize(4)/2 - 300, 800, 600];
figure('Position', pos);

subplot(2,2,1)
plot(real(z),-imag(z),'-o');
hold on
plot(real(Zmodel),-imag(Zmodel),'-o')
legend({'exp-data','KK-fit'});grid on;
xlabel('Real(ohm)');ylabel('-Imag(ohm)')

subplot(2,2,2)
semilogx(f,real(residuals)*100,'-o',DisplayName='Real(%)');
hold on
semilogx(f,imag(residuals)*100,'-o',DisplayName='Imag(%)');
legend();xlabel('f(Hz)'),ylabel('residuals(%)');grid on;
title('(Z_{fit}-Z_{exp})/|Z|')

subplot(2,2,3)
semilogx(f,-imag(z),'-o',DisplayName='Real(exp)')
hold on
semilogx(f,-imag(Zmodel),'-o',DisplayName='Real(KK-fit)')
legend();xlabel('f(Hz)'),ylabel('-Imag(ohm)');grid on;

subplot(2,2,4)
semilogx(f,real(z),'-o',DisplayName='Real(exp)')
hold on
semilogx(f,real(Zmodel),'-o',DisplayName='Real(KK-fit)')
legend();xlabel('f(Hz)'),ylabel('Real(ohm)');grid on;

end
