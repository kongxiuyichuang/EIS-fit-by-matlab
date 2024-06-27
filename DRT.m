function output= DRT(Z,options)
% input:
%Z:the experiment impedance data,size(:,3);
%lambda: the regularization parameter,default is 0.001;
%FWHM: is the full width at half maximum of gaussian function,defualt is 0.5;
%----------------------------
%output:
%output is structure with fields:
%f: frequency [Hz],size(:,1);
%Z_exp:
%tau:the calculated tau range [s],size(:,1)
%Z_DRT:the Reconstructed impedance spectrum by DRT,a complex matrix,size(:,1);
%gtau:the gtau is distribution function of resistances at the continuous
%relaxation times,[tau_fine,gtau],size(10*length(htau),2)
%htau:the R value of RC,size is double of frequency,[tau,htau],size(2*length(f),2];
%----------------------------
% the DRT caculation method is based on presented in:
%1. Danzer, Michael A. "Generalized distribution of relaxation times analysis
% for the characterization of impedance spectra." Batteries 5.3 (2019): 53.
%2.Hahn, Markus, et al. "Optimized process parameters for a reproducible 
% distribution of relaxation times analysis of electrochemical systems."
% Batteries 5.2 (2019): 43.
%3.Wildfeuer, Leo, Philipp Gieler, and Alexander Karger. "Combining the 
% distribution of relaxation times from EIS and time-domain data for 
% parameterizing equivalent circuit models of lithium-ion batteries." 
% Batteries 7.3 (2021): 52.

arguments
    Z (:,3) double
    options.lambda (1,1) =0.001;
    options.FWHM (1,1) =0.5;
end

Z_exp=Z;
if Z_exp(1,1)<Z_exp(end,1)
    Z_exp=flip(Z_exp,1);
else
end

f=Z_exp(:,1);
z=Z_exp(:,2)+1j*Z_exp(:,3);

tauL =floor(log10(1/2/pi/max(f)))-0.5;
tauH = ceil(log10(1/2/pi/min(f)))+0.5;
% tauH=ceil(max(log10(1./f)))+1;
% tauL=floor(min(log10(1./f)))-1;
M=2*length(f);
tau = logspace(tauL,tauH,M);
tau = tau(:)';       % Force row vector.

z = z(:); f = f(:);  % Force column vector.
magZ = abs(z);
omega = 2*pi*f;
s = 1j*omega;
M = length(tau);
coeffRn = 1./(1+s*tau);

% Build measurement vector (y).
%  y = [real(z)./magZ; imag(z)./magZ];
y = [real(z); imag(z)];
% Build measurement matrix (H). Columns correspond to elements
Hrealz = zeros(length(z),3+M);   Himagz = zeros(length(z),3+M);
Hrealz(:,1) = 1;                 Himagz(:,1) = 0;           % R0
Hrealz(:,2) = 0;                 Himagz(:,2) = omega;       % L0
Hrealz(:,3) = 0;                 Himagz(:,3) = -1./omega;   % 1/C0
Hrealz(:,4:end) = real(coeffRn); Himagz(:,4:end) = imag(coeffRn); % R1...RM
%  H = [Hrealz./magZ; Himagz./magZ];
H = [Hrealz; Himagz];

%% Find parameter vector by least-squares regression.
lambda=options.lambda;
L=zeros(length(tau),1);
L=L+diag(ones(1,length(tau)))+diag(0.5*ones(1,length(tau)-1),1)+diag(0.5*ones(1,length(tau)-1),-1);
L=[zeros(length(tau),3),L];
B=[H;lambda*L];

U=[y(:);zeros(length(tau),1)];
x=lsqnonneg(B,U);

%% output matrix
elem.R0 = x(1);
elem.L0 = x(2);
elem.C0inv = x(3);
elem.Rn = x(4:end);
Zmodel = elem.R0 + s*elem.L0 + elem.C0inv./s + coeffRn*elem.Rn;
Rn = elem.Rn;
output.f=f;
output.Z_exp=Z;
output.tau=tau';
output.Z_DRT=Zmodel;

%% plot
figure("Units","centimeters","Position",[20,10,16*2,12])
subplot(1,2,1)
plot(real(z),-imag(z),'o',DisplayName='exp',LineWidth=2)
hold on;
plot(real(Zmodel),-imag(Zmodel),'-',LineWidth=2,DisplayName='DRT');
grid on;legend("Location",'northwest');xlabel('Real(Z)(\Omega)');ylabel('-Imag(Z)(\Omega)')

% calculate RBF of weight
[peaks,x_index]=findpeaks(Rn);
tau_peak=log10(tau(x_index));

FWHM=options.FWHM;
sigma=FWHM/(2*sqrt(2*log(2)));
% gaussian=@(x,xc) exp(-(x-xc).^2/(2*sigma^2));
% x_sample=log10(tau);y_sample=Rn
num=length(peaks);
Phi=zeros(num,num);
for i=1:num
    Phi(i,:)= exp(-(tau_peak-tau_peak(i)).^2./(2*sigma^2));
end
% w=Phi\peaks;
weight=lsqnonneg(Phi,peaks);
y_fit=Phi*weight;

% calculate RBF
tau_fine=logspace(tauL,tauH,M*10);
y_fine=zeros(length(tau_fine),1);
for i=1:length(tau_fine)
    for k=1:length(weight)
        y_fine(i)=y_fine(i)+weight(k)*exp(-(log10(tau_fine(i))-tau_peak(k)).^2./(2*sigma^2));
    end
end
subplot(1,2,2)
semilogx(tau_fine,y_fine,LineWidth=2),xlabel('\tau(s)');ylabel('g(h(k))(\Omega)');grid on;

% [a,b]=findpeaks(y_fine);
% sum(Rn)*a(2)/sum(a)

%plot kernels
y_kernels=zeros(length(tau_peak),length(tau_fine));
for i=1:length(tau_peak)
    for k=1:length(tau_fine)
        y_kernels(i,k)=weight(i)*exp(-(log10(tau_fine(k))-tau_peak(i)).^2./(2*sigma^2));

    end
end
hold on
semilogx(tau_fine,y_kernels),xlabel('\tau(s)');ylabel('h(k)(\Omega)');grid on;

output.gtau=[tau_fine',y_fine];
output.htau=[tau',Rn];

end