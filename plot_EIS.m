function plot_EIS(w,Z,options)
%plot the EIS (Nyquist or Bode)
%Inputs
%---------------
% w:= Angular frequency [1/s]
% Z:the impedance data, a complex matrix size:(1*length(w))
% options:is a structure with fields:
% Nyquist:Nyquist Graphic Switch,the default is 'on'
% Bode:Bode Graphic Switch,the default is 'off',
%   'on': plot Real of impedance & -Imag vs frequency
%   'real':plot Real of impedance vs frequency
%   'module':plot module of impedance vs frequency
% 'xlim_switch':Limit axis range Switch,the default is 'off'
arguments
    w (1,:) double
    Z (1,:) double
    options.Nyquist (1,1) string='on'
    options.Bode (1,1) string='off'
    options.xlim_switch(1,1) string ='off';
end
f=w./(2*pi);
xmax=ceil(max(max(real(Z)),max(-imag(Z))));
ymin=floor(min(-imag(Z)));

Nyquist=options.Nyquist;
Bode=options.Bode;
xlim_switch=options.xlim_switch;
switch  Nyquist
    case 'on'
        figure
        plot(real(Z),-imag(Z),'-o',LineWidth=2);
        grid on;
        switch xlim_switch
            case 'on'
                xlim([0,xmax]);ylim([ymin,xmax]);
            otherwise
        end

        xlabel('Real(ohm)'),ylabel('-Imag(ohm)')
    otherwise

end
switch Bode
    case 'on'
        figure
        semilogx(f,real(Z),'-o',LineWidth=2);
        hold on
        semilogx(f,-imag(Z),'-o',LineWidth=2);
        grid on; hold off
        xlabel('f(Hz)'),ylabel('Real,-Imag (ohm)');
        legend({'Real','-Imag'})
    case 'real'
        figure
        semilogx(f,real(Z),'-o',LineWidth=2);
        grid on;
        xlabel('f(Hz)'),ylabel('Real(ohm)');
        legend({'Real'})

    case 'imag'
        figure
        semilogx(f,-imag(Z),'-o',LineWidth=2);
        grid on;
        xlabel('f(Hz)'),ylabel('-Imag (ohm)');
        legend({'-Imag'})

    case 'module'
        figure
        semilogx(f,abs(Z),'-o',LineWidth=2);
        grid on;
        xlabel('f(Hz)'),ylabel('|Z| (ohm)');
        legend({'module'})
    case 'off'

    otherwise
        fprintf('plot Bode is not surpport %s\n',Bode)

end

end