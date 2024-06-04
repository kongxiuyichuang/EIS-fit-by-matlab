function plot_fit_result(w,exp_data,sim_circuit,x,options)
%plot the EIS (Nyquist or Bode) of experament and fit
%Inputs
%---------------
% w:= Angular frequency [1/s]
% exp_data:the impedance data, a complex matrix size:(1*length(w))
% sim_circuit:object circuit function handle
% x:circuit parameter of slover,size:(1:N)
% options:is a structure with fields:
% Nyquist:Nyquist Graphic Switch,the default is 'on'
% Bode:Bode Graphic Switch,the default is 'off',
%   'on': plot Real of impedance & -Imag vs frequency
%   'real':plot Real of impedance vs frequency
%   'module':plot module of impedance vs frequency
% 'xlim_switch':Limit axis range Switch,the default is 'off'
% 'legend':the legend of data,type is cell,size:(1,2),the default is {'exp-data','sim-data'}
% 'exp_LineStyle': the LineStyle of experament data,the default is '-o'
% 'exp_LineWidth': the LineWidth of experament data,the default is 2
% 'sim_LineStyle': the LineStyle of fit data,the default is '-o'
% 'sim_LineWidth': the LineWidth of fit data,the default is 1
arguments
    w (1,:) double;
    exp_data double;
    sim_circuit;
    x(1,:) double;
    options.Nyquist (1,1) string='on'
    options.Bode (1,1) string='off'
    options.legend ={'exp-data','sim-data'}
    options.exp_LineStyle (1,:) string = '-o';
    options.exp_LineWidth (1,1) {mustBeNumeric} = 2;
    options.sim_LineStyle (1,:) string = '-o';
    options.sim_LineWidth (1,1) {mustBeNumeric} = 1;
    options.xlim_switch(1,1) string ='off';
end

Nyquist=options.Nyquist;
Bode=options.Bode;
legend_text=options.legend;
exp_LineStyle=options.exp_LineStyle;
exp_LineWidth=  options.exp_LineWidth;
sim_LineStyle=  options.sim_LineStyle;
sim_LineWidth=  options.sim_LineWidth;
xlim_switch=options.xlim_switch;
sim_data=sim_circuit(x);
f=w./(2*pi);

switch  Nyquist
    case 'on'
        figure
        plot(real(exp_data),-imag(exp_data),exp_LineStyle,LineWidth=exp_LineWidth);
        hold on
        plot(real(sim_data),-imag(sim_data),sim_LineStyle,LineWidth=sim_LineWidth)
        grid on;
        xmax=ceil(max([max(real(exp_data)),max(-imag(exp_data)),max(real(sim_data)),max(-imag(exp_data))]));
        switch xlim_switch
            case 'on'
                xlim([0,xmax]);ylim([ymin,xmax]);
            otherwise
        end
        xlabel('Real(ohm)'),ylabel('-Imag(ohm)'),legend(legend_text);
    case 'off'

    otherwise
        fprintf('plot Nyquist is not surpport %s\n',Nyquist)

end

switch Bode
    case 'on'
        figure
        semilogx(f,real(exp_data),exp_LineStyle,LineWidth=exp_LineWidth);
        hold on
        semilogx(f,-imag(exp_data),exp_LineStyle,LineWidth=exp_LineWidth);

        semilogx(f,real(sim_data),sim_LineStyle,LineWidth=sim_LineWidth);
        semilogx(f,-imag(sim_data),sim_LineStyle,LineWidth=sim_LineWidth);
        grid on; hold off
        xlabel('f(Hz)'),ylabel('Real,-Imag (ohm)');
        legend({['Real-',legend_text{1}],['Imag-',legend_text{1}],['Real-',legend_text{2}],['Imag-',legend_text{2}]})
    case 'real'
        figure
        semilogx(f,real(exp_data),exp_LineStyle,LineWidth=exp_LineWidth);
        hold on
        semilogx(f,real(sim_data),sim_LineStyle,LineWidth=sim_LineWidth);
        grid on;
        xlabel('f(Hz)'),ylabel('Real(ohm)');
        legend(legend_text)

    case 'imag'
        figure
        semilogx(f,-imag(exp_data),exp_LineStyle,LineWidth=exp_LineWidth);
        hold on
        semilogx(f,-imag(sim_data),sim_LineStyle,LineWidth=sim_LineWidth);
        grid on;
        xlabel('f(Hz)'),ylabel('-Imag (ohm)');
        legend(legend_text)

    case 'module'
        figure
        semilogx(f,abs(exp_data),exp_LineStyle,LineWidth=exp_LineWidth);
        hold on
        semilogx(f,abs(sim_data),sim_LineStyle,LineWidth=sim_LineWidth);
        grid on;
        xlabel('f(Hz)'),ylabel('|Z| (ohm)');
        legend({'module'})
    case 'off'

    otherwise
        fprintf('plot Bode is not surpport %s\n',Bode)

end


end