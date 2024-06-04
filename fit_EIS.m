function [x,fav]=fit_EIS(w,exp_data,sim_circuit,param_0,options)
%Return the circuit paramter value and erros
%Inputs
%---------------
% w:= Angular frequency [1/s]
% exp_data:a complex matrix size:(1*length(w))
% sim_circuit:object circuit function handle
% param_0: the inital of circuit parameter,size:(1:N)
% options:is a structure with fields:
% method:slover method,the default is 'fminunc';Include:'lsqnonlin','fminunc','fmincon','pso'
% ub:the upper of circuit paramter, size:(1:N)
% lb:the lower of circuit paramter, size:(1:N)
% error_type:The type of error function to calculate,the default is 'Chi-2';Include:'Chi-2','RMSE'

arguments
    w (1,:) double;
    exp_data double;
    sim_circuit ;
    param_0 (1,:) double;
    options.method string ='fminunc';
    options.ub (1,:) double=[];
    options.lb (1,:) double=[];
    options.error_type string='Chi-2';

end
method=options.method;
ub=options.ub;
lb=options.lb;
error_type=options.error_type;

object_fun=@(param) error_cacl(exp_data,sim_circuit,param,error_type);
switch method
    case 'lsqnonlin'
        %         object_fun=@(param) error_cacl(exp_data,sim_circuit,param);
        object_fun=@(param) error_cacl_lsq(exp_data,sim_circuit,param,error_type);
        [x,fav]=lsqnonlin(object_fun,param_0,lb,ub);
    case 'fminunc'

        [x,fav]=fminunc(object_fun,param_0);
    case 'fmincon'
        [x,fav]=fmincon(object_fun,param_0,[],[],[],[],lb,ub);
    case 'pso'
        pso_options= optimoptions('particleswarm','SwarmSize',300,'HybridFcn',@fmincon,'MaxIterations',200);
        %          pso_options= optimoptions('particleswarm','SwarmSize',100,'MaxIterations',100);
        [x,fav]= particleswarm(object_fun,length(lb),lb,ub,pso_options);
    otherwise
        fprintf('solver method is not match %s\n',method)

end

    function error= error_cacl(exp_data,sim_circuit,param,error_type)

        sim_data=sim_circuit(param);
        %Chi-2
        Chi_2=sum(1*((real(sim_data)-real(exp_data)).^2)./abs(exp_data)+1*((imag(sim_data)-imag(exp_data)).^2)./abs(exp_data));
        %残差均方根
        RMSE=mean(((real(sim_data)-real(exp_data)).^2+(imag(sim_data)-imag(exp_data)).^2).^0.5./(abs(exp_data)));
        %平均百分比误差
        %         PE=sum(0.5*(real(sim_data)-real(exp_data))./abs(exp_data)+0.5*(imag(sim_data)-imag(exp_data))./abs(exp_data));
        switch error_type
            case 'Chi-2'
                error=Chi_2;
            case 'RMSE'
                error=RMSE;
                %             case 'PE'
                %                 error=PE;

        end

    end
    function error=error_cacl_lsq(exp_data,sim_circuit,param,error_type)
        sim_data=sim_circuit(param);
        Chi_2=  ((real(sim_data)-real(exp_data)).^2)+((imag(sim_data)-imag(exp_data)).^2);
        error=Chi_2;
        RMSE=((real(sim_data)-real(exp_data)).^2+(imag(sim_data)-imag(exp_data)).^2).^0.5./(abs(exp_data));

        switch error_type
            case 'Chi-2'
                error=Chi_2;
            case 'RMSE'
                error=RMSE;


        end

    end


end