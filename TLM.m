function Z=TLM(w,options)
%Return the impedance of TLM
% complex size(1:length(w))
%Inputs
%---------------
%w = Angular frequency [1/s]
%Ri:ion resistance of electrode;
%Ri:electronic resistance of electrode [ohm/cm]
%Zt:impedance for the interfacial of electrode [ohm*cm]
%L:length of electrode [cm]
%Ref:-----------------
% J. Bisquert et al. / Electrochemistry Communications 1 (1999) 429–435
% ECS Transactions, 69 (18) 71-80 (2015)
% https://doi.org/10.1021/acsmeasuresciau.2c00070
% https://mp.weixin.qq.com/s/6FoIglc22il9k_nTO7U4kA

arguments
    w (1,:) double
    options.Ri (1,:) double=[];
    options.Re (1,:) double=[];
    options.Zt (1,:) double=[];
    options.L (1,1) double=[];

end
Ri=options.Ri;
Re=options.Re;
Zt=options.Zt;
L=options.L;
if isempty(Re)
    Z=(Ri.*Zt).^0.5.*coth(L.*(Ri./Zt).^0.5);
else
    A=Re.*Ri./(Re+Ri);
    B=(Re.^2+Ri.^2)./(Re+Ri);
    k=(Zt./(Re+Ri)).^0.5;
    Z=A.*(L+2*k./sinh(L./k))+k.*B.*coth(L./k);
end

end