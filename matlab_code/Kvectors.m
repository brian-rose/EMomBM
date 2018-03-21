function [K1 K2 K3] = Kvectors(K,phiu,phiq,param)

cosK = cos(phiu).*K;  
K1  = cosK(1:end-1)./cos(phiq)*param.delt/param.a^2/param.dphi^2;
K3 = cosK(2:end)./cos(phiq)*param.delt/param.a^2/param.dphi^2;
K2 = [0; K1(2:end)]+ [K3(1:end-1);0];