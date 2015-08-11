function [ sl,sh ] = lowpass_nl(s,L, mu)
% testor for doing a "lowpass" via nonlocal windows
% sl = argmin_y (1/2)*||s-y||^2 + (mu/2)*y'Ly
% s is 2d image


sl = zeros(size(s));
for i = 1:length(L)
    Ltemp = L{i};
    if ~isfield(Ltemp,'Maxeig')
        Meig = 1;
    else
        Meig = Ltemp.Maxeig;
    end        
    I1 = Ltemp.ind1(1):Ltemp.ind2(1);
    I2 = Ltemp.ind1(2):Ltemp.ind2(2);
    sp = reshape(s(I1,I2),length(I1)*length(I2),1);
    sp_hat = (Ltemp.phi)'*sp;
    sp_par = Ltemp.phi*sp_hat;
    sp_perp = sp-sp_par;    
    sl_patch = sp_perp/(1+mu*Meig) + Ltemp.phi*(sp_hat./(1+Ltemp.E));
    sl_patch = reshape(sl_patch,length(I1),length(I2));
    sl(I1,I2) = sl_patch;
end
sh = s-sl;



end

