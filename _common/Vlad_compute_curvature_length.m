function [k,L]=Vlad_compute_curvature_length(px,py,pz)
    [L,~,K]=curvature([px, py, pz]);
    ks=sqrt(K(:,1).^2 + K(:,2).^2 + K(:,3).^2);
    k=(1/L(end))*(sum(ks(2:end).*(L(2:end)-L(1:end-1)),'omitnan'));
    L=L(end);
end