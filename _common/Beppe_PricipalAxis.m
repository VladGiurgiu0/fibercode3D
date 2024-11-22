function [PrincAxis] = Beppe_PricipalAxis(vol)

%generate points from the volume
[Xa, Ya, Za] = meshgrid(1:size(vol,2), 1:size(vol,1), 1:size(vol,3));
points = [Xa(:) Ya(:) Za(:) vol(:)];
nonzero = points(:,4)>0; %threshold to retain points
filtered = points(nonzero,1:3);

%   j=boundary(filtered(1,:)',filtered(2,:)',filtered(3,:)',0.1);
%generate an orientation vector from points
cen = mean(filtered,1);
[~,~,PrincAxis] = svd((filtered-cen),'econ'); %singular value decomposition, econ for speed
end