function [cf,df] = plane_coef(plane)
% Calculates the coefficients and constant terms of the equation of a
% plane in implicit form
%
% This function is taken from https://www.mathworks.com/matlabcentral/fileexchange/48509-computational-geometry-toolbox
%

[~,m]=size(plane);
pdiff=diff(plane);
cf=zeros(m,1);
sign=1;
r=1:m;
for i=r
    cf(i)=sign*det(pdiff(:,r~=i));
    sign=-sign;
end
cf=cf/norm(cf);
df=-plane(1,:)*cf;