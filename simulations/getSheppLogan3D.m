function n = getSheppLogan3D(simP)

% Nphantom = simP.Lphantom/simP.dx;
% Nphantom = Nphantom - mod(Nphantom,2);
% [~,E] = phantom3d('Modified Shepp-Logan',Nphantom);
% %arbitrary, but at most max param.dn
% E(:,1) = simP.dn*E(:,1);
% E(2,1) = -0.8*simP.dn;%.0.08 for simP.dn = 0.1;
% E(3,1) = -(simP.dn + E(2,1)); E(4,1) = E(3,1);
% E(5:end,1) = 0.3*simP.dn;% similar
% n = phantom3d(E, Nphantom);
% n = padarray(n,max((simP.Nroi - Nphantom)/2,0),0,'both') + simP.n0;

Nphantom = simP.Lphantom/simP.dx;
Nphantom = Nphantom - mod(Nphantom,2);
[~,E] = phantom3d('Modified Shepp-Logan',Nphantom);
%arbitrary, but at most max param.dn. The value of the ellipses are added
E(:,1) = 0.5*simP.dn*E(:,1);%largest ellipse "background"
%E(2,1) = -0.8*simP.dn;%.0.08 for simP.dn = 0.1;%It was the outer layer
E(3,1) = -0.1*simP.dn; E(4,1) = E(3,1);%inner ellipse ("the eyes")
E(5:end,1) = 0.5*simP.dn;% similar, inner ellipse, "mouth" and small craps above
E(end+1,:) = E(end,:);
E(end,5:7) = [0,0.01,-0.6]; E(end,2:4) = [0.05,0.1,0.05];%adds a bright spot at z=-0.7
E(end+1,:) = E(end,:); E(end,7) = 0.65;%adds a bright spot at z= 0.7
E = E([1,3:end],:);
n = phantom3d(E, Nphantom);
n = padarray(n,max((simP.Nroi - Nphantom)/2,0),0,'both') + simP.n0;

end
