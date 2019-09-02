function stereoGrid(dPhi,dPsi,Color)
if nargin==0
    dPhi = 10;
    dPsi = 10;
elseif nargin==1
    dPsi = dPhi;
end
if nargin<3
    Color = [0.7,0.7,0.7];
end
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    zMax = max(get(gca,'ZLim'));
    for PsiGrid=0:dPsi:90
        Xs = cosd(0:360)*sind(PsiGrid);
        Ys = sind(0:360)*sind(PsiGrid);
        Zs = cosd(PsiGrid)*ones(size(0:360));
        Xp = Xs.*(1-Zs./(Zs+1));
        Yp = Ys.*(1-Zs./(Zs+1));
        Zp = ones(size(Xp))*zMax;
        plot3(Xp,Yp,Zp,'Color',Color)
    end
    
    for PhiGrid=0:dPhi:350
        X = cosd(PhiGrid)*sind(0:90);
        Y = sind(PhiGrid)*sind(0:90);
        Z = cosd(0:90);
        Xp = X.*(1-Z./(Z+1));
        Yp = Y.*(1-Z./(Z+1));
        Zp = ones(size(Xp))*zMax;
        plot3(Xp,Yp,Zp,'Color',Color)
    end
    
end

