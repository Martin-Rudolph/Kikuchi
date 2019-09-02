function [Xp,Yp,X,Y,Z] = stereoProj(Phi,Psi)
    % Resize scalar values
    % This allows for example a Psi vector with constant Phi and vise versa.
    % Not necessary -> user friendly option.
    % Otherwise the user has to pass input vectors (Phi,Psi) of the same size.
    if numel(Psi)==1
        Psi = Psi*ones(size(Phi));
    end
    if numel(Phi)==1
        Phi = Phi*ones(size(Psi));
    end
    % calculate the cartesian from the spherical coordinates (r=1)
    X = cosd(Phi).*sind(Psi);
    Y = sind(Phi).*sind(Psi);
    Z = cosd(Psi);
    % calculate the stereographic projection
    P = 1-Z./(1+Z);
    Xp = X.*P;
    Yp = Y.*P;
end

