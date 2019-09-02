clear all; close all; addpath('BKD');
%% Input
% Plot type
PlotType = 'provisional'; % all,result,provisional
% Reconstruction type
Recon = 'overlayed'; % cut, overlayed, average, points
% Resolution of reconstruction (not used for points)
GridRes = 300;
% Projection type
Projection = 'stereo'; % sphere, stereo
% Stereo grid
GridColor = 'none';
GridStep = [30,10];
% Define view direction of Kikuchi spehere
VD = [0,90];
% Selection radius (relative)
Rsel = 0.95;
% Pattern list
% {Image, Relative distance (Rrel) between specimen and screen with respect to screen diameter (D), Relative position of pattern center [x,y] -> Origin upper-left, Rotation angle (Phi0) in degree}
PatternList = {'Mart0.bmp', 0.6070, [0.44945,0.3101], 0, 70;
%                'Mart30.bmp', 0.6099, [0.4726,0.3182], -30, 70; % gut passend sind 0° und 30°
%               'Mart60.bmp', 0.6115, [0.4910,0.3070], -60, 70;
%                'Mart90.bmp', 0.6028, [0.4902,0.3096], -90, 70;
%                'Mart150.bmp', 0.5998 [0.4925,0.2949], -150, 70;
%                'Mart_180.bmp', 0.60542, [0.4890,0.2764], -180, 70;
%                'Mart240.bmp', 0.6042, [0.4683,0.2907], -260, 70 % pattern list Name, Center (x,y) [relative position], angle [degree]
};
% PatternList = {'Test_0.bmp', 0.6, [0.5, 0.25], 0;
%                'Test_30.bmp', 0.6, [0.5, 0.25], 30;
%                'Test_60.bmp', 0.6, [0.5, 0.25], 60;
%                'Test_0.bmp', 0.6, [0.5, 0.25], 90;
%                'Test_0.bmp', 0.6, [0.5, 0.25], 180;
%                'Test_0.bmp', 0.6, [0.5, 0.25], 270;
%                };

%% Actual Code do not edit

% Gray Values
MaxGrayVal = 255;

% Rotation matrices
Rz = @(Phi) [ cosd(Phi) sind(Phi) 0;
             -sind(Phi) cosd(Phi) 0;
                      0         0 1];

Ry = @(Phi) [ cosd(Phi) 0 -sind(Phi);
                      0 1          0;
              sind(Phi) 0  cosd(Phi)]; 

Rx = @(Phi) [1          0         0;
             0  cosd(Phi) sind(Phi);
             0 -sind(Phi) cosd(Phi)];
         
% Radius of Kikuchi sphere
R = 1;
Part = 0;

XspAll = [];
YspAll = [];
ZspAll = [];
IspAll = [];
for i = 1:size(PatternList,1)
    % Open pattern and load information
    Pattern = imread(PatternList{i,1});
    Rrel = PatternList{i,2};
    PCn  = PatternList{i,3};
    Phi0 = PatternList{i,4};
    Psi0 = PatternList{i,5};
    % Calculate screen radius
    D = R/Rrel;
    % Define Kikuchi sphere 
    [Kxx,Kyy,Kzz] = sphere(20); 
    Kxx = Kxx*R;
    Kyy = Kyy*R;
    Kzz = Kzz*R;
    if size(Pattern,3)>1
        Pattern = rgb2gray(Pattern);
    end
    % Calculate pattern center for image
    PCpix = PCn.*fliplr(size(Pattern));
    
    % Open figure and plot pattern
    if ~strcmpi(PlotType,'result')
        Fig = figure(i);
        set(Fig,'Units','normalized', 'Position',[0, 0.1, 1, 0.8])
        subplot(2,2,1)
        imshow(Pattern)
        hold on 
            plot(PCpix(1),PCpix(2),'ro','MarkerFaceColor','r')
            axis on
            set(gca,'PlotBoxAspectRatio',[1 1 1], 'DataAspectRatio',[1 1 1])
        hold off
    end
    
    % Define laboratory and specimen coordinate axes
    L  = R*[1 0 0;
           0 1 0
           0 0 1];
    S  = L*Rz(Phi0)*Ry(Psi0);
    PB = [0 0 1.2*R];
    
    % Set correct viewing direction of pattern
    Pattern = flipud(Pattern);
    Intensity = double(Pattern(:))./MaxGrayVal;
    % Create pattern mesh
    Y = linspace(0,D,size(Pattern,2));
    Z = linspace(0,D,size(Pattern,1));
    [YY,ZZ] = meshgrid(Y,Z);
    Y = YY(:);
    Z = ZZ(:);
    I = Intensity;
    X = R*ones(size(Y));
    
    %Calculate screen center
    Cs  = [R, 0.5*D, 0.5*D];
    % Remove edges
    Idx = find((Y-Cs(2)).^2+(Z-Cs(3)).^2>(Rsel*D/2).^2);
    Y(Idx) = []; Z(Idx) = []; X(Idx) = []; I(Idx) = [];  
    
    % Calculate pattern center
    PCn = [1, PCn(1), 1-PCn(2)];
    PCs = PCn.*[R, D, D];
    % Shift pattern center to L1
    Cs(2:3) = Cs(2:3)-PCs(2:3);
    Y = Y-PCs(2);
    Z = Z-PCs(3);
    PCs = [R, 0, 0];

    
    % Plot provisional result
    if ~strcmpi(PlotType,'result')
        subplot(2,2,2)
        surf(Kxx,Kyy,Kzz,'EdgeColor','none','FaceColor','g','FaceAlpha',0.1)
        hold on
            plot3([0,PB(1)], [0,PB(2)], [0,PB(3)], 'y','LineWidth',3)
            plotAxes(L,[0 0.5 0],{'L_1^{PC}','L_2','L_3^{PB}'})
            plotAxes(0.5*S,[0.75 0 0],{'S_1^{RD}','S_2^{TD}','S_3^{ND}'})
            scatter3(X,Y,Z,1,I*ones(1,3))
            plot3(PCs(1),PCs(2),PCs(3),'ro','MarkerFaceColor','r')
            plot3(Cs(1),Cs(2),Cs(3),'bo','MarkerFaceColor','b')
            set(gca,'PlotBoxAspectRatio',[1 1 1], 'DataAspectRatio',[1 1 1])
            view(VD(1),VD(2))
        hold off
    end
    
    % Use the specimen coordinate system as reference system
    S = R*[1 0 0;
           0 1 0;
           0 0 1];
    L = S*Ry(-Psi0)*Rz(-Phi0);
    PB = PB*Ry(-Psi0)*Rz(-Phi0);
    Cs = Cs*Ry(-Psi0)*Rz(-Phi0);
    PCs = PCs*Ry(-Psi0)*Rz(-Phi0);
    XYZ = [X,Y,Z]*Ry(-Psi0)*Rz(-Phi0);
    X = XYZ(:,1);
    Y = XYZ(:,2);
    Z = XYZ(:,3);
    
    % Calculate circle
    if strcmpi(PlotType,'cut')
        x0 = Cs(1);
        y0 = Cs(2);
        z0 = Cs(3);
        a  = PCs(1);
        b  = PCs(2);
        c  = PCs(3);
        v1 = (Cs-PCs);
        v1 = v1/norm(v1);
        v2 = cross(v1,PCs/norm(PCs));
        p  = Cs;
        t  = (0:360)';
        r  = Rsel*D/2;
        XYZcircle = ones(size(t))*p + r*cosd(t)*v1 + r*sind(t)*v2;
        Xc = XYZcircle(:,1);
        Yc = XYZcircle(:,2);
        Zc = XYZcircle(:,3);
    end
    
    
    % Plot provisional result
    if ~strcmpi(PlotType,'result')
        subplot(2,2,3)
        surf(Kxx,Kyy,Kzz,'EdgeColor','none','FaceColor','g','FaceAlpha',0.1)
        hold on
            plot3([0,PB(1)], [0,PB(2)], [0,PB(3)], 'y','LineWidth',3)
            plotAxes(L,[0 0.5 0],{'L_1^{PC}','L_2','L_3^{PB}'})
            plotAxes(0.5*S,[0.75 0 0],{'S_1^{RD}','S_2^{TD}','S_3^{ND}'})
            scatter3(X,Y,Z,1,I*ones(1,3))
            plot3(PCs(1),PCs(2),PCs(3),'ro','MarkerFaceColor','r')
            plot3(Cs(1),Cs(2),Cs(3),'bo','MarkerFaceColor','b')
            if strcmpi(PlotType,'cut')
                for T = 0:0.1:1.2
                    plot3(T*Xc,T*Yc,T*Zc,'b-','LineWidth',2)
                end
            end
            set(gca,'PlotBoxAspectRatio',[1 1 1], 'DataAspectRatio',[1 1 1])
            view(VD(1),VD(2))
        hold off
    end
    
    [Theta,Phi] = cart2sph(X,Y,Z);
    [X, Y, Z] = sph2cart(Theta,Phi,R*ones(size(Theta)));
    if strcmpi(PlotType,'cut')
        try    
            [Theta,Phi] = cart2sph(Xc,Yc,Zc);
            [Xc, Yc, Zc] = sph2cart(Theta,Phi,R*ones(size(Theta)));
            pp(:,i) = EllipseDirectFit([Xc,Yc]);
        %     Idx = find(pp(1,i).*X.^2+pp(2,i).*X.*Y+pp(3,i).*Y.^2+pp(4,i).*X+pp(5,i).*Y+pp(6,i)>0);
        %     Y(Idx) = []; Z(Idx) = []; X(Idx) = []; I(Idx) = [];
            if strcmpi(Recon,'cut')
                for j=i-1:-1:1
                    Idx = find(pp(1,j).*X.^2+pp(2,j).*X.*Y+pp(3,j).*Y.^2+pp(4,j).*X+pp(5,j).*Y+pp(6,j)<0);
                    Y(Idx) = []; Z(Idx) = []; X(Idx) = []; I(Idx) = [];
                end
            end
            warning('Cut does not work properly please use overlayed instead. A further development is unlikely due to the complexity of the deviation of the elliptical equation. See cone.mlx!')
        catch
            warning('Cut-option is under development and currently not available.')
        end
    end
    
    
    % Plot provisional result
    if ~strcmpi(PlotType,'result')
        subplot(2,2,4)
        surf(Kxx,Kyy,Kzz,'EdgeColor','none','FaceColor','g','FaceAlpha',0.1)
        hold on
            plot3([0,PB(1)], [0,PB(2)], [0,PB(3)], 'y','LineWidth',3)
            plotAxes(L,[0 0.5 0],{'L_1^{PC}','L_2','L_3^{PB}'})
            plotAxes(0.5*S,[0.75 0 0],{'S_1^{RD}','S_2^{TD}','S_3^{ND}'})
            scatter3(X,Y,Z,1,I*ones(1,3))
            plot3(PCs(1),PCs(2),PCs(3),'ro','MarkerFaceColor','r')
            plot3(Cs(1),Cs(2),Cs(3),'bo','MarkerFaceColor','b')
            plot3(Xc,Yc,Zc,'b-','LineWidth',2)
            plot3(Xc,Yc,zeros(size(Xc)),'c--')
%             pp = fit_ellipse(Xc,Yc,gca)
            set(gca,'PlotBoxAspectRatio',[1 1 1], 'DataAspectRatio',[1 1 1])
            view(VD(1),VD(2))
        hold off
    end
    
    Part(i) = numel(X);
    XspAll  = [XspAll; X];
    YspAll  = [YspAll; Y];
    ZspAll  = [ZspAll; Z];
    IspAll  = [IspAll; I];
end

S = R*[1 0 0;
       0 1 0;
       0 0 1];
   
if strcmpi(Projection,'stereo')
%     [Theta,Phi] = cart2sph(XspAll,YspAll,ZspAll);
%     Psi = Theta*180/pi-90;
%     Phi = Phi*180/pi; 
    Psi = acosd(ZspAll);
    Phi = (1-abs(sign(YspAll))+sign(YspAll)).*acosd(XspAll./sind(Psi));
    [XspAll,YspAll] = stereoProj(Phi,Psi);
    ZspAll = zeros(size(XspAll));
end

if ~strcmpi(PlotType,'provisional')
    figure;
    if strcmpi(Projection,'sphere')
        surf(Kxx,Kyy,Kzz,'EdgeColor','none','FaceColor','g','FaceAlpha',0.1)
        hold on
            plotAxes(0.5*S,[0.75 0 0],{'S_1^{RD}','S_2^{TD}','S_3^{ND}'})
        hold off
    end
    switch Recon
        case 'points'
            hold on
                scatter3(XspAll,YspAll,ZspAll,1,IspAll*ones(1,3))
                set(gca,'PlotBoxAspectRatio',[1 1 1], 'DataAspectRatio',[1 1 1])
                view(VD(1),VD(2))
            hold off
        case {'cut','average'}
            Xg = linspace(min(XspAll),max(XspAll),GridRes);
            Yg = linspace(min(YspAll),max(YspAll),GridRes);
            [XXg,YYg] = meshgrid(Xg,Yg);
            if strcmpi(Projection,'stereo')
                ZZg = zeros(size(XXg));
            else
                ZZg = real(sqrt(R^2-XXg.^2-YYg.^2));
            end
            IIg = griddata(XspAll,YspAll,IspAll,XXg,YYg);
            hold on
                surf(XXg,YYg,ZZg,IIg,'EdgeColor','none')
                colormap 'gray'
                caxis([0,1])
                set(gca,'PlotBoxAspectRatio',[1 1 1], 'DataAspectRatio',[1 1 1])
                view(VD(1),VD(2))
            hold off
        case 'overlayed'
            Part = [0, cumsum(Part)];
            for i=numel(Part):-1:2
                Xsp = XspAll(Part(i-1)+1:Part(i));
                Ysp = YspAll(Part(i-1)+1:Part(i));
                Zsp = ZspAll(Part(i-1)+1:Part(i));
                Isp = IspAll(Part(i-1)+1:Part(i));
                Xg = linspace(min(Xsp),max(Xsp),GridRes);
                Yg = linspace(min(Ysp),max(Ysp),GridRes);
                [XXg,YYg] = meshgrid(Xg,Yg);
                if strcmpi(Projection,'stereo')
                    ZZg = zeros(size(XXg))+(i-2)*0.001;
                else
                    ZZg = real(sqrt((R+(i-2)*0.001)^2-XXg.^2-YYg.^2));
                end
                IIg = griddata(Xsp,Ysp,Isp,XXg,YYg);
                hold on
                    surf(XXg,YYg,ZZg,IIg,'EdgeColor','none')
                    colormap 'gray'
                    caxis([0,1])
                    set(gca,'PlotBoxAspectRatio',[1 1 1], 'DataAspectRatio',[1 1 1])
                    view(VD(1),VD(2))
                hold off
            end
    end
    if strcmpi(Projection,'stereo')
        hold on
            stereoGrid(GridStep(1),GridStep(2),GridColor)
        hold off
    end
end