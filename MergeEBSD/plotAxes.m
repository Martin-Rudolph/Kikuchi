function plotAxes(AB,Color,Label,uh,vk,wl)
%PLOTAXES Summary of this function goes here
%   Detailed explanation goes here
    Hold = ishold;
    Points = false;
    if nargin==1
        Color = 'b';
    end
    if ischar(Color) && strcmpi(Color,'auto')
        Color = lines(size(AB,1));
    end
    if ischar(Color)
        Color = {Color};
    end
    if iscell(Color) && numel(Color)==1
        Color = repmat(Color,size(AB,1),1);
    elseif isnumeric(Color) && size(Color,1)==1
        Color = repmat(Color,size(AB,1),1);
    end
    if nargin<3
        Label = {'a','b','c'};
    end
    if ischar(Label)
        switch(lower(Label))
            case {'reciprocal','rec','rep','rp','re','r'}
                Label = {'a^*','b^*','c^*'};
            case {'lattice', 'latt', 'lat', 'l'}
                N = (diff(uh)+1)*diff(vk+1)+diff(wl+1)+1;
                Label = zeros(N,3);
                i = 1;
                for u=uh(1):uh(2)
                    for v=vk(1):vk(2)
                        for w=wl(1):wl(2)
                            i=i+1;
                            Label(i,:) = [u,v,w];
                        end
                    end
                end
                AB = Label*AB;
                Points = true;
                if iscell(Color)
                    Color = Color{1};
                else
                    Color = Color(1,:);
                end
        end
    end
    if isnumeric(Label)
        Label = strtrim(cellstr(num2str(Label)));
        Label = regexprep(Label,'\s+', ' ');
    end
    if Points == true
        plot3(AB(:,1),AB(:,2),AB(:,3),'.','Color',Color)
        hold on
        text(AB(:,1),AB(:,2),AB(:,3),Label, 'HorizontalAlignment','right', 'VerticalAlignment','bottom', 'FontSize',8, 'Color',Color)
    else
        for i=1:size(AB,1)
            if iscell(Color)
                Col = Color{i};
            else
                Col = Color(i,:);
            end
            quiver3(0,0,0,AB(i,1),AB(i,2),AB(i,3),'Color',Col,'AutoScale','off')
            if i==1
                hold on
            end
            if i<=numel(Label)
                text(AB(i,1),AB(i,2),AB(i,3),[' ',Label{i}], 'Color',Col, 'FontWeight','bold', 'FontSize',12)
            end
        end
    end
    set(gca,'DataAspectRatio',[1 1 1], 'PlotBoxAspectRatio',[1 1 1])
    xlabel('x')
    ylabel('y')
    zlabel('z')
    view(3)
    if ~Hold
        hold off
    end
end
