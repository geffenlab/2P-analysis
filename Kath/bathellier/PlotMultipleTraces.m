function h=PlotMultipleTraces(X,dT,offT,offY,h)

% Sizes
[m n p]=size(X);

colors = [   
    150 50 255
    233 217 143
    255 168 211
    250 120 0
    0   0   0
    0   0   255
    0   100 160
    0   255 0
    0   138 0
    191 0   191
    255 255 0
    200 0 0
    255 0   0
    ]/255;

% Size graph

t=repmat(permute((1:p)*dT,[1 3 2]),[m n 1]);
for i=1:n
    t(:,i,:)=t(:,i,:)+(p*dT+offT)*(i-1);
end
for i=1:m
    X(i,:,:)=X(i,:,:)+offY*(i-1);
end

if(nargin<5)
    figure; h=subplot(1,1,1);
    if(m==13)
        for i=1:m
            plot(squeeze(t(i,:,:))',squeeze(X(i,:,:))','Parent',h,'Color',colors(i,:))
            hold on
        end
        hold off
    else
        plot(reshape(t,[m*n p])',reshape(X,[m*n p])','Parent',h)
    end
else
    hold(h,'on')
    plot(reshape(t,[m*n p])',reshape(X,[m*n p])','-k','Parent',h)
    hold(h,'off')
end

xlim(h,[0 max(t(:))]);
