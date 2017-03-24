function SetFigure(x,y)
box off
%sets up a figure of dimensions x,y
set(gcf, 'Color', [1 1 1])                  %Set background color
% try
%  xx=get(gca,'Title');
%  set(gcf, 'Name', xx)        %Name the chart
% catch errmsg='no title'
% end
% set(gcf, 'PaperSize', [20 20])                 %Set paper size
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
[a] =  get(gcf, 'PaperPosition');
set(gcf, 'PaperPosition', [a(1) a(2) x y]);       
%set(gcf,'PaperUnits','centimeters', 'Paperposition',[1 1 16 16]); %set the paper position in cm
set(gcf,'units','centimeters','position',get(gcf,'paperposition')); %this command shifts the 4 corners the number of units you tell it from where it is now.  the 3rd and 4th will actually alter the size 0 0 0 0 implements what you said in fig size.
