function [c] = LogCAxis(bounds)
caxis([log10(bounds(1)) log10(bounds(2))])
c = colorbar; shading flat
global fontsize
if size(fontsize,1) > 0
   set(c,'fontsize',fontsize)
end
if max(abs(bounds)) > 1
   n = floor(log10(bounds(1)));
   set(c,'YTickLabel',Num2CellStr(roundn(10.^(get(c,'YTick')),n)))
else
   set(c,'YTickLabel',Num2CellStr(10.^(get(c,'YTick'))))
end