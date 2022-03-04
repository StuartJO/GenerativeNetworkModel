function Lbl = AddPlotLabel(axis,Label,FontSize,Xscale,Yscale,ScaleOrTrans)

if nargin < 6
    ScaleOrTrans = 1;
end

if nargin < 4
    if ScaleOrTrans
   Xscale = 1; 
    elseif ScaleOrTrans == 2
   Xscale = 0;      
    end
end

if nargin < 5
    if ScaleOrTrans
   Yscale = 1; 
    elseif ScaleOrTrans == 2
   Yscale = 0;      
    end
end

Pos = get(axis,'Position');

Lblsize = [0.0588 0.0588];

if ScaleOrTrans == 1

Lbl = annotation('textbox',[Pos(1)-(Lblsize(1)/Xscale) Pos(2)+(Pos(4)/Yscale) Lblsize],'String',Label,'EdgeColor','none','FontSize',FontSize);

else
    
Lbl = annotation('textbox',[Pos(1)+Xscale Pos(2)+Yscale Lblsize],'String',Label,'EdgeColor','none','FontSize',FontSize);
    
end