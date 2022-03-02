function ax = my_subplot(m,n,p,PosType)

if nargin < 4
    PosType = 'OuterPosition';
end

nplots = n*m;

if p > nplots
   error('Wrong plot index!') 
end

row_spacing = 1/m;
col_spacing = 1/n;

row_index = ceil(p/(nplots/m));
col_index = mod(p,(nplots/m));

if col_index == 0
    col_index = n;
end

row_start = 1-(row_index*row_spacing);
col_start = (col_index*col_spacing)-col_spacing;

%ax = axes('OuterPosition',[col_start row_start col_spacing row_spacing]);
ax = axes(PosType,[col_start row_start col_spacing row_spacing]);