function alpha_cmap = findAlphaColor(cmap,alpha)

% Sometimes you want to find the color of a colormap when an alpha value is
% applied. This code does just that. Basically the alpha value shifts the
% color towards being more white. This code calculates the difference
% between the current color and white (i.e. an RGB value of [1 1 1]),
% applies the alpha value to that difference and then adds it to the
% current color.

alpha_cmap = ((1-cmap)*(1-alpha))+cmap;