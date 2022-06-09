function squarezoom(IMAGE, zoom, width,x,y)
% function that zoom in on a square with width and cordinate (x,y)
n = size(IMAGE,1);

if nargin < 4
    x = n/2-15;
    y = n/2-15;
elseif nargin < 3
    width = 29;
    x = n/2-15;
    y = n/2-15;
end


im     = IMAGE;

imzoom = imresize(im(x:x+width,y:y+width),zoom);
X      = n-zoom*(width+1);
Y      = n-zoom*(width+1);
WIDTH  = zoom*(width+1);

im(n-zoom*(width+1)+1:end,n-zoom*(width+1)+1:end) = imzoom;

imagesc(im)
hold on
rectangle('Position',[y,x,width,width],...
  'EdgeColor', 'white',...
  'LineWidth', 2,...
  'LineStyle','-')
rectangle('Position',[Y,X,WIDTH,WIDTH],...
  'EdgeColor', 'white',...
  'LineWidth', 3,...
  'LineStyle','-')
colormap gray
end