x = 0:pi/100:2*pi;
x_c = 0:pi/20:2*pi;
plot(x_c,sin(x_c))

hold on
y2 = sin(x_c) + 0.5;
plot(x_c,y2,'.')

y3 = sin(x_c) - 0.5;
plot(x_c,y3,'.')


%x = linspace(0,2 * pi,1000);
%y = linspace(-1.5,1.5,1000);
%[X,Y] = meshgrid(x,y);
%Z = (0.5 + Y + sin(X)) * (0.5 - sin(X) - Y);

%contour(X,Y,Z,50)
