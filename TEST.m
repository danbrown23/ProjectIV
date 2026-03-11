% Plot a parabola and a marker
x = -10:0.5:10;
y = x.^2;
p = plot(x,y,"-o",MarkerFaceColor="red");

% Move the marker along the parabola and capture frames in a loop
for i=1:41
    p.MarkerIndices = i;
    exportgraphics(gca,"parabola.gif",Append=true)
end