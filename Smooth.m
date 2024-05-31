% simple kernel smoothing for group state densities
% plot them all together
hold on
for i = 1:grid
    x(i) = i/grid;
end
for j = 16:19 % which group densities to plot
    n = ceil(100*rand);
    smooth = zeros(1,grid); 
    for i = 3:grid-2
        smooth(i) = (state(n,i-2)+state(n,i-1)+state(n,i)+state(n,i+1)+state(n,i+2))/5;
    end
    plot(x,smooth)
end