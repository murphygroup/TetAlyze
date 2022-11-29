clc
clear all

figname = 'E:\BB Project\generative model images\2- Cell cycle generative model\N=1\010622 N=1 B1868 Poc1-mCh dividing cells_3-sld - 1005\Alignment.fig';

h1=openfig(figname, 'invisible');

ax = gca; 
h = findobj(gca,'Type','line'); 

x = [];
y = [];
z = [];
figure();
for i=1:length(h)   
    x = h(i).XData;
    y = h(i).YData;
    z = h(i).ZData;
    

    z_ = 2*(z - min(z))/(max(z) - min(z)) - 1;
    theta = atan(y./x);
    theta = mean(theta)/pi*180-90;
    
    plot3(x', y', z', 'o-');
    axis equal
    axis off
%     hold on
    xlim([-18 18])
    ylim([-18 18])
    zlim([-3 55])

    
end
