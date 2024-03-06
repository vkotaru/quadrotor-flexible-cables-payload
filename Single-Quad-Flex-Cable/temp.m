clear;

load('export_freq3.mat');
close all;
t = [0:0.1:1]';

for i = 1:length(export)
    
    output = export{i,1};
    xL = output.xL;
    xQ = output.xQ;
    l = sqrt(sum((xL(:,1:3)-xQ(:,1:3)).^2,2));
    
    figure(1); hold on;
    plot3(xL(:,1),xL(:,2),xL(:,3));
    plot3(xQ(:,1),xQ(:,2),xQ(:,3));
    xlabel('x-axis');ylabel('y-axis');zlabel('z-axis');
    grid on;
    
    figure(2); hold on;
    plot(t,l);
    grid on;
    
    output2 = export{i,2};
    xL2 = output2.xL;
    xQ2 = output2.xQ;
    l2 = sqrt(sum((xL2(:,1:3)-xQ2(:,1:3)).^2,2));
    
    figure(3); hold on;
    plot3(xL2(:,1),xL2(:,2),xL2(:,3));
    plot3(xQ2(:,1),xQ2(:,2),xQ2(:,3));
    xlabel('x-axis');ylabel('y-axis');zlabel('z-axis');
    grid on;
    
    figure(4); hold on;
    plot(t,l2);
    grid on;
    
end

