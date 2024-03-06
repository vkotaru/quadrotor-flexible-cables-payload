 function animate_3dquad_flexcable(torig, xorig,params)
% function to animate Quadrotor with suspended flexible cable
%
% Last Updated: 23-May-2016
% ===========================================================

% Extracting parameters
% =====================
n = params.n; % No. of links
m = params.m; % Mass of links
l = params.l; % Length of links
mQ = params.mQ; % Quad-rotor mass
g = params.g ; 
e1 = params.e1 ;
e2 = params.e2 ;
e3 = params.e3 ;
J = params.J;


RATE = 25 * 1;
%==========================================
% initialize the animation figure and axes
%==========================================
    figure_x_limits = [-200 200];
    figure_y_limits = [-200 200];
    figure_z_limits = [-200 200] ;
    fig1 = figure;

    set(0,'Units','pixels')
    scnsize = get(0,'ScreenSize');

    screen_width = scnsize(3);
    screen_height = scnsize(4);

    % find the minimum scaling factor
    figure_x_size = figure_x_limits(2) - figure_x_limits(1);
    figure_y_size = figure_y_limits(2) - figure_y_limits(1);

    xfactor = screen_width/figure_x_size;
    yfactor = screen_height/figure_y_size;

    if (xfactor < yfactor)
      screen_factor = 0.5*xfactor;
    else
      screen_factor = 0.5*yfactor;
    end

    % calculate screen offsets
    screen_x_offset = (screen_width - screen_factor*figure_x_size)/2;
    screen_y_offset = (screen_height - screen_factor*figure_y_size)/2;

    % draw figure and axes
    set(fig1,'Position', [screen_x_offset screen_y_offset screen_factor*figure_x_size screen_factor*figure_y_size]);
    set(fig1,'MenuBar', 'none');
    axes1 = axes;
    set(axes1,'XLim',figure_x_limits,'YLim',figure_y_limits);
%     set(axes1,'Position',[0 0 1 1]);
%     set(axes1,'Color','w');
    %set(axes1,'TickDir','out');
    axis equal ;
    
    box on;

    x = xorig ;
%     xLdorig = xLd ;
    [t, x] = even_sample(torig, x, RATE);
    t = t+torig(1) ;
    
hist = 250 ;
MAKE_MOVIE = 0;

% aviobj = avifile('sample2.avi','compression','None');
    for i=1:length(t)
        %set(axes1,'XLim',figure_x_limits+pH(i,1)) ;
        drawone(axes1, x(i,:)',params);
        plot3(x(max(1,i-hist):i, 1), x(max(1,i-hist):i, 2), x(max(1,i-hist):i, 3), 'b') ;

%         s = sprintf('Running\n t = %1.2fs \n 1/%d realtime speed',t(i), RATE/25);
%         text(-1.3,2.4,s,'FontAngle','italic','FontWeight','bold');
        drawnow;
        set(axes1,'XLim',figure_x_limits,'YLim',figure_y_limits,'ZLim',figure_z_limits);

    end
%     aviobj = close(aviobj);
end

function drawone(parent, x,params)
    tem = get(parent,'Children');
    delete(tem);
    
    s.L = 0.175; %length of quadrotor boom
    s.R = 0.1; %radius of propeller prop
    
% Extracting parameters
% =====================
n = params.n; % No. of links
m = params.m; % Mass of links
l = params.l; % Length of links
mQ = params.mQ; % Quad-rotor mass
g = params.g ; 
e1 = params.e1 ;
e2 = params.e2 ;
e3 = params.e3 ;
J = params.J;

    
    % Extract state x
    % ===============
    xL = x(1:3);
    vL = x(4:6);
    R = reshape(x(7:15),3,3);
    Omega = x(16:18);
    q = reshape(x(19:3*n+18),3,n);
    omega = reshape(x(3*n+19:6*n+18),3,n);

%     Omega = x(22:24) ;
    xQ = xL - sum(repmat(l(1:n),3,1).*q,2);
    xChain = xQ;
    for i = 1:n
       xChain(:,i+1) = xChain(:,i) + l(i)*q(:,i);
    end
%     plot3(xQ(1), xQ(2), xQ(3), 'r.') ;

    BRW = R' ;

point1 = BRW'*[s.L,0,0]';
point2 = BRW'*[0,s.L,0]';
point3 = BRW'*[-s.L,0,0]';
point4 = BRW'*[0,-s.L,0]';

nprop = 40;
propangs = linspace(0,2*pi,nprop);
proppts = s.R*BRW'*[cos(propangs);sin(propangs);zeros(1,nprop)];

wp = xQ ;
wp1 = wp + point1;
wp2 = wp + point2;
wp3 = wp + point3;
wp4 = wp + point4;
wp_cable_attach = wp + BRW'*[0;0;0] ; %[0;0;-params.cable_attach_d];

prop1 = proppts + wp1*ones(1,nprop);
prop2 = proppts + wp2*ones(1,nprop);
prop3 = proppts + wp3*ones(1,nprop);
prop4 = proppts + wp4*ones(1,nprop);

lwp = 2 ;
lw1 = 1 ;
lwc = 2 ;
lwl = 2 ;

s.qhandle1 = line([wp1(1),wp3(1)],[wp1(2),wp3(2)],[wp1(3),wp3(3)]); hold on ;
s.qhandle2 = line([wp2(1),wp4(1)],[wp2(2),wp4(2)],[wp2(3),wp4(3)]);
set(s.qhandle1,'Color','k', 'LineWidth',lw1);
set(s.qhandle2,'Color','k', 'LineWidth',lw1);

s.hprop1 = plot3(prop1(1,:),prop1(2,:),prop1(3,:),'r-', 'LineWidth',lwp);
s.hprop2 = plot3(prop2(1,:),prop2(2,:),prop2(3,:),'b-', 'LineWidth',lwp);
s.hprop3 = plot3(prop3(1,:),prop3(2,:),prop3(3,:),'b-', 'LineWidth',lwp);
s.hprop4 = plot3(prop4(1,:),prop4(2,:),prop4(3,:),'b-', 'LineWidth',lwp);

s.hload = plot3(xL(1,1), xL(2,1), xL(3,1), 'ko', 'LineWidth',lwl) ;
s.cable = plot3([wp_cable_attach(1) xChain(1,:)], [wp_cable_attach(2) xChain(2,:)], [wp_cable_attach(3) xChain(3,:)], '-k.') ;
s.cable_attach = plot3([wp(1) wp_cable_attach(1)], [wp(2) wp_cable_attach(2)], [wp(3) wp_cable_attach(3)], 'r') ;

offset = 0.5;

% axis equal
% set(handles.axes10,'Xlim',[max(min(s.xs(s.goodv))-offset,-20), min(max(s.xs(s.goodv))+offset,20)]);
% set(handles.axes10,'Ylim',[max(min(s.ys(s.goodv))-offset,-20), min(max(s.ys(s.goodv))+offset,20)]);
% set(handles.axes10,'Zlim',[max(min(s.zs(s.goodv))-offset,-20), min(max(s.zs(s.goodv))+offset,20)]);
grid on
view([-1,-1,1])
xlabel('X')
ylabel('Y')
zlabel('Z')


end