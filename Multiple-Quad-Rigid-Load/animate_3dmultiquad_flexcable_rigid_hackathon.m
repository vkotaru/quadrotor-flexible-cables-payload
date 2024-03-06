function animate_3dmultiquad_flexcable_rigid_hackathon(torig, xorig,xLcom, data)
% function to animate Quadrotor with suspended flexible cable
%
% Last Updated: 23-May-2016
% ===========================================================
params = data.params;
global vertices
vertices = data.vert;

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
    figure_x_limits = [-5 5];
    figure_y_limits = [-5 5];
    figure_z_limits = [-5 5] ;
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
    
%     sp1 = subplot(2,2,1,'Tag','f1');
%     sp2 = subplot(2,2,2,'Tag','f2');
%     sp3 = subplot(2,2,3,'Tag','f3');
%     sp4 = subplot(2,2,4,'Tag','f4');


    axes1 = axes;
    set(axes1,'XLim',figure_x_limits,'YLim',figure_y_limits);
%     set(axes1,'Position',[0 0 1 1]);
%     set(axes1,'Color','w');
    %set(axes1,'TickDir','out');
    axis equal ;
%     
    box on;

%     xLdorig = xLd ;
    for i = 1:size(xorig,3)
        [t, tempx] = even_sample(torig, xorig(:,:,i), RATE);
         t = t+torig(1) ;
         x(:,:,i) = tempx;
         [t, tempxLcom] = even_sample(torig, xLcom(:,:,i), RATE);
         xcom(:,:,i) = tempxLcom;
    end
    
hist = 250 ;
MAKE_MOVIE = 0;

% aviobj = avifile('sample2.avi','compression','None');
    for i=1:length(t)
        %set(axes1,'XLim',figure_x_limits+pH(i,1)) ;
%         ax1 = findobj(fig1,'Type','axes','-and','Tag','f1');
        drawone(axes1, x(i,:,:),xcom(i,:,:),params);
%         ax2 = findobj(fig1,'Type','axes','-and','Tag','f2');
%         drawone(ax2, x(i,:,:),xcom(i,:,:),params,2);
%         ax3 = findobj(fig1,'Type','axes','-and','Tag','f3');
%         drawone(ax3, x(i,:,:),xcom(i,:,:),params,3);
%         ax4 = findobj(fig1,'Type','axes','-and','Tag','f4');
%         drawone(ax4, x(i,:,:),xcom(i,:,:),params,4);
%         
%         subplot(2,2,1);    
        plot3(xcom(max(1,i-hist):i, 1,1), xcom(max(1,i-hist):i, 2,1), xcom(max(1,i-hist):i, 3,1), 'b') ;
%         subplot(2,2,2);    
%         plot3(xcom(max(1,i-hist):i, 1,1), xcom(max(1,i-hist):i, 2,1), xcom(max(1,i-hist):i, 3,1), 'b') ;
%         subplot(2,2,3);    
%         plot3(xcom(max(1,i-hist):i, 1,1), xcom(max(1,i-hist):i, 2,1), xcom(max(1,i-hist):i, 3,1), 'b') ;
%         subplot(2,2,4);    
%         plot3(xcom(max(1,i-hist):i, 1,1), xcom(max(1,i-hist):i, 2,1), xcom(max(1,i-hist):i, 3,1), 'b') ;

%         s = sprintf('Running\n t = %1.2fs \n 1/%d realtime speed',t(i), RATE/25);
%         text(-1.3,2.4,s,'FontAngle','italic','FontWeight','bold');
        drawnow;
        set(axes1,'XLim',figure_x_limits,'YLim',figure_y_limits,'ZLim',figure_z_limits);
%         set(ax2,'XLim',figure_x_limits,'YLim',figure_y_limits,'ZLim',figure_z_limits);
%         set(ax3,'XLim',figure_x_limits,'YLim',figure_y_limits,'ZLim',figure_z_limits);
%         set(ax4,'XLim',figure_x_limits,'YLim',figure_y_limits,'ZLim',figure_z_limits);

    end
%     aviobj = close(aviobj);
end

function[parent]= drawone(parent, xd,xLcom,params)
tem = get(parent,'Children');
delete(tem);

points = [];
for j = 1:size(xd,3)   
    
    x = xd(1,:,j)';
 
        s(j).L = 0.175; %length of quadrotor boom
        s(j).R = 0.1; %radius of propeller prop

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
        
        points = [points;xL'];

    %     Omega = x(22:24) ;
        xQ = xL - sum(repmat(l(1:n),3,1).*q,2);
        xChain = xQ;
        for i = 1:n
           xChain(:,i+1) = xChain(:,i) + l(i)*q(:,i);
        end
    %     plot3(xQ(1), xQ(2), xQ(3), 'r.') ;

        BRW = R' ;

    point1 = BRW'*[s(j).L,0,0]';
    point2 = BRW'*[0,s(j).L,0]';
    point3 = BRW'*[-s(j).L,0,0]';
    point4 = BRW'*[0,-s(j).L,0]';

    nprop = 40;
    propangs = linspace(0,2*pi,nprop);
    proppts = s(j).R*BRW'*[cos(propangs);sin(propangs);zeros(1,nprop)];

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

    s(j).qhandle1 = line([wp1(1),wp3(1)],[wp1(2),wp3(2)],[wp1(3),wp3(3)]); hold on ;
    s(j).qhandle2 = line([wp2(1),wp4(1)],[wp2(2),wp4(2)],[wp2(3),wp4(3)]);
    set(s(j).qhandle1,'Color','k', 'LineWidth',lw1);
    set(s(j).qhandle2,'Color','k', 'LineWidth',lw1);

    s(j).hprop1 = plot3(prop1(1,:),prop1(2,:),prop1(3,:),'r-', 'LineWidth',lwp);
    s(j).hprop2 = plot3(prop2(1,:),prop2(2,:),prop2(3,:),'b-', 'LineWidth',lwp);
    s(j).hprop3 = plot3(prop3(1,:),prop3(2,:),prop3(3,:),'b-', 'LineWidth',lwp);
    s(j).hprop4 = plot3(prop4(1,:),prop4(2,:),prop4(3,:),'b-', 'LineWidth',lwp);

    s(j).hload = plot3(xL(1,1), xL(2,1), xL(3,1), 'k.', 'LineWidth',lwl) ;
    s(j).xLcom = plot3(xLcom(1,1,j), xLcom(1,2,j), xLcom(1,3,j), 'ko', 'LineWidth',lwl) ;
    s(j).cable = plot3([wp_cable_attach(1) xChain(1,:)], [wp_cable_attach(2) xChain(2,:)], [wp_cable_attach(3) xChain(3,:)], '-k.') ;
    s(j).cable_attach = plot3([wp(1) wp_cable_attach(1)], [wp(2) wp_cable_attach(2)], [wp(3) wp_cable_attach(3)], 'r') ;

    offset = 0.5;

    % axis equal
    % set(handles.axes10,'Xlim',[max(min(s.xs(s.goodv))-offset,-20), min(max(s.xs(s.goodv))+offset,20)]);
    % set(handles.axes10,'Ylim',[max(min(s.ys(s.goodv))-offset,-20), min(max(s.ys(s.goodv))+offset,20)]);
    % set(handles.axes10,'Zlim',[max(min(s.zs(s.goodv))-offset,-20), min(max(s.zs(s.goodv))+offset,20)]);
    grid on
%     if v == 1
%        view(0,0); 
%     end
%     if v == 2
%        view(90,0); 
%     end
%     if v == 3
%        view(0,90); 
%     end
%     if v == 4
%        view([-1,-1,1]);
%     end
    view([-1,-1,1])
%     view(0,0);
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
end
global vertices

xLnew = mean(points)
R = Polyhedron('V',xLnew+vertices);
R.plot('Color','r','alpha',0.5);
% h = fill3(points(:,1),points(:,2),points(:,3),'y');
% set(h,'facealpha',.5)

end