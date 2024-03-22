% clear
% for f = 1:10000
     clc;
     clear
    %% -----parameters------
    BoxL = 16.1875; %x length of box
    BoxW = 8.375;  % y length of box
    alphar = 2; %self propulsion coef
    alphaf = 2; %self propulsion coef
    Dr = 1;   %random motion coef
    Dtw = 1.1;   %toe twitch coef
    dt = 0.1; %time step
    T =  500; %number of timesteps
    beta = 1.5; %frog motivation
    vr0 = 1.7; %constant roach velocity 
    vf0 = 1.5; %constant frog velocity
    tong_rad = 2; %Tongue Radius
    
    M = 1; %number of roaches
    N = 1; %number of frogs
    %% -----initial position, orientation,velocity frog------
    for j = 1:N
        xf(1,j) = rand(1) * BoxL;%4.2 for calibration check
        yf(1,j) = rand(1) * BoxW;%6.3 for calibration check
        thetad = rand * 2 * pi;
        dxf(j) = cos(thetad);
        dyf(j) = sin(thetad);
        theta2 = rand * 2 * pi; 
        vfx(1,j) = 0;
        vfy(1,j) = 0;
        sight(j) = 0;
    end
    %% -----inital position, orientation,velocity roach-----
    for k = 1:M
        xr(1,k) = rand(1) * BoxL; %8.8 for calibration check
        yr(1,k) = rand(1) * BoxW; %1.6 for calibration check
        theta1 = rand * 2 * pi; 
        vrx(1,k) = vr0*cos(theta1);
        vry(1,k) = vr0*sin(theta1);
        roach_alive(k) = 1.0; % 1 means alive, 0 means dead
    end
    
    
    %% -----Roach Dynamics-----
    for i = 1:T 
        Time = (i-1)*dt;
        %-compute if frog in roach field of vision-
        for k = 1:M
        for j = 1:N
        if roach_alive(k) == 1
           theta(k,j) = acos(((xr(i,k)-xf(i,j))*dxf(j) + (yr(i,k)-yf(i,j))*dyf(j))/(sqrt((xr(i,k)-xf(i,j))^2 + (yr(i,k)-yf(i,j))^2))); 
        if abs(theta(k,j)) < pi/2
            phi(k,j) = 1.0;
        else
            phi(k,j) = 0.0;
        end
        if phi(k,j) == 1.0
            distance(k,j) = sqrt((xr(i,k)-xf(i,j))^2 + (yr(i,k)-yf(i,j))^2);
        else
            distance(k,j) = BoxL*100;
    
        end
        elseif roach_alive(k) == 0
                phi(k,j) = 0.0;
                distance(k,j) = BoxL*100;
        end
        end
        end
    
       [Q,I] = min(distance); %Q(j) = min distance for frog j
    
        %-new roach velocities- 
        for k = 1:M  %trying to add all the twitching terms before computing new velocity.
            if (roach_alive(k) == 1)
        twitchingx(k) = 0;
        twitchingy(k) = 0;
        for j = 1:N % Dont want normal random motion in this sum
            twitchingx(k) = twitchingx(k) + dt*(1-max(phi(k,j)))*Dtw*(xf(i,j)-xr(i,k))/sqrt((xf(i,j)-xr(i,k))^2+(yf(i,j)-yr(i,k))^2);
            twitchingy(k) = twitchingy(k) + dt*(1-max(phi(k,j)))*Dtw*(yf(i,j)-yr(i,k))/sqrt((xf(i,j)-xr(i,k))^2+(yf(i,j)-yr(i,k))^2);
        end
        end
        end
        for k = 1:M
            if roach_alive(k) == 1
        vrx(i+1,k) = vrx(i,k) + dt*(alphar*vrx(i,k)*(vr0*vr0-(vrx(i,k)^2+vry(i,k)^2))) + Dr*normrnd(0,sqrt(dt)) + twitchingx(k);
        vry(i+1,k) = vry(i,k) + dt*(alphar*vry(i,k)*(vr0*vr0-(vrx(i,k)^2+vry(i,k)^2))) + Dr*normrnd(0,sqrt(dt)) + twitchingy(k);
        %-Update bug position-
        xr(i+1,k) = xr(i,k) + dt*vrx(i+1,k);
        yr(i+1,k) = yr(i,k) + dt*vry(i+1,k);
            end
        end
    %% -----Frog Dynamics----- 
    for j = 1:N
        vfx(i+1,j) = max(phi(:,j))*vfx(i,j) + dt*(max(phi(:,j))*alphaf*vfx(i,j)*(vf0^2-(vfx(i,j)^2+vfy(i,j)^2))) + max(phi(:,j))*beta*(xr(i,I(j))-xf(i,j))/(sqrt((xr(i,I(j))-xf(i,j))^2 + (yr(i,I(j))-yf(i,j))^2)); 
        vfy(i+1,j) = max(phi(:,j))*vfy(i,j) + dt*(max(phi(:,j))*alphaf*vfy(i,j)*(vf0^2-(vfx(i,j)^2+vfy(i,j)^2))) + max(phi(:,j))*beta*(yr(i,I(j))-yf(i,j))/(sqrt((xr(i,I(j))-xf(i,j))^2 + (yr(i,I(j))-yf(i,j))^2));
        %-update frog position-
        xf(i+1,j) = xf(i,j) + dt*vfx(i+1,j);
        yf(i+1,j) = yf(i,j) + dt*vfy(i+1,j);
    end
    %% -----boundary conditions-----
    for k = 1:M
        if roach_alive(k) == 1
        if xr(i+1,k) < 0
            xr(i+1,k)=0.1;
            vrx(i+1,k) = -vrx(i+1,k);
        elseif xr(i+1,k) > BoxL
            xr(i+1,k)=BoxL-0.1;
            vrx(i+1,k) = -vrx(i+1,k);
        end
         if yr(i+1,k) < 0
            yr(i+1,k)=0.1;
            vry(i+1,k) = -vry(i+1,k);
        elseif yr(i+1,k) > BoxW
            yr(i+1,k)=BoxW-0.1;
            vry(i+1,k) = -vry(i+1,k);
         end
        end
    end
    for j = 1:N
         if xf(i+1,j) < 0
            xf(i+1,j)=0.1;
            vfx(i+1,j) = -vfx(i+1,j);
         elseif xf(i+1,j) > BoxL
            xf(i+1,j)=BoxL-0.1;
            vfx(i+1,j) = -vfx(i+1,j);
        end
         if yf(i+1,j) < 0
            yf(i+1,j)=0.1;
            vfy(i+1,j) = -vfy(i+1,j);
         elseif yf(i+1,j) > BoxW
            yf(i+1,j)=BoxW-0.1;
            vfy(i+1,j) = -vfy(i+1,j);
         end
    
        if (sqrt(vfx(i+1,j)^2+vfy(i+1,j)^2)>0.0)
            sight(j) = 1;
        end
        if (sight(j) == 1) %if the frog is moving, update unit vector for velocity
           % sight=0;
            thetad1(j) = atan2(vfy(i+1,j),vfx(i+1,j));
            dxf(j) = dxf(j) + dt*cos(thetad1(j));
            dyf(j) = dyf(j) + dt*sin(thetad1(j));
            thetad(j) = atan2(dyf(j),dxf(j));
            dxf(j) = cos(thetad(j));
            dyf(j) = sin(thetad(j));
        end
    end
    
    %% -----plot motion-----
       clf()
        hold on
       for k = 1:M
           if roach_alive(k) == 1
               plot(xr(i+1,k),yr(i+1,k),'ko','MarkerSize',10,"MarkerFaceColor",'k')
           end
       end
       for j = 1:N
        plot(xf(i+1,j),yf(i+1,j),'go','MarkerSize',75,"MarkerFaceColor","g") 
       end
       axis([0 BoxL 0 BoxW])
       xt = [BoxL-1.0, BoxL-.65, BoxL-.25];
       yt = [BoxW-7.75, BoxW-7.75,BoxW-7.75];
       str={'t = ',num2str(Time),'s'};
       text(xt,yt,str,'FontSize',13,'Color','blue' )
        
        for j = 1:N
            if (Q(j)<tong_rad)
    
                    'Frog Caught It, Yeaaaah Boy'
                    if ( sqrt((xf(i+1,j)-xr(i+1,I(j)))^2+(yf(i+1,j)-yr(i+1,I(j)))^2) < BoxL)
                   plot([xf(i+1,j) xr(i+1,I(j))], [yf(i+1,j) yr(i+1,I(j))], 'r-')
                    end
                    
                    roach_alive(I(j)) = 0.0; % Now Dead if in tong rad
                    roach_death_time(I(j)) = Time;
            end
           
        end
            hold off
            drawnow
    
            if max(roach_alive) == 0
                break
            end
%     vision(i) = sight(j);
    end
%     time(f,2) = vision(1);
%     time(f,1) =  Time;
% end
% writematrix(time,'Simulation2.12.24.csv');