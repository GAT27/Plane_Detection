function [P,R,C,O,plane] = PlaneFiltering_final(I,aspect,k_max,l,eta,S,epsilon,a_in)
%Takes in a depth map I
%aspect     =   Aspect ratio
%k_max      =   Maximum number of iterations
%l          =   Number of local samples (should be greater than 3)
%eta        =   Neighborhood for global samples (in pixels)
%S          =   Plane size in world space for local samples
%epsilon    =   Maximum plane offset error for inliers
%a_in       =   Minimum fraction of inliers to accept local sample
%all the above parameters need general tweeking to get best plane tracking

%Horizontal field of view and vertical field of view of kinect
fh=57;
fv=43;
%Setup depth map and aquire height and width
imagesc(I);
hold on;
I = double(I);
[h,w] = size(I);
w = w/3;
%List of 3D points
P = [];
%List of normals
R = [];
%List of polygons
C = [];
%List of outliers
O = [];
%Found plane
plane = [];

%Setup height and width of plane
plane_h = S * sqrt((aspect(2)^2) / ((aspect(1)^2) + (aspect(2)^2)));
plane_w = plane_h * (aspect(1) / aspect(2));

%Continue sampling until either max sampling reached
n = 0;
k = 0;

while (k < k_max)
    %Display iteration count and place a center marker
    plot(w/2,h/2,'rx','MarkerSize',20);
    k = k + 1;
    display(k);
    display('||||||||||||||||||');
    
    %Choose 3 random points on depth map
    d0 = [randi([-eta,eta]) + w/2,randi([-eta,eta]) + h/2];
    d0 = round(d0);
    %display(d0);
    if ((d0(1)>w) || (d0(1)<=0))
        display('skip d0');
        continue;
    end
    if (d0(2)>h) || (d0(2)<=0)
        display('skip d0');
        continue;
    end
    
    d1 = d0 + [randi([-eta,eta]),randi([-eta,eta])];
    d1 = round(d1);
    %display(d1);
    if (d1(1)>w) || (d1(1)<=0)
        display('skip d1');
        continue;
    end
    if (d1(2)>h) || (d1(2)<=0)
        display('skip d1');
        continue;
    end
    
    d2 = d0 + [randi([-eta,eta]),randi([-eta,eta])];
    d2 = round(d2);
    %display(d2);
    if (d2(1)>w) || (d2(1)<=0)
        display('skip d2');
        continue;
    end
    if (d2(2)>h) || (d2(2)<=0)
        display('skip d2');
        continue;
    end
    
    %Recreate the points into 3D points, 0 is at the center
    p0 = [I(d0(2),d0(1)) * ((d0(1)/(w-1))-0.5) * tand(fh/2) ...
          I(d0(2),d0(1)) * ((d0(2)/(h-1))-0.5) * tand(fv/2) ...
          I(d0(2),d0(1))];
    p0 = (p0 .* [w/2,h/2,1]) + [w/2,h/2,0];
      
    p1 = [I(d1(2),d1(1)) * ((d1(1)/(w-1))-0.5) * tand(fh/2) ...
          I(d1(2),d1(1)) * ((d1(2)/(h-1))-0.5) * tand(fv/2) ...
          I(d1(2),d1(1))];
    p1 = (p1 .* [w/2,h/2,1]) + [w/2,h/2,0];
      
    p2 = [I(d2(2),d2(1)) * ((d2(1)/(w-1))-0.5) * tand(fh/2) ...
          I(d2(2),d2(1)) * ((d2(2)/(h-1))-0.5) * tand(fv/2) ...
          I(d2(2),d2(1))];
    p2 = (p2 .* [w/2,h/2,1]) + [w/2,h/2,0];
    
    %Compute plane normals and setup to find relative points, but skip if
    %the depth average is 0
    rc = cross(p1-p0,p2-p0);
    r = rc ./ abs(rc);
    r(isnan(r)) = 0;
    z_avr = (p0(3) + p1(3) + p2(3)) / 3;
    if z_avr==0
        display('skip z');
        continue;
    end
    %display(z_avr);
    
    %Projection parameter
    plane_d = 31.5 + 126 * (1 - (z_avr/255));
    wr = ((1 * sind(90-fh)) / sind(fh));
    hr = ((1 * sind(90-fv)) / sind(fv));
    w_pri = 2 * (hr/wr * plane_w * plane_d) / 31.5;
    h_pri = 2 * (hr/wr * plane_h * plane_d) / 31.5;
    %rectangle('Position', ...
    %          [(w - w_pri)/2,(h - h_pri)/2,w_pri,h_pri], ...
    %          'LineWidth',0.2,'LineStyle','--');
    plane = [(w - w_pri)/2,(h - h_pri)/2,w_pri,h_pri];
    
    %Setup for number of inliers
    inliers = 0;
    P_hat = [];
    R_hat = [];
    c_hat = [];
    
    %Search for relative points around the staring point that was choosen
    %and grow
    dj = d0;
    rx = round(w_pri/2);
    ry = round(h_pri/2);
    for j=3:l
        e = 0;
        %dp = dj;
        
        %Choose a relative random point and convert it into 3D, making sure
        %it is within bounds
        dj = d0 + [randi([-rx,rx]),randi([-ry,ry])];
        %display(dj);
        if (dj(1)>w) || (dj(1)<=0)
            %display('skip dj');
            %dj = dp;
            continue;
        end
        if (dj(2)>h) || (dj(2)<=0)
            %display('skip dj');
            %dj = dp;
            continue;
        end
        
        pj = [I(dj(2),dj(1)) * ((dj(1)/(w-1))-0.5) * tand(fh/2) ...
              I(dj(2),dj(1)) * ((dj(2)/(h-1))-0.5) * tand(fv/2) ...
              I(dj(2),dj(1))];
        pj = (pj .* [(w*wr)/(3*hr),(h*wr)/(3*hr),1]) + [w/2,h/2,0];
        
        if (pj(1)>((w/2)+rx)) || (pj(1)<=((w/2)-rx))
            %display('skip pj');
            %dj = dp;
            continue;
        end
        if (pj(2)>((h/2)+ry)) || (pj(2)<=((h/2)-ry))
            %display('skip pj');
            %dj = dp;
            continue;
        end
        
        %Plane fit error checking to see if the point is an inlier
        e = abs(r * (pj - p0).');
        display(e);
        if e < epsilon
            %if it is, then temporary record the point and its normal
            P_hat = [P_hat;pj];
            R_hat = [R_hat;r];
            inliers = inliers + 1;
        end
        
    end
    
    %Scenario to skip if all chosen points are in the middle
    if e==0
        display('skip e');
        continue;
    end
    
    %Percentage of inliers to be accepted into the main list
    if (inliers > (a_in * l))
        display('+++++++++++++++++++++++++++++++++++++++++++pass');
        P_hat(:,1) = P_hat(:,1);
        P_hat(:,2) = P_hat(:,2);
        P = [P;P_hat];
        R = [R;R_hat];
        
        %Planes are stored, with points gathered and planes drawn onto the
        %depth map
        n = n + inliers;
    else
        %inliers that were not accepted are classed as outliers
        display('--------------------------------------------fail');
        O = [O;P_hat];
    end
    
end

%Draw found frames
if n>0
    xx = P(:,1);
    yy = P(:,2);
    c_hat = convhull(double(xx),double(yy));
    C = [C;c_hat];
    plot(xx,yy,'b.');
    plot(xx(c_hat),yy(c_hat),'r-');
end

hold off;

end
