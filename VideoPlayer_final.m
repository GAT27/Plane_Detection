function [P,R,C,O] = VideoPlayer_final(feed,n_max,k_max,l,eta,S,epsilon,a_in)
%List of 3D points
P = [];
%List of normals
R = [];
%List of polygons
C = [];
%List of outliers
O = [];
%Plane checker
pchk = 0;

%Video data
frames = feed.NumberOfFrames;
vid_height = feed.Height;
vid_width = feed.Width;

%Crop video
h1 = (vid_height - 480) / 2;
h2 = vid_height - (480/2);
w1 = (vid_width - 640) / 2;
w2 = vid_width - (640/2);

%Skip frames for specific testing cases 
for k=300:frames-400
    I = read(feed,k);
    I = I(h1:h2,w1:w2,:);
    
    %Filter a frame, display planes, and store data with seperation buffer
    [p,r,c,o,plane] = PlaneFiltering_final(I,n_max,k_max,l,eta,S,epsilon,a_in);
    P = [P;p];
    R = [R;r];
    C = [C;c];
    O = [O;o];
    
    %If plane was found, then add to checker
    if isempty(p)==0
        pchk = pchk + 1;
        if pchk > 5
            pchk = 5;
        end
    else
        %otherwise, subtract from checker
        pchk = pchk - 2;
        if pchk < 0
            pchk = 0;
        end
    end
    
    %If plane is consistenly found, then project
    if pchk > 2
        rectangle('Position', ...
                  [plane(1),plane(2),plane(3),plane(4)], ...
                  'FaceColor','g');
    end
    
end

end