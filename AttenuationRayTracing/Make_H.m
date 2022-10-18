function SM=Make_H(scan,start_x, start_y, Nx, Ny, FOVcenter, FOVradius,Pix)
dbstop if error
H = sparse( 0, 0 ); % 按角度顺序横着排列
% 反投影backprojection（这里采用pixel-driven）
grid.x=linspace(start_x,Nx+start_x,Nx+1); % Nx+1表示，这些都是点数，Nx+1个点
grid.y=linspace(start_y,Ny+start_y,Ny+1);
scan.FanAngle=scan.FanAngle*pi/180;
dgamma=scan.FanAngle/scan.Ndet;
gamma=-(scan.FanAngle-dgamma)/2:dgamma:(scan.FanAngle-dgamma)/2;
if length(scan.theta) ==0
    theta = (180:-1:-179)*pi/180;
    cosTheta = cos(theta);
    sinTheta = sin(theta);
else
    theta=scan.theta*pi/180;
    cosTheta = cos(theta);
    sinTheta = sin(theta);
end
srcPosition0=[-scan.os;0];
detPosition0=[(scan.os+scan.od).*cos(gamma)-scan.os;(scan.os+scan.od).*sin(gamma)]; % 都是以网格点为基准表达长度的，但是计算的却是单位中心的坐标
%deltaBeta = abs(theta(2) - theta(1));

    
for Angle = 1 : scan.Nangle % column
    disp(Angle)
    H_temp = zeros(scan.Ndet, Nx*Ny);
    szb.x=srcPosition0(1)*cosTheta(Angle)+srcPosition0(2)*sinTheta(Angle);%源坐标
    szb.y=-srcPosition0(1)*sinTheta(Angle)+srcPosition0(2)*cosTheta(Angle);
    for k=1:scan.Ndet%探测器循环
        xx0 = detPosition0(1,k);
        yy0 = detPosition0(2,k);
        dzb.x= xx0*cosTheta(Angle)+yy0*sinTheta(Angle);%探测器坐标
        dzb.y= yy0*cosTheta(Angle)-xx0*sinTheta(Angle);

        sects = get_sect(szb,dzb,grid);%获得交点
        if length(sects)<2
            H_temp(k,:)=0;
        else
            SM = get_li_xyx(sects,grid,FOVcenter, FOVradius);
            H_temp(k,SM(1,:))=SM(2,:);
        end
        
    end          
    H = [H;H_temp];
end
Hsys=sparse(H)*Pix;
save Hsys Hsys
end