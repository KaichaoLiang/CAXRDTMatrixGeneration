function lid = get_li_xyx(sects,grid)
% sects: 交点
% grid: 网格
% FOVcenter: 2 dimentional vector, [x, y]
% FOVradious: 扫描区域的半径
if length(sects)<2
    lid=[];
else
    Nx = length(grid.x)-1;
    Ny = length(grid.y)-1;
    nsect = length(sects);
    midx = (sects(1,1:nsect-1)+sects(1,2:nsect))/2;  %物理坐标
    midy = (sects(2,1:nsect-1)+sects(2,2:nsect))/2;  %物理坐标
    l = sqrt((sects(2,1:nsect-1)-sects(2,2:nsect)).^2+(sects(1,1:nsect-1)-sects(1,2:nsect)).^2); %所有交线长度，是一个向量
    pixelSize = grid.x(2)-grid.x(1);%pixel size,按照以像素的单位的方法，pixelSize永远是1
    
    % 对在FOV圆外的像素，令其对应的系统矩阵系数为0
%     dist2FOVcenter = (midx-FOVcenter(1)).*(midx-FOVcenter(1))+(midy-FOVcenter(2)).*(midy-FOVcenter(2));
%     ii = find( dist2FOVcenter > (FOVradius+pixelSize)^2);  %留一个像素的余量
%     l(ii)=0;
    
    nx = floor((midx-grid.x(1))/pixelSize)+1;
    ny = floor((midy-grid.y(1))/pixelSize)+1;
    id=(nx-1)*Ny+ny;
    id2 = id(find(l>0));  %去除FOV外的区域,实际上，除了相切这种极小概率，是不会等于0的
    l2 = l(find(l>0));    %去除FOV外的区域
    lid = [id2;l2];
    %lid = [id;l];
    
end