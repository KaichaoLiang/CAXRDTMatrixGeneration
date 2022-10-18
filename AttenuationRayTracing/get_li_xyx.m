function lid = get_li_xyx(sects,grid)
% sects: ����
% grid: ����
% FOVcenter: 2 dimentional vector, [x, y]
% FOVradious: ɨ������İ뾶
if length(sects)<2
    lid=[];
else
    Nx = length(grid.x)-1;
    Ny = length(grid.y)-1;
    nsect = length(sects);
    midx = (sects(1,1:nsect-1)+sects(1,2:nsect))/2;  %��������
    midy = (sects(2,1:nsect-1)+sects(2,2:nsect))/2;  %��������
    l = sqrt((sects(2,1:nsect-1)-sects(2,2:nsect)).^2+(sects(1,1:nsect-1)-sects(1,2:nsect)).^2); %���н��߳��ȣ���һ������
    pixelSize = grid.x(2)-grid.x(1);%pixel size,���������صĵ�λ�ķ�����pixelSize��Զ��1
    
    % ����FOVԲ������أ������Ӧ��ϵͳ����ϵ��Ϊ0
%     dist2FOVcenter = (midx-FOVcenter(1)).*(midx-FOVcenter(1))+(midy-FOVcenter(2)).*(midy-FOVcenter(2));
%     ii = find( dist2FOVcenter > (FOVradius+pixelSize)^2);  %��һ�����ص�����
%     l(ii)=0;
    
    nx = floor((midx-grid.x(1))/pixelSize)+1;
    ny = floor((midy-grid.y(1))/pixelSize)+1;
    id=(nx-1)*Ny+ny;
    id2 = id(find(l>0));  %ȥ��FOV�������,ʵ���ϣ������������ּ�С���ʣ��ǲ������0��
    l2 = l(find(l>0));    %ȥ��FOV�������
    lid = [id2;l2];
    %lid = [id;l];
    
end