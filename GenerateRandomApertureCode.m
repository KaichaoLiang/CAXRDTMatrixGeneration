ApertureRangeY = 90;
ApertureRangeZ = 90;
codeSize = 1;
codeNum = 1800;
GridPosY=[];
GridPosZ=[];
for code =1:1800
    if(numel(GridPosY)==0)
        GridPosY(end+1)=(rand(1)-0.5)*ApertureRangeY;
        GridPosZ(end+1)=(rand(1)-0.5)*ApertureRangeZ;
        continue;
    end
    posY=GridPosY(end);
    posZ=GridPosZ(end);
    while(CheckOverlap(posY,posZ, GridPosY,GridPosZ,codeSize))
        posY=(rand(1)-0.5)*ApertureRangeY;
        posZ=(rand(1)-0.5)*ApertureRangeZ;
    end
    GridPosY(end+1)=posY;
    GridPosZ(end+1)=posZ;
end
save GridPosYUn GridPosY
save GridPosZUn GridPosZ

function res = CheckOverlap(posY,posZ, GridPosY,GridPosZ,codeSize)
res=false;
for index = 1:numel(GridPosY)
    cmpY=GridPosY(index);
    cmpZ=GridPosZ(index);
    if(abs(posY-cmpY)<=codeSize && abs(posZ-cmpZ)<=codeSize)
        res=true;
        break;
    end
end
end