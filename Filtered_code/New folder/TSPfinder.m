
function [Lstridetime,Lstridelength,Lsteptime,Lsteplength, Rstridetime,Rstridelength,Rsteptime,Rsteplength, Rspeed,Lspeed]=TSPfinder(data,data2)
time=(1:size(data,1))/100;
N(1).trdata=zeros(size(data));
N(1).trdata(1300:2390,:)=data(1300:2390,:);
N(2).trdata=zeros(size(data2));
N(2).trdata(1300:2390,:)=data2(1300:2390,:);
for i=1:2
    for j=1:3
N(i).freeinGCS(:,j)=N(i).trdata(:,j);
    end
end
for i=1:2
for j=1:3
N(i).Vel(:,j)=cumtrapz(time,N(i).freeinGCS(:,j));
end
end

for i=1:2
for j=1:3
N(i).Slope2(j)=(N(i).Vel(2390,j)-N(i).Vel(1300,j))/1090;
end
end
for i=1:2
    for j=1:3
        for k=1:size(data,1)
        if k>=1300 && k<=2390
N(i).ModVel(k,j)= N(i).Vel(k,j)-N(i).Slope2(j)*(k-1300)-N(i).Vel(1300,j);
        else
N(i).ModVel(k,j)= 0;            
        end
        end
    end
end


for i=1:2
for j=1:3
N(i).Position(:,j)=cumtrapz(time,N(i).ModVel(:,j));
end
end
% 

%% Temporal Spatial Parameters
Vel=N(1).ModVel;
Vel2=N(2).ModVel;
Position=N(1).Position;
Position2=N(2).Position;
[PKSR,LOCSR]= findpeaks(-Vel(:,3),'MinPeakDistance',100);
[PKSL,LOCSL]= findpeaks(-Vel2(:,3),'MinPeakDistance',100);

for k=1:size(LOCSL)-1
    
TS.Lstridetime(k)=(LOCSL(k+1)-LOCSL(k))/100;
TS.Lstridelength(k)=Position2(LOCSL(k+1),1)-Position2(LOCSL(k),1);
TS.Lstridespeed(k)=TS.Lstridelength(k)/TS.Lstridetime(k);

end
for k=1:size(LOCSR)-1
TS.Rstridetime(k)=(LOCSR(k+1)-LOCSR(k))/100;
TS.Rstridelength(k)=Position(LOCSR(k+1),1)-Position(LOCSR(k),1);
TS.Rstridespeed(k)=TS.Rstridelength(k)/TS.Rstridetime(k);
end
for k=1:size(LOCSR)
    for z=1:size(LOCSL)-1
if LOCSR(k)> LOCSL(z) && LOCSR(k)< LOCSL(z+1)
TS.Lsteptime(k)=(LOCSL(z+1)-LOCSR(k))/100;
TS.Rsteptime(k)=(-LOCSL(z)+LOCSR(k))/100;
TS.Lsteplength(k)=Position2(LOCSL(z+1),1)-Position(LOCSR(k),1);
TS.Rsteplength(k)=-Position2(LOCSL(z),1)+Position(LOCSR(k),1);
TS.Lstepspeed(k)=TS.Lsteplength(k)/TS.Lsteptime(k);
TS.Rstepspeed(k)=TS.Rsteplength(k)/TS.Rsteptime(k);
end
    end
end

[PKSR2,LOCSR2]= findpeaks(Vel(:,3),'MinPeakDistance',100);
[PKSL2,LOCSL2]= findpeaks(Vel2(:,3),'MinPeakDistance',100);
for k=2:size(LOCSL2)
  for z=1:100
   if Vel(LOCSL2(k)-z,3)-Vel(LOCSL2(k)-z-1,3)<0
   break
   end
TS.Lskatetime(k-1)=(LOCSL(k)-z-LOCSL(k-1))/100;   
  end
end

for k=2:size(LOCSR2)-1
  for z=1:100
   if Vel(LOCSR2(k)-z,3)-Vel(LOCSR2(k)-z-1,3)<0
   break
   end
TS.Rskatetime(k-1)=(LOCSR(k)-z-LOCSR(k-1))/100;   
  end
end
Lstridetime=TS.Lstridetime;
Lstridelength=TS.Lstridelength;
Rstridetime=TS.Rstridetime;
Rstridelength=TS.Rstridelength;
Lsteptime=TS.Lsteptime;
Lsteplength=TS.Lsteplength;
Rsteptime=TS.Rsteptime;
Rsteplength=TS.Rsteplength;
Rspeed=TS.Rstridespeed;
Lspeed=TS.Lstridespeed;
end

