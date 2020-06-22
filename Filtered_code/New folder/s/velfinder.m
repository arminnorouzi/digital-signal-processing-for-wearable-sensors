function [time,vel]=velfinder(data,data2,n)
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
vel=[N(1).ModVel(:,n),N(2).ModVel(:,n)];
end