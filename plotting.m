clc
clear all
close all
[A.num,A.txt,A.raw]=xlsread('Trial1.xlsx');
Mtrajectories = A.num(5:end,:);
ylabeling={'x [mm]','y [mm]','z [mm]'};
Mtrajectories(:,1)=Mtrajectories(:,1)/100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('units','normalized','outerposition',[0 0 1 1])
sgtitle ('Right Foot');
for i=1:4
for j=1:3
    subplot(3,1,j)
    plot(Mtrajectories(:,1),Mtrajectories(:,3*i+j-1)); hold on;
    xlabel('Time [s]');
    ylabel(ylabeling(j));

end

end
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('units','normalized','outerposition',[0 0 1 1])
sgtitle ('Right Shank');
for i=5:9
for j=1:3
    subplot(3,1,j)
    plot(Mtrajectories(:,1),Mtrajectories(:,3*i+j-1)); hold on;
    xlabel('Time [s]');
    ylabel(ylabeling(j));
end
end
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('units','normalized','outerposition',[0 0 1 1])
sgtitle ('Left Shank');
for i=10:14
for j=1:3
    subplot(3,1,j)
    plot(Mtrajectories(:,1),Mtrajectories(:,3*i+j-1)); hold on;
    xlabel('Time [s]');
    ylabel(ylabeling(j));
end
end
hold off;
figure('units','normalized','outerposition',[0 0 1 1])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sgtitle ('Left Foot');
for i=15:18
for j=1:3
    
    subplot(3,1,j)
    plot(Mtrajectories(:,1),Mtrajectories(:,3*i+j-1)); hold on;
    xlabel('Time [s]');
    ylabel(ylabeling(j));
end
end
hold off;