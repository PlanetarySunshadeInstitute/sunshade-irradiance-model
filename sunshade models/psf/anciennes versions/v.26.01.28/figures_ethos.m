%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	"Shading map"    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


figure;
imagesc(x,y,shading_total); axis equal; set(gca,'YDir','normal');
title({sprintf('Cumulative Penumbral Shading (%%) – %d Objects',num_objs), ...
       sprintf('Avg. sunlight reduction on Earth disk: %.2e %%',avg_reduction)});
xlabel('X (km)'); ylabel('Y (km)');
colormap(jet); colorbar; clim([0,max(shading_total(:))]);
hold on; theta = linspace(0,2*pi,360);
plot(EarthRadius*cos(theta),EarthRadius*sin(theta),'k','LineWidth',2); hold off;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	"Placement plot(s)"    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


figure;

subplot(1,2,1)
scatter(D_mill_km,offset_km,50,'filled');
xlabel('Distance (10^6 km)'); ylabel('Δy (km)');
title('Offset vs Distance'); grid on;

subplot(1,2,2)
%   X-axis : distance,  Y-axis : Δx (along-track),  Z-axis : Δy (cross-track ↑)
scatter3(D_mill_km, track_km, offset_km, 50, offset_km, 'filled');
xlabel('Distance (10^6 km)');
ylabel('Δx (km)');
zlabel('Δy (km)');
title('3-D Distribution (Distance, Δx, Δy)');
grid on; view(30,25);     