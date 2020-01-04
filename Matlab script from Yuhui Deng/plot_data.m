function f=plot_data
k=10; % scale for the size of arrows
centre = [0,0,0];
two_dim_norm = [0,0,-1]; % norm direction of 2D plot
data = csvread('Tdtomato_Dermal_df_E11_5.csv',1,2);
[r,l]=size(data);
theta = zeros(1,r);
figure(1);
plot(0,0,'r*');
for i = 1:r
    x_start = data(i,4:6)-centre;
    x_end = data(i,1:3)-centre;
    theta(i)=acos((-x_start)*(x_end-x_start)'/(normest(-x_start)*normest(x_end-x_start)));    
    if x_start(3) < 0
       two_dim_norm = [0,0,1];
    end
    data_norm=cross(x_start,x_end);
    cos_angle = two_dim_norm*data_norm'/(normest(two_dim_norm)*normest(data_norm));
    rotation_axis = cross(two_dim_norm,data_norm);
    rotation_axis = rotation_axis/normest(rotation_axis);
    x_start_new = rotate(x_start,rotation_axis,cos_angle);
    x_end_new = rotate(x_end,rotation_axis,cos_angle);
    direction = x_end_new - x_start_new;
    direction = direction / normest(direction);
    hold on;
    quiver(x_start_new(1),x_start_new(2),k*direction(1),k*direction(2),'MaxHeadSize',1);
end
figure(2);
edges=[0 90 180];
histogram(theta*180/pi,edges);
figure(3);
alpha = 0:0.01:2*pi;
plot(cos(alpha),sin(alpha),'k-');
axis([-1.5 1.5 -1.5 1.5]);
text(1.1,0,'0^o');
text(0,1.1,'90^o');
text(-1.2,0,'180^o');
text(0,-1.1,'270^o');
text(1,-1,['n=',num2str(r)]);
for i = 1:r
    hold on;
    quiver(0,0,cos(theta(i)),sin(theta(i)),1,'k');
end
end


function f=rotate(x,axis,angle)
f=(x-(x*axis')*axis)*angle+cross(x,axis)*sqrt(1-angle^2)+(x*axis')*axis;
end



