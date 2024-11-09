clc;clear;
%所有物理量都是国际单位
%mdh参数，不是dh参数
dh(1,:) = [    0, -0.285,      0,    pi, 0, 0];
dh(2,:) = [    0,      0,      0,  pi/2, 0, pi/2];
dh(3,:) = [    0,      0, -0.398,    pi, 0, 0];
dh(4,:) = [    0, -0.194,  -0.36,    pi, 0, -pi/2];
dh(5,:) = [    0,  -0.17,      0, -pi/2, 0, 0];
dh(6,:) = [    0,  -0.16,      0,  pi/2, 0, 0];

% 记第i个连杆由MDH建模方法决定的随体坐标系为{i}
% 记第i个连杆的质心坐标系为{ic}
% {ic}与{i}方向相同
Inertia = [[0.017, 0.466, 0.017, 0, -0.001, 0],[0, 0.043, 0.005], 5.675;
           [0.014, 0.124, 0.122, 0, 0, -0.009],[-0.298, 0, -0.126], 5.466;
           [0.01, 0.081, 0.076, 0, 0, -0.01],[-0.274, 0, 0.061],4.294;
           [0.005, 0.003, 0.005, 0, 0, 0],[0, -0.055, 0.001],2.732;
           [0.005, 0.003, 0.005, 0, 0, 0],[0, 0.055, 0.001],2.732;
           [0.0001, 0.0001, 0.0001, 0, 0, 0],[0, 0, 0.017], 0.28];      
%Inertia(i,:) = [[Ixx,Iyy,Izz,Ixy,Iyz,Ixz],[rx,ry,rz],m]
%[Ixx,Iyy,Izz,Ixy,Iyz,Ixz]:第i个连杆在{i}下的惯量
%[rx,ry,rz]：第i个连杆质心在{i}下的坐标
%m：连杆质量       
for i = 1:6
L(i) = Link(dh(i,:),'modified');
I = zeros(3,3); I(1,2:3) = Inertia(i,[4,6]);I(2,3) = Inertia(i,5);
I = I+I'; I = I + diag(Inertia(i,1:3)) ;
L(i).I = I;
L(i).r = Inertia(i,7:9); %L(i).r=[rx,ry,rz]
L(i).m = Inertia(i,10);
L(i).I =L(i).I + L(i).m*skew(L(i).r)^2; %平行轴定理,L(i).I是第i个连杆在{ic}下的惯量
L(i).Jm = 0;
end
robot = SerialLink(L,'name', 'qinzhe');

syms t;
time_ub=1.5;
theta0=[0.5 0.7 0.1 -0.3 0.2 0.12];
theta_test1=[sin(t) sin(t^2) 1-cos(t) sin(t^3) 2-2*cos(t) 1-cos(t^2)];

time=0;
theta1=[];
theta2=[];
theta3=[];
theta4=[];
theta5=[];
theta6=[];

theta1_ik=[];
theta2_ik=[];
theta3_ik=[];
theta4_ik=[];
theta5_ik=[];
theta6_ik=[];


for time=0:0.01:time_ub
    theta_test=[sin(time) sin(time^2) 1-cos(time) sin(time^3) 2-2*cos(time) 1-cos(time^2)];
    
    t1=theta_test(1)+theta0(1);
    t2=theta_test(2)+theta0(2);
    t3=theta_test(3)+theta0(3);
    t4=theta_test(4)+theta0(4);
    t5=theta_test(5)+theta0(5);
    t6=theta_test(6)+theta0(6);
    T01=[1 0 0 0;0 -1 0 0;0 0 -1 0;0 0 0 1]*[cos(t1) -sin(t1) 0 0;sin(t1) cos(t1) 0 0;0 0 1 -0.285;0 0 0 1];
    T12=[1 0 0 0;0 0 -1 0;0 1 0 0;0 0 0 1]*[-sin(t2) -cos(t2) 0 0;cos(t2) -sin(t2) 0 0;0 0 1 0;0 0 0 1];
    T23=[1 0 0 -0.398;0 -1 0 0;0 0 -1 0;0 0 0 1]*[cos(t3) -sin(t3) 0 0;sin(t3) cos(t3) 0 0;0 0 1 0;0 0 0 1];
    T34=[1 0 0 -0.36;0 -1 0 0;0 0 -1 0;0 0 0 1]*[sin(t4) cos(t4) 0 0;-cos(t4) sin(t4) 0 0;0 0 1 -0.194;0 0 0 1];
    T45=[1 0 0 0;0 0 1 0;0 -1 0 0;0 0 0 1]*[cos(t5) -sin(t5) 0 0;sin(t5) cos(t5) 0 0;0 0 1 -0.17;0 0 0 1];
    T56=[1 0 0 0;0 0 -1 0;0 1 0 0;0 0 0 1]*[cos(t6) -sin(t6) 0 0;sin(t6) cos(t6) 0 0;0 0 1 -0.16;0 0 0 1];
    
    T_test=T01*T12*T23*T34*T45*T56;
    t_test=ARM_ik(theta0,T_test);
    
    
    % theta_ik=robot.ikine(T_test);
    % theta1_ik(end+1)=theta_ik(1)+theta0(1);
    % theta2_ik(end+1)=theta_ik(2)+theta0(2);
    % theta3_ik(end+1)=theta_ik(3)+theta0(3);
    % theta4_ik(end+1)=theta_ik(4)+theta0(4);
    % theta5_ik(end+1)=theta_ik(5)+theta0(5);
    % theta6_ik(end+1)=theta_ik(6)+theta0(6);


    theta1_ik(end+1)=theta_test(1);
    theta2_ik(end+1)=theta_test(2);
    theta3_ik(end+1)=theta_test(3);
    theta4_ik(end+1)=theta_test(4);
    theta5_ik(end+1)=theta_test(5);
    theta6_ik(end+1)=theta_test(6);

    t_test=t_test(1,:);
    robot.fkine(t_test+theta0);
    theta1(end+1)=t_test(1,1);
    theta2(end+1)=t_test(1,2);
    theta3(end+1)=t_test(1,3);
    theta4(end+1)=t_test(1,4);
    theta5(end+1)=t_test(1,5);
    theta6(end+1)=t_test(1,6);
end

x=0:0.01:time_ub;
figure;
subplot(2,3,1);
plot(x,theta1,'-r',LineWidth=4);
hold on
plot(x,theta1_ik,'-b',LineWidth=2);
hold off
legend('解析解','生成轨迹');
xlabel('Time/s');
ylabel('θ1/rad');

subplot(2,3,2);
plot(x,theta2,'-r',LineWidth=4);
hold on
plot(x,theta2_ik,'-b',LineWidth=2);
hold off
legend('解析解','生成轨迹');
xlabel('Time/s');
ylabel('θ2/rad');

subplot(2,3,3);
plot(x,theta3,'-r',LineWidth=4);
hold on
plot(x,theta3_ik,'-b',LineWidth=2);
hold off
legend('解析解','生成轨迹');
xlabel('Time/s');
ylabel('θ3/rad');

subplot(2,3,4);
plot(x,theta4,'-r',LineWidth=4);
hold on
plot(x,theta4_ik,'-b',LineWidth=2);
hold off
legend('解析解','生成轨迹');
xlabel('Time/s');
ylabel('θ4/rad');

subplot(2,3,5);
plot(x,theta5,'-r',LineWidth=4);
hold on
plot(x,theta5_ik,'-b',LineWidth=2);
hold off
legend('解析解','生成轨迹');
xlabel('Time/s');
ylabel('θ5/rad');

subplot(2,3,6);
plot(x,theta6,'-r',LineWidth=4);
hold on
plot(x,theta6_ik,'-b',LineWidth=2);
hold off
legend('解析解','生成轨迹');
xlabel('Time/s');
ylabel('θ6/rad');


function theta = ARM_ik(theta0,T)
    theta=[];
    P06=T(1:3,4);
    d4=0.194;
    d6=0.16;
    a3=0.398;
    a4=0.36;
    
    t1=zeros(2,1);
    t2=zeros(4,1);
    t3=zeros(4,1);
    t4=zeros(4,1);
    t5=zeros(4,1);
    t6=zeros(4,1);

    %求解θ1
    P05=T*[0;0;d6;1];
    phi1=-atan2(P05(2),P05(1));
    t1(1)=phi1+acos(d4/sqrt(P05(1)^2+P05(2)^2))-pi/2;
    % t1(2)=phi1+acos(d4/sqrt(P05(1)^2+P05(2)^2))-pi/2-pi;


    %求解θ5
    T01_1=[1 0 0 0;0 -1 0 0;0 0 -1 0;0 0 0 1]*[cos(t1(1)) -sin(t1(1)) 0 0;sin(t1(1)) cos(t1(1)) 0 0;0 0 1 -0.285;0 0 0 1];
    R01=T01_1(1:3,1:3);
    P16=transpose(R01)*P06;
    t5(1)=acos((P16(2)-d4)/d6);
    t5(2)=-acos((P16(2)-d4)/d6);


    
    %求解θ6
    T16=inv(T01_1)*T;
    R16=T16(1:3,1:3);
    R61=transpose(R16);
    Y61=R61(:,2);
    t6(1)=atan2((Y61(2))/(-sin(t5(1))),Y61(1)/sin(t5(1)));
    t6(2)=atan2((Y61(2))/(-sin(t5(2))),Y61(1)/sin(t5(2)));



    %求解θ3
    T45=[1 0 0 0;0 0 1 0;0 -1 0 0;0 0 0 1]*[cos(t5(1)) -sin(t5(1)) 0 0;sin(t5(1)) cos(t5(1)) 0 0;0 0 1 -0.17;0 0 0 1];
    T56=[1 0 0 0;0 0 -1 0;0 1 0 0;0 0 0 1]*[cos(t6(1)) -sin(t6(1)) 0 0;sin(t6(1)) cos(t6(1)) 0 0;0 0 1 -0.16;0 0 0 1];
    T14_1=inv(T01_1)*T*inv(T56)*inv(T45);
    t3(1)=acos((T14_1(1,4)^2+T14_1(3,4)^2-a3^2-a4^2)/(2*a3*a4));

    %求解θ2
    phi1=atan2(T14_1(1,4),-T14_1(3,4));
    phi2=asin((a4*sin(t3(1)))/sqrt(T14_1(1,4)^2+T14_1(3,4)^2));
    t2(1)=phi1+phi2;

    %求解θ4
    T12=[1 0 0 0;0 0 -1 0;0 1 0 0;0 0 0 1]*[-sin(t2(1)) -cos(t2(1)) 0 0;cos(t2(1)) -sin(t2(1)) 0 0;0 0 1 0;0 0 0 1];
    T23=[1 0 0 -0.398;0 -1 0 0;0 0 -1 0;0 0 0 1]*[cos(t3(1)) -sin(t3(1)) 0 0;sin(t3(1)) cos(t3(1)) 0 0;0 0 1 0;0 0 0 1];
    T34=inv(T23)*inv(T12)*T14_1;
    t4(1)=atan2(T34(1,1),T34(2,1));

    theta(end+1,:)=[t1(1)-theta0(1) t2(1)-theta0(2) t3(1)-theta0(3) t4(1)-theta0(4) t5(1)-theta0(5) t6(1)-theta0(6)];




    % T45=[1 0 0 0;0 0 1 0;0 -1 0 0;0 0 0 1]*[cos(t5(2)) -sin(t5(2)) 0 0;sin(t5(2)) cos(t5(2)) 0 0;0 0 1 -0.17;0 0 0 1];
    % T56=[1 0 0 0;0 0 -1 0;0 1 0 0;0 0 0 1]*[cos(t6(2)) -sin(t6(2)) 0 0;sin(t6(2)) cos(t6(2)) 0 0;0 0 1 -0.16;0 0 0 1];
    % T14_2=inv(T01_1)*T*inv(T56)*inv(T45);
    % if (T14_2(1,4)^2+T14_2(3,4)^2-a3^2-a4^2)/(2*a3*a4)<1
    %     t3(2)=acos((T14_2(1,4)^2+T14_2(3,4)^2-a3^2-a4^2)/(2*a3*a4));
    % 
    %     phi1=atan2(T14_2(1,4),-T14_2(3,4));
    %     phi2=asin((a4*sin(t3(2)))/sqrt(T14_2(1,4)^2+T14_2(3,4)^2));
    %     t2(2)=phi1+phi2;
    % 
    %     T12=[1 0 0 0;0 0 -1 0;0 1 0 0;0 0 0 1]*[-sin(t2(2)) -cos(t2(2)) 0 0;cos(t2(2)) -sin(t2(2)) 0 0;0 0 1 0;0 0 0 1];
    %     T23=[1 0 0 -0.398;0 -1 0 0;0 0 -1 0;0 0 0 1]*[cos(t3(2)) -sin(t3(2)) 0 0;sin(t3(2)) cos(t3(2)) 0 0;0 0 1 0;0 0 0 1];
    %     T34=inv(T23)*inv(T12)*T14_2;
    %     t4(2)=atan2(T34(1,1),T34(2,1));
    % 
    %     theta(end+1)=[t1(1)-theta0(1) t2(2)-theta0(2) t3(2)-theta0(3) t4(2)-theta0(4) t5(2)-theta0(5) t6(2)-theta0(6)];
    % end

    [row,col]=size(theta);

    for i=1:row
        for j=1:col
            while theta(i,j)<-2*pi
                theta(i,j)=theta(i,j)+2*pi;
            end
            while theta(i,j)>2*pi
                theta(i,j)=theta(i,j)-2*pi;
            end
        end
    end

end