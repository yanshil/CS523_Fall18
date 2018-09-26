d = [-1 1];
[x,y,z]= meshgrid(d,d,d);
points = [x(:),y(:),z(:)];

r = points';
r0 = points';

Ibody = 0;
for i = 1:8
    Ibody = Ibody + 1/8 *(r0(:,i)' * r0(:,i) * eye(3) - r0(:,i) * r0(:,i)');
end

Ibodyinv = inv(Ibody);

%% Rigid Body

mass = 1;
v = [0,0,9.8]';
omega = [0.05,0.02,0.01]';

%% Initialization
x = [0,0,0]';
F = [0,0,-9.8]' * mass;
q = [0,0,0,1];
q = quatnormalize(q);
P = v * mass;
R = quat2rotm(q);
I = R * Ibody * R';
L = I * omega;

torque = getTorque(r,x,F);

faces = [1 2 6 5;
    1 5 7 3;
    5 6 8 7;
    7 3 4 8;
    4 8 6 2;
    1 3 4 2];

%% Iterate

h = 0.008; % timestep
frame = 0;
while x >= -1
    
    frame = frame + 1;
    
    v = P / mass;
    x = x + v * h;

    q  = q + 1/2 * quatmultiply([0;omega]', q);
    q = quatnormalize(q);

    R = quat2rotm(q);

    P = P + F * h;
    
    I = R * Ibody * R'; % Keeps the same
    
    L = L + torque * h; % torque = 0
    
    omega = inv(I) * L; % Keeps the same
    
    for j = 1:8
        r(:,j) = R * r0(:,j) + x;
    end
    
    torque = getTorque(r,x,F);
    
    newplot
    patch('Faces',faces,'Vertices',r', 'FaceVertexCData',hsv(8),'FaceColor','interp')
    title(sprintf('v = [%.2f,%.2f,%.2f], omega = [%.2f,%.2f,%.2f]', v(1),v(2),v(3),omega(1),omega(2),omega(3)))
    view(3)
    xlim([-4 4])
    ylim([-4 4])
    zlim([-2 8])
    axis vis3d
    alpha(0.3) 
    
    saveas(gcf,sprintf('./png/%04d.png', frame))
    
 end





