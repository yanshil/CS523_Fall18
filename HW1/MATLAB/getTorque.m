function [torque] = getTorque(r,x,F)
torque = 0;
for i =1:8
    torque = torque + cross((r(:,i) - x), F/8);
end

end

