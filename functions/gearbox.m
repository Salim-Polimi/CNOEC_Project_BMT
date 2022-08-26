function [tau] = gearbox(vref)

gear_vector = [2.14, 2.68, 3.84, 5.86, 9.97];

if vref < 10
    tau = gear_vector(1);
elseif vref < 20
    tau = gear_vector(2);
elseif vref < 40
    tau = gear_vector(3);
elseif vref < 70
    tau = gear_vector(4);
else
    tau = gear_vector(5);
end


