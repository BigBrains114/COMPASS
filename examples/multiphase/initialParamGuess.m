function params = initialParamGuess(obj)

t = obj.time;
kp = [1; obj.auxdata.kp(:)];
np = obj.np;

params = zeros(np,1);
for k = 1:np
    params(k) = t(kp(k)+1) - t(kp(k));
end