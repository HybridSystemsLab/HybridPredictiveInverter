function [value, isterminal, direction ] = odeJumpEvent( ~, y )

value = single(~D_inverter(y));
isterminal = 1;
direction = -1;
end

