function [value,isterminal,direction] = EventFcn_1(~,x)


value(1) = x(5) - 0.5;
direction(1) = 0;
isterminal(1) = 1;