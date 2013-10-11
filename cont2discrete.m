function [ A, B ] = cont2discrete( F, G, Ts )
%cont2discrete creates the discrete time version of the continous state
%space model, F and G with sampling time Ts

A = expm(F*Ts); % calculate Discrete time free response
B = (F*expm(F*Ts))*G; % calculate discrete time forced response


end

