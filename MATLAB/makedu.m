%% makedu
% Compute the model mean as a function of position at the peak model
% parameters T
%
% define shape function 
% x = postions
% T = shape params
function du = makedu( x, T )

% extract shape param
xp= T.x0;
I = T.I0;
L = T.L0;
D = T.D0;
a = T.a0;
Lg= T.Lg;

% use periodic BC
Lx = numel( x );
xx = [x-Lg,x,x+Lg];


Lcut = 16;

if 0
%%
du =  I./(1+(abs((xx-xp)/L)).^D).^(a/D);
du(abs((xx-xp)/L)>16)=0;
%%
else
du = I.*(  (1+(abs((xx-xp)/L)).^D).^(-a/D)-(1+Lcut.^D).^(-a/D) )/(1-(1+Lcut.^D).^(-a/D) ) ;
du(du<0) = 0;
end


% reduce back to the right size.
du = du( (1:Lx)+Lx );
end