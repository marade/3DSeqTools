%% genF2( N, model )
% this function generated simulated null hypothesis for computing the test
% statistic using monte carlo.
function t = genF2( N, model )


p  = model.p0;
b1 = model.k1;
b2 = model.k2;
x0 = model.R0;


F = rand( [1,N] );

t = zeros([1,N]);
R = exprnd( 1, [1,N] );
t2 = -R/b1+x0;
t3 = R/b2+x0;

p1 = p;
p2 = p+(1-p)*b2/(b1+b2);

%flag0 = (F<p1);
flag2 = and(F>p1,F<p2);
flag3 = F>p2;


t(flag2) = t2(flag2);
t(flag3) = t3(flag3);

end