function [tab, tab2] = makeTab( data )

width = drill(data.Tvec, '.L0')';
height = drill(data.Tvec, '.I0')';
position = drill(data.Tvec, '.x0')';
log10pval = drill(data.Tvec, '.lpval')';
position_err = drill(data.Tvec, '.fisherI(1,1)')';
a0 = drill(data.Tvec, '.a0')';

tab = table( log10pval, position, position_err, height, width, a0 );

pos = data.xdec';
model     = data.u_l';

tab2 = table( pos, model );

end