function [data] = decimateData3D( data )


Lg = data.Lg;
dL = data.dL;
x  = data.x;
DD = data.D;

nbin = ceil(Lg/dL) 

xdec = nan( [1, nbin] );
Ddec = nan( [1, nbin] );
wdec = nan( [1, nbin] );

h = waitbar( 0, 'Decimating data' );

for ii = 1:nbin;
    
    waitbar(ii/nbin, h)
    
    xmin = (ii-1)*dL;
    xmax = (ii)*dL;

    ind = and( x>xmin, x<xmax);
       
    xdec(ii) = sum(DD(ind).*x(ind))/sum(DD(ind));
    
    if isnan(xdec(ii))
        xdec(ii) = mean(x(ind));
    end
    
    Ddec(ii) = mean(DD(ind));
    wdec(ii) = sum( ind );
end

close(h);

flag = and(~isnan(Ddec),~isnan(xdec));
xdec = xdec(flag);
Ddec = Ddec(flag);
wdec = wdec(flag);


data.xdec = xdec;
data.Ddec = Ddec;
data.wdec = wdec;


% lets glue the mean to zero 
u0 = mean( Ddec );
data.T0.u0 = u0;

% normalize data to same total enzyme activity
s = std( Ddec );
data.T0.s = s;

if ~isfield( data, 'salt_and_pep_flag' )
    data.salt_and_pep_flag = false;
end


if data.salt_and_pep_flag
    %% De salt and pepper
    Ddecp = medfilt1( Ddec, 3);
    
    dd = (Ddec-Ddecp)/s;
    
    flagger = dd > data.s_p_cut;
    
    
    Ddecp(~flagger) = Ddec(~flagger);
    
    data.xdec = xdec;
    data.Ddec = Ddecp;
    data.wdec = wdec;
end

end