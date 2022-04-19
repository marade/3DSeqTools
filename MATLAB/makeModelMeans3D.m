function data = makeModelMeans3D( data )

xx = data.xdec;

%% Make the model means
u0 = data.T0.u0;

if isfield( data, 'Tvec' )
    du = u0;
    du_h = u0;
    for ii = 1:numel( data.Tvec )
        That = data.Tvec{ii};
        
        pval = That.pval;
        lp = makelp3D( That.lam, data.datap );
        
        umax = That.umax;
        xmax = That.xmax;
        
        % If the peak is significant, include it
        if pval < data.p_sig;
            
            du_ = makedu(xx, That);
            du_h_ = makedu(data.x, That);
            
            du   = du +   du_;
            du_h = du_h + du_h_;
        end
    end
end

data.u_h = du_h;
data.u_l = du;

end