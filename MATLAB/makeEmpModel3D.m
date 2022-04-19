function model_ = makeEmpModel3D( Rval )




kk = 1/std(Rval(Rval>0));


% guess where the transition point is
[yy,rr] = hist( Rval, 50 );
[~,ind] = max(yy(2:end));

yy = yy/sum(yy)/(rr(2)-rr(1));

R0_ = rr(ind+1);

T0 = [.2,kk*.9,kk,R0_];

hot_ind = 1:3;
T_ = fminsearch( @fitter, T0(hot_ind) );
T0(hot_ind) = T_;

hot_ind = 1:4;
T_ = fminsearch( @fitter, T0(hot_ind) );

    function Y = fitter( T )
        if T(1)>1;
            T(1) = 1;
        elseif T(1)<0;
            T(1) = 0;
        end
        
        if T(2)<0;
            T(2) = 0;
        end
        
        if T(3)<0;
            T(3) = 0;
        end
        
        if numel(hot_ind) > 3
            if T(4)<0;
                T(4) = 0;
            end
            model.R0 = T(4);
        else
            model.R0 = T0(4);
        end
        
        
        model.p0 = T(1);
        model.k1  = T(2);
        model.k2  = T(3);
        
        %T-T0
        Y = -sum(log( empProbFun( Rval, model )));
        
        figure(5);
        clf;
        semilogy(  rr, yy );
        hold on;
        semilogy(  rr, empProbFun( rr, model ) );
        drawnow;
        
        
    end

T0
T_

model_.p0  = T_(1);
model_.k1  = T_(2);
model_.R0  = T_(4);
model_.k2  = T_(3);

end


