%%  doGumbelFit
% fit the gumbel parameters
function model = doGumbelFit2( lam )

    lam = sort(lam);

    N = numel( lam );
    p = 1-((1:N)-1)/(N-1);
    
    ind = (N-floor( N/10 )):N;
    
    % make fit guesses
   
    
    
    T0 = [mean(lam),std(lam)];

   Tfit = lsqnonlin( @fitFun, T0 );
   %Tfit = fminsearch( @fitFun, T0 );

    model.u = Tfit(1);
    model.s = Tfit(2);

    function dy = fitFun( T )
        
        u = T(1);
        s = T(2);
        
        theory = 1-exp(-exp( -(lam(ind)-u)/s ));
        dy = (theory-p(ind))/p(ind(1));
        
         figure(5);
         clf;
         
         loglog( lam, p, 'b.-' );
         hold on;
         loglog( lam, 1-exp(-exp( -(lam-u)/s)), 'r.-' );
 drawnow;
        
        
    end

end