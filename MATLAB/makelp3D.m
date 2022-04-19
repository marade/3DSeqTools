function [lp] = makelp3D( lam, datap )


lp = -(lam-datap.ull)/datap.sll/log(10);


end