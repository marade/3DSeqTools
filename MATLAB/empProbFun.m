function p = empProbFun( R, model)

p = zeros(size(R));

p(R==0) = model.p0;

ranger = and(R>0,R<model.R0);

p(ranger) =  (model.k1^-1+model.k2^-1)^-1*(1-model.p0)*exp( (R(ranger)-model.R0)*model.k1);

ranger = R>model.R0;

p(ranger) =  (model.k1^-1+model.k2^-1)^-1*(1-model.p0)*exp( -(R(ranger)-model.R0)*model.k2);

end

