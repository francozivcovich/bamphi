function [ u, info ] = cuschss68_bamphi( u, k, t, A, g, opts, info )

  mu = g( NaN, NaN );

  % persistent store_info

  % u = exp( ( - 1i * mu * k / 2 ) * ( u .* conj( u ) ) ) .* u;
  % if length( store_info ) > 0,
  %   info = store_info{ 1 };
  % end
  % [ u, info ] = bamphi( 1i * k, A, [], u, opts, info );
  % store_info{ 1 } = info;
  % u = exp( ( - 1i * mu * k / 2 ) * ( u .* conj( u ) ) ) .* u;


  % [ u, info ] = bamphi( 1i * k / 2, A, [], u, opts, info );
  % u = exp( ( - 1i * mu * k ) * ( u .* conj( u ) ) ) .* u;
  % [ u, info ] = bamphi( 1i * k / 2, A, [], u, opts, info );

  u = exp( ( - 1i * mu * k / 2 ) * ( u .* conj( u ) ) ) .* u;
  [ u, info ] = bamphi( 1i * k, A, [], u, opts, info );
  u = exp( ( - 1i * mu * k / 2 ) * ( u .* conj( u ) ) ) .* u;

end
