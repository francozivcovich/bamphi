function [ u, info ] = exprk4s6_bamphi( u, k, t, A, g, opts, info, B )
% exprk4s6_bamphi:
%
% There are 3 main ways to carry on this exponential integration with bamphi:
%
% 1. compute 1 (IOM)-Arnoldi and recycle info for each call and each expint step;
% 2. (default) compute 1 (IOM)-Arnoldi PER CALL inside the FIRST expint step and
%    recycle info;
% 3. compute (IOM)-Arnoldi each time bamphi is called.
%
% You should run 1. if you are really against (IOM)-Arnoldi (mat-vec product
% very cheap but huge matrix sizes make storage a problem) and 3. if (IOM)-Arnoldi
% doesn't bother you a lot (mat-vec really expensive but modest matrix sizes).
%
% (If in 1. you may also consider the low_storage option (see bamphi help).)
%
% The three strategies are coded in the following: only the 2nd one is left
% uncommented as it is the suggested one for this kind of task (mat-vec
% relatively cheap, matrix sizes relatively large).
%
% There are other strategies: e.g. say you have a matrix B similar to A but for
% some reason (see sine-Gordon example in BAMPHI paper) you prefer computing
% Ritz's values with B and run calculations with A. Then you just give B as
% input to exprk4s6_bamphi and any strategy is overruled by this.
%
% Hybrid strategies, such as the 2nd but recomputing (IOM)-Arnoldi every P expint
% steps are possible but have to be handled in concert with the outside of this
% function.
%

  gn = g( t, u );
  Fn = A( u ) + gn;
  z = zeros( length( u ), 1 );
  U = zeros( length( u ), 2 );

  % % 1st STRATEGY:
  %
  % if isempty( info ) && ( nargin == 8 )
  %   aux = bamphi_fov( length( u ), B ); % B is a matrix we rather use for computing A's FoV
  %   info = aux;
  % end
  % [ U(:,1), info ] = bamphi( 1/2 * k, A, [], [ z, Fn ], opts, info );
  % U(:,1) = u + U(:,1);
  %
  % [ U, info ] = bamphi( k * [ 1/3, 1/2 ], A, [], [z,Fn,2/k*(g(t+1/2*k,U(:,1)) - gn)], opts, info );
  % U(:,2) = g( t + 1/2 * k, u + U(:,2) ) - gn;
  % U(:,1) = g( t + 1/3 * k, u + U(:,1) ) - gn;
  %
  % [ U,info ] = bamphi( k * [ 1/3, 5/6 ], A, [], [z,Fn,6/k*(3/2*U(:,1)-2/3*U(:,2)),12/(k^2)*(2*U(:,2)-3*U(:,1))], opts, info );
  % U(:,2) = g( t + 5/6 * k, u + U(:,2) ) - gn;
  % U(:,1) = g( t + 1/3 * k, u + U(:,1) ) - gn;
  %
  % [ z,info ] = bamphi( k, A, [], [z,Fn,2/k*(5/2*U(:,1)-2/5*U(:,2)),4/(k^2)*(6/5*U(:,2)-3*U(:,1))], opts, info );
  % u = u + z;

  % 2nd STRATEGY:
  persistent store_info

  if length( store_info ) > 0
    info = store_info{ 1 };
  else
    if ( nargin == 8 )
      aux = bamphi_fov( length( u ), B ); % B is a matrix we rather use for computing A's FoV
      info = aux;
    end
  end
  [ U(:,1), info ] = bamphi( 1/2 * k, A, [], [ z, Fn ], opts, info );
  store_info{ 1 } = info;
  info = [];
  U(:,1) = u + U(:,1);

  if length( store_info ) > 1
    info = store_info{ 2 };
  else
    if ( nargin == 8 )
      info = aux;
    end
  end
  [ U, info ] = bamphi( k * [ 1/3, 1/2 ], A, [], [z,Fn,2/k*(g(t+1/2*k,U(:,1)) - gn)], opts, info );
  store_info{ 2 } = info;
  info = [];
  U(:,2) = g( t + 1/2 * k, u + U(:,2) ) - gn;
  U(:,1) = g( t + 1/3 * k, u + U(:,1) ) - gn;

  if length( store_info ) > 2
    info = store_info{ 3 };
  else
    if ( nargin == 8 )
      info = aux;
    end
  end
  [ U,info ] = bamphi( k * [ 1/3, 5/6 ], A, [], [z,Fn,6/k*(3/2*U(:,1)-2/3*U(:,2)),12/(k^2)*(2*U(:,2)-3*U(:,1))], opts, info );
  store_info{ 3 } = info;
  info = [];

  U(:,2) = g( t + 5/6 * k, u + U(:,2) ) - gn;
  U(:,1) = g( t + 1/3 * k, u + U(:,1) ) - gn;

  if length( store_info ) > 3
    info = store_info{ 4 };
  else
    if ( nargin == 8 )
      info = aux;
    end
  end
  [ z,info ] = bamphi( k, A, [], [z,Fn,2/k*(5/2*U(:,1)-2/5*U(:,2)),4/(k^2)*(6/5*U(:,2)-3*U(:,1))], opts, info );
  store_info{ 4 } = info;
  info = [];
  u = u + z;

  % % 3rd STRATEGY:
  % persistent store_info suggested_arnoldi_size
  %
  % if ( length( suggested_arnoldi_size ) > 0 )
  %   opts.r_arn = suggested_arnoldi_size( 1 );
  % end
  % if length( store_info ) > 0
  %   info = store_info{ 1 };
  % else
  %   if ( nargin == 8 )
  %     aux = bamphi_fov( length( u ), B ); % B is a matrix we rather use for computing A's FoV
  %     info = aux;
  %   end
  % end
  % [ U(:,1), info ] = bamphi( 1/2 * k, A, [], [ z, Fn ], opts, info );
  % suggested_arnoldi_size( 1 ) = bamphi_suggest_arnoldi_size( opts, info );
  % info = rmfield( info, 'A' );
  % store_info{ 1 } = info;
  % info = [];
  % U(:,1) = u + U(:,1);
  %
  % if ( length( suggested_arnoldi_size ) > 1 )
  %   opts.r_arn = suggested_arnoldi_size( 2 );
  % end
  % if length( store_info ) > 1
  %   info = store_info{ 2 };
  % else
  %   if ( nargin == 8 )
  %     info = aux;
  %   end
  % end
  % [ U, info ] = bamphi( k * [ 1/3, 1/2 ], A, [], [z,Fn,2/k*(g(t+1/2*k,U(:,1)) - gn)], opts, info );
  % suggested_arnoldi_size( 2 ) = bamphi_suggest_arnoldi_size( opts, info );
  % info = rmfield( info, 'A' );
  % store_info{ 2 } = info;
  % info = [];
  % U(:,2) = g( t + 1/2 * k, u + U(:,2) ) - gn;
  % U(:,1) = g( t + 1/3 * k, u + U(:,1) ) - gn;
  %
  % if ( length( suggested_arnoldi_size ) > 2 )
  %   opts.r_arn = suggested_arnoldi_size( 3 );
  % end
  % if length( store_info ) > 2
  %   info = store_info{ 3 };
  % else
  %   if ( nargin == 8 )
  %     info = aux;
  %   end
  % end
  % [ U,info ] = bamphi( k * [ 1/3, 5/6 ], A, [], [z,Fn,6/k*(3/2*U(:,1)-2/3*U(:,2)),12/(k^2)*(2*U(:,2)-3*U(:,1))], opts, info );
  % suggested_arnoldi_size( 3 ) = bamphi_suggest_arnoldi_size( opts, info );
  % info = rmfield( info, 'A' );
  % store_info{ 3 } = info;
  % info = [];
  % U(:,2) = g( t + 5/6 * k, u + U(:,2) ) - gn;
  % U(:,1) = g( t + 1/3 * k, u + U(:,1) ) - gn;
  %
  % if ( length( suggested_arnoldi_size ) > 3 )
  %   opts.r_arn = suggested_arnoldi_size( 4 );
  % end
  % if length( store_info ) > 3
  %   info = store_info{ 4 };
  % else
  %   if ( nargin == 8 )
  %     info = aux;
  %   end
  % end
  % [ z,info ] = bamphi( k, A, [], [z,Fn,2/k*(5/2*U(:,1)-2/5*U(:,2)),4/(k^2)*(6/5*U(:,2)-3*U(:,1))], opts, info );
  % suggested_arnoldi_size( 4 ) = bamphi_suggest_arnoldi_size( opts, info );
  % info = rmfield( info, 'A' );
  % store_info{ 4 } = info;
  % info = [];
  % u = u + z;


end
