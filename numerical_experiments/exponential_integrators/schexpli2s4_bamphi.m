function [ u, info ] = schexpli2s4_bamphi( u, k, t, A, g, opts, info, B )
% schexpli2s4_bamphi
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
% input to schexpli2s4_bamphi and any strategy is overruled by this.
%
% Hybrid strategies, such as the 2nd but recomputing (IOM)-Arnoldi every P expint
% steps are possible but have to be handled in concert with the outside of this
% function.
%

  mu = g( NaN, NaN );
  ut = u;

  % % 1st STRATEGY:
  %
  % if isempty( info ) && ( nargin == 8 ),
  %   bamphi_fov( length( u ), B ); % B is a matrix we rather use for computing A's FoV
  %   info = aux;
  % end
  % [ U, info ] = bamphi( [ -1i * k, 1i * k ], A, [], ut, opts, info );
  %
  % [ u, info ] = bamphi( -2i * k, A, [], [ zeros( size( ut ) ), conj( U(:,1) ) / ( -2i * k ), conj( U( :,1 ) ) / ( 4 * k^2 ) ], opts, info );
  % U( :,2 ) = ( -1i * mu * k ) * ( U( :,2 ).^2 .* u );
  %
  % [ u, info ] = bamphi( -2i * k, A, [], [ zeros( length( ut ), 2 ), -conj( ut ) / ( 4 * k^2 ) ], opts, info );
  % u = ut - ( 1i * mu * k ) * ut.^2 .* u - ( mu^2 * k^2 * (1/2) ) * ( ut .* conj( ut ) ).^2 .* ut;
  %
  % [ u, info ] = bamphi( 1i * k, A, [], u, opts, info );
  % u = u + U( :,2 );

  % 2nd STRATEGY:
  persistent store_info

  if length( store_info ) > 0,
    info = store_info{ 1 };
  end
  [ U, info ] = bamphi( [ -1i * k, 1i * k ], A, [], ut, opts, info );
  store_info{ 1 } = info;

  if length( store_info ) > 1,
    info = store_info{ 2 };
  end
  [ u, info ] = bamphi( -2i * k, A, [], [ zeros( size( ut ) ), conj( U(:,1) ) / ( -2i * k ), conj( U( :,1 ) ) / ( 4 * k^2 ) ], opts, info );
  store_info{ 2 } = info;

  U( :,2 ) = ( -1i * mu * k ) * ( U( :,2 ).^2 .* u );

  if length( store_info ) > 2,
    info = store_info{ 3 };
  end
  [ u, info ] = bamphi( -2i * k, A, [], [ zeros( length( ut ), 2 ), -conj( ut ) / ( 4 * k^2 ) ], opts, info );
  store_info{ 3 } = info;

  u = ut - ( 1i * mu * k ) * ut.^2 .* u - ( mu^2 * k^2 * (1/2) ) * ( ut .* conj( ut ) ).^2 .* ut;

  if length( store_info ) > 3,
    info = store_info{ 4 };
  end
  [ u, info ] = bamphi( 1i * k, A, [], u, opts, info );
  store_info{ 4 } = info;

  u = u + U( :,2 );

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
  % [ U, info ] = bamphi( [ -1i * k, 1i * k ], A, [], ut, opts, info );
  % suggested_arnoldi_size( 1 ) = bamphi_suggest_arnoldi_size( opts, info );
  % info = rmfield( info, 'A' );
  % store_info{ 1 } = info;
  % info = [];
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
  % [ u, info ] = bamphi( -2i * k, A, [], [ zeros( size( ut ) ), conj( U(:,1) ) / ( -2i * k ), conj( U( :,1 ) ) / ( 4 * k^2 ) ], opts, info );
  % suggested_arnoldi_size( 2 ) = bamphi_suggest_arnoldi_size( opts, info );
  % info = rmfield( info, 'A' );
  % store_info{ 2 } = info;
  % info = [];
  %
  % U( :,2 ) = ( -1i * mu * k ) * ( U( :,2 ).^2 .* u );
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
  % [ u, info ] = bamphi( -2i * k, A, [], [ zeros( length( ut ), 2 ), -conj( ut ) / ( 4 * k^2 ) ], opts, info );
  % suggested_arnoldi_size( 3 ) = bamphi_suggest_arnoldi_size( opts, info );
  % info = rmfield( info, 'A' );
  % store_info{ 3 } = info;
  % info = [];
  %
  % u = ut - ( 1i * mu * k ) * ut.^2 .* u - ( mu^2 * k^2 * (1/2) ) * ( ut .* conj( ut ) ).^2 .* ut;
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
  % [ u, info ] = bamphi( 1i * k, A, [], u, opts, info );
  % suggested_arnoldi_size( 4 ) = bamphi_suggest_arnoldi_size( opts, info );
  % info = rmfield( info, 'A' );
  % store_info{ 4 } = info;
  % info = [];
  %
  % u = u + U( :,2 );

end
