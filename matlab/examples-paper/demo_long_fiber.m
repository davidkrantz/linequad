% Evaluate Stokes slender body potential from closed-loop squiggle curve.
% Computes using interpolatory quadrature, and compares to semi-smart adaptive
% quadrature, which is reasonably fast and seems to compute the right thing.
% Slow computations of reference using integral() is also available.
%
% Matlab code attempts to be fast by avoiding inner loop vector ops, so that
% we can compute meaningful running times.

% Alex variant also includes all targs close to curve test - see targtype.

% Alex added (6) case for 3-digit slice test at nquad=4

function demo_long_fiber(varargin) % Function makes Matlab JIT compile better

SAVEPLOTS = false;

if nargin==0
    test_no = 2;
else
    test_no = varargin{1};
end

%% Default setup
slender_eps = 1e-3;
%WHICH_REFERENCE = 'integral';
WHICH_REFERENCE = 'adaptive';

UPSAMPLE = true; % Upsample before applying interpolatory quadrature
nquad = 16;
rho = 3; % Interpolatory quadrature limit rule
Hlim = 1; % Adaptive quadrature distance rule
%rho = 1 % Disable specquad

% When running near points
dist = 1e-3;            % the dist of all pts from the curve

%% Case setup
switch test_no
  case 0
    % Play
    targtype = 'n';   % 's' for original slice, or 'n' for all targs near curve at random locs.
    tol = 1e-6;
  case 1
    % Paper plot for field
    targtype = 's';
    tol = 1e-10;      
  case 2
    % Comparison for nearby points
    targtype = 'n';    
    tol = 1e-6;
    dist = 1e-2;
  case 3
    targtype = 'n';    
    tol = 1e-10;
    dist = 1e-2;
  case 4
    targtype = 'n';    
    tol = 1e-6;
    dist = 1e-4;
  case 5
    targtype = 'n';    
    tol = 1e-10;
    dist = 1e-4;
  case 6                    % try lower nquad for lower acc.
   nquad = 4;  % #digits + 1
   targtype = 's';    
   tol = 1e-3;
  otherwise
    error('Unknown test no')
end

%% Run

% Setup fiber
[x, y, z, xp, yp, zp] = squiggle();
s = @(t) sqrt(xp(t).^2 + yp(t).^2 + zp(t).^2);

% Artificial density (trivial)
f1 = @(t) x(t);
f2 = @(t) y(t);
f3 = @(t) z(t);

% Discretize
fprintf('* Discretizing, tol=%.1e\n', tol)
[tj, wj, npan] = adaptive_panelization(s, nquad, tol);
fprintf('nquad=%d, npan=%d\n', nquad, npan);
xj = x(tj); yj = y(tj); zj = z(tj);
xpj = xp(tj); ypj = yp(tj); zpj = zp(tj);
sj = s(tj);
f1j = f1(tj); f2j = f2(tj); f3j = f3(tj);


if targtype=='s'     % Evaluation surface - slice
  Neval = 100;    % in each direction
  if SAVEPLOTS
      Neval = 200;
  end
  [X, Z, Y] = meshgrid( linspace(-1.4,1.4, Neval), ...                      
                        linspace(-1.4,1.4, Neval), ...
                        0.25);
elseif targtype=='n'  % all targs near curve, in random normal directions a fixed dist away
  Ne = 5e3;  % # targs
  t = rand(1,Ne); 
  v = randn(3,Ne); utang = [xp(t);yp(t);zp(t)]./s(t); % sloppy unit tangents
  vdotutang = sum(v.*utang,1); v = v - utang.*vdotutang;  % orthog v against the tangent
  v = v./sqrt(sum(v.*v,1));    % normalize all the v vecs
  X = x(t) + dist*v(1,:);   % displace by v vecs from pts on curve
  Y = y(t) + dist*v(2,:);
  Z = z(t) + dist*v(3,:);
end


% Compute ref solution
[uref1, uref2, uref3] = deal(zeros(size(X)));
switch WHICH_REFERENCE
  case 'integral'
    % Use integral() reference (slooooow, so only do one component)
    uref1 = zeros(size(X));
    disp('* Reference: integral()')
    tic
    parfor i=1:numel(X)
        uref1(i) = integral(@(t) s(t) .* ...
                            slender_body_kernel(x(t)-X(i), y(t)-Y(i), z(t)-Z(i), f1(t), f2(t), f3(t), slender_eps), ...
                            0, 1, 'abstol', 1e-15, 'reltol', 1e-15);
    end
    toc
  case 'adaptive'
    % Adaptive quadrature
    disp('* Reference: Adaptive')
    % Make sure that we get a different discretization for reference computations,
    % otherwise errors get artifically small.
    nquad_ref = nquad+2;
    tol_ref = 5e-14;
    if nquad<16, tol_ref = 1e-2*tol;  end    % low-acc test cases
    Hlim_ref = 1.0;
    [tj_ref, wj_ref, npan_ref] = adaptive_panelization(s, nquad_ref, tol_ref);
    fprintf('Discretization: nquad=%d, npan=%d\n', nquad_ref, npan_ref);
    [uref1, uref2, uref3] = adaptive_quadrature(x(tj_ref), y(tj_ref), z(tj_ref), s(tj_ref), ...
                                                wj_ref, ...
                                                f1(tj_ref), f2(tj_ref), f3(tj_ref), ...
                                                X, Y, Z, nquad_ref, slender_eps, Hlim_ref);
end

% Compute adaptive quadrature
disp(' ')
disp('* Adaptive quadrature')
[adquad1, adquad2, adquad3, adstats] = adaptive_quadrature(xj, yj, zj, sj, wj, f1j, f2j, f3j, ...
                                             X, Y, Z, nquad, slender_eps, Hlim);

% Compute special quadrature
disp(' ')
disp('* Interpolatory quadrature')
[specquad1,specquad2,specquad3, specstats] = interpolatory_quadrature(...
    xj, yj, zj, sj, wj, f1j, f2j, f3j, X, Y, Z, ...
    nquad, rho, UPSAMPLE, slender_eps);

% Compute errors
if strcmp(WHICH_REFERENCE, 'integral')
    % We only have reference for 1st component
    adquad2 = uref2;
    adquad3 = uref3;
    specquspec2 = uref2;
    specquspec3 = uref3;
end

adquad_errmax = compute_error(uref1, uref2, uref3, adquad1, adquad2, adquad3);
adquad_errmaxmax = max(adquad_errmax(:))

specquad_errmax = compute_error(uref1, uref2, uref3, specquad1, specquad2, specquad3);
specquad_errmaxmax = max(specquad_errmax(:))


% Plot
if SAVEPLOTS
    pubfig = @publication_fig;
else
    pubfig = @() false;
end

sfigure(1);
clf; pubfig();
plot3(xj, yj, zj, '.-k')
if targtype=='n', hold on; plot3(X,Y,Z,'r.'); end
axis equal tight vis3d
grid on
box on

sfigure(2);
clf; pubfig();
if targtype=='s'
  surf(X, Y, Z, log10(specquad_errmax))                  
  shading flat
  xlabel(colorbar(), 'log_{10} E_{rel}')  
  caxis auto
  hold on
  plot3(xj, yj, zj, '.-k')
  axis equal tight vis3d
  grid off
  box on
else
  semilogy(t,specquad_errmax,'.'); xlabel('t'); ylabel('err')
end
if ~SAVEPLOTS && targtype=='s'
    title('Interpolatory quadrature')
end

  
sfigure(3);
clf; pubfig();
if targtype=='s'
  surf(X, Y, Z, log10(adquad_errmax))                  
  shading flat
  colorbar
  caxis auto
  hold on
  plot3(xj, yj, zj, '.-k')
  axis equal tight vis3d
  grid off
  box on
else
  semilogy(t,adquad_errmax,'.'); xlabel('t'); ylabel('err')
end
title('Adaptive quadrature')

% Test-specific plot saving
if SAVEPLOTS && test_no==1
    sfigure(2);
    axis off
    caxis([-16,-13])
    text(-1,0,-2, ['$\max E_{rel}=$', ...
                        sprintf('%.1e',specquad_errmaxmax)],'interpreter','latex')    
    saveplot(1, 'long_fiber_geo.png')
    saveplot(2, 'long_fiber_field.png')
end

% table line
adstats
fprintf(' %.1e &  %.1e & %.1e & %.2f & %.2f & %.1e & %.1e & %.2f & %.2f & %.1e \\\\ ', dist, tol, ...
        adstats.kerevals_near, adstats.time_ker_near, adstats.time_interp, adquad_errmaxmax,...
        specstats.kerevals_near, specstats.time_ker_near, specstats.time_weights, specquad_errmaxmax)
disp(['% case ' num2str(test_no)])
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function saveplot(num, filename)
    sfigure(num);
    publication_fig();
    path = ['../../linequadpaper/fig/', filename];
    print('-dpng', '-r200', path)
    system(['mogrify -trim ' path])    
    disp(['Wrote ' path])
end

%%%%% INTERPOLATORY QUADRATURE

function [specquad1,specquad2,specquad3,stats] = interpolatory_quadrature(...
    xj, yj, zj, sj, wj, f1j, f2j, f3j, X, Y, Z, nquad, rho, UPSAMPLE, slender_eps)
    [specquad1,specquad2,specquad3] = deal(zeros(size(X)));
    npan = numel(xj)/nquad;    
    [tgl, wgl] = legendre.gauss(nquad);
    if UPSAMPLE
        nquad2 = 2*nquad;
        [tgl2, wgl2] = legendre.gauss(nquad2);
        B = bclag_interp_matrix(tgl, tgl2);
        %rho = sqrt(rho);
    else
        nquad2 = nquad;
        B = 1;
        tgl2 = tgl;
        wgl2 = wgl;
    end
    time_weights = 0;
    kerevals_near = 0;
    maintic = tic();
    time_ker_near = 0;
    time_far = 0;
    for j=1:npan
        % Load panel
        idx = (1:nquad) + nquad*(j-1);
        xjpan = xj(idx);
        yjpan = yj(idx);
        zjpan = zj(idx);
        sjpan = sj(idx);
        wjpan = wj(idx);
        f1jpan = f1j(idx);
        f2jpan = f2j(idx);
        f3jpan = f3j(idx);    
        % Upsample panel
        xjpan_up  = xjpan *B';
        yjpan_up  = yjpan *B';
        zjpan_up  = zjpan *B';
        sjpan_up  = sjpan *B';
        f1jpan_up = f1jpan*B';
        f2jpan_up = f2jpan*B';
        f3jpan_up = f3jpan*B';
        % Compute quadrature weights
        atic = tic();
        [all_w1, all_w3, all_w5, specquad_needed] = line3_near_weights(tgl2, wgl2, xjpan_up, yjpan_up, zjpan_up, ...
                                                          X, Y, Z, rho);
        % Evaluation count
        time_weights = time_weights + toc(atic);
        kerevals_near = kerevals_near + nquad2*sum(specquad_needed(:));
        % Evaluate each panel-to-point pair
        for i=1:numel(X)    
            Xi = X(i);
            Yi = Y(i);
            Zi = Z(i);
            q1 = 0; q2 = 0; q3 = 0;
            if specquad_needed(i)
                atic = tic();
                for k=1:nquad2
                    r1 = xjpan_up(k)-Xi;
                    r2 = yjpan_up(k)-Yi;
                    r3 = zjpan_up(k)-Zi;
                    [u1R1, u1R3, u1R5, u2R1, u2R3, u2R5, u3R1, u3R3, u3R5] ...
                        = slender_body_kernel_split(r1, r2, r3, f1jpan_up(k), f2jpan_up(k), f3jpan_up(k), ...
                                                    slender_eps);
                    q1 = q1 + ...
                         all_w1(k,i)*sjpan_up(k)*u1R1 + ...
                         all_w3(k,i)*sjpan_up(k)*u1R3 + ...
                         all_w5(k,i)*sjpan_up(k)*u1R5;                
                    q2 = q2 + ...
                         all_w1(k,i)*sjpan_up(k)*u2R1 + ...
                         all_w3(k,i)*sjpan_up(k)*u2R3 + ...
                         all_w5(k,i)*sjpan_up(k)*u2R5;                
                    q3 = q3 + ...
                         all_w1(k,i)*sjpan_up(k)*u3R1 + ...
                         all_w3(k,i)*sjpan_up(k)*u3R3 + ...
                         all_w5(k,i)*sjpan_up(k)*u3R5;                                
                end      
                % Rescale (weights are for [-1,1])
                q1 = q1/2*sum(wjpan);
                q2 = q2/2*sum(wjpan);
                q3 = q3/2*sum(wjpan);            
                time_ker_near =  time_ker_near + toc(atic);
            else            
                atic = tic();
                [q1, q2, q3] = quadsum(xjpan, yjpan, zjpan, sjpan, wjpan, f1jpan, f2jpan, f3jpan, ...
                                       Xi, Yi, Zi, nquad, slender_eps);
                time_far = time_far + toc(atic);
            end
            specquad1(i) = specquad1(i) + q1;        
            specquad2(i) = specquad2(i) + q2;        
            specquad3(i) = specquad3(i) + q3;                
        end
    end
    toc(maintic)
    fprintf('Near field kernel evals: %e\n', kerevals_near);
    fprintf('Near field kernel evals time: %f\n', time_ker_near)
    fprintf('Total time line3_near_weights: %f\n', time_weights)
    fprintf('Far field time: %f\n', time_far)
    stats.kerevals_near = kerevals_near;
    stats.time_weights = time_weights;
    stats.time_ker_near = time_ker_near;
end    
