function [z,GammaDCTy,s,exitflag,dz] = smoothn(varargin)

%SMOOTHN Fast smoothing of 1-D to N-D data.
%   Z = SMOOTHN(Y) automatically smoothes the uniformly-sampled array Y. Y
%   can be any N-D noisy array (time series, images, 3D data,...). Non
%   finite data (NaN or Inf) are treated as missing values.
%
%   Z = SMOOTHN(Y,S) smoothes the array Y using the smoothing parameter S.
%   S must be a real positive scalar. The larger S is, the smoother the
%   output will be. If the smoothing parameter S is omitted (see previous
%   option) or empty (i.e. S = []), it is automatically determined using
%   the generalized cross-validation (GCV) method.
%
%   Z = SMOOTHN(Y,W) or Z = SMOOTHN(Y,W,S) specifies a weighting array W of
%   real positive values, that must have the same size as Y. Note that a
%   nil weight corresponds to a missing value.
%
%   Robust smoothing
%   ----------------
%   Z = SMOOTHN(...,'robust') carries out a robust smoothing that minimizes
%   the influence of outlying data.
%
%   Upsampling
%   ----------
%   Z = SMOOTHN(...,'upsample',N) increases sampling rate of Z by N along
%   non singleton dimensions. If N is omitted, a factor N = 2 is used, i.e
%   SMOOTHN(...,'upsample') is the same as SMOOTHN(...,'upsample',2);
%
%   [Z,S] = SMOOTHN(...) also returns the calculated value for S so that
%   you can fine-tune the smoothing subsequently if needed.
%
%   An iteration process is used in the presence of weighted and/or missing
%   values. Z = SMOOTHN(...,OPTION_NAME,OPTION_VALUE) smoothes with the
%   termination parameters specified by OPTION_NAME and OPTION_VALUE. They
%   can contain two criteria:
%       -----------------
%       TolZ:       Termination tolerance on Z (default = 1e-3)
%                   TolZ must be in ]0,1[
%       MaxIter:    Maximum number of iterations allowed (default = 30)
%       -----------------
%   Syntax: [Z,...] = SMOOTHN(...,'MaxIter',50,'TolZ',1e-4)
%
%   [Z,S,EXITFLAG] = SMOOTHN(...) returns a boolean value EXITFLAG that
%   describes the exit condition of SMOOTHN:
%       1       SMOOTHN converged.
%       0       Maximum number of iterations was reached.
%
%   Gradient
%   --------
%   [Z,S,EXIFLAG,DZ] = SMOOTHN(...) returns the numerical gradient of Z.
%   The gradient option is still in a Beta version. It works with 1D and 2D
%   arrays only and the image processing & partial differential equation
%   toolboxes are required. If Z is a matrix then DZ is a cell: DZ{1}
%   corresponds to dZ/dx, the differences in x (horizontal) direction, and
%   DZ{2} corresponds to dZ/dy, the differences in y (vertical) direction.
%   The spacing between points in each direction is assumed to be one.
%
%   Over- & under-smoothing
%   -----------------------
%   [...] = SMOOTHN(...,'over') or [...] = SMOOTHN(...,'under') forces
%   over- or under-smoothing. The smoothing parameter S is first determined
%   automatically then multiplied (over-smoothing) or divided (under-
%   smoothing) by 100.
%
%   Class Support
%   -------------
%   Input array can be numeric or logical. The returned array is of class
%   double.
%
%   Notes
%   -----
%   The N-D (inverse) discrete cosine transform functions <a
%   href="matlab:web('http://www.biomecardio.com/matlab/dctn.html')"
%   >DCTN</a> and <a
%   href="matlab:web('http://www.biomecardio.com/matlab/idctn.html')"
%   >IDCTN</a> are required.
%
%   Reference
%   --------- 
%   Garcia D, Fast and robust automatic smoothing of one-dimensional to
%   multidimensional data. Dealing with occurrence of missing values.
%   Computational Statistics & Data Analysis, 2009, under review.
%
%   Examples:
%   --------
%   % 1-D example
%   x = linspace(0,100,2^8);
%   y = cos(x/10)+(x/50).^2 + randn(size(x))/10;
%   y([70 75 80]) = [5.5 5 6];
%   z = smoothn(y); % Regular smoothing
%   zr = smoothn(y,'robust'); % Robust smoothing
%   subplot(121), plot(x,y,'r.',x,z,'k','LineWidth',2)
%   axis square, title('Regular smoothing')
%   subplot(122), plot(x,y,'r.',x,zr,'k','LineWidth',2)
%   axis square, title('Robust smoothing')
%
%   % 2-D example
%   xp = 0:.02:1;
%   [x,y] = meshgrid(xp);
%   f = exp(x+y) + sin((x-2*y)*3);
%   fn = f + randn(size(f))*0.5;
%   fs = smoothn(fn);
%   subplot(121), surf(xp,xp,fn), zlim([0 8]), axis square
%   subplot(122), surf(xp,xp,fs), zlim([0 8]), axis square
%
%   % 2-D example with missing data
%   n = 300;
%   y0 = peaks(n);
%   y = y0 + rand(size(y0))*3;
%   I = randperm(n^2);
%   y(I(1:n^2*0.5)) = NaN; % lose 1/2 of data
%   y(50:100,200:250) = NaN; % create a hole
%   z = smoothn(y); % smooth data
%   subplot(2,2,1:2), imagesc(y), axis equal off
%   title('Noisy corrupted data')
%   subplot(223), imagesc(z), axis equal off
%   title('Recovered data ...')
%   subplot(224), imagesc(y0), axis equal off
%   title('... compared with original data')
%
%   % 3-D example
%   [x,y,z] = meshgrid(-2:.2:2);
%   xslice = [-0.8,1]; yslice = 2; zslice = [-2,0];
%   vn = x.*exp(-x.^2-y.^2-z.^2) + randn(size(x))*0.06;
%   subplot(121), slice(x,y,z,vn,xslice,yslice,zslice,'cubic')
%   title('Noisy data')
%   v = smoothn(vn);
%   subplot(122), slice(x,y,z,v,xslice,yslice,zslice,'cubic')
%   title('Smoothed data')
%
%   % Cardioid
%   t = linspace(0,2*pi,1000);
%   x = 2*cos(t).*(1-cos(t)) + randn(size(t))*0.1;
%   y = 2*sin(t).*(1-cos(t)) + randn(size(t))*0.1;
%   z = smoothn(complex(x,y));
%   plot(x,y,'r.',real(z),imag(z),'k','linewidth',2)
%   axis equal tight
%
%   % Cellular vortical flow
%   [x,y] = meshgrid(linspace(0,1,24));
%   Vx = cos(2*pi*x+pi/2).*cos(2*pi*y);
%   Vy = sin(2*pi*x+pi/2).*sin(2*pi*y);
%   Vx = Vx + sqrt(0.05)*randn(24,24); % adding Gaussian noise
%   Vy = Vy + sqrt(0.05)*randn(24,24); % adding Gaussian noise
%   I = randperm(numel(Vx));
%   Vx(I(1:30)) = (rand(30,1)-0.5)*5; % adding outliers
%   Vy(I(1:30)) = (rand(30,1)-0.5)*5; % adding outliers
%   Vx(I(31:60)) = NaN; % missing values
%   Vy(I(31:60)) = NaN; % missing values
%   Vs = smoothn(complex(Vx,Vy),'robust'); % automatic smoothing
%   subplot(121), quiver(x,y,Vx,Vy,2.5), axis square
%   title('Noisy velocity field')
%   subplot(122), quiver(x,y,real(Vs),imag(Vs)), axis square
%   title('Smoothed velocity field')
%
%   See also SMOOTH, DCTN, IDCTN.
%
%   -- Damien Garcia -- 2009/03, revised 2009/09

% Check input arguments
error(nargchk(1,10,nargin));

%% Test & prepare the variables
%---
k = 0;
while k<nargin && ~ischar(varargin{k+1}), k = k+1; end
%---
% y = array to be smoothed
y = double(varargin{1});
sizy = size(y);
noe = prod(sizy); % number of elements
if noe<2, z = y; return, end
%---
% Smoothness parameter and weights
W = ones(sizy);
s = [];
if k==2
    if isempty(varargin{2}) || isscalar(varargin{2}) % smoothn(y,s)
        s = varargin{2}; % smoothness parameter
    else % smoothn(y,W)
        W = varargin{2}; % weight array
    end
elseif k==3 % smoothn(y,W,s)
        W = varargin{2}; % weight array
        s = varargin{3}; % smoothness parameter
end
if ~isequal(size(W),sizy)
        error('MATLAB:smoothn:SizeMismatch',...
            'Arrays for data and weights must have same size.')
elseif ~isempty(s) && (~isscalar(s) || s<0)
    error('MATLAB:smoothn:IncorrectSmoothingParameter',...
        'The smoothing parameter must be a scalar >=0')
end
%---
% "Maximal number of iterations" criterion
I = find(strcmpi(varargin,'MaxIter'),1);
if isempty(I)
    MaxIter = 30; % default value for MaxIter
else
    try
        MaxIter = varargin{I+1};
    catch
        error('MATLAB:smoothn:IncorrectMaxIter',...
            'MaxIter must be an integer >=1')
    end
    if ~isscalar(MaxIter) || MaxIter<1 || MaxIter~=round(MaxIter)
        error('MATLAB:smoothn:IncorrectMaxIter',...
            'MaxIter must be an integer >=1')        
    end    
end
%---
% "Tolerance on smoothed output" criterion
I = find(strcmpi(varargin,'TolZ'),1);
if isempty(I)
    TolZ = 1e-3; % default value for TolZ
else
    try
        TolZ = varargin{I+1};
    catch
        error('MATLAB:smoothn:IncorrectTolZ',...
            'TolZ must be in ]0,1[')
    end
    if ~isscalar(TolZ) || TolZ<=0 || TolZ>=1 
        error('MATLAB:smoothn:IncorrectTolZ',...
            'TolZ must be in ]0,1[')
    end    
end
%---
% "Upsample" criterion
I = find(strcmpi(varargin,'upsample'),1);
if isempty(I)
    nupsample = 1;
else
    try
        nupsample = varargin{I+1};
    catch
        nupsample = 2;
    end
%     if ~isnumeric(nupsample) || ~isscalar(nupsample)
%         nupsample = 2;
%     elseif nupsample<1 
%         error('MATLAB:smoothn:IncorrectUpSample',...
%             'Upsampling factor must be >= 1')
%     end    
end
%---
% Weights. Zero weights are assigned to not finite values (Inf or NaN),
% (Inf/NaN values = missing data).
IsFinite = isfinite(y);
nof = nnz(IsFinite); % number of finite elements
W = W.*IsFinite;
if any(W<0)
    error('MATLAB:smoothn:NegativeWeights',...
        'Weights must all be >=0')
else 
    W = W/max(W(:));
end
%---
% Weighted or missing data?
isweighted = true;
if all(W(:)==1), isweighted = false; end
%---
% Robust smoothing?
isrobust = false;
if any(strcmpi(varargin,'robust')), isrobust = true; end
%---
% Automatic smoothing?
isauto = false;
if isempty(s), isauto = true; end
%---
% Over- or under-smoothing?
isover = any(strcmpi(varargin,'over'));
isunder = any(strcmpi(varargin,'under'));
if isover&&isunder
    isover = false; isunder = false;
elseif isover&&~isauto
    error('MATLAB:smoothn:OverOption',...
        '''over'' option requires automatic smoothing.')
elseif isunder&&~isauto
    error('MATLAB:smoothn:UnderOption',...
        '''under'' option requires automatic smoothing.')    
end
%---
% DCTN and IDCTN are required
test4DCTNandIDCTN

%% Creation of the Lambda tensor
%---
% Lambda contains the eingenvalues of the difference matrix used in this
% penalized least squares process.
d = ndims(y);
Lambda = zeros(sizy);
for i = 1:d
    siz0 = ones(1,d);
    siz0(i) = sizy(i);
    Lambda = bsxfun(@plus,Lambda,...
        cos(pi*(reshape(1:sizy(i),siz0)-1)/sizy(i)));
end
Lambda = -2*(d-Lambda);

%% Upper and lower bound for the smoothness parameter
% The average leverage (h) is by definition in [0 1]. Weak smoothing occurs
% if h is close to 1, while over-smoothing appears when h is near 0. Upper
% and lower bounds for h are given to avoid under- or over-smoothing. See
% equation relating h to the smoothness parameter.
N = sum(sizy~=1); % tensor rank of y-array
hMin = 1e-6; hMax = 0.99;
sMinBnd = (((1+sqrt(1+8*hMax.^(2/N)))/4./hMax.^(2/N)).^2-1)/16;
sMaxBnd = (((1+sqrt(1+8*hMin.^(2/N)))/4./hMin.^(2/N)).^2-1)/16;

%% Initialize before iterating
%---
Wtot = W;
y(~IsFinite) = 0; % arbitrary values for missing y-data
z = y;
if exist('inpaint_nans','file') && N<3 && isreal(z)
    z(~IsFinite) = NaN;
    z = inpaint_nans(z);
elseif exist('inpaint_nans3','file') && N==3 && isreal(z)
    z(~IsFinite) = NaN;
    z = inpaint_nans3(z);
else
    z(~IsFinite) = mean(y(IsFinite));
end


z0 = z;
tol = Inf;
RobustIterativeProcess = true;
RobustStep = 0;
nit = 0;
%--- Error on p. Smoothness parameter s = 10^p
errp = 0.05;
opt = optimset('TolX',errp);
%---

%% Main iterative process
%---
while RobustIterativeProcess
    %--- "amount" of weights (see the function GCVscore)
    aow = sum(Wtot(:))/noe; % 0 < aow <= 1
    %---
    while tol>TolZ && nit<MaxIter
        nit = nit+1;
        if ~isauto
            if ~exist('Gamma','var')
                Gamma = 1./(1+s*Lambda.^2);
            end
            GammaDCTy = Gamma.*dctn(Wtot.*(y-z)+z);
            z = idctn(GammaDCTy);
        else
            %---
            % The generalized cross-validation (GCV) method is used.
            % We seek the smoothing parameter s that minimizes the GCV
            % score i.e. s = Argmin(GCVscore)
            %---
            DCTy = dctn(Wtot.*(y-z)+z);
            fminbnd(@gcv,log10(sMinBnd),log10(sMaxBnd),opt);
            GammaDCTy = Gamma.*DCTy;
            z = idctn(GammaDCTy);
        end
        if ~isweighted
            tol = 0; % no iteration if no weighted/missing data
        else
            tol = norm(z0(:)-z(:))/norm(z(:));
            z0 = z;
        end
    end
    exitflag = nit<MaxIter;
    
    %------
    % After one robust iterative loop, the GCV method is stopped and the
    % last smoothness parameter (s) is used. It makes the algorithm faster
    % without altering the final results significantly.
    if RobustStep==1, isauto = false; end
    %------

    if isrobust % Robust Smoothing: iteratively re-weighted process
        %--- average leverage
        h = sqrt(1+16*s);
        h = sqrt(1+h)/sqrt(2)/h;
        h = h^N;
        %--- robust weights
        Wrobust = RobustWeights(y-z,IsFinite,h);
        Wtot = W.*Wrobust;
        %--- re-initialize for another iterative weighted process
        isweighted = true; tol = Inf; nit = 0; 
        %---
        RobustStep = RobustStep+1;
        RobustIterativeProcess = RobustStep<5; % 5 robust steps are enough.
    else
        RobustIterativeProcess = false; % stop the whole process
    end
end

%% Over- / Under-smoothing
%---
if isunder || isover
    s = isunder*s/100 + isover*s*100;
    Gamma = 1./(1+s*Lambda.^2);
    GammaDCTy = Gamma.*DCTy;
    z = idctn(GammaDCTy);
end

%% Warning messages
%---
if abs(log10(s)-log10(sMinBnd))<errp
    warning('MATLAB:smoothn:SLowerBound',...
            ['s = ' num2str(s,'%.3e') ': the lower bound for s ',...
            'has been reached. Put s as an input variable if required.'])
elseif abs(log10(s)-log10(sMaxBnd))<errp
    warning('MATLAB:smoothn:SUpperBound',...
            ['s = ' num2str(s,'%.3e') ': the upper bound for s ',...
            'has been reached. Put s as an input variable if required.'])    
end
if nargout<3 && ~exitflag
    warning('MATLAB:smoothn:MaxIter',...
            ['Maximum number of iterations (' int2str(MaxIter) ') has ',...
            'been exceeded. Increase MaxIter option or decrease TolZ value.'])    
end

%% Upsampling
%---
if nupsample>1
    upsiz = round(nupsample.*sizy);
    upsiz(sizy==1) = 1;
    z = zeros(upsiz);
    I = ['1:' int2str(sizy(1))];
    for i = 2:length(sizy)
        I = strcat(I,[',1:',int2str(sizy(i))]);
    end
    eval(['z(' I ') = GammaDCTy;'])
    z = idctn(z*sqrt(prod(upsiz)/noe));
end

%% Gradient
%---
% Beta version:
% 1) For vectors and matrices only
% 2) Requires IDCT and DST (Image Processing & Partial Differential
%    Equation Toolboxes)
%---
% To do: DCTD and IDSTD which will allow DCT and IDST for N-dimensional
% array along dimension D.
%---
test = license('test','image_toolbox') & license('test','pde_toolbox');
if nargout==4    
    if isvector(z) && test
       Z = GammaDCTy(:);
       Z(1) = [];
       nz = length(z);
       Z = -Z*sqrt(2/noe)/noe*pi;
       Z = Z.*(1:length(Z))';
       dz = dst(Z,2*nz-1);
       dz = reshape(dz(1:2:end,:),size(z));
       dz = dz/length(z)*length(y);
    elseif N==2 && test
        dz = cell(1,2);
        for k = 1:2
            GammaDCTy = shiftdim(GammaDCTy,1);
            Z = GammaDCTy;
            nz = numel(z)/size(z,k);
            nZ = size(Z,1);
            Z(1,:) = [];
            Z = -Z*sqrt(2/nZ)/nZ*pi;
            Z = bsxfun(@times,Z,(1:nZ-1)');
            dz{k} = dst(Z,2*nz-1);
            dz{k} = dz{k}(1:2:end,:);
            dz{k} = shiftdim(dz{k},1);
            dz{k} = idct(dz{k},numel(z)/nz);
            dz{k} = dz{k}*sqrt(nZ)/sqrt(nz);
        end
        dz{2} = -dz{2}.'; 
    else
        dz = [];
        warning('MATLAB:smoothn:BetaVersion',...
            ['The gradient option is still in a Beta version. ',...
            'It works with 1D and 2D arrays only. In addition, both ',...
            'image processing & PDE toolboxes are required.'])
    end
end


%% GCV score
%---
function GCVscore = gcv(p)
    % Search the smoothing parameter s that minimizes the GCV score
    %---
    s = 10^p;
    Gamma = 1./(1+s*Lambda.^2);
    %--- RSS = Residual sum-of-squares
    % if ~isweighted
    if aow>0.90 % aow = 1 means that all of the data are equally weighted
        % very much faster: does not require any inverse DCT
        RSS = norm(DCTy(:).*(Gamma(:)-1))^2;
    else
        % take account of the weights to calculate RSS:
        yhat = idctn(Gamma.*DCTy);
        RSS = norm(sqrt(Wtot(IsFinite)).*(y(IsFinite)-yhat(IsFinite)))^2;
    end
    %---
    TrH = sum(Gamma(:));
    GCVscore = RSS/nof/(1-TrH/noe)^2;
end

end

%% Robust weights
function W = RobustWeights(r,I,h)
    % weights for robust smoothing.
    MAD = median(abs(r(I)-median(r(I)))); % median absolute deviation
    u = abs(r/(1.4826*MAD)/sqrt(1-h)); % studentized residuals
    c = 4.685; W = (1-(u/c).^2).^2.*((u/c)<1); % bisquare weights
    % c = 2.385; W = 1./(1+(u/c).^2); % Cauchy weights
    W(isnan(W)) = 0; 
end

%% Test for DCTN and IDCTN
function test4DCTNandIDCTN
    if ~exist('dctn','file')
        error('MATLAB:smoothn:MissingFunction',...
            ['DCTN and IDCTN are required. Download <a href="matlab:web(''',...
            'http://www.biomecardio.com/matlab/dctn.html'')">DCTN</a>.'])
    elseif ~exist('idctn','file')
        error('MATLAB:smoothn:MissingFunction',...
            ['DCTN and IDCTN are required. Download <a href="matlab:web(''',...
            'http://www.biomecardio.com/matlab/idctn.html'')">IDCTN</a>.'])
    end
end
