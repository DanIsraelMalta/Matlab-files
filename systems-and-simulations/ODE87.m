%{
% An eight order dorman prince explicit integrator.
% Interface is equivalent to any matlab ODE solver.
%
% usage example:

dydt = @(y) [3 - 0.3 * y(2)^0.7; 0.5 * (1-y(1)^0.22)*y(2)-sqrt(y(1))];
[t1, y1] = ODE87(@vdp1,[0 20],[3 7]);
[t2, y2] = ode45(@vdp1,[0 20],[3 7]);
plot(t1, y1(:,1), 'b', t2, y2(:, 1), 'r');
legend('ODE87', 'ode45');
grid on;

%
% Dan I. Malta (2010)
%}
function [tout,xout] = ODE87(odefun,tspan,x0,options,varargin)
% coefficients
c_i =  [ 1/18, 1/12, 1/8, 5/16, 3/8, 59/400, 93/200, 5490023248/9719169821, 13/20, 1201146811/1299019798, 1, 1]';
a_i_j = [ 1/18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                   1/48, 1/16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                    1/32, 0, 3/32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                    5/16, 0, -75/64, 75/64, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                    3/80, 0, 0, 3/16, 3/20, 0, 0, 0, 0, 0, 0, 0, 0;
                    29443841/614563906, 0, 0, 77736538/692538347, -28693883/1125000000, 23124283/1800000000, 0, 0, 0, 0, 0, 0, 0;
                    16016141/946692911, 0, 0, 61564180/158732637, 22789713/633445777, 545815736/2771057229, -180193667/1043307555, 0, 0, 0, 0, 0, 0;
                    39632708/573591083, 0, 0, -433636366/683701615, -421739975/2616292301, 100302831/723423059, 790204164/839813087, 800635310/3783071287, 0, 0, 0, 0, 0;
                    246121993/1340847787, 0, 0, -37695042795/15268766246, -309121744/1061227803, -12992083/490766935, 6005943493/2108947869, 393006217/1396673457, 123872331/1001029789, 0, 0, 0, 0;
                    -1028468189/846180014, 0, 0, 8478235783/508512852, 1311729495/1432422823, -10304129995/1701304382, -48777925059/3047939560, 15336726248/1032824649, -45442868181/3398467696, 3065993473/597172653, 0, 0, 0;
                    185892177/718116043, 0, 0, -3185094517/667107341, -477755414/1098053517, -703635378/230739211, 5731566787/1027545527, 5232866602/850066563, -4093664535/808688257, 3962137247/1805957418, 65686358/487910083, 0, 0;
                    403863854/491063109, 0, 0, -5068492393/434740067, -411421997/543043805, 652783627/914296604, 11173962825/925320556, -13158990841/6184727034, 3936647629/1978049680, -160528059/685178525, 248638103/1413531060, 0, 0]';
b_8 = [ 14005451/335480064, 0, 0, 0, 0, -59238493/1068277825, 181606767/758867731,   561292985/797845732,   -1041891430/1371343529,  760417239/1151165299, 118820643/751138087, -528747749/2220607170,  1/4]';
b_7 = [ 13451932/455176623, 0, 0, 0, 0, -808719846/976000145, 1757004468/5645159321, 656045339/265891186,   -3867574721/1518517206,   465885868/322736535,  53011238/667516719, 2/45, 0]';

% step control power
stepPow = 1/8;

% input arguments check
if nargin < 5
    varargin={};
end
if nargin < 4
    options = [];
    if nargin < 3
        error('Not enough input arguments.  See ODE87.'); %#ok<ERTAG>
    end
end

% output parameter checks
haveoutfun = 1;
outfun = odeget(options,'OutputFcn',[],'fast');
if isempty(outfun)
    haveoutfun = 0;
end

%  relative error tolerance determination
tol=odeget(options,'RelTol');
if isempty(tol)
    tol = 1.e-6;
end

% housekeeping
nstep = 0; %#ok<NASGU>
t0 = tspan(1);
tfinal = tspan(2); %#ok<NASGU>
dt = tfinal - t0;
t = t0;
hmin = 16*eps*abs(t);   % Minimal step size
reject = 0;         % constant for step rejection

% maximal step size determination
hmax=odeget(options,'MaxStep');
if isempty(hmax)
    hmax =dt / 2.5;
end

% initial step size determination
h=odeget(options,'InitialStep');
if isempty(h)
    h = dt / 50;
    if h>0.1
        h=0.1;
    end
    if h>hmax
        h = hmax;
    end
end

% start point
x = x0(:);
f = x * zeros(1,13);  % right hand side calculations
tout = t;
xout = x.';
tau = tol * max(norm(x,'inf'), 1);  %#ok<NASGU> % accuracy

% Initial output
if haveoutfun
    feval(outfun,t,x,'init',varargin{:});
end

% iterative solution
while (t < tfinal) && (h >= hmin)
    if (t + h) > tfinal
        h = tfinal - t;
    end

    % Compute right-hand-side for step of method
    f(:,1) = feval(odefun,t,x, varargin{:});
    for j = 1: 12,
        f(:,j+1) = feval(odefun, t+c_i(j)*h, x+h*f*a_i_j(:,j), varargin{:});
    end

    % Two solution
    sol2=x+h*f*b_8;
    sol1=x+h*f*b_7;

    % Truncation error
    error_1 = norm(sol1-sol2);

    % estimate acceptable error
    Error_step = norm(error_1,'inf');
    tau = tol*max(norm(x,'inf'),1.0);

    % solution update only when error is acceptable
    if Error_step <= tau
        t = t + h;
        x = sol2;
        tout = [tout; t];
        xout = [xout; x.'];
        if haveoutfun
            status=feval(outfun,t,x,'',varargin{:});
            if status == 1
                return;
            end
        else
            t; %#ok<VUNUS>
            x.'; %#ok<VUNUS>
        end
        reject = reject - 1;
    else
        reject = 1;
    end

    % Step control
    if Error_step == 0.0
        Error_step = eps*10.0;
    end
    h = min(hmax, 0.9*h*(tau/Error_step)^stepPow);
    if (abs(h) <= eps)
        if reject == 0
            disp('Warning!!! ode87. Step is very small!!!');
            h = eps * 100; %#ok<NASGU>
            return
        else
            disp('Error in ode87. Step is too small.');
            return;
        end
    end
end

if (t < tfinal)
    disp('Error in ODE87...')
    t %#ok<NOPRT>
end

if haveoutfun
    feval(outfun,t,x,'done',varargin{:});
end
end
