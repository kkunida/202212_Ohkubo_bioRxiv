function x1 = discretizeStateFcn(x0,u,Ts,STATEFCN)

% Integrates the continuous model for one control interval using forward
% Euler formula

Nstep = 100;
xnew = x0;
dt = Ts/Nstep;
for i = 1:Nstep
    xnew = xnew + dt*STATEFCN(xnew,u);
end
x1 = xnew;

end

