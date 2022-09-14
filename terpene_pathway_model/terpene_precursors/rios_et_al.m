
tspan = [0 200];
x0 = 2*ones(22,1);
fun = @supplementary;
[t, c] = ode23s(fun, tspan, x0);
figure();
plot(t, c)