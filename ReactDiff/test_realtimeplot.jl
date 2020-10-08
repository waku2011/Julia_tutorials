using Plots
gr()

t  = 1:0.25:25
x1 = sin.(t./pi)
x2 = rand(25) .+ 2

p1 = scatter(t, x1)
p2 = scatter(t, x2)

display(plot(p1, p2, layout=(2, 1)))

for i in 1:100
    
    global x1 = sin.((t.+i)/pi)
    global x2 = rand(25) .+ 2

    global p1 = scatter(t, x1)
    global p2 = scatter(t, x2)

    display(plot(p1, p2, layout=(2, 1)))
    
    sleep(0.1)
end

readline()

