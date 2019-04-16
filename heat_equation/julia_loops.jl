using LinearAlgebra

function heat(t,D,v,Nx,Nt)
    T = typeof(D)
    dt = t / Nt
    dx = 1 / Nx
    x = collect(0:1/(Nx-1):1)
    u = zeros(T, Nx)
    u_next = zeros(T, Nx)
    f(x) = sin(2pi*x) - sin(4pi*x)
    u .= f.(x)
    for it=1:Nt
        for ix=2:Nx-1
            u_next[ix] = (
                  u[ix]
                + D * dt / (dx * dx) * (u[ix+1] - 2u[ix] + u[ix-1])
                - v * dt / dx * (u[ix+1] - u[ix])
            )
        end
        u .= u_next
    end
    return u
end

function main()
    Nx = 1000
    Nt = 1000000
    D = 1.14;
    v = .85;
    t = 0.1;
    heat(t, D, v, Nx, Nt)
end
