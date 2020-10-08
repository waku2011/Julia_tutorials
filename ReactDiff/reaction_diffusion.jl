#
# 2-D simple reaction-diffusion system using cyclic conditions
# FitzHugh-Nagumo equation
#

using Plots
gr()

# computation area
lx = 1.0
ly = 1.0

# grid size for x,y-dimensions
mx = 64
my = 64 

# grid spacings
dx = lx/mx
dy = ly/mx
dv = dx*dy

x = zeros(Float64, mx)
y = zeros(Float64, my)
for i=1:mx
  x[i] = 0.5*dx + (i-1)*dx 
end
for j=1:my
  y[j] = 0.5*dy + (j-1)*dy
end

# time vars.
Time    = 0.0
dTime   = 0.01
timeMax =10.0
Nsteps=Int(timeMax/dTime)

# main vars.
u = zeros(Float64, mx, my)
unew = copy(u)
v = zeros(Float64, mx, my)
vnew = copy(v)

dflux_ux = zeros(Float64, mx+1, my)
dflux_uy = zeros(Float64, mx, my+1)

dflux_vx = zeros(Float64, mx+1, my)
dflux_vy = zeros(Float64, mx, my+1)

RHSU = zeros(Float64, mx, my)
RHSV = zeros(Float64, mx, my)

# parameters
Difu = 0.006
Difv = 0.07
a = 1.2 
b = 0.5
c = 0.0

# initializations
for i=1:mx, j=1:my
   u[i,j] = 2*rand()-1
   v[i,j] = 2*rand()-1
end

# functions
function set_dflux()
   for i=2:mx, j=1:my
      global dflux_ux[i,j] = -Difu*(u[i,j]-u[i-1,j])/dx
      global dflux_vx[i,j] = -Difv*(v[i,j]-v[i-1,j])/dx
   end
   for i=1:mx, j=2:my
      global dflux_uy[i,j] = -Difu*(u[i,j]-u[i,j-1])/dy
      global dflux_vy[i,j] = -Difv*(v[i,j]-v[i,j-1])/dy
   end
   # boundaries for x (W/E)
   for j=1:my
      global dflux_ux[1,j] = -Difu*(u[1,j]-u[mx,j])/dx
      global dflux_vx[1,j] = -Difv*(v[1,j]-v[mx,j])/dx
      global dflux_ux[mx+1,j] = dflux_ux[1,j]
      global dflux_vx[mx+1,j] = dflux_vx[1,j]
   end
   # boundaries for y (S/N)
   for i=1:mx
      global dflux_uy[i,1] = -Difu*(u[i,1]-u[i,my])/dy
      global dflux_vy[i,1] = -Difv*(v[i,1]-v[i,my])/dy
      global dflux_uy[i,my+1] = dflux_uy[i,1]
      global dflux_vy[i,my+1] = dflux_vy[i,1]
   end
end 

function set_RHS()
   for i=1:mx, j=1:my
     global RHSU[i,j] =-( dflux_ux[i+1,j]*dy
                  -dflux_ux[i  ,j]*dy
                  +dflux_uy[i,j+1]*dx
                  -dflux_uy[i  ,j]*dx)/dv
                  + u[i,j]-u[i,j]^3-v[i,j]
     global RHSV[i,j] =-( dflux_vx[i+1,j]*dy
                  -dflux_vx[i  ,j]*dy
                  +dflux_vy[i,j+1]*dx
                  -dflux_vy[i  ,j]*dx)/dv
                  +a*(u[i,j]-b*v[i,j]-c)
   end
end

# main
for iTime in 1:Nsteps
   global Time = Time + dTime
   #println(Time)
   
   set_dflux()
   
   set_RHS()
   
   global unew[:,:] = u[:,:] + RHSU[:,:]*dTime
   global vnew[:,:] = v[:,:] + RHSV[:,:]*dTime
   global u[:,:] = unew[:,:]
   global v[:,:] = vnew[:,:]
   
   p1=contour(x, y, u, title="U: camera 30 30 (default)", size=[800,400] );
   p2=contour(x, y, v, title="V: camera 30 30 (default)", size=[800,400] );
   plot(p1,p2,layout=2)
   
   sleep(2)
end

#savefig( "test.png" )

