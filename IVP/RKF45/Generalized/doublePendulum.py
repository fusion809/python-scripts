from matplotlib import pyplot as plt
import matplotlib
import numpy as np
import time
from RKF45 import RKF45
from matplotlib.animation import FuncAnimation
from scipy.interpolate import interp1d

class paramObj:
    def __init__(self, g, r1, r2, m1r, m1b, m2r, m2b, b1r, b1b, b2r, b2b,
                 c1r, c1b, c2r, c2b):
        self.g = g
        self.r1 = r1
        self.r2 = r2
        self.m1r = m1r
        self.m1b = m1b
        self.m2r = m2r
        self.m2b = m2b
        self.b1r = b1r
        self.b1b = b1b
        self.b2r = b2r
        self.b2b = b2b
        self.c1r = c1r
        self.c1b = c1b
        self.c2r = c2r
        self.c2b = c2b

def doubPen(params, t, vars, dt):
    g, r1, r2, m1r = params.g, params.r1, params.r2, params.m1r
    m1b, m2r, m2b, b1r = params.m1b, params.m2r, params.m2b, params.b1r
    b1b, b2r, b2b, c1r = params.b1b, params.b2r, params.b2b, params.c1r
    c1b, c2r, c2b = params.c1b, params.c2r, params.c2b
    theta1, thetaDot1, theta2, thetaDot2 = vars[0], vars[1], vars[2], vars[3]
    outercoef = 1/((m2r/12+m2b)*r2 - (m2b**2*r2*np.cos(theta1-theta2)**2)/(m1r/12 + m1b+m2b))
    mass = m1r/2 + m2r + m1b + m2b
    v1b = r1*thetaDot1
    v1r = v1b/2
    v2b = np.sqrt(r1**2*thetaDot1**2+r2**2*thetaDot2**2+2*r1*r2*thetaDot1*thetaDot2*np.cos(theta1-theta2))
    v2r = np.sqrt(r1**2*thetaDot1**2+(r2**2*thetaDot2**2)/4 + r2*thetaDot1*thetaDot2*np.cos(theta1-theta2))
    drag1b = (b1b + c1b*v1b)*v1b
    drag1r = (b1r + c1r*v1r)*v1r/2
    drag2b = (b2b + c2b*v2b)*(r1*thetaDot1 + r2*thetaDot2*np.cos(theta1-theta2))
    drag2b2 = (b2b + c2b*v2b)*(r1*thetaDot1*np.cos(theta1-theta2) + r2*thetaDot2)
    drag2r = (b2r + c2r*v2r)*(r1*thetaDot1 + r2*thetaDot2*np.cos(theta1-theta2)/2)
    drag2r2 = 1/4*(b2r + c2r*v2r)*(2*r1*thetaDot1*np.cos(theta1-theta2) + r2*thetaDot2)
    innercoef = -m2b*np.cos(theta1-theta2)/(m1r/12+m1b+m2b)
    inner = innercoef*(-m2b*r2*thetaDot2**2*np.sin(theta1-theta2) - g*np.cos(theta1)*(m1r/2+m2r+m1b+m2b)-drag1b-drag2b-drag1r-drag2r)
    extra = m2b*(r1*thetaDot1**2*np.sin(theta1-theta2)-g*np.cos(theta2))-drag2r2-drag2b2
    thetaDDot2 = outercoef*(inner+extra)
    outercoef1 = 1/((m1r/12+m1b+m2b)*r1)
    inner11 = -m2b*r2*(thetaDDot2*np.cos(theta1-theta2)+thetaDot2**2*np.sin(theta1-theta2))
    inner12 = -g*np.cos(theta1)*mass - drag1b - drag2b - drag1r - drag2r
    thetaDDot1 = outercoef1*(inner11 + inner12)
    return [dt * thetaDot1, dt * thetaDDot1, dt*thetaDot2, dt*thetaDDot2]

def animate_doubPen(params, theta1, theta2, t, N):
    tuni = np.linspace(t[0], t[len(t)-1], num=N)
    theta1int = interp1d(t, theta1)
    theta2int = interp1d(t, theta2)
    theta1uni = theta1int(tuni)
    theta2uni = theta2int(tuni)
    x1 = params.r1 * np.cos(theta1uni)
    x2 = x1 + params.r2 * np.cos(theta2uni)
    y1 = params.r1 * np.sin(theta1uni)
    y2 = y1 + params.r2 * np.sin(theta2uni)

    # Set up figure
    fig, ax = plt.subplots()
    ax.set_xlim(-2.2, 2.2)
    ax.set_ylim(-2.2, 2.2)
    ax.set_aspect('equal')
    ax.axis('off')

    line1, = ax.plot([], [], 'o-', lw=2, color='blue')
    line2, = ax.plot([], [], 'o-', lw=2, color='red')
    time_text = ax.text(0.05, 0.95, '', transform=ax.transAxes)

    # Init function
    def init():
        line1.set_data([], [])
        line2.set_data([], [])
        time_text.set_text('')
        return line1, line2, time_text

    # Animation function
    def update(i):
        # First rod: origin to (x1, y1)
        line1.set_data([0, x1[i]], [0, y1[i]])
        # Second rod: (x1, y1) to (x2, y2)
        line2.set_data([x1[i], x2[i]], [y1[i], y2[i]])
        time_text.set_text(f"t = {tuni[i]:.2f} s")
        return line1, line2, time_text

    # Create the animation
    ani = FuncAnimation(fig, update, frames=len(tuni), init_func=init, blit=False, interval=tuni[1]-tuni[0])
    plt.show()

def main():
    # Problem parameters
    t0 = 0.0
    tf = 100
    epsilon = 1e-8
    theta10, thetaDot10, theta20, thetaDot20 = 0, 0, 0, 0
    dtInitial = 1e-2       
    conds = np.array([[theta10, thetaDot10, theta20, thetaDot20]])    
    params = paramObj(g=9.81, r1=1.0, r2=1.0, m1r=1.0, m1b=1.0, m2r=1.0,
                      m2b=1.0, b1r=0.0, b1b=0.0, b2r=0.0, b2b=0.0, c1r=0.0,
                      c1b=0.0, c2r=0.0, c2b=0.0)
    t, vars, dt = RKF45(doubPen, t0, tf, dtInitial, params, conds, epsilon)
    dt = np.array(dt)
    theta1, thetaDot1, theta2, thetaDot2 = vars[:,0], vars[:,1], vars[:,2], vars[:,3]

    animate_doubPen(params, theta1, theta2, t, 1000)

if __name__ == "__main__":
    main()