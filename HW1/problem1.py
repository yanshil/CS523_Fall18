## CS523.HW1.problem1
## Yanshi Luo

## Code from jupyter notebook, plots are displayed in the pdf

## Backward Euler
def backwardEuler(u,v,h):
    for i in range(1,100):
    # Reference: https://stackoverflow.com/questions/22742951/solve-an-equation-using-a-python-numerical-solver-in-numpy
        func = lambda tau : tau - u[i-1] - h * u[i-1] * ( (v[i-1] / (1- h * (1-tau))) - 2)

        tau = np.linspace(-10, 10, 201)
        tau_initial_guess = 0.5

        u_new = fsolve(func, tau_initial_guess)

        v_new = v[i-1] / (1 - h * (1 - u_new))

        u.append(u_new)
        v.append(v_new)

    plt.plot(u,v,'ro')
    plt.axis([0,6,0,10])
    plt.title('u0=%s, v0=%s, h=%s' % (u[0], v[0], h))
    plt.show()
    return;

backwardEuler([4],[8],.12)


## Sympletic Euler
def sympleticEuler(u,v,h):
    for i in range(1,100):
        v_new = v[i-1] + h * v[i-1] * (1.0 -u[i-1])
        u_new = u[i-1] + h * u[i-1] * (v_new - 2)

        u.append(u_new)
        v.append(v_new)

    plt.plot(u,v,'ro')
    plt.axis([0,6,0,10])
    plt.title('u0=%s, v0=%s, h=%s' % (u[0], v[0], h))
    plt.show()
    return;

sympleticEuler([4],[2],.12)
sympleticEuler([6],[2],.12)

