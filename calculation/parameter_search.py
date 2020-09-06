from thepopulationdynamicsofmaternaleffectselfishgenes import *
import pdb

def search():
    gen = 1000
    _ = plt.rc('axes.spines', **{'bottom': True, 'left': True, 'right': False, 'top': False})
    fig = plt.figure(figsize=(10, 8))

    axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    xs = []
    ys = []
    cs = []
    for s in np.arange(-0.02,1.0001,0.03):
        for k in np.arange(0.01,1.0001,0.03):
            freq = one_simulation(mm=0, mo=0.1, oo=0.9, k=k, s=s, h=0.5, t=1, gen=gen)
            # if len(freq) < gen and freq[-1] > 0 and freq[-1] < 1:
            #     cs.append('blue')
            if freq[-1] > freq[0] and freq[-1] < 0.9999:
                cs.append('blue')
            else:
                cs.append('red')
            xs.append(s)
            ys.append(k)

    _ = plt.scatter(xs,ys,c=cs,s=15)
    _ = plt.xlabel('fitness cost s')
    _ = plt.ylabel('outcrossing rate k')
    # plt.legend()
    plt.ylim(0,1.02)
    plt.show()
    pdb.set_trace()
    # plt.savefig('parameter_space_balancing_selection_2.svg', dpi=300)


search()