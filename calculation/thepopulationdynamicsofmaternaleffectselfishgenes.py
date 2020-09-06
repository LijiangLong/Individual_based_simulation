import numpy as np
import pdb
import matplotlib.pyplot as plt
def three_genetypes_onegeneration(mm,mo,oo,s=0.1,h=1,t=1):
    new_mm,new_mo,new_oo =  0,0,0
    # family 1
    new_mm += (1-s)*mm*mm
    # family 2
    new_mm += 0.5 * (1-s) * mo * mm
    new_mo += 0.5 * (1-s) * mo * mm
    # family 3
    new_mo += (1-s) * oo * mm
    # family 4
    new_mm += 0.5 * (1-h*s) * mm * mo
    new_mo += 0.5 * (1-h*s) * mm * mo
    # family 5
    new_mm += 0.25 * (1-h*s) * mo * mo
    new_mo += 0.5 * (1-h*s) * mo * mo
    new_oo += 0.25 * (1-h*s) * mo * mo * (1-t)
    # family 6
    new_mo += 0.5 * (1-h*s) * oo * mo
    new_oo += 0.5 * (1-h*s) * oo * mo * (1-t)
    # family 7
    new_mo += mm * oo
    # family 8
    new_mo += 0.5 * mo * oo
    new_oo += 0.5 * mo * oo
    # family 9
    new_oo += oo * oo

    new_total = new_mm + new_mo + new_oo
    return new_mm/new_total,new_mo/new_total,new_oo/new_total

def partial_selfing_maternal(mm,mo,oo,k=0.5,s=0.1,h=1,t=1):
    new_mm, new_mo, new_oo = 0, 0, 0
    # mm selfing
    new_mm += (1-k) * (1-s) * mm
    # mo selfing
    new_mm += (1-k) * (1-h*s) * mo * 0.25
    new_mo += (1-k) * (1-h*s) * mo * 0.5
    # oo selfing
    new_oo += (1-k) * oo
    # family 1
    new_mm += (1-s)*mm*mm * k
    # family 2
    new_mm += 0.5 * (1-s) * mo * mm* k
    new_mo += 0.5 * (1-s) * mo * mm* k
    # family 3
    new_mo += (1-s) * oo * mm* k
    # family 4
    new_mm += 0.5 * (1-h*s) * mm * mo* k
    new_mo += 0.5 * (1-h*s) * mm * mo* k
    # family 5
    new_mm += 0.25 * (1-h*s) * mo * mo* k
    new_mo += 0.5 * (1-h*s) * mo * mo* k
    new_oo += 0.25 * (1-h*s) * mo * mo * (1-t)* k
    # family 6
    new_mo += 0.5 * (1-h*s) * oo * mo* k
    new_oo += 0.5 * (1-h*s) * oo * mo * (1-t)* k
    # family 7
    new_mo += mm * oo* k
    # family 8
    new_mo += 0.5 * mo * oo* k
    new_oo += 0.5 * mo * oo* k
    # family 9
    new_oo += oo * oo* k

    new_total = new_mm + new_mo + new_oo
    return new_mm/new_total,new_mo/new_total,new_oo/new_total


def partial_selfing_paternal(mm,mo,oo,k=0.5,s=0.1,h=1,t=1):
    new_mm, new_mo, new_oo = 0, 0, 0
    # mm selfing
    new_mm += (1-k) * (1-s) * mm
    # mo selfing
    new_mm += (1-k) * (1-h*s) * mo * 0.25
    new_mo += (1-k) * (1-h*s) * mo * 0.5
    # oo selfing
    new_oo += (1-k) * oo
    # family 1
    new_mm += (1-s)*mm*mm * k
    # family 2
    new_mm += 0.5 * (1-s) * mo * mm* k
    new_mo += 0.5 * (1-s) * mo * mm* k
    # family 3
    new_mo += (1-s) * oo * mm* k
    # family 4
    new_mm += 0.5 * (1-h*s) * mm * mo* k
    new_mo += 0.5 * (1-h*s) * mm * mo* k
    # family 5
    new_mm += 0.25 * (1-h*s) * mo * mo* k
    new_mo += 0.5 * (1-h*s) * mo * mo* k
    new_oo += 0.25 * (1-h*s) * mo * mo * (1-t)* k
    # family 6
    new_mo += 0.5 * (1-h*s) * oo * mo* k
    new_oo += 0.5 * (1-h*s) * oo * mo * k
    # family 7
    new_mo += mm * oo* k
    # family 8
    new_mo += 0.5 * mo * oo* k
    new_oo += 0.5 * mo * oo* k* (1-t)
    # family 9
    new_oo += oo * oo* k

    new_total = new_mm + new_mo + new_oo
    return new_mm/new_total,new_mo/new_total,new_oo/new_total
def allele_freq(mm,mo,oo):
    return mm+mo/2

def one_simulation(mm=0,mo=0.1,oo=0.9,k=0.9, s=0.75, h=0, t=1,gen=100):
    freq = []

    for i in range(gen):
        mm,mo,oo = partial_selfing_paternal(mm, mo, oo, k, s, h, t)
        af = allele_freq(mm,mo,oo)
        freq.append(af)
        if af == 1 or af == 0:
            return freq
        try:
            if freq[-1] == freq[-2] and freq[-2] == freq[-3] and freq[-3] == freq[-4]:
                return freq
        except:
            pass
    return freq
# def main():
#     gen = 1000
#     _ = plt.rc('axes.spines', **{'bottom': True, 'left': True, 'right': False, 'top': False})
#     fig = plt.figure(figsize=(10, 8))
#
#     axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])
#
#     # for k in np.arange(0,1.1,0.1):
#     for k in [1,0.7,0.4,0.1]:
#         freq = one_simulation(mm=0, mo=0.1, oo=0.9, k=k, s=0.05, h=0, t=1, gen=gen)
#         axes.plot(range(len(freq)),freq, label='outcrossing rate k={k:.2f}'.format(k=k))
#
#     plt.legend()
#     plt.ylim(0,1)
#     plt.xlim(0,gen)
#     plt.plot((0, gen), (1, 1), linestyle='-.')
#     plt.show()
#     pdb.set_trace()

    # plt.savefig('vark_2.svg', dpi=300)

# def main():
#     gen = 2000
#     _ = plt.rc('axes.spines', **{'bottom': True, 'left': True, 'right': False, 'top': False})
#     fig = plt.figure(figsize=(11, 7))
#
#     axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])
#
#     for s in np.arange(0,0.1,0.01):
#         freq = one_simulation(mm=0, mo=0.1, oo=0.9, k=1, s=s, h=0, t=1, gen=gen)
#         axes.plot(range(len(freq)),freq, label='female fecundity loss s={s:.2f}'.format(s=s))
#
#     plt.legend()
#     plt.plot((0,gen),(1,1),linestyle='-.')
#     plt.ylim(0.9,1)
#     plt.xlim(0,gen)
#     plt.show()
#     pdb.set_trace()
    # plt.savefig('vars_3.svg', dpi=300)

def main():
    gen = 200
    _ = plt.rc('axes.spines', **{'bottom': True, 'left': True, 'right': False, 'top': False})
    fig = plt.figure(figsize=(11, 7))

    axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])

    for k in [1, 0.5]:
        for s in [0,-0.2]:
            freq = one_simulation(mm=0, mo=0.1, oo=0.9, k=k, s=s, h=0, t=1, gen=gen)
            axes.plot(range(len(freq)),freq, label='s={s:.1f},k={k:.1f}'.format(s=s,k=k))

    plt.legend()
    plt.plot((0,gen),(1,1),linestyle='-.')
    plt.ylim(0,1)
    plt.xlim(0,gen)
    # plt.show()
    # pdb.set_trace()
    plt.savefig('fixeds,vark.svg', dpi=300)


# def main():
#     gen = 100
#     _ = plt.rc('axes.spines', **{'bottom': True, 'left': True, 'right': False, 'top': False})
#     fig = plt.figure(figsize=(11, 7))
#
#     axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])
#
#     freq = one_simulation(mm=0, mo=0.1, oo=0.9, k=0.9, s=0, h=0, t=1, gen=gen)
#     axes.plot(range(len(freq)),freq)
#
#     # plt.legend()
#     plt.ylim(0,1)
#     # pdb.set_trace()
#     plt.savefig('s=0.svg', dpi=300)
# def main():
#     gen = 500
#     _ = plt.rc('axes.spines', **{'bottom': True, 'left': True, 'right': False, 'top': False})
#     fig = plt.figure(figsize=(10, 8))
#
#     axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])
#
#     for s in np.arange(0,1.1,0.4):
#         for k in [0.1,0.9]:
#             freq = one_simulation(mm=0, mo=0.1, oo=0.9, k=k, s=s, h=0, t=1, gen=gen)
#             axes.plot(range(len(freq)),freq, label='k={k:.2f},s={s:.2f}'.format(k=k,s=s))
#
#     plt.legend()
#     plt.ylim(0,1)
#     plt.plot((0, gen), (1, 1), linestyle='-.')
#     plt.show()
#     pdb.set_trace()

if __name__ == "__main__":
    main()



