import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, Normalize
import matplotlib.pyplot as plt


class Network:
    def __init__(self, prefix):
        self.prefix = prefix
        self.read_aux_file()
        self.read_crd_file()
        self.read_net_file()
        self.read_dual_file()
        self.init_ring_colours()

    def read_aux_file(self):
        f = open(self.prefix + "_aux.dat")
        self.n = int(f.readline())
        self.geom = str(f.readline())
        f.readline()
        self.pb = np.array([float(a) for a in f.readline().split()])
        self.rpb = np.array([float(a) for a in f.readline().split()])
        f.close()

    def read_crd_file(self):
        self.crds = np.genfromtxt(self.prefix + "_crds.dat", dtype=float)

    def read_net_file(self):
        f = open(self.prefix + "_net.dat")
        self.netCnxs = []
        for line in f:
            self.netCnxs.append(np.array([int(a) for a in line.split()]))
        self.netCnxs = np.array(self.netCnxs)
        f.close()

    def read_dual_file(self):
        f = open(self.prefix + "_dual.dat")
        self.dualCnxs = []
        for line in f:
            self.dualCnxs.append(np.array([int(a) for a in line.split()]))
        self.dualCnxs = np.array(self.dualCnxs)
        self.avCnxs = np.average(np.array([r.size for r in self.dualCnxs]))
        if self.avCnxs > 4.5:
            self.avCnxs = 6
        elif self.avCnxs > 3.5:
            self.avCnxs = 4
        else:
            self.avCnxs = 3
        f.close()

    def generate_ring_crds(self, rings, av_ring_size):
        self.init_ring_colours2(av_ring_size)
        self.ring_crds = []
        for ring in rings:
            ring_crds = []
            for i, node in enumerate(ring):
                ring_crds.append(self.crds[node, :])
            self.ring_crds.append(np.array(ring_crds))

    def plot_rings(self, ax, lw=0.1, zorder=1):
        # p=[]
        # p.append(self.crds[0,:])
        # p.append(self.crds[1,:])
        # p.append(self.crds[9,:])
        # p.append(self.crds[10,:])
        # print p
        # ax.add_collection3d(Poly3DCollection([self.ring_crds[0]]))
        #
        # ring_colours = []
        lim = 0.0
        for i, ring in enumerate(self.ring_crds):
            if np.max(self.ring_crds[i]) > lim:
                lim = np.max(self.ring_crds[i])

            crds_3d = np.zeros_like(self.ring_crds[i])
            crds_3dp = np.zeros_like(self.ring_crds[i])
            crds_3d[:, :] = self.ring_crds[i][:, :]
            xy = crds_3d[:, 0] ** 2 + crds_3d[:, 1] ** 2
            crds_3dp[:, 0] = np.sqrt(xy + crds_3d[:, 2] ** 2)
            crds_3dp[:, 1] = np.arctan2(np.sqrt(xy), crds_3d[:, 2])
            crds_3dp[:, 2] = np.arctan2(crds_3d[:, 1], crds_3d[:, 0])
            delta = (np.max(crds_3dp[:, 2]) + np.pi - np.pi * -0.0) % (2 * np.pi)
            # if delta<np.pi:
            #    alpha=0
            # else:
            alpha = (delta) / (2 * np.pi)
            alpha = 1
            ax.add_collection3d(
                Poly3DCollection(
                    [self.ring_crds[i]],
                    alpha=alpha,
                    facecolor=self.ring_colours[ring[:, 0].size],
                    edgecolor="k",
                    linewidth=lw,
                )
            )
        #             ring_colours.append(self.ring_colours[ring[:, 0].size])
        #         ax.add_collection(PatchCollection(patches, facecolor=ring_colours, edgecolor='k', linewidths=lw, alpha=alpha, zorder=zorder))
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_zlim(-lim, lim)
        return ax

    def plot_connections(self, ax, lw=1):
        ax.scatter(self.crds[:, 0], self.crds[:, 1], self.crds[:, 2], zorder=2, c="k")
        for i, cnxs in enumerate(self.netCnxs):
            ax.text(self.crds[i, 0], self.crds[i, 1], self.crds[i, 2], i, zorder=2)
            for j in cnxs:
                ax.plot(
                    xs=(self.crds[i, 0], self.crds[j, 0]),
                    ys=(self.crds[i, 1], self.crds[j, 1]),
                    zs=(self.crds[i, 2], self.crds[j, 2]),
                    c="k",
                    lw=lw,
                )
        #         ax.set_xlim(-self.pb[0]*0.2,self.pb[0]*1.2)
        #         ax.set_ylim(-self.pb[0]*0.2,self.pb[0]*1.2)
        return ax

    def init_ring_colours2(self, av_ring_size):
        map_lower = cm.get_cmap("Blues_r", 128)
        map_upper = cm.get_cmap("Reds", 128)
        map_mean = cm.get_cmap("Greys")
        map_lower = ListedColormap(map_lower(np.arange(30, 100)))
        map_upper = ListedColormap(map_upper(np.arange(30, 100)))

        norm_lower = Normalize(vmin=av_ring_size - 3, vmax=av_ring_size)
        norm_upper = Normalize(vmin=av_ring_size, vmax=av_ring_size + 6)
        colour_mean = map_mean(50)

        # colour_mean='ivory'
        # map = np.vstack((map_upper(np.linspace(0, 1, 128)),map_lower(np.linspace(0, 1, 128))))
        # colormap = ListedColormap(map, name='custom')
        # av_ring_size=6
        self.ring_colours = []
        for i in range(30):
            if i < 3:
                self.ring_colours.append("white")
            elif np.abs(i - av_ring_size) < 1e-6:
                self.ring_colours.append(colour_mean)
            elif i < av_ring_size:
                self.ring_colours.append(map_lower(norm_lower(i)))
            else:
                self.ring_colours.append(map_upper(norm_upper(i)))

    def init_ring_colours(self):

        colormap_turquoise = plt.cm.get_cmap("GnBu")
        colormap_greens = plt.cm.get_cmap("Greens")
        colormap_blues = plt.cm.get_cmap("Blues")
        colormap_greys = plt.cm.get_cmap("Greys")
        colormap_reds = plt.cm.get_cmap("Reds")
        colormap_oranges = plt.cm.get_cmap("YlOrBr")
        colormap_purples = plt.cm.get_cmap("PuRd")
        colormap_pinks = plt.cm.get_cmap("RdPu")

        turquoise = colormap_turquoise(140)
        green = colormap_greens(100)
        blue = colormap_blues(150)
        grey = colormap_greys(90)
        red = colormap_reds(105)
        orange = colormap_oranges(100)
        purple = colormap_purples(100)
        pink = colormap_pinks(80)

        self.ring_colours = []
        for i in range(3):
            self.ring_colours.append("white")
        self.ring_colours.append(turquoise)
        self.ring_colours.append(green)
        self.ring_colours.append(blue)
        self.ring_colours.append(grey)
        self.ring_colours.append(red)
        self.ring_colours.append(orange)
        self.ring_colours.append(purple)
        self.ring_colours.append(pink)
        for i in range(20):
            self.ring_colours.append("black")


def options(flags):

    options_a = []
    options_b = []
    options_s = []

    if len(flags) == 2:
        options_a.append("r")

    if len(flags) >= 3:
        if "a" in flags[2]:
            options_a = flags[2]
        if "b" in flags[2]:
            options_b = flags[2]
        if "s" in flags[2]:
            options_s.append("s")
        if "S" in flags[2]:
            options_s.append("S")
    if len(flags) == 4:
        if "a" in flags[3]:
            options_a = flags[3]
        if "b" in flags[3]:
            options_b = flags[3]
    return options_a, options_b, options_s


def main():

    prefix = sys.argv[1]
    options_a, options_b, options_s = options(sys.argv)

    network_a = Network(prefix + "_A")
    network_b = Network(prefix + "_B")

    network_a.generate_ring_crds(network_b.dualCnxs, network_b.avCnxs)
    network_b.generate_ring_crds(network_a.dualCnxs, network_a.avCnxs)

    params = {"figure.figsize": (6, 6)}
    pylab.rcParams.update(params)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.set_axis_off()

    if "r" in options_a:
        ax = network_a.plot_rings(ax, lw=1.0)
    if "n" in options_a:
        ax = network_a.plot_connections(ax)

    if "r" in options_b:
        ax = network_b.plot_rings(ax, lw=0.5)
    if "n" in options_b:
        ax = network_b.plot_connections(ax)

    ax.view_init(14, 170)
    if "s" in options_s:
        plt.savefig("./{:}.png".format(prefix), dpi=400)
    if "S" in options_s:
        plt.savefig("./{:}.pdf".format(prefix))
    plt.show()


main()
