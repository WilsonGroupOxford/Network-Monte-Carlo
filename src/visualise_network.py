import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, Normalize
import matplotlib.pylab as pylab


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
        self.crds[:, 0] -= self.pb[0] * np.round(self.crds[:, 0] * self.rpb[0])
        self.crds[:, 1] -= self.pb[1] * np.round(self.crds[:, 1] * self.rpb[1])
        if self.crds.shape[1] == 3:
            crds_3d = np.zeros_like(self.crds)
            crds_3dp = np.zeros_like(self.crds)
            crds_3d[:, :] = self.crds[:, :]
            angx = 0.0
            angy = 0.0
            angz = 0.0
            for i, c in enumerate(crds_3d):
                cc = np.zeros(3)
                cc[0] = c[0]
                cc[1] = np.cos(angz) * c[1] - np.sin(angz) * c[2]
                cc[2] = np.sin(angz) * c[1] + np.cos(angz) * c[2]
                crds_3d[i, :] = cc[:]
            for i, c in enumerate(crds_3d):
                cc = np.zeros(3)
                cc[1] = c[1]
                cc[0] = np.cos(angz) * c[0] + np.sin(angz) * c[2]
                cc[2] = -np.sin(angz) * c[0] + np.cos(angz) * c[2]
                crds_3d[i, :] = cc[:]
            for i, c in enumerate(crds_3d):
                cc = np.zeros(3)
                cc[0] = np.cos(angz) * c[0] - np.sin(angz) * c[1]
                cc[1] = np.sin(angz) * c[0] + np.cos(angz) * c[1]
                cc[2] = c[2]
                crds_3d[i, :] = cc[:]
            xy = crds_3d[:, 0] ** 2 + crds_3d[:, 1] ** 2
            crds_3dp[:, 0] = np.sqrt(xy + crds_3d[:, 2] ** 2)
            crds_3dp[:, 1] = np.arctan2(np.sqrt(xy), crds_3d[:, 2])
            crds_3dp[:, 2] = np.arctan2(crds_3d[:, 1], crds_3d[:, 0])
            print(np.min(crds_3dp[:, 1]), np.max(crds_3dp[:, 1]))
            self.crds = np.zeros((crds_3d.shape[0], 2))
            self.crds[:, 0] = (1 + np.abs(crds_3dp[:, 1])) ** 2 * np.cos(crds_3dp[:, 2])
            self.crds[:, 1] = (1 + np.abs(crds_3dp[:, 1])) ** 2 * np.sin(crds_3dp[:, 2])

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
        # self.avCnxs=5
        f.close()

    def generate_ring_crds(self, rings, av_ring_size):
        self.init_ring_colours2(av_ring_size)
        self.ring_crds = []
        x_cut = 0.5 * self.pb[0]
        y_cut = 0.5 * self.pb[1]
        for ring in rings:
            x = np.zeros(ring.size)
            y = np.zeros(ring.size)
            for i, node in enumerate(ring):
                x[i] = self.crds[node, 0]
                y[i] = self.crds[node, 1]
            origin_x = x[np.argmin(np.abs(x - x_cut))]  # +np.abs(y-y_cut))]
            origin_y = y[np.argmin(np.abs(y - y_cut))]  # +np.abs(y-y_cut))]
            x -= origin_x
            y -= origin_y
            x[x > x_cut] -= self.pb[0]
            y[y > y_cut] -= self.pb[1]
            x[x < -x_cut] += self.pb[0]
            y[y < -y_cut] += self.pb[1]
            x += origin_x
            y += origin_y
            self.ring_crds.append(np.array(list(zip(x, y))))
        self.ring_crds = np.array(self.ring_crds)

    def generate_ring_clusters(self, r_network):

        clst_type = 5
        min_clst_cnxs = 3
        min_aux_cnxs = 2
        active_clst = np.zeros(r_network.n, dtype=int)
        active_aux = np.zeros(r_network.n, dtype=int)
        for i, cnxs in enumerate(r_network.netCnxs):
            if len(cnxs) == clst_type:
                n = 0
                for j in cnxs:
                    if len(r_network.netCnxs[j]) == clst_type:
                        n += 1
                if n >= min_clst_cnxs:
                    active_clst[i] = 1
                    active_aux[i] = 0
                elif n >= min_aux_cnxs:
                    active_clst[i] = 0
                    active_aux[i] = 1
                else:
                    active_clst[i] = 0
                    active_aux[i] = 0
            else:
                active_clst[i] = 0
                active_aux[i] = 0
        clusters = []
        for i in range(r_network.n):
            if active_clst[i]:
                clst = set([])
                search0 = set([])
                search1 = set([])
                search1a = set([])
                clst.add(i)
                search0.add(i)
                while True:
                    for j in search0:
                        ids = r_network.netCnxs[j]
                        for k in ids:
                            if active_clst[k]:
                                search1.add(k)
                            elif active_aux[k]:
                                search1a.add(k)
                    if len(search1) == 0 and len(search1a) == 0:
                        break
                    for j in clst:
                        if j in search1:
                            search1.remove(j)
                        if j in search1a:
                            search1a.remove(j)
                    search0 = search1.copy()
                    for j in search1:
                        clst.add(j)
                    for j in search1a:
                        clst.add(j)
                    search1 = set([])
                    search1a = set([])
                    if len(search0) == 0:
                        break
                for j in clst:
                    active_clst[j] = 0
                    active_aux[j] = 0
                clusters.append(np.array(list(clst)))
        self.clusters = clusters
        self.ring_cluster = np.zeros(r_network.n, dtype=int)
        self.ring_cluster[:] = -1
        for i, clst in enumerate(self.clusters):
            for j in clst:
                self.ring_cluster[j] = i
        print(np.max(np.array([len(x) for x in self.clusters])))

    def plot_rings(self, ax, lw=0.1, alpha=1.0, zorder=1, x_shift=0, y_shift=0):

        patches = []
        ring_colours = []
        x_shift *= self.pb[0]
        y_shift *= self.pb[1]
        order = np.argsort(np.array([np.max(np.abs(c)) for c in self.ring_crds]))[::-1]
        for i, ring in enumerate(self.ring_crds[order]):
            ring[:, 0] += x_shift
            ring[:, 1] += y_shift
            plt.scatter(ring[:, 0], ring[:, 1], zorder=2, c="k", alpha=alpha, s=1)
            patches.append(Polygon(np.array(ring), True))
            ring_colours.append(self.ring_colours[ring[:, 0].size])
            ring[:, 0] -= x_shift
            ring[:, 1] -= y_shift
        ax.add_collection(
            PatchCollection(
                patches,
                facecolor=ring_colours,
                edgecolor="k",
                linewidths=lw,
                alpha=alpha,
                zorder=zorder,
            )
        )
        ax.set_xlim(-self.pb[0] * 0.5, self.pb[0] * 1.5)
        ax.set_ylim(-self.pb[0] * 0.5, self.pb[0] * 1.5)
        return ax

    def plot_connections(self, ax, lw=1, alpha=1, x_shift=0, y_shift=0):

        x_shift *= self.pb[0]
        y_shift *= self.pb[1]
        x_cut = 0.5 * self.pb[0]
        y_cut = 0.5 * self.pb[1]
        self.crds[:, 0] += x_shift
        self.crds[:, 1] += y_shift
        # mask_3 = np.array([self.netCnxs[i].size==3 for i in range(self.netCnxs.size)])
        # mask_4 = np.array([self.netCnxs[i].size==4 for i in range(self.netCnxs.size)])
        plt.scatter(
            self.crds[:, 0], self.crds[:, 1], zorder=2, c="k", alpha=alpha, s=10
        )
        # plt.scatter(self.crds[mask_3,0],self.crds[mask_3,1],zorder=2,c="darkorange",alpha=alpha,s=10)
        # plt.scatter(self.crds[mask_4,0],self.crds[mask_4,1],zorder=2,c="forestgreen",alpha=alpha,s=10)
        for i, cnxs in enumerate(self.netCnxs):
            # plt.text(self.crds[i,0],self.crds[i,1],i)
            for j in cnxs:
                x = self.crds[j, 0] - self.crds[i, 0]
                y = self.crds[j, 1] - self.crds[i, 1]
                if x > x_cut:
                    x -= self.pb[0]
                elif x < -x_cut:
                    x += self.pb[0]
                if y > y_cut:
                    y -= self.pb[1]
                elif y < -y_cut:
                    y += self.pb[1]
                plt.plot(
                    (self.crds[i, 0], self.crds[i, 0] + x),
                    (self.crds[i, 1], self.crds[i, 1] + y),
                    c="k",
                    lw=lw,
                    alpha=alpha,
                    zorder=1,
                )
        self.crds[:, 0] -= x_shift
        self.crds[:, 1] -= y_shift
        return ax

    def plot_clusters(self, ax, lw=0.1, alpha=1.0, zorder=1, x_shift=0, y_shift=0):

        patches = []
        rainbow = cm.get_cmap("rainbow")
        random_generator = np.random.RandomState(0)
        cluster_colours = [
            rainbow(random_generator.uniform(0, 1)) for x in range(len(self.clusters))
        ]
        cluster_colours.append("white")
        ring_colours = []
        x_shift *= self.pb[0]
        y_shift *= self.pb[1]
        order = np.argsort(np.array([np.max(np.abs(c)) for c in self.ring_crds]))[::-1]
        for i, ring in enumerate(self.ring_crds):
            ring[:, 0] += x_shift
            ring[:, 1] += y_shift
            patches.append(Polygon(np.array(ring), True))
            ring_colours.append(cluster_colours[self.ring_cluster[i]])
            ring[:, 0] -= x_shift
            ring[:, 1] -= y_shift
        ax.add_collection(
            PatchCollection(
                patches,
                facecolor=ring_colours,
                edgecolor="k",
                linewidths=lw,
                alpha=alpha,
                zorder=zorder,
            )
        )

        return ax

    def init_ring_colours2(self, av_ring_size):
        map_lower = cm.get_cmap("Blues_r", 128)
        map_upper = cm.get_cmap("Reds", 128)
        map_mean = cm.get_cmap("Greys")
        map_lower = ListedColormap(map_lower(np.arange(20, 100)))
        map_upper = ListedColormap(map_upper(np.arange(20, 100)))

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

        colormap_turquoise = cm.get_cmap("GnBu")
        colormap_greens = cm.get_cmap("Greens")
        colormap_blues = cm.get_cmap("Blues")
        colormap_greys = cm.get_cmap("Greys")
        colormap_reds = cm.get_cmap("Reds")
        colormap_oranges = cm.get_cmap("YlOrBr")
        colormap_purples = cm.get_cmap("PuRd")
        colormap_pinks = cm.get_cmap("RdPu")

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
        if "r" in flags[2]:
            options_a.append("r")
        if "n" in flags[2]:
            options_a.append("n")
        if "p" in flags[2]:
            options_a.append("p")
        if "a" in flags[2]:
            options_a.append("a")
        if "b" in flags[2]:
            options_b.append("b")
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
    ax = fig.add_subplot(111)
    ax.set_axis_off()

    if "r" in options_a:
        ax = network_a.plot_rings(ax, lw=0.5)
        if "p" in options_a:
            for x in [-1, 0, 1]:
                for y in [-1, 0, 1]:
                    if abs(x) + abs(y) > 0:
                        ax = network_a.plot_rings(
                            ax, lw=0.5, x_shift=x, y_shift=y, alpha=1.0
                        )
    if "n" in options_a:
        ax = network_a.plot_connections(ax)
        if "p" in options_a:
            for x in [-1, 0, 1]:
                for y in [-1, 0, 1]:
                    if abs(x) + abs(y) > 0:
                        ax = network_a.plot_connections(
                            ax, x_shift=x, y_shift=y, alpha=0.5
                        )
    if "c" in options_a:
        network_a.generate_ring_clusters(network_b)
        ax = network_a.plot_clusters(ax, lw=1.0)
        # if "p" in options_a:
        #     for x in [-1,0,1]:
        #         for y in [-1,0,1]:
        #             if abs(x)+abs(y)>0:
        #                 ax=network_a.plot_rings(ax,lw=1.0,x_shift=x,y_shift=y,alpha=0.5)

    if "r" in options_b:
        ax = network_b.plot_rings(ax, lw=0.5)
        if "p" in options_b:
            for x in [-1, 0, 1]:
                for y in [-1, 0, 1]:
                    if abs(x) + abs(y) > 0:
                        ax = network_b.plot_rings(
                            ax, lw=0.5, x_shift=x, y_shift=y, alpha=0.5
                        )
    if "n" in options_b:
        ax = network_b.plot_connections(ax, alpha=1.0)
        if "p" in options_b:
            for x in [-1, 0, 1]:
                for y in [-1, 0, 1]:
                    if abs(x) + abs(y) > 0:
                        ax = network_b.plot_connections(
                            ax, x_shift=x, y_shift=y, alpha=0.5
                        )

    ax.set_xlim(np.min(network_a.crds[:, 0]), np.max(network_a.crds[:, 0]))
    ax.set_ylim(np.min(network_a.crds[:, 1]), np.min(network_a.crds[:, 1]))

    if "s" in options_s:
        plt.savefig("{:}.png".format(prefix), dpi=400)
    if "S" in options_s:
        plt.savefig("{:}.pdf".format(prefix))

    plt.show()


main()
