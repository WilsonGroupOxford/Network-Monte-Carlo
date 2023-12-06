import numpy as np
import os
import argparse
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.colors import ListedColormap, Normalize
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
        with open(f"{self.prefix}_aux.dat", "r") as aux_file:
            self.n = int(aux_file.readline())
            self.geom = str(aux_file.readline())
            aux_file.readline()
            self.pb = np.array([float(a) for a in aux_file.readline().split()])
            self.rpb = np.array([float(a) for a in aux_file.readline().split()])

    def read_crd_file(self):
        self.crds = np.genfromtxt(f"{self.prefix}_crds.dat", dtype = float)
        self.crds[:, 0] -= self.pb[0] * np.round(self.crds[:, 0] * self.rpb[0])
        self.crds[:, 1] -= self.pb[1] * np.round(self.crds[:, 1] * self.rpb[1])

    def read_net_file(self):
        with open(f"{self.prefix}_net.dat", "r") as net_file:
            self.netCnxs = []
            for line in net_file:
                self.netCnxs.append(np.array([int(a) for a in line.split()]))

    def read_dual_file(self):
        with open(f"{self.prefix}_dual.dat", "r") as dual_file:
            self.dualCnxs = []
            for line in dual_file:
                self.dualCnxs.append(np.array([int(a) for a in line.split()]))
            self.avCnxs = np.average(np.array([r.size for r in self.dualCnxs]))

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
            origin_x = x[np.argmin(np.abs(x - x_cut))]
            origin_y = y[np.argmin(np.abs(y - y_cut))]
            x -= origin_x
            y -= origin_y
            x[x > x_cut] -= self.pb[0]
            y[y > y_cut] -= self.pb[1]
            x[x < -x_cut] += self.pb[0]
            y[y < -y_cut] += self.pb[1]
            x += origin_x
            y += origin_y

            self.ring_crds.append(np.asarray(list(zip(x, y))))

    def generate_ring_clusters(self, r_network):
        clst_type = 5
        min_clst_cnxs = 3
        min_aux_cnxs = 2
        active_clst = np.zeros(r_network.n, dtype = int)
        active_aux = np.zeros(r_network.n, dtype = int)
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

    def plot_rings(self, ax, lw=0.1, alpha=1.0, zorder=1, x_shift=0, y_shift=0):
        patches = []
        ring_colours = []
        x_shift *= self.pb[0]
        y_shift *= self.pb[1]
        for i in range(len(self.ring_crds)):
            ring = np.asarray(self.ring_crds[i])
            ring[:,0] += x_shift
            ring[:,1] += y_shift
            plt.scatter(ring[:,0], ring[:,1], zorder = 2, c = "k", alpha = alpha, s = 1)
            patches.append(Polygon(np.array(ring), closed = True))
            ring_colours.append(self.ring_colours[ring.shape[0]])
            ring[:, 0] -= x_shift
            ring[:, 1] -= y_shift
        ax.add_collection(PatchCollection(patches,
                                          facecolor = ring_colours,
                                          edgecolor = "k",
                                          linewidths = lw,
                                          alpha = alpha,
                                          zorder = zorder))
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
        plt.scatter(self.crds[:, 0], self.crds[:, 1], zorder=2, c="k", alpha=alpha, s=10)
        for i in range(len(self.netCnxs)):
            cnxs = self.netCnxs[i]
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
                plt.plot((self.crds[i, 0], self.crds[i, 0] + x),
                         (self.crds[i, 1], self.crds[i, 1] + y),
                         c="k",
                         lw=lw,
                         alpha=alpha,
                         zorder=1)
        self.crds[:, 0] -= x_shift
        self.crds[:, 1] -= y_shift
        return ax

    def plot_clusters(self, ax, lw=0.1, alpha=1.0, zorder=1, x_shift=0, y_shift=0):
        patches = []
        rainbow = matplotlib.colormaps.get_cmap("rainbow")
        random_generator = np.random.RandomState(0)
        cluster_colours = [rainbow(random_generator.uniform(0, 1)) for x in range(len(self.clusters))]
        cluster_colours.append("white")
        ring_colours = []
        x_shift *= self.pb[0]
        y_shift *= self.pb[1]
        for i in range(len(self.ring_crds)):
            ring = self.ring_crds[i]
            ring[:, 0] += x_shift
            ring[:, 1] += y_shift
            patches.append(Polygon(np.array(ring), True))
            ring_colours.append(cluster_colours[self.ring_cluster[i]])
            ring[:, 0] -= x_shift
            ring[:, 1] -= y_shift
        ax.add_collection(PatchCollection(patches,
                                          facecolor=ring_colours,
                                          edgecolor="k",
                                          linewidths=lw,
                                          alpha=alpha,
                                          zorder=zorder))
        return ax

    def init_ring_colours2(self, av_ring_size):
        map_lower = matplotlib.colormaps.get_cmap("Blues_r")
        map_upper = matplotlib.colormaps.get_cmap("Reds")
        map_mean  = matplotlib.colormaps.get_cmap("Greys")
        map_lower = ListedColormap(map_lower(np.arange(20, 100)))
        map_upper = ListedColormap(map_upper(np.arange(20, 100)))

        norm_lower = Normalize(vmin=av_ring_size - 3, vmax=av_ring_size)
        norm_upper = Normalize(vmin=av_ring_size, vmax=av_ring_size + 6)
        colour_mean = map_mean(50)
        self.ring_colours = []
        for i in range(340):
            if i < 3:
                self.ring_colours.append("white")
            elif np.abs(i - av_ring_size) < 1e-6:
                self.ring_colours.append(colour_mean)
            elif i < av_ring_size:
                self.ring_colours.append(map_lower(norm_lower(i)))
            else:
                self.ring_colours.append(map_upper(norm_upper(i)))

    def init_ring_colours(self):
        colormap_turquoise 	= matplotlib.colormaps.get_cmap("GnBu")
        colormap_greens 	    = matplotlib.colormaps.get_cmap("Greens")
        colormap_blues 		= matplotlib.colormaps.get_cmap("Blues")
        colormap_greys 		= matplotlib.colormaps.get_cmap("Greys")
        colormap_reds 		= matplotlib.colormaps.get_cmap("Reds")
        colormap_oranges 	= matplotlib.colormaps.get_cmap("YlOrBr")
        colormap_purples 	= matplotlib.colormaps.get_cmap("PuRd")
        colormap_pinks 		= matplotlib.colormaps.get_cmap("RdPu")

        turquoise = colormap_turquoise(140)
        green     = colormap_greens(100)
        blue      = colormap_blues(150)
        grey      = colormap_greys(90)
        red       = colormap_reds(105)
        orange    = colormap_oranges(100)
        purple    = colormap_purples(100)
        pink      = colormap_pinks(80)

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
        for i in range(30):
            self.ring_colours.append("black")
            

def export_plot(prefix, potential_path, file_format):
        if type(potential_path) == bool:
            save_path = prefix
        elif os.path.isdir(potential_path):
            save_path = os.path.join(potential_path, prefix)
        plt.savefig(f"{save_path}.{file_format}", format = file_format)

def main():
    parser = argparse.ArgumentParser(description = "Visualise NetMC Output Files")
    parser.add_argument("prefix", type = str, help = "Prefix of NetMC output files")
    parser.add_argument("-a", "--network_a", type = str, nargs = "*", choices = ["rings", "nodes", "periodic"],
                        metavar = "rings nodes periodic", help = "What to plot for network a", default = False)
    parser.add_argument("-b", "--network_b", type = str, nargs = "*", choices = ["rings", "nodes", "periodic"],
                        metavar = "rings nodes periodic", help = "What to plot for network b", default = False)
    save_arg_parser = parser.add_mutually_exclusive_group(required = False)
    save_arg_parser.add_argument("-s", "--save_as_png", type = str, metavar = "directory", help = "Save as png",
                                 required = False, nargs = "?", default = None, const = True)
    save_arg_parser.add_argument("-S", "--save_as_pdf", type = str, metavar = "directory", help = "Save as pdf",
                                 required = False, nargs = "?", default = None, const = True)
    args = parser.parse_args()
    
    network_parser = argparse.ArgumentParser(description = "Parse network arguments")
    network_parser.add_argument("-rings",    action = "store_true", default = False, required = False)
    network_parser.add_argument("-nodes",    action = "store_true", default = False, required = False)
    network_parser.add_argument("-periodic", action = "store_true", default = False, required = False)
    
    if args.network_a:
        args.network_a = [f"-{arg}" for arg in args.network_a]
        network_a_args = network_parser.parse_args(args.network_a)
    else:
        network_a_args = network_parser.parse_args([])
    network_a_rings    = network_a_args.rings
    network_a_nodes    = network_a_args.nodes
    network_a_periodic = network_a_args.periodic
    if args.network_b:
        args.network_b = [f"-{arg}" for arg in args.network_b]
        network_b_args = network_parser.parse_args(args.network_b)
    else:
        network_b_args = network_parser.parse_args([])
    network_b_rings    = network_b_args.rings
    network_b_nodes    = network_b_args.nodes
    network_b_periodic = network_b_args.periodic
        
    prefix = args.prefix

    network_a = Network(prefix + "_A")
    network_b = Network(prefix + "_B")

    network_a.generate_ring_crds(network_b.dualCnxs, network_b.avCnxs)
    network_b.generate_ring_crds(network_a.dualCnxs, network_a.avCnxs)

    params = {"figure.figsize": (6, 6)}
    pylab.rcParams.update(params)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_axis_off()

    if network_a_rings:
        ax = network_a.plot_rings(ax, lw=0.5)
        if network_a_periodic:
            for x in [-1, 0, 1]:
                for y in [-1, 0, 1]:
                    if abs(x) + abs(y) > 0:
                        ax = network_a.plot_rings(ax, lw=0.5, x_shift=x, y_shift=y, alpha=1.0)
    if network_a_nodes:
        ax = network_a.plot_connections(ax)
        if network_a_periodic:
            for x in [-1, 0, 1]:
                for y in [-1, 0, 1]:
                    if abs(x) + abs(y) > 0:
                        ax = network_a.plot_connections(ax, x_shift=x, y_shift=y, alpha=0.5)
# =============================================================================
#     if "c" in options_a:
#         network_a.generate_ring_clusters(network_b)
#         ax = network_a.plot_clusters(ax, lw=1.0)
#         if "p" in options_a:
#             for x in [-1,0,1]:
#                  for y in [-1,0,1]:
#                      if abs(x)+abs(y)>0:
#                          ax=network_a.plot_rings(ax,lw=1.0,x_shift=x,y_shift=y,alpha=0.5)
# =============================================================================
    if network_b_rings:
        ax = network_b.plot_rings(ax, lw=0.5)
        if network_b_periodic:
            for x in [-1, 0, 1]:
                for y in [-1, 0, 1]:
                    if abs(x) + abs(y) > 0:
                        ax = network_b.plot_rings(ax, lw=0.5, x_shift=x, y_shift=y, alpha=0.5)
    if network_b_nodes:
        ax = network_b.plot_connections(ax, alpha=1.0)
        if network_b_periodic:
            for x in [-1, 0, 1]:
                for y in [-1, 0, 1]:
                    if abs(x) + abs(y) > 0:
                        ax = network_b.plot_connections(ax, x_shift=x, y_shift=y, alpha=0.5)

    ax.set_xlim(np.min(network_a.crds[:, 0]), np.max(network_a.crds[:, 0]))
    ax.set_ylim(np.min(network_a.crds[:, 1]), np.max(network_a.crds[:, 1]))


    if args.save_as_png:
        export_plot(prefix, args.save_as_png, "png")
    elif args.save_as_pdf:
        export_plot(prefix, args.save_as_pdf, "pdf")
    plt.show()

if __name__ == "__main__":
    main()