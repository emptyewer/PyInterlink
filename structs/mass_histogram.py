import matplotlib as mpl
mpl.use("Qt4Agg")
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np
from PyQt4 import QtGui

class MplCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi, tight_layout=True)
        self.fig.patch.set_facecolor('white')
        self.fig.patch.set_alpha = 1.0
        self.axes = self.fig.add_subplot(111)
        self.scatter = None
        # We want the axes cleared every time plot() is called
        self.axes.hold(True)
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.setup_plot()

    def plot_figure(self, read):
        pass

    def setup_plot(self):
        pass

class View(MplCanvas):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        super(View, self).__init__(parent)
        self.annotations = []
        self.points = []
        self.visible_annotation = None

    def get_nearest_position(self, value):
        try:
            idx = (np.abs(np.array(self.points - value))).argmin()
            return idx
        except (TypeError, ValueError):
            return None

    def on_plot_hover(self, event):
        nearest = self.get_nearest_position(event.xdata)
        if nearest:
            annot = self.annotations[nearest]
            if annot != self.visible_annotation:
                if self.visible_annotation:
                    self.visible_annotation.set_visible(False)
                annot.set_visible(True)
                self.visible_annotation = annot
                self.fig.canvas.draw()

    def setup_plot(self):
        self.axes.tick_params(axis=u'both', which=u'both', length=0)
        self.axes.set_xlabel('Mass')
        self.axes.set_ylabel('Intensity')
        ticks_font = mpl.font_manager.FontProperties(family='times new roman', style='normal', size=14, weight='normal',
                                                     stretch='normal')

        labels = [self.axes.title, self.axes.xaxis.label, self.axes.yaxis.label]
        labels += self.axes.get_xticklabels() + self.axes.get_yticklabels()
        for item in labels:
            item.set_fontproperties(ticks_font)
            item.set_fontsize(16)
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_plot_hover)

    def plot_figure(self, read):
        self.axes.clear()
        self.annotations = []
        self.points = []
        masses = [x[0] for x in read.mass_int]
        intensities = [x[1] for x in read.mass_int]
        self.axes.bar(masses, intensities, width=1.0, linewidth=0)
        for i in range(len(read.mass_int)):
            annotation = self.axes.annotate("Mass: %.3f" % read.mass_int[i][0],
                                     xy=(read.mass_int[i][0], read.mass_int[i][1]), xycoords='data',
                                     xytext=(read.mass_int[i][0], read.mass_int[i][1]), textcoords='data',
                                     horizontalalignment="center",
                                     bbox=dict(boxstyle="round", facecolor="w",
                                               edgecolor="0.5", alpha=0.7)
                                     )
            annotation.set_visible(False)
            self.annotations.append(annotation)
            self.points.append(read.mass_int[i][0])
        self.scatter = None
        self.setup_plot()
        self.fig.canvas.draw()

    def plot_fragment(self, fragments):
        if self.scatter:
            self.scatter.remove()
        y_lim = self.axes.get_ylim()[1]
        x = [f.mass for f in fragments]
        y = [y_lim/2] * len(x)
        # print x, y
        self.scatter = self.axes.bar(x, y, width=1.0, linewidth=0, facecolor='r')
        self.axes.relim()
        self.axes.autoscale_view()
        self.fig.canvas.draw()


# class View(QtWidgets.QWidget):
#     def __init__(self, parent=None):
#         super(View, self).__init__(parent)
#         self.dpi = 300
#         self.fig = Figure((4.5, 4.5), dpi=self.dpi)
#         self.axes = self.fig.add_subplot(111)
#         self.canvas = FigureCanvas(self.fig)
#         self.canvas.setParent(self)
#         self.canvas.mpl_connect('button_press_event', self._onpick)
#         self.layout = QtWidgets.QVBoxLayout()
#         self.layout.addWidget(self.canvas)
#         self.layout.setStretchFactor(self.canvas, 1)
#         self.setLayout(self.layout)
#         # if self.pmap.use_ca:
#         #     self.xcoor = self.pmap.residue_numbers_ca[self.pmap.parent.current_model]
#         # else:
#         #     self.xcoor = self.pmap.residue_numbers_cb[self.pmap.parent.current_model]
#         # self.ycoor = self.pmap.histogram_maps[self.pmap.parent.current_model]
#         self.bar = self.axes.bar(1, 1, width=1.0, linewidth=0)
#         self.canvas.show()
#         self.set_parameters()
#
#     def _get_clicked_residues(self, event):
#         xmin, xmax = self.axes.get_xlim()
#         return(int(math.ceil(event.xdata - xmin))-1)
#
#     def _onpick(self, event):
#         if self.pmap.use_ca:
#             atom = 'CA'
#             chains = self.pmap.chain_names_ca[self.pmap.parent.current_model]
#             residue_numbers = self.pmap.residue_numbers_ca[self.pmap.parent.current_model]
#         else:
#             atom = 'CB'
#             chains = self.pmap.chain_names_cb[self.pmap.parent.current_model]
#             residue_numbers = self.pmap.residue_numbers_cb[self.pmap.parent.current_model]
#
#         index = self._get_clicked_residues(event)
#         self.pymol.select_density(residue_numbers[index], atom, self.pmap.cutoff, chains[index])
#
    def set_parameters(self):
        self.axes.tick_params(axis=u'both', which=u'both', length=0)
        # self.axes.set_xlim(min(self.xcoor), max(self.xcoor))
        # self.axes.set_ylim(0, max(self.ycoor) + 1)
        self.axes.set_xlabel('Mass')
        self.axes.set_ylabel('Intensity')
        # fractions = self.ycoor / max(self.ycoor)
        # normalized_colors = colors.Normalize(fractions.min(), fractions.max())
        # count = 0
        # for rect in self.bar:
        #     c = cm.jet(normalized_colors(fractions[count]))
        #     rect.set_facecolor(c)
        #     count += 1
        # self.fig.patch.set_facecolor((0.886, 0.886, 0.886))
        ticks_font = mpl.font_manager.FontProperties(family='times new roman', style='normal', size=12, weight='normal',
                                                     stretch='normal')
        labels = [self.axes.title, self.axes.xaxis.label, self.axes.yaxis.label]
        labels += self.axes.get_xticklabels() + self.axes.get_yticklabels()
        for item in labels:
            item.set_fontproperties(ticks_font)
            item.set_fontsize(5)

    def update_graph(self):
        self.axes.clear()
        if self.pmap.use_ca:
            self.xcoor = self.pmap.residue_numbers_ca[self.pmap.parent.current_model]
        else:
            self.xcoor = self.pmap.residue_numbers_cb[self.pmap.parent.current_model]
        self.ycoor = self.pmap.histogram_maps[self.pmap.parent.current_model]
        self.bar = self.axes.bar(self.xcoor, self.ycoor, width=1.0, linewidth=0)
        self.set_parameters()
        self.canvas.draw()
