from PyQt4 import QtGui, uic, QtCore
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar2QT
from settings import *
from matplotlib import lines
from matplotlib.widgets import SpanSelector
from scipy.signal import savgol_filter
import numpy as np
import csv

class PointBrowser(object):
    def __init__(self,parent=None, toolbar=None):
        self.cur_line = None
        self.marker = None
        self.parent=parent
        self.toolbar=toolbar
        
        self.configurator=None

        self.baseline=[]

        self.baseline_regions_x=[]
        self.baseline_regions_y=[]
        self.span=None

    def resize(self, event):
        try:
            if self.parent.figure.get_axes():
                self.parent.figure.tight_layout()
                event.canvas.draw()
        except:
            pass
            
    def onpick(self, event):
        ax=event.artist.axes
        ax.hold(True)
        if event.mouseevent.button == 1:
            if event.mouseevent.dblclick:
                if isinstance(event.artist, lines.Line2D):
                    if self.marker is not None:
                        self.configurator=Graph_Configurator(self, ax)
            else:
                if isinstance(event.artist, lines.Line2D):
                    if self.marker is not None and self.cur_line is event.artist:
                        self.marker.remove()
                        self.marker = None
                    elif self.marker is not None and self.cur_line is not event.artist and self.marker in ax.lines:
                        self.marker.remove()
                        self.marker = None
                        self.cur_line = event.artist
                        ax.autoscale_view(False, False, False)
                        self.marker = ax.plot(event.artist.get_xdata(), event.artist.get_ydata(), ms=12, linewidth=3, alpha=0.4,
                                             color='yellow')[0]
                        ax.autoscale_view(True, True, False)
                    else:
                        self.cur_line = event.artist
                        ax.autoscale_view(False, False, False)
                        self.marker = ax.plot(event.artist.get_xdata(), event.artist.get_ydata(), ms=12, linewidth=3, alpha=0.4,
                                             color='yellow')[0]
                        ax.autoscale_view(True, True, False)
                    event.artist.get_figure().canvas.draw()
        elif event.mouseevent.button == 3:
            if isinstance(event.artist, lines.Line2D):
                if self.marker is None:
                    self.cur_line = event.artist
                    ax.autoscale_view(False, False, False)
                    self.marker = ax.plot(event.artist.get_xdata(), event.artist.get_ydata(), ms=12, linewidth=3, alpha=0.4,
                                         color='yellow')[0]
                    ax.autoscale_view(True, True, False)
                elif self.cur_line is not event.artist and self.marker in ax.lines:
                    self.marker.remove()
                    self.marker = None
                    self.cur_line = event.artist
                    ax.autoscale_view(False, False, False)
                    self.marker = ax.plot(event.artist.get_xdata(), event.artist.get_ydata(), ms=12, linewidth=3, alpha=0.4,
                                         color='yellow')[0]
                    ax.autoscale_view(True, True, False)
                event.artist.get_figure().canvas.draw()
                if self.cur_line is event.artist:
                    self.popup_menu = QtGui.QMenu()
                    def _delete_line():
                        figure=event.artist.get_figure()
                        self.marker.remove()
                        self.marker = None
                        legend_labels=ax.legend_.get_texts()
                        for i in range(len(legend_labels)):
                            if legend_labels[i].get_text() == self.cur_line.get_label():
                                del legend_labels[i]
                                break
                        self.cur_line.remove()
                        self.cur_line = None
                        if len(legend_labels):
                            ax.legend(framealpha=config["legend_framealpha"], loc=config["legend_loc"])
                        else:
                            ax.legend_.remove()
                        ax.relim()
                        ax.autoscale_view(True, True, True)
                        figure.canvas.draw()
                    self.popup_menu.addAction("Delete line", _delete_line)
                    self.popup_menu.popup(QtGui.QCursor.pos())
        ax.hold(self.toolbar._actions["hold_figure"].isChecked())

    def onpick_tdms(self, event):
        ax=event.artist.axes
        ax.hold(True)
        if event.mouseevent.button == 1:
            if event.mouseevent.dblclick:
                if isinstance(event.artist, lines.Line2D):
                    if self.marker is not None:
                        self.configurator=Graph_Configurator(self, ax)
            else:
                if isinstance(event.artist, lines.Line2D):
                    def _select_action():
                        def _onselect(xmin, xmax):
                            if xmin != xmax:
                                indmin, indmax = np.searchsorted(self.baseline[0], (xmin, xmax))
                                indmax = min(len(self.baseline[0])-1, indmax)
                                
                                self.baseline_regions_x+=self.baseline[0][indmin:indmax]
                                self.baseline_regions_y+=self.baseline[1][indmin:indmax]

                                self.parent.label_regions.setText((self.parent.label_regions.text()+" ["+str("%10.4f"%xmin).strip()+" - "+str("%10.4f"%xmax).strip()+"]").replace("Please select regions for baselinefit", ""))
                        self.cur_line = event.artist
                        ax.autoscale_view(False, False, False)
                        self.marker = ax.plot(event.artist.get_xdata(), event.artist.get_ydata(), ms=12, linewidth=3, alpha=0.4,
                                             color='yellow')[0]
                        ax.autoscale_view(True, True, False)
                        if self.parent.toolButton_set_baseline.isChecked():
                            x_data, y_data=self.cur_line.get_data()
                            self.baseline=[list(x_data), list(y_data)]
                            self.parent.label_baseline.setText(self.cur_line.get_label().replace("$",""))
                            if not config["automatic_baseline"]:
                                self.parent.label_regions.setText("Please select regions for baselinefit")
                                self.span = SpanSelector(ax, _onselect, 'horizontal', useblit=True,
                                                    rectprops=dict(alpha=0.5, facecolor='red'))
                            else:
                                end_index=len(self.baseline[0])-1
                                self.baseline_regions_x+=self.baseline[0][5:15]+self.baseline[0][end_index-15:end_index]
                                self.baseline_regions_y+=self.baseline[1][5:15]+self.baseline[1][end_index-15:end_index]
                                self.parent.label_regions.setText("["+str("%10.4f"%self.baseline[0][5]).strip()+" - "+str("%10.4f"%self.baseline[0][15]).strip()+"] ["+str("%10.4f"%self.baseline[0][5])-strip()+" - "+str("%10.4f"%self.baseline[0][15])-strip()+"]")
                                self.parent.toolButton_set_baseline.setChecked(0)
                    if self.marker is not None and self.cur_line is event.artist:
                        self.marker.remove()
                        self.marker = None
                    elif self.marker is not None and self.cur_line is not event.artist and self.marker in ax.lines:
                        self.marker.remove()
                        self.marker = None
                        _select_action()
                                

                    else:
                        _select_action()
                    event.artist.get_figure().canvas.draw()
        elif event.mouseevent.button == 3:
            if isinstance(event.artist, lines.Line2D):
                if self.marker is None:
                    self.cur_line = event.artist
                    ax.autoscale_view(False, False, False)
                    self.marker = ax.plot(event.artist.get_xdata(), event.artist.get_ydata(), ms=12, linewidth=3, alpha=0.4,
                                         color='yellow')[0]
                    ax.autoscale_view(True, True, False)
                elif self.cur_line is not event.artist and self.marker in ax.lines:
                    self.marker.remove()
                    self.marker = None
                    self.cur_line = event.artist
                    ax.autoscale_view(False, False, False)
                    self.marker = ax.plot(event.artist.get_xdata(), event.artist.get_ydata(), ms=12, linewidth=3, alpha=0.4,
                                         color='yellow')[0]
                    ax.autoscale_view(True, True, False)
                event.artist.get_figure().canvas.draw()
                if self.cur_line is event.artist:
                    self.popup_menu = QtGui.QMenu()
                    def _delete_line():
                        figure=event.artist.get_figure()
                        self.marker.remove()
                        self.marker = None
                        legend_labels=ax.legend_.get_texts()
                        for i in range(len(legend_labels)):
                            if legend_labels[i].get_text() == self.cur_line.get_label():
                                del legend_labels[i]
                                break
                        self.cur_line.remove()
                        self.cur_line = None
                        if len(legend_labels):
                            ax.legend(framealpha=config["legend_framealpha"], loc=config["legend_loc"])
                        else:
                            ax.legend_.remove()
                        ax.relim()
                        ax.autoscale_view(True, True, True)
                        figure.canvas.draw()
                    def _savgol_filter():
                        figure=event.artist.get_figure()
                        self.marker.remove()
                        self.marker = None
                        x_data, y_data = self.cur_line.get_data()
                        ax.autoscale_view(False, False, False)
                        ax.plot(x_data, savgol_filter(y_data, config["savgol_window_length"], config["savgol_polyorder"]), label=self.cur_line.get_label()+"$\;Savgol\;filter$", picker=5)
                        ax.legend(framealpha=config["legend_framealpha"], loc=config["legend_loc"])
                        ax.autoscale_view(True, True, False)
                        figure.canvas.draw()
                        
                    self.popup_menu.addAction("Savgol Golay filter", _savgol_filter)
                    self.popup_menu.addAction("Delete line", _delete_line)
                    self.popup_menu.popup(QtGui.QCursor.pos())
        ax.hold(self.toolbar._actions["hold_figure"].isChecked())
        
class NavigationToolbar(NavigationToolbar2QT):
    def __init__(self,canvas,parent,figure,remove=[]):
        self.canvas = canvas
        self.figure = figure
        self.span = {}
        self.toolitems = (
            ('Home', 'Reset original view', 'home', 'home'),
            ('Back', 'Back to  previous view', 'back', 'back'),
            ('Forward', 'Forward to next view', 'forward', 'forward'),
            (None, None, None, None),
            ('Pan', 'Pan axes with left mouse, zoom with right', 'move', 'pan'),
            ('Zoom', 'Zoom to rectangle', 'zoom_to_rect', 'zoom'),
            (None, None, None, None),
            ('Subplots', 'Configure subplots', 'subplots', 'configure_subplots'),
            ('Save', 'Save the figure', 'filesave', 'save_figure'),
            ('Csv Save', 'Save data to CSV File', '', 'csv_save'),
            ('Legend', 'Toggle legend', '', 'legend_switch'),
            ('Hold', 'Hold graph', '', 'hold_figure'),
            ('Cut', 'Cut figure', '', 'cut_figure'),
            
            )
        NavigationToolbar2QT.__init__(self,canvas,parent)
        actions = self.findChildren(QtGui.QAction)
        self._actions["csv_save"].setIcon(QtGui.QIcon("images/CSV.png"))
        self._actions["legend_switch"].setIcon(QtGui.QIcon("images/legend.png"))
        self._actions["hold_figure"].setIcon(QtGui.QIcon("images/toggle.png"))
        self._actions["cut_figure"].setIcon(QtGui.QIcon("images/Cut.png"))
        self._actions["cut_figure"].setCheckable(1)
        self._actions["legend_switch"].setCheckable(1)
        self._actions["hold_figure"].setCheckable(1)
        for a in actions:
            if a.text() == 'Customize':
                self.removeAction(a)
                break

        for i in remove:
            self._actions[i].hide()
        self._actions["hold_figure"].setChecked(config["hold_graph_on_start"])

    def cut_figure(self):
        if self._actions["cut_figure"].isChecked():
            def onselect(xmin, xmax):
                ax=self.figure.gca()
                if xmin != xmax:
                    if self.parent.browser.marker is not None:
                        self.parent.browser.marker.remove()
                        self.parent.browser.marker=None
                    for i in ax.lines.copy():
                        x_data, y_data=i.get_data()
                        indmin, indmax = np.searchsorted(x_data, (xmin, xmax))
                        indmax = min(len(x_data)-1, indmax)
                        label=i.get_label()
                        i.remove()
                        ax.plot(x_data[indmin:indmax], y_data[indmin:indmax], label=label, picker=5)
                    ax.legend(framealpha=config["legend_framealpha"], loc=config["legend_loc"])
                    ax.grid(True)
                    ax.relim()
                    ax.autoscale_view(True, True, True)
                    self.span={}
                    self._actions["cut_figure"].setChecked(0)
                    self.canvas.draw()
            ax=self.figure.get_axes()
            for i in ax:
                self.span[i]=SpanSelector(i, onselect, 'horizontal', useblit=True,
                                                    rectprops=dict(alpha=0.5, facecolor='red'))
            self.canvas.draw()
        else:
            self.span={}
                        

                    
    def hold_figure(self):
        ax=self.figure.get_axes()
        for i in ax:
            if self._actions["hold_figure"].isChecked():
                i.hold(True)
            else:
                i.hold(False)

                    
    def save_figure(self):
        if self.figure.get_axes():
            filename=self.figure.get_axes()[0].get_title().replace("$","").replace("\;", " ")
            path = QtGui.QFileDialog.getSaveFileName(self, "Save Figure", filename, "Portable Network Graphics (*.png);;Encapsulated Postscript (*.eps)\
                                                                                ;;Portable Document Format (*.pdf);;PGF code for LaTeX (*.pgf);;\
                                                                                Postscript (*.ps);;Raw RGBA bitmap (*.raw,*.rgba);;\
                                                                                 Scalable Vector Graphics (*.svg,*svgz)")
            if path:
                fig_size=self.figure.get_size_inches()
                self.figure.set_size_inches(9,9)
##                self.figure.savefig(path, dpi=200, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format=None, transparent=False, bbox_inches="tight", pad_inches=0.1, frameon=None)
                self.figure.savefig(path, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format=None, transparent=False, bbox_inches="tight", pad_inches=0.1, frameon=None)
                self.figure.set_size_inches(fig_size[0],fig_size[1])
                self.canvas.draw()
        

                
    def csv_save(self):
        if self.figure.get_axes():
            filename=self.figure.get_axes()[0].get_title().replace("$","").replace("\;", " ")
            path = QtGui.QFileDialog.getSaveFileName(self, "Speicherort auswählen", filename, "CSV-Datei (*.csv)")
            if path:
                ax=self.figure.get_axes()
                with open(path, "w") as f:
                    writer=csv.writer(f, delimiter=";", lineterminator="\n")
                    for j in ax:
                        for i in j.lines:
                            writer.writerow([i.get_label()])
                            writer.writerow(["x","y"])
                            x_werte=i.get_xdata()
                            y_werte=i.get_ydata()
                            for j in range(len(x_werte)):
                                writer.writerow([x_werte[j],y_werte[j]])
                            writer.writerow([""])
                f.close()
        
    def legend_switch(self):
        ax=self.figure.get_axes()
        for i in ax:
            if i.legend_:
                if self._actions["legend_switch"].isChecked():
                    i.legend_.set_visible(True)
                else:
                    i.legend_.set_visible(False)
        self.canvas.draw()

class Graph_Configurator(QtGui.QWidget):
    def __init__(self, parent=None, axis=None):
        super(Graph_Configurator, self).__init__()
        uic.loadUi("ui_files/graph_configurator.ui", self)
        self.parent = parent
        self.ax = axis


        self.cur_line_label=self.parent.cur_line.get_label()
        self.x_label=self.ax.get_xlabel()
        self.y_label=self.ax.get_ylabel()
        self.title=self.ax.get_title()

        
        self.lineEdit_line_label.setText(self.parent.cur_line.get_label())
        self.lineEdit_x_label.setText(self.ax.get_xlabel())
        self.lineEdit_y_label.setText(self.ax.get_ylabel())
        self.lineEdit_title.setText(self.ax.get_title())
        
        
        self.pushButton_save.clicked.connect(self.save)
        self.pushButton_cancel.clicked.connect(self.close)
        self.setWindowModality(QtCore.Qt.ApplicationModal)
        self.show()

    def save(self):
        try:
            self.ax.set_title(self.lineEdit_title.text())
            self.ax.set_xlabel(self.lineEdit_x_label.text())
            self.ax.set_ylabel(self.lineEdit_y_label.text())
            legend_labels=self.parent.parent.figure.gca().legend_.get_texts()
            for i in range(len(legend_labels)):
                if legend_labels[i].get_text() == self.parent.cur_line.get_label():
                    self.parent.cur_line.set_label(self.lineEdit_line_label.text())
                    legend_labels[i].set_text(self.lineEdit_line_label.text())
                    if self.parent.marker:
                        self.parent.marker.remove()
                        self.parent.marker = None
                    self.parent.parent.canvas.draw()
            self.close()
        except ValueError as e:
            QtGui.QMessageBox.critical(self, "Error", str(e), QtGui.QMessageBox.Ok)
            self.ax.set_title(self.title)
            self.ax.set_xlabel(self.x_label)
            self.ax.set_ylabel(self.y_label)
            self.parent.cur_line.set_label(self.cur_line_label)
            self.ax.legend(framealpha=config["legend_framealpha"], loc=config["legend_loc"])
        
            self.parent.parent.canvas.draw()


            
