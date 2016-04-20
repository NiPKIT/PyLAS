import sys
from PyQt4 import QtGui, uic, QtCore
import hapi
import matplotlib as mpl
mpl.use('QT4Agg')
import matplotlib.pyplot as plt
from functools import partial
from scipy.interpolate import interp1d
from nptdms import TdmsFile
import os
import csv
import json
import gc
import datetime
import numpy as np
from NavigationToolBar import *
from settings import *



class MainWindow(QtGui.QMainWindow):
##    MainWindow with all Tabs and connections to actionmenu
    def __init__(self):
        super(MainWindow, self).__init__()
        uic.loadUi("ui_files/root.ui", self)
        self.show()

        self.unit_blocker=False
        if not config["unit_wave"]:
            self.tableWidget_lineselection.horizontalHeaderItem(3).setText("Wave length")
            self.label_wave.setText("Wave length [nm]")
            self.doubleSpinBox_vmin.setMinimum(1)
        if not config["unit_pressure"]:
            self.label_pressure.setText("Pressure [bar]")
        if not config["unit_temperature"]:
            self.doubleSpinBox_temperature.setRange(70-273.15, 3000-273.15)
            self.label_temperature.setText("Temperature [°C]")

        self.doubleSpinBox_vmin.setValue(config["vmin"][config["unit_wave"]])
        self.doubleSpinBox_vmax.setValue(config["vmax"][config["unit_wave"]])
        self.doubleSpinBox_pressure.setValue(config["pressure"][config["unit_pressure"]])
        self.doubleSpinBox_temperature.setValue(config["temperature"][config["unit_temperature"]])
        if config["broadening"] == "gamma_self":
            broadening_text="self"
        elif config["broadening"] == "gamma_air":
            broadening_text="air"
        self.broadening_label.setText("Broadening: "+broadening_text+" ["+config["profile"]+"]")
        self.path_length_label.setText("Path length: "+str(config["path_length"])+" cm")


        self.doubleSpinBox_vmin.valueChanged.connect(self.set_wave_min_value)
        self.doubleSpinBox_vmax.valueChanged.connect(self.set_wave_max_value)
        self.doubleSpinBox_pressure.valueChanged.connect(self.set_pressure_value)
        self.doubleSpinBox_temperature.valueChanged.connect(self.set_temperature_value)

        self.hitran_tab=Spectrum_Tab(self,self.grid_layout_hitran)

        self.lineselection_tab=Lineselection_Tab(self)

        self.tdlas_tab=TDMS_File_Browser(self, self.tdlas_layout)

        self.actionEinstellungen.triggered.connect(self.preferences)
        self.actionAbout.triggered.connect(self.about)

        mpl.rc('font', family='serif')
##        mpl.rc('font', size=15.0)
##        mpl.rc('legend', fontsize=15.0)
##        mpl.rc('font', weight='normal')
##        mpl.rc('text', usetex=True)


    def preferences(self):
        self.prefWindow = prefWindow(self)
        self.prefWindow.setWindowModality(QtCore.Qt.ApplicationModal)
        self.prefWindow.show()


    def about(self):
        QtGui.QMessageBox.about(self, "About", "PyLAS uses the HITRAN-2012/HITEMP-2010 databases and a modified version of the HITRAN Application Programming Interface to process the line-by-line data given by HITRAN/HITEMP.<br><a href=\"https://www.cfa.harvard.edu/hitran/\">HITRAN Page</a><br><br>All icons from <a href=\"https://icons8.com/\">Icons8.com</a>")


    def closeEvent(self, event):
        json.dump(config, open("txt/config.txt", "w"))
        event.accept()

##    write changed parameters into config dictionary
    def set_pressure_value(self, value):
        if config["unit_pressure"]:
            config["pressure"][1] = value
            config["pressure"][0] = value*1.01325
        else:
            config["pressure"][1] = value/1.01325
            config["pressure"][0] = value

    def set_temperature_value(self, value):
        if config["unit_temperature"]:
            config["temperature"][1] = value
            config["temperature"][0] = value - 273.15
        else:
            config["temperature"][1] = value + 273.15
            config["temperature"][0] = value

    def set_wave_min_value(self, value):
        if config["unit_wave"]:
            config["vmin"][1]=value
            config["vmax"][0]=1e7/value
        else:
            config["vmin"][0]=value
            config["vmax"][1]=1e7/value
##        if not self.unit_blocker:
##            self.doubleSpinBox_vmin.setMaximum(self.doubleSpinBox_vmax.value()-0.01)

    def set_wave_max_value(self, value):
        if config["unit_wave"]:
            config["vmax"][1]=value
            config["vmin"][0]=1e7/value
        else:
            config["vmax"][0]=value
            config["vmin"][1]=1e7/value
##        if not self.unit_blocker:
##            self.doubleSpinBox_vmax.setMinimum(self.doubleSpinBox_vmin.value()+0.01)


class TDMS_Data_Tab(QtGui.QWidget):
    def __init__(self, parent=None):
        super(TDMS_Data_Tab, self).__init__()
        uic.loadUi("ui_files/tdms_data_tab.ui", self)
        self.parent = parent

##TDMS GUI for plotted channels of tdms file
class TDMS_Graph(QtGui.QWidget):
    def __init__(self, parent=None, full_x=None):
##        main GUI buttons and figure
        super(TDMS_Graph, self).__init__()
        uic.loadUi("ui_files/tdms_graph.ui", self)
        self.parent = parent
        self.is_baseline_plotted=False
        self.full_x=full_x

        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self, self.figure)
        self.graph_layout.addWidget(self.toolbar,0,0)
        self.graph_layout.addWidget(self.canvas,1,0)

        self.browser = PointBrowser(self, self.toolbar)
        self.canvas.mpl_connect('pick_event', self.browser.onpick_tdms)
        self.canvas.mpl_connect('resize_event', self.browser.resize)

        self.toolButton_set_baseline.clicked.connect(self.set_baseline)
        self.pushButton_optical_density.clicked.connect(self.plot_optical_density)
        self.pushButton_plot_spectrum.clicked.connect(partial(self.load_data, "absorptioncoefficient", self.plot_spectrum))

        self.baseline_graphs={}
        self.fitted_baseline=None
        self.optical_density=None
        self.baseline_graph_counter=1


        self.file_loader = QtFileLoader(self)

        self.statusBar = QtGui.QStatusBar(self)
        self.statusBar.setMaximumHeight(40)
        self.statusBar.setSizeGripEnabled(False)
        self.progress_label = QtGui.QLabel("")
        self.progress_label.setMaximumHeight(40)
        self.progressBar = QtGui.QProgressBar(self)
        self.progressBar.setMaximumHeight(40)
        self.terminate_button = QtGui.QToolButton()
        self.terminate_button.setMaximumHeight(40)
        self.terminate_button.setIcon(QtGui.QIcon("images/cancel.png"))
        self.terminate_button.setAutoRaise(1)
        self.terminate_button.hide()
        self.progressBar.hide()
        self.progress_label.hide()
        self.progressBar.setMaximumSize(QtCore.QSize(200,20))
        self.statusBar.addPermanentWidget(self.progress_label)
        self.statusBar.addPermanentWidget(self.progressBar)
        self.statusBar.addPermanentWidget(self.terminate_button)
        self.statusbar_layout.addWidget(self.statusBar)

        self.terminate_button.clicked.connect(self.terminate)

        self.toolButton_concentration.hide()
        self.label_concentration.hide()
        self.toolButton_concentration.clicked.connect(self.calculate_concentration)

        self.final_measurement=[]
        self.theor_spectrum=[]

        self.span_conc=None
        
##    file loader specific methods
    def error(self, message):
        self.box("Error", message)

    def log_print(self, message, IDS_to_uncheck=[]):
        if config["uncheck_missing"]:
            for i in IDS_to_uncheck:
                self.parent.parent.parent.tdmsmoltree.selected_ISO_ID[i].setCheckState(0, 0)
        with open("txt/logfile.txt", "a") as f:
            f.write(str(datetime.datetime.now())+" "+message+"\n")
        f.close()

    def line_counter(self, RowID, nline):
        self.progress_label.setText("processing data")
        self.progressBar.setRange(0,nline)
        self.progressBar.setValue(RowID)

    def terminate(self):
        self.progress_label.setText("cancel process")
        self.file_loader.stop()

    def box(self, title, text):
        QtGui.QMessageBox.critical(self, title, text, QtGui.QMessageBox.Ok)

##    begin to load HITRAN/HITEMP data in initialized QThread
    def load_data(self, spectrum, function):
        if not self.parent.parent.parent.tdmsmoltree.selected_ISO_ID:
            self.box("Error", "Please select a molecule.")
        elif self.file_loader.isRunning():
            self.box("Warning", "Wait until the calculation has finished.")
        elif config["vmin"][1] >= config["vmax"][1]:
            self.box("Warning", "Please check wave range.")
        else:
            ISO_IDS = list(self.parent.parent.parent.tdmsmoltree.selected_ISO_ID.keys())
            self.file_loader.update("tdms_tab", ISO_IDS, self.parent.parent.parent.tdmsmoltree.concentrations, spectrum, config["profile"])
            if self.file_loader.new_data:
                self.progress_label.setText("loading data")
                self.progressBar.setRange(0,0)
                self.progress_label.show()
                self.progressBar.show()
                self.terminate_button.show()
                self.file_loader.finished.connect(function)
                self.file_loader.start()


##    plot theoretical spectrum in front of measured / cross correlation
    def plot_spectrum(self):
        self.file_loader.finished.disconnect()
        self.progress_label.hide()
        self.progressBar.hide()
        self.terminate_button.hide()
        if self.file_loader.processed:
            ax=self.figure.gca()
            density_line=ax.lines[0].get_ydata()
            x_spec=self.file_loader.x
            y_spec=np.array(self.file_loader.y)*config["path_length"]
            if not config["unit_wave"]:
                x_spec = [1e7/i for i in reversed(x_spec)]
                y_spec = [i for i in reversed(y_spec)]

##            Messbereich in cm^-1
##            Measurerange in cm^-1
            cutted_line = len(density_line)/len(self.parent.full_x)*config["measure_range"]
            x_measure_range=np.arange(0,cutted_line,cutted_line/len(density_line))

##            Werteabstände Messbereich mit Werteabstände theoretischen Spektren gleichsetzen
##            equalize distances of theoretical and measure values
            interpolation=interp1d(x_measure_range, density_line)
            x_measure_range=np.arange(0, x_measure_range[len(x_measure_range)-2], config["omega_step"])
            final_measurement=interpolation(x_measure_range)
            
##            Timeshift / Wellenshift durch Kreuzkorrelation ermitteln
##            Cross correlation to get the wave shift of the measurement
            correlation=np.correlate(y_spec, final_measurement, "full")
            shift=len(x_measure_range)-1-np.argmax(correlation)
            
##            korrelierter Wellenbereich für Messung
##            set wave range for measured data
            density_x=x_spec[0]+x_measure_range-(shift*config["omega_step"])

            if not config["unit_wave"]:
                x_spec              = [1e7/i for i in reversed(x_spec)]
                density_x           = [1e7/i for i in reversed(density_x)]
                y_spec              = [i for i in reversed(y_spec)]
                final_measurement   = [i for i in reversed(final_measurement)]

            self.theor_spectrum=[x_spec, y_spec]
            self.final_measurement=[density_x, final_measurement]
            
            ax.hold(False)
            ax.plot(x_spec, np.array(y_spec)*(max(final_measurement)/max(y_spec)), label=", ".join([self.parent.parent.parent.tdmsmoltree.selected_tex_molecules[i] for i in self.file_loader.molecule_names]), picker=5)
            ax.hold(True)
            ax.plot(density_x, final_measurement, label="$Measurement$", picker=5)
            ax.legend(framealpha=config["legend_framealpha"], loc=config["legend_loc"])
            ax.grid(True)
            self.figure.tight_layout()
            ax.autoscale_view(True, True, False)
            ax.get_xaxis().get_major_formatter().set_useOffset(False)
            ax.hold(self.toolbar._actions["hold_figure"].isChecked())
            self.canvas.draw()
            self.toolButton_concentration.show()
            self.label_concentration.show()
##            self.pushButton_plot_spectrum.hide()
            

##    use selected range, calculate areas of measured line and theoretical spectrum and then calculate the concentration
    def calculate_concentration(self, clicked):
        if clicked:
            ax=self.figure.gca()
            def _onselect(xmin, xmax):
                if xmin != xmax:
                    indmin, indmax = np.searchsorted(self.final_measurement[0], (xmin, xmax))
                    indmax = min(len(self.final_measurement[0])-1, indmax)
                    final_y=self.final_measurement[1][indmin:indmax]
                    
                    indmin, indmax = np.searchsorted(self.theor_spectrum[0], (xmin, xmax))
                    indmax = min(len(self.theor_spectrum[0])-1, indmax)
                    meas_y=self.theor_spectrum[1][indmin:indmax]
                    
                    conc=np.trapz(final_y, dx=0.01)/np.trapz(meas_y, dx=0.01)*100
                    
                    self.span_conc=None
                    self.toolButton_concentration.setChecked(0)
                    self.label_concentration.setText("Concentration: "+"%10.2f"%conc+" %")
                
            self.span_conc = SpanSelector(ax, _onselect, 'horizontal', useblit=True,
                                                        rectprops=dict(alpha=0.5, facecolor='red'))
            self.canvas.draw()
            self.label_concentration.setText("Please select region for concentration calculation")
        else:
            self.span_conc=None


##    use fitted baseline and measurement to calculate the optical density
    def plot_optical_density(self):
        if not self.toolButton_set_baseline.isChecked() and self.is_baseline_plotted:
            base_graph=TDMS_Graph(self)
            tab_title="Graph: Baseline "+str(self.baseline_graph_counter)
            
            base_graph.toolButton_set_baseline.hide()
            base_graph.label_baseline.hide()
            base_graph.label_regions.hide()
            base_graph.pushButton_optical_density.hide()
            ax=base_graph.figure.gca()
            self.optical_density=np.log(self.fitted_baseline/np.array(self.browser.baseline[1]))
            ax.plot(self.browser.baseline[0], self.optical_density, label="$Optical\;density$", picker=5)
            ax.set_xlabel(self.figure.gca().get_xlabel())
            ax.set_ylabel("$Optical\;density\;[-ln(I/I_{0})]$")
            base_graph.figure.tight_layout()
            ax.grid(True)
            ax.legend(framealpha=config["legend_framealpha"], loc=config["legend_loc"])
            base_graph.toolbar._actions["legend_switch"].setChecked(True)
            self.baseline_graph_counter+=1
            self.parent.file_tab.addTab(base_graph, tab_title)
            self.baseline_graphs[tab_title]=base_graph

##    define regions to use for baselinefitting
    def set_baseline(self, clicked):
        if clicked:
            self.is_baseline_plotted=False
            self.label_baseline.setText("")
            self.label_regions.setText("Please select line for baselinefit")
            self.browser.baseline=[]
            self.browser.baseline_regions_x=[]
            self.browser.baseline_regions_y=[]
        else:
            self.label_baseline.setText("")
            self.label_regions.setText("")
            self.browser.span=None
            if self.browser.baseline_regions_x:
                self.is_baseline_plotted=True
                self.fitted_baseline=np.polyval(np.polyfit(self.browser.baseline_regions_x,self.browser.baseline_regions_y, config["baseline_polyorder"]), self.browser.baseline[0])
                ax=self.figure.gca()
                ax.plot(self.browser.baseline[0], self.fitted_baseline, label="$Baseline$", picker=5)
                ax.legend(framealpha=config["legend_framealpha"], loc=config["legend_loc"])
                self.canvas.draw()
        
        
##TDMS GUI for loaded tdms file with main header info 
class TDMS_File_Tab(QtGui.QWidget):
    def __init__(self, parent=None, path=None):
        super(TDMS_File_Tab, self).__init__()
        uic.loadUi("ui_files/tdms_file.ui", self)
        self.parent = parent
        self.path=path

        self.needed_properties=["wf_increment",
                                "wf_samples",
                                "NI_UnitDescription",
                                "NI_ChannelName",
                                "NI_ExpStartTimeStamp",
                                "wf_start_offset",
                                "wf_xname",
                                "wf_xunit_string",
                                "wf_time_pref",
                                "frequency"]

        self.show_properties={"NI_ExpStartTimeStamp":"Start time",
                              "wf_increment":"Interval",
                              "data_sets":"Data sets",
                              "flame_noise":"Flame noise",
                              "frequency":"Frequency",
                              "wf_samples":"Samples",
                              "NI_ChannelName":"Channelname",
                              "NI_UnitDescription":"Unit",
                              "wf_xname":"X name",
                              "wf_xunit_string":"X unit",
                              "wf_start_offset":"Start offset",
                              "wf_time_pref":"Time format",}


        self.properties={}
        self.data={}

        self.groups={}
        self.channels={}
        self.averaged=0

        self.x_axis={}
        self.y_axes={}

        self.noise_flame={}
        self.data_sets={}

        self.header={}
        self.header_labels={}

        self.tabs={}
        self.hold_graph={}
        self.graphs={}
        self.graph_count={}
        self.holded_graph=None
        self.done = True

        self.toolButton_csv_save.clicked.connect(self.save_tab)
        self.toolButton_x_axis.clicked.connect(self.set_x_axis)
        self.toolButton_y_axes.clicked.connect(self.set_y_axes)
        self.file_tab.currentChanged.connect(self.change_axe_label)
        self.file_tab.tabCloseRequested.connect(self.close_tab)
        self.pushButton_plot.clicked.connect(self.plot_graph)
        self.toolButton_hold_graph.clicked.connect(self.select_hold_graph)

        self.open_file(path)

##    functions for column selection and plotting selected columns
    def set_x_axis(self, clicked):
        if clicked:
            if self.toolButton_y_axes.isChecked():
                self.toolButton_y_axes.setChecked(0)
            self.axes_hint_label.setText("<font color=\"red\">Please select X - axis</font>")
        else:
            self.axes_hint_label.setText("")
        
    def set_y_axes(self, clicked):
        current_Tab=self.file_tab.tabText(self.file_tab.currentIndex())
        if clicked:
            self.y_axes[current_Tab]=[]
            self.y_axes_label.setText("")
            if self.toolButton_x_axis.isChecked():
                self.toolButton_x_axis.setChecked(0)
            self.axes_hint_label.setText("<font color=\"red\">Please select Y - axes</font>")
        else:
            if not self.toolButton_x_axis.isChecked():
                self.axes_hint_label.setText("")
        
    def change_axe_label(self, index):
        current_Tab=self.file_tab.tabText(index)
        if current_Tab in self.groups and current_Tab in self.x_axis:
            self.x_axis_label.setText(self.x_axis[current_Tab])
            self.y_axes_label.setText(", ".join(self.y_axes[current_Tab]))
    
    def select_axes(self, item):
        current_Tab=self.file_tab.tabText(self.file_tab.currentIndex())
        clicked_Header=self.header_labels[current_Tab][item.column()]
        if clicked_Header in self.header_labels[current_Tab] and current_Tab in self.x_axis:
            if self.toolButton_x_axis.isChecked():
                self.x_axis[current_Tab]=clicked_Header
                self.x_axis_label.setText(self.x_axis[current_Tab])
                self.toolButton_x_axis.setChecked(0)
                self.axes_hint_label.setText("")
            elif self.toolButton_y_axes.isChecked():
                if clicked_Header not in self.y_axes[current_Tab]:
                    self.y_axes[current_Tab].append(clicked_Header)
                    self.y_axes_label.setText(", ".join(self.y_axes[current_Tab]))

    def select_hold_graph(self, clicked):
        current_Tab=self.file_tab.tabText(self.file_tab.currentIndex())
        current_Group=None
        for i in self.groups:
            if i in current_Tab:
                current_Group=i
        if clicked and current_Group:
            items=[i for i in self.graphs[current_Group]]
            if not items:
                self.toolButton_hold_graph.setChecked(0)
            else:
                item, ok = QtGui.QInputDialog.getItem(self, "Select Graph", "Select Graph to plot in",
                                                              items, editable=False)
                if item and ok:
                    self.holded_graph=item
                    self.label_hold_graph.setText(item)
                else:
                    self.toolButton_hold_graph.setChecked(0)
        else:
            self.holded_graph=None
            self.label_hold_graph.setText("")
            

    def plot_graph(self):
        current_Tab=self.file_tab.tabText(self.file_tab.currentIndex())
        if current_Tab in self.groups:
            if self.x_axis[current_Tab] and self.y_axes[current_Tab]:
                self.toolButton_x_axis.setChecked(0)
                self.toolButton_y_axes.setChecked(0)
                self.axes_hint_label.setText("")

                if not self.toolButton_hold_graph.isChecked() or not self.graph_count[current_Tab] or not self.holded_graph:
                    new_graph=TDMS_Graph(self, full_x=self.data[current_Tab][self.x_axis[current_Tab]])
                    new_graph.pushButton_plot_spectrum.hide()
                    if self.graph_count[current_Tab]:
                        tab_title="Graph "+str(self.graph_count[current_Tab])+": "+current_Tab
                    else:
                        tab_title="Graph: "+current_Tab
                    self.graphs[current_Tab][tab_title]=new_graph
                    self.graph_count[current_Tab]+=1
                    self.file_tab.addTab(new_graph, tab_title)
                    ax=new_graph.figure.gca()
                    for i in self.y_axes[current_Tab]:
                        ax.plot(self.data[current_Tab][self.x_axis[current_Tab]], self.data[current_Tab][i], label="$"+i+"$", picker=5)
                    ax.set_xlabel("$"+self.x_axis[current_Tab]+"$")
                    ax.set_ylabel("$"+"/".join(self.y_axes[current_Tab])+"$")
                    new_graph.figure.tight_layout()
                    ax.grid(True)
                    ax.legend(framealpha=config["legend_framealpha"], loc=config["legend_loc"])
                    new_graph.toolbar._actions["legend_switch"].setChecked(True)
                    ax.autoscale_view(True, True, False)
                    new_graph.canvas.draw()
                else:
                    current_graph=self.graphs[current_Tab][self.holded_graph]
                    ax=current_graph.figure.gca()
                    for i in self.y_axes[current_Tab]:
                        ax.plot(self.data[current_Tab][self.x_axis[current_Tab]], self.data[current_Tab][i], label="$"+i+"$", picker=5)
                    ax.set_xlabel("$"+self.x_axis[current_Tab]+"$")
                    ax.set_ylabel("$"+"/".join(self.y_axes[current_Tab])+"$")
                    current_graph.figure.tight_layout()
                    ax.grid(True)
                    ax.legend(framealpha=config["legend_framealpha"], loc=config["legend_loc"])
                    current_graph.toolbar._actions["legend_switch"].setChecked(True)
                    ax.autoscale_view(True, True, False)
                    current_graph.canvas.draw()
                
            else:
                QtGui.QMessageBox.critical(self, "Error", "Please select axes.", QtGui.QMessageBox.Ok)


##    function for reading tdms file and filtering main data
    def open_file(self, path):

        if path.endswith(".tdms"):
            tdms_file = TdmsFile(path)
            self.groups=tdms_file.groups()
            tab_counter=0
            
            for group in self.groups:
                self.hold_graph[group]=False
                self.graphs[group]={}
                self.graph_count[group]=0
                
                currentWidget=TDMS_Data_Tab(self)
                currentWidget.file_table.itemClicked.connect(self.select_axes)
                self.file_tab.addTab(currentWidget, group)
                self.tabs[group]=currentWidget
                items=[]
                channels=[i for i in tdms_file.group_channels(group) if type(i.data) is np.ndarray]
                if not channels:
                    QtGui.QMessageBox.critical(self, "Error", "No data available in Group \""+group+"\"", QtGui.QMessageBox.Ok)
                    continue
                self.channels[group]=channels
                self.properties[group]={}
                for channel in channels:
                    channel_str=os.path.basename(channel.path).replace("'","")
                    items.append(channel_str)
                    self.properties[group][channel_str]={}
                    for i in self.needed_properties:
                        self.properties[group][channel_str][i]=""
                        if i in channel.properties:
                            self.properties[group][channel_str][i]=channel.properties[i]
                
                if not config["saved_tdms_column"] or config["saved_tdms_column"]:
                    item, ok = QtGui.QInputDialog.getItem(self, "Average", "Please choose channel for average calculation.<br>Channels for group \""+str(group)+"\"",
                                                          items, editable=False)
                else:
                    ok=True
                    item=config["saved_tdms_column"]
                if ok and item:
                    tdms_object=tdms_file.object(group, item)
                    data=tdms_object.data
                    time_track=tdms_object.time_track()
                    if config["start_index_manually"]:
                        value, bool_manual_start = QtGui.QInputDialog.getInt(self, "Start index", "Set start index for average calculation", 0, 1)
                        if bool_manual_start:
                            start_index_manual=value-1
                        else:
                            self.done=False
                            return
                    if not config["save_frequency"]:
                        value, ok = QtGui.QInputDialog.getItem(self, "Set Frequency", "Set the frequency", config["frequency"], 1)
                        if ok:
                            config["frequency"]=value
                    if ok:
                        self.data[group]={}
                        data_set_length=1/(config["frequency"]*self.properties[group][item]["wf_increment"])
                        if config["start_index_manually"]:
                            start_index=start_index_manual
                        else:
                            start_index=np.where(data[0:data_set_length] == data[0:data_set_length].min())[0][-1]+1
                        data=np.array_split(data, np.arange(start_index, len(data),data_set_length))[1:-1]
                        self.data_sets[group]=len(data)
                        data=np.mean(data, 0)
                        self.noise_flame[group]=data.min()
                        self.data[group][item]=data
                        for channel in channels:
                            channel_str=os.path.basename(channel.path).replace("'","")
                            if channel_str != item:
                                data=tdms_file.channel_data(group, channel_str)
                                self.data[group][channel_str]=np.mean(np.array_split(data, np.arange(start_index, len(data),data_set_length))[1:-1], 0)

                        current_file_table=currentWidget.file_table
                        current_file_table.setRowCount(data_set_length)

                        column_counter=0
                        temp_time_track=[]
                        double_header_labels=[]
                        self.header_labels[group]=[]
                        data_plot={}
                        x_axes_counter=0
                        self.y_axes[group]=[]
                        for channel in channels:
                            channel_str=os.path.basename(channel.path).replace("'","")
                            if set(temp_time_track) != set(tdms_file.object(group, channel_str).time_track()):
                                current_file_table.setColumnCount(column_counter+1)
                                if "wf_xname" in self.properties[group][channel_str]:
                                    current_header=self.properties[group][channel_str]["wf_xname"]+" ["+self.properties[group][channel_str]["wf_xunit_string"]+"]"
                                    if current_header in double_header_labels:
                                        double_header_labels.append(current_header)
                                        current_header=current_header+" "+double_header_labels.count(current_header)
                                    current_file_table.setHorizontalHeaderItem(column_counter, QtGui.QTableWidgetItem(current_header))
                                else:
                                    current_header=channel_str
                                    if current_header in double_header_labels:
                                        double_header_labels.append(current_header)
                                        current_header=current_header+" "+double_header_labels.count(current_header)
                                    current_file_table.setHorizontalHeaderItem(column_counter, QtGui.QTableWidgetItem(channel_str))
                                if x_axes_counter == 0:
                                    self.x_axis[group]=current_header
                                    self.x_axis_label.setText(current_header)
                                x_axes_counter+=1
                                time=tdms_file.object(group, channel_str).time_track()[0:int(data_set_length)]
                                data_plot[current_header]=time
                                for i in range(int(data_set_length)):
                                    current_file_table.setItem(i,column_counter, QtGui.QTableWidgetItem(str(time[i])))
                                self.header_labels[group].append(current_header)
                                column_counter+=1
                            current_file_table.setColumnCount(column_counter+1)
                            if "NI_UnitDescription" in self.properties[group][channel_str]:
                                current_header=channel_str+" ["+self.properties[group][channel_str]["NI_UnitDescription"]+"]"
                                if current_header in double_header_labels:
                                    double_header_labels.append(current_header)
                                    current_header=current_header+" "+double_header_labels.count(current_header)
                                current_file_table.setHorizontalHeaderItem(column_counter, QtGui.QTableWidgetItem(current_header))
                            else:
                                current_header=channel_str
                                if current_header in double_header_labels:
                                    double_header_labels.append(current_header)
                                    current_header=current_header+" "+double_header_labels.count(current_header)
                                current_file_table.setHorizontalHeaderItem(column_counter, QtGui.QTableWidgetItem(current_header))
                            data_plot[current_header]=self.data[group][channel_str]
                            for i in range(int(data_set_length)):
                                current_file_table.setItem(i,column_counter, QtGui.QTableWidgetItem(str(self.data[group][channel_str][i])))
                            temp_time_track=tdms_file.object(group, channel_str).time_track()
                            column_counter+=1
                            self.header_labels[group].append(current_header)
                        self.data[group]=data_plot


                        groupname=QtGui.QTreeWidgetItem(self.file_tree)
                        groupname.setText(0, group)
                        for channel, property_items in self.properties[group].items():
                            channelname = QtGui.QTreeWidgetItem(groupname)
                            channelname.setText(0, channel)
                            for i in self.needed_properties:
                                if i in property_items:
                                    if i == "NI_ExpStartTimeStamp":
                                        time_stamp=property_items[i].strftime("%d.%m.%Y %H:%M:%S")
                                        propertyname = QtGui.QTreeWidgetItem(channelname)
                                        propertyname.setText(0, self.show_properties[i])
                                        propertyname.setText(1, time_stamp)
                                    elif i == "wf_increment":
                                        propertyname = QtGui.QTreeWidgetItem(channelname)
                                        propertyname.setText(0, self.show_properties[i])
                                        propertyname.setText(1, str("%10.2E"%property_items[i]).strip())
                                    elif i == "frequency":
                                        propertyname = QtGui.QTreeWidgetItem(channelname)
                                        propertyname.setText(0, self.show_properties[i])
                                        propertyname.setText(1, str(config["frequency"])+" Hz")
                                    else:
                                        propertyname = QtGui.QTreeWidgetItem(channelname)
                                        propertyname.setText(0, self.show_properties[i])
                                        propertyname.setText(1, str(property_items[i]))
                        noise_of_group = QtGui.QTreeWidgetItem(groupname)
                        noise_of_group.setText(0, "Flame noise")
                        noise_of_group.setText(1, str("%10.2E"%self.noise_flame[group]).strip())
                        data_sets = QtGui.QTreeWidgetItem(groupname)
                        data_sets.setText(0, "Data Sets")
                        data_sets.setText(1, str(self.data_sets[group]))
                        self.header[group]=self.properties[group]
                        self.header[group]["flame_noise"]=self.data_sets[group]
                        self.header[group]["data_sets"]=self.noise_flame[group]
                    else:
                        self.done=False
                        return
                else:
                    self.done=False
        
        
    def close_tab(self, index):
        current_Tab=self.file_tab.tabText(index)
        current_Group=None
        if current_Tab not in self.groups:
            if "Baseline" in current_Tab and current_Tab in self.file_tab.widget(index).parent.baseline_graphs:
##                self.file_tab.widget(index).parent.baseline_graph_counter-=1
                del self.file_tab.widget(index).parent.baseline_graphs[current_Tab]
                self.file_tab.removeTab(index)
            else:
                for i in self.groups:
                    if i in current_Tab:
                        current_Group=i
                if current_Tab in self.graphs[current_Group]:
                    if current_Tab == self.holded_graph:
                        self.toolButton_hold_graph.setChecked(0)
                        self.holded_graph=None
                        self.label_hold_graph.setText("")
                    del self.graphs[current_Group][current_Tab]
##                    self.graph_count[current_Group]-=1
                    self.file_tab.removeTab(index)

##    save averaged tdms file to csv file
    def save_tab(self, index):
        filename=self.parent.tdms_tab_widget.tabText(self.parent.tdms_tab_widget.currentIndex())
        path = QtGui.QFileDialog.getSaveFileName(self, "Save File", os.path.splitext(filename)[0], "CSV-File (*.csv)")
        if path:
            try:
                with open(path, "w") as f:
                    writer=csv.writer(f, delimiter=";", lineterminator="\n")
                    for group in self.data:
                        writer.writerow([group])
                        writer.writerow([""])
                        writer.writerow(["Header"])
                        writer.writerow([self.header[group]])
                        writer.writerow([""])
                        writer.writerow(self.header_labels[group])
                        writer.writerows(np.transpose([self.data[group][i] for i in self.header_labels[group]]))
                        writer.writerow([""])
                f.close()
            except Exception as e:
                QtGui.QMessageBox.critical(self, "Error", str(e), QtGui.QMessageBox.Ok)
                pass
                

##TDMS Tab in MainGUI
class TDMS_File_Browser(QtGui.QWidget):
    def __init__(self, parent=None, tabLayout=None):
        super(TDMS_File_Browser, self).__init__()
##        uic.loadUi("ui_files/tdms_file.ui", self)
        self.parent = parent
        self.tabLayout=tabLayout

        self.tdmsmoltree=MoleculeTree(self, select="tdms")
        self.tabLayout.addWidget(self.tdmsmoltree.tree, 0, 0)

        self.tdms_tab_widget=QtGui.QTabWidget(self)
        self.tdms_tab_widget.setTabsClosable(True)
        self.tdms_tab_widget.setMovable(True)
        self.tdms_tab_widget.tabCloseRequested.connect(self.close_tab)

        self.tabLayout.addWidget(self.tdms_tab_widget, 0, 1)

        self.parent.actionoeffnen.triggered.connect(self.new_tab)

        self.tabname={}

    def close_tab(self, index):
        del self.tabname[self.tdms_tab_widget.tabText(index)]
        self.tdms_tab_widget.removeTab(index)


    def new_tab(self):
        path=QtGui.QFileDialog.getOpenFileName(self, "Open TDMS File", "", "Labview (*.tdms)")
        if path:
            if os.path.basename(path) in self.tabname:
                self.tdms_tab_widget.setCurrentIndex(self.tdms_tab_widget.IndexOf(self.tabname[os.path.basename(path)]))
            else:
                new_tab=TDMS_File_Tab(self, path)
                self.tabname[os.path.basename(path)]=new_tab
                if new_tab.done:
                    self.tdms_tab_widget.addTab(new_tab, os.path.basename(path))
                    self.parent.tabWidget.setCurrentIndex(2)
                else:
                    new_tab=None


class Lineselection_Tab(QtGui.QWidget):
    def __init__(self, parent=None):
        super(Lineselection_Tab, self).__init__()
##        uic.loadUi("ui_files/lineselection.ui", self)
        self.parent = parent

        self.statusBar = QtGui.QStatusBar(self)
        self.statusBar.setMaximumHeight(40)
        self.statusBar.setSizeGripEnabled(False)
        self.progress_label = QtGui.QLabel("")
        self.progressBar = QtGui.QProgressBar(self)
        self.terminate_button = QtGui.QToolButton()
        self.terminate_button.setIcon(QtGui.QIcon("images/cancel.png"))
        self.terminate_button.setAutoRaise(1)
        self.terminate_button.hide()
        self.progressBar.hide()
        self.progress_label.hide()
        self.progressBar.setMaximumSize(QtCore.QSize(200,20))
        self.statusBar.addPermanentWidget(self.progress_label)
        self.statusBar.addPermanentWidget(self.progressBar)
        self.statusBar.addPermanentWidget(self.terminate_button)
        self.parent.gridLayout_lineselection_status.addWidget(self.statusBar, 2, 0)
        
        
        self.moltree=MoleculeTree(self, select="lineselection")
        self.moltreeWidget=self.moltree.tree
        self.parent.molecule_layout.addWidget(QtGui.QLabel("Spectrum"),0,0)
        self.parent.molecule_layout.addWidget(self.moltreeWidget,1,0)

        self.moltree_selected=MoleculeTree(self, select="lineselection")
        self.moltree_selected_Widget=self.moltree_selected.tree
        self.parent.molecule_layout.addWidget(QtGui.QLabel("Moleculeselection"),2,0)
        self.parent.molecule_layout.addWidget(self.moltree_selected_Widget,3,0)

        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self, self.figure)
        self.spectrum_layout = QtGui.QVBoxLayout()
        self.spectrum_layout.addWidget(self.toolbar)
        self.spectrum_layout.addWidget(self.canvas)
        
        self.browser = PointBrowser(self, self.toolbar)
        self.canvas.mpl_connect('pick_event', self.browser.onpick)
        self.canvas.mpl_connect('resize_event', self.browser.resize)

        self.parent.gridLayout_lineselection.addLayout(self.spectrum_layout,0,0)


        self.file_loader_selected = QtFileLoader(self)
        
        self.terminate_button.clicked.connect(self.terminate)

        self.connect(self.file_loader_selected, self.file_loader_selected.line_counter_signal, self.line_counter)
        self.connect(self.file_loader_selected, self.file_loader_selected.error_signal, self.error)
        self.connect(self.file_loader_selected, self.file_loader_selected.log_signal, self.log_print_select)

        self.file_loader = QtFileLoader(self)

        self.connect(self.file_loader, self.file_loader.line_counter_signal, self.line_counter)
        self.connect(self.file_loader, self.file_loader.error_signal, self.error)
        self.connect(self.file_loader, self.file_loader.log_signal, self.log_print_main)

        self.parent.pushButtonstart_lineselection.clicked.connect(self.load_data_selected)

        self.lineselection_processor = QtLineprocessor(self)

        self.connect(self.lineselection_processor, self.lineselection_processor.line_counter_signal, self.line_counter_selection)
        self.connect(self.lineselection_processor, self.lineselection_processor.error_signal, self.error)
        self.connect(self.lineselection_processor, self.lineselection_processor.line_full_overlap_signal, self.selection_full_overlap)

    def line_counter_selection(self, line, full, text):
        self.progress_label.setText(text)
        self.progressBar.setRange(0,full)
        self.progressBar.setValue(line)

##    plot lineselection to view if full spectrum is overlayed
    def selection_full_overlap(self):

        self.parent.tableWidget_lineselection.clearContents()
        self.parent.tableWidget_lineselection.setRowCount(0)
            
        self.progress_label.hide()
        self.progressBar.hide()
        self.terminate_button.hide()
        QtGui.QMessageBox.critical(self, "Warning", "The complete spectrum is overlayed or no lines are available for given parameters", QtGui.QMessageBox.Ok)
        self.lineselection_processor.finished.disconnect()
        ax=self.figure.gca()
        if not config["unit_wave"]:
            self.lineselection_processor.x_main=[1e7/i for i in reversed(self.lineselection_processor.x_main)]
            self.lineselection_processor.Y_main=np.array([i for i in reversed(self.lineselection_processor.y_main)])
##            self.lineselection_processor.x_select=[1e7/i for i in reversed(self.lineselection_processor.x_select)]
##            self.lineselection_processor.Y_select=np.array([i for i in reversed(self.lineselection_processor.y_select)])
        ax.plot(self.lineselection_processor.x_main, self.lineselection_processor.y_main, label=", ".join([self.moltree.selected_tex_molecules[i] for i in self.file_loader.molecule_names]), picker=5)
        ax.hold(True)
        ax.plot(self.lineselection_processor.x_select, self.lineselection_processor.y_select, label=", ".join([self.moltree_selected.selected_tex_molecules[i] for i in self.file_loader_selected.molecule_names]), picker=5)
        if config["unit_wave"]:
            ax.set_xlabel(r"$Wavenumber\;\nu\;[cm^{-1}]$")
        else:
            ax.set_xlabel(r"$Wavelength\;\lambda\;[nm]$")
        if config["unit"]:
            ax.set_ylabel(r"$Absorption\;coefficient\;[cm{^2}/molecule]$")
        else:
            ax.set_ylabel(r"$Absorption\;coefficient\;[cm^{-1}]$")
        if config["unit_pressure"]:
            pressure_str=str("%10.2f"%(config["pressure"][1])).strip()+" atm"
        else:
            pressure_str=str("%10.2f"%(config["pressure"][0])).strip() + " bar"
        if config["unit_temperature"]:
            temperature_str=str(config["temperature"][1]) + " K"
        else:
            temperature_str=str("%10.2f"%(config["temperature"][0])).strip() + " °C"
        ax.set_title(r"$Lineselection\;"+", ".join([pressure_str, temperature_str]).replace(" ","\;")+"$")
        ax.legend(framealpha=config["legend_framealpha"], loc=config["legend_loc"])
        self.toolbar._actions["legend_switch"].setChecked(True)
        self.figure.tight_layout()
        ax.grid(True)
        ax.autoscale_view(True, True, False)
        ax.get_xaxis().get_major_formatter().set_useOffset(False)
        self.canvas.draw()
        ax.hold(self.toolbar._actions["hold_figure"].isChecked())

##    file loading specific functions
    def error(self, message):
        self.box("Error", message)
        
    def log_print_select(self, message, IDS_to_uncheck=[]):
        if config["uncheck_missing"]:
            for i in IDS_to_uncheck:
                self.moltree_selected.selected_ISO_ID[i].setCheckState(0, 0)
        with open("txt/logfile.txt", "a") as f:
            f.write(str(datetime.datetime.now())+" "+message+"\n")
        f.close()

    def log_print_main(self, message, IDS_to_uncheck=[]):
        if config["uncheck_missing"]:
            for i in IDS_to_uncheck:
                self.moltree.selected_ISO_ID[i].setCheckState(0, 0)
        with open("txt/logfile.txt", "a") as f:
            f.write(str(datetime.datetime.now())+" "+message+"\n")
        f.close()
        
    def line_counter(self, RowID, nline):
        self.progress_label.setText("processing data")
        self.progressBar.setRange(0,nline)
        self.progressBar.setValue(RowID)

    def terminate(self):
        self.progress_label.setText("cancel process")
        self.file_loader.stop()
        self.file_loader_selected.stop()
        self.lineselection_processor.stop()

    def closeEvent(self, event):
        json.dump(config, open("txt/config.txt", "w"))
        event.accept()

    def box(self, title, text):
        QtGui.QMessageBox.critical(self, title, text, QtGui.QMessageBox.Ok)

##    load HITRAN/HITEMP data for selected molecule
    def load_data_selected(self):
        self.file_loader_selected.table=None
        self.file_loader.table=None
        if len(self.moltree_selected.selected_tex_molecules) > 1:
            self.box("Error", "Please select only one molecule for Lineselection.")
        elif not self.moltree.selected_ISO_ID:
            self.box("Error", "Please select a molecule.")
        elif self.file_loader_selected.isRunning() or self.file_loader.isRunning() or self.lineselection_processor.isRunning():
            self.box("Warning", "Wait until the calculation has finished.")
        elif config["vmin"][1] >= config["vmax"][1]:
            self.box("Warning", "Please check wave range.")
        else:
            ISO_IDS = list(self.moltree_selected.selected_ISO_ID.keys())
##            print(ISO_IDS)
            self.file_loader_selected.update("lineselection_selected", ISO_IDS, self.moltree_selected.concentrations, "absorptioncoefficient", config["profile"])
            self.progress_label.setText("loading data")
            self.progressBar.setRange(0,0)
            self.progress_label.show()
            self.progressBar.show()
            self.terminate_button.show()
            self.file_loader_selected.finished.connect(self.load_spectrum_data)
            self.file_loader_selected.start()

##    load HITRAN/HITEMP data for main molecules               
    def load_spectrum_data(self):
        self.file_loader_selected.finished.disconnect()
        ISO_IDS = list(self.moltree.selected_ISO_ID.keys())
##        print(ISO_IDS)
        self.file_loader.update("lineselection_main", ISO_IDS, self.moltree.concentrations, "absorptioncoefficient", config["profile"])
        self.progress_label.setText("loading data")
        self.progressBar.setRange(0,0)
        self.progress_label.show()
        self.progressBar.show()
        self.terminate_button.show()
        self.file_loader.finished.connect(self.process_lineselection)
        self.file_loader.start()

##    start lineprocessing thread
    def process_lineselection(self):
        self.progress_label.hide()
        self.progressBar.hide()
        self.terminate_button.hide()
        self.file_loader.finished.disconnect()
        if self.file_loader.runs:
            self.progress_label.setText("loading data")
            self.progressBar.setRange(0,0)
            self.progress_label.show()
            self.progressBar.show()
            self.terminate_button.show()
            
            self.lineselection_processor.x_main=self.file_loader.x
            self.lineselection_processor.y_main=self.file_loader.y
            self.lineselection_processor.x_select=self.file_loader_selected.x
            self.lineselection_processor.y_select=self.file_loader_selected.y
            self.lineselection_processor.concentrations=self.file_loader_selected.concentrations
            self.lineselection_processor.IDS=self.file_loader_selected.IDS

            self.lineselection_processor.finished.connect(self.plot_lineselection)
            self.lineselection_processor.start()

##    plotting finished lineselection and fill table with calculated parameters
    def plot_lineselection(self):
        self.progress_label.hide()
        self.progressBar.hide()
        self.terminate_button.hide()
        if self.lineselection_processor.runs:
            self.lineselection_processor.finished.disconnect()
            ax=self.figure.gca()
            if not config["unit_wave"]:
                self.lineselection_processor.x_main=[1e7/i for i in reversed(self.lineselection_processor.x_main)]
                self.lineselection_processor.y_main=np.array([i for i in reversed(self.lineselection_processor.y_main)])
            ax.plot(self.lineselection_processor.x_main, self.lineselection_processor.y_main, label=", ".join([self.moltree.selected_tex_molecules[i] for i in self.file_loader.molecule_names]), picker=5)
            ax.hold(True)
            plot_x=[]
            plot_y=[]
            for i, j in self.lineselection_processor.lines_selected.items():
                if not config["unit_wave"]:
                    x=[1e7/i for i in reversed(j[0])]
                    y=np.array([i for i in reversed(j[1])])
                else:
                    x=j[0]
                    y=j[1]
##                ax.plot(x, y, color="k", picker=5)
                plot_x+=x
                plot_y+=list(y)
            for i, j in self.lineselection_processor.lines_for_overlap.items():
                if not config["unit_wave"]:
                    x=[1e7/i for i in reversed(j[0])]
                    y=np.array([i for i in reversed(j[1])])
                else:
                    x=j[0]
                    y=j[1]
##                ax.plot(x, y, color="k", picker=5)
                plot_x+=x
                plot_y+=list(y)
            ax.plot(plot_x, plot_y, color="k", label=", ".join([self.moltree_selected.selected_tex_molecules[i] for i in self.file_loader_selected.molecule_names]), picker=5)
            ax.grid(True)
            if config["unit_wave"]:
                ax.set_xlabel(r"$Wavenumber\;\nu\;[cm^{-1}]$")
            else:
                ax.set_xlabel(r"$Wavelength\;\lambda\;[nm]$")
            if config["unit"]:
                ax.set_ylabel(r"$Absorption\;coefficient\;[cm{^2}/molecule]$")
            else:
                ax.set_ylabel(r"$Absorption\;coefficient\;[cm^{-1}]$")
            if config["unit_pressure"]:
                pressure_str=str("%10.2f"%(config["pressure"][1])).strip()+" atm"
            else:
                pressure_str=str("%10.2f"%(config["pressure"][0])).strip() + " bar"
            if config["unit_temperature"]:
                temperature_str=str(config["temperature"][1]) + " K"
            else:
                temperature_str=str("%10.2f"%(config["temperature"][0])).strip() + " °C"
            ax.set_title(r"$Lineselection\;"+", ".join([pressure_str, temperature_str]).replace(" ","\;")+"$")
            ax.legend(framealpha=config["legend_framealpha"], loc=config["legend_loc"])
            self.toolbar._actions["legend_switch"].setChecked(True)
            self.figure.tight_layout()
            ax.autoscale_view(True, True, False)
            ax.get_xaxis().get_major_formatter().set_useOffset(False)
            self.canvas.draw()
            ax.hold(self.toolbar._actions["hold_figure"].isChecked())

            self.parent.tableWidget_lineselection.clearContents()
            self.parent.tableWidget_lineselection.setRowCount(0)
            self.parent.tableWidget_lineselection.setSortingEnabled(0)

##            QTableWidget with selected lines
            for i,j in reversed(sorted(self.lineselection_processor.lines_selected.items())):
                if not config["overlap_filter"][0] or j[9] <= config["overlap_filter"][1]/100:
                    row = self.parent.tableWidget_lineselection.rowCount()
                    self.parent.tableWidget_lineselection.setRowCount(row+1)
                    self.parent.tableWidget_lineselection.setCellWidget(row, 0, QtGui.QLabel( MoleculeTree.transform_html(self.file_loader_selected, self.file_loader_selected.molecule_names[0]) ))
                    self.parent.tableWidget_lineselection.setItem(row, 1, QtGui.QTableWidgetItem( str(j[3]).strip()))
                    self.parent.tableWidget_lineselection.setItem(row, 2, QCustomTableWidgetItem( "%-10.4G"%j[5]))
                    if not config["unit_wave"]:
                        rounded = str("%10.2f"%(1e7/j[4])).strip()
##                        plot_clicked[rounded] = j[4]
                    else:
                        rounded = str("%10.2f"%j[4]).strip()
##                        plot_clicked[rounded] = j[4]
                    self.parent.tableWidget_lineselection.setItem(row, 3, QtGui.QTableWidgetItem( rounded))
                    self.parent.tableWidget_lineselection.setItem(row, 4, QCustomTableWidgetItem( "%-10.4G"%j[6]))
                    self.parent.tableWidget_lineselection.setItem(row, 5, QCustomTableWidgetItem( str("%10.2f"%(j[8] * 100)).strip()))
                    self.parent.tableWidget_lineselection.setItem(row, 6, QCustomTableWidgetItem( str("%10.2f"%(j[9] * 100)).strip()))
                    
            self.parent.tableWidget_lineselection.setSortingEnabled(1)
##Lineprocessing thread to calculate overlay and overlap
class QtLineprocessor(QtCore.QThread):
    def __init__(self,parent=None):
        QtCore.QThread.__init__(self,parent)
        self.parent=parent

        self.x_main=None
        self.y_main=None
        self.x_select=None
        self.y_select=None
        self.concentrations=None
        self.IDS=None
        self.lines_selected=None
        self.lines_for_overlap=None

        self.runs = None

        self.error_signal = QtCore.SIGNAL("error")
        self.line_counter_signal = QtCore.SIGNAL("line_counter")
        self.line_full_overlap_signal = QtCore.SIGNAL("full_overlap")

    def absorptioncoef(self):
        if config["profile"]=="Voigt":
            return hapi.absorptionCoefficient_Voigt(Concentrations=self.concentrations, IDS=self.IDS, SourceTables='lineselection_selected_line', OmegaRange=(config["vmin"][1],config["vmax"][1]), OmegaStep=config["omega_step"], HITRAN_units=config["unit"], GammaL=config["broadening"], Environment={'p':config["pressure"][1],'T':config["temperature"][1]})
        elif config["profile"]=="Lorentz":
            return hapi.absorptionCoefficient_Lorentz(Concentrations=self.concentrations, IDS=self.IDS, SourceTables='lineselection_selected_line', OmegaRange=(config["vmin"][1],config["vmax"][1]), OmegaStep=config["omega_step"], HITRAN_units=config["unit"], GammaL=config["broadening"], Environment={'p':config["pressure"][1],'T':config["temperature"][1]})
        elif config["profile"]=="Gauss":
            return hapi.absorptionCoefficient_Doppler(Concentrations=self.concentrations, IDS=self.IDS, SourceTables='lineselection_selected_line', OmegaRange=(config["vmin"][1],config["vmax"][1]), OmegaStep=config["omega_step"]/10, HITRAN_units=config["unit"], GammaL=config["broadening"], Environment={'p':config["pressure"][1],'T':config["temperature"][1]})



    def run(self):
        try:
            self.runs=True
            self.lines_selected ={}
            self.lines_for_overlap={}
            line_temp=hapi.getColumn('lineselection_selected',"nu")
            counter=0
            if not config["unit_wave"]:
                self.x_main=[1e7/i for i in reversed(self.x_main)]
                self.y_main=np.array([i for i in reversed(self.y_main)])
##            calculate parameters for single lines
            for j in line_temp:
                if not self.runs:
                    return
                hapi.select('lineselection_selected', 'lineselection_selected_line', Conditions=('between','nu',j,j))
                nu_line, coef_line = self.absorptioncoef()
                nu_line=list(nu_line)
                max_coef_line = max(coef_line)
                if not config["spectrum_filter"][0] or max_coef_line > config["spectrum_filter"][1]/100 * coef_interp( nu_line[coef_line.index(max_coef_line)]):
                    coef_line_interp = interp1d(nu_line, coef_line)
                    integ_line=np.trapz(coef_line, nu_line)
                    if any([1 for i in coef_line-self.y_main if i >= 0]):
##                        calculate area of the line which is overlayed by main spectrum
                        counter_2=0
                        coef_up=[]
                        for u in coef_line:
                            if u > self.y_main[counter_2]:
                                coef_up.append(self.y_main[counter_2])
                            else:
                                coef_up.append(u)
                            counter_2+=1
##                        if hapi.getColumn('lineselection_selected_line',"local_lower_quanta")[0] =="     P  9e     ":
##                            self.parent.figure.gca().fill_between([1e7/i for i in reversed(self.x_main)], 0, [i for i in reversed(coef_up)])
##                            self.parent.figure.gca().fill_between(self.x_main, 0, coef_up)
##                            
                        integ_up=np.trapz(coef_up, dx=0.01)
                        area=integ_line-integ_up
                        if config["noise_filter"][0] and integ_up/integ_line > config["noise_filter"][1]/100:
##                            dictionary:
##                                {overlayfree area:[x-values,
##                                                   y-values,
##                                                   interpolationfunction,
##                                                   linename,
##                                                   x-value of wave peak,
##                                                   lineintensity (peak),
##                                                   overlayfree area,
##                                                   full area of line,
##                                                   overlay,
##                                                   overlap,
##                                                   moleculename]}
                            self.lines_for_overlap[area]=[nu_line, coef_line, coef_line_interp, hapi.getColumn('lineselection_selected_line',"local_lower_quanta")[0], nu_line[list(coef_line).index(max_coef_line)], max_coef_line, area, integ_line, integ_up/integ_line, 100]
                        else:
                            self.lines_selected[area]=[nu_line, coef_line, coef_line_interp, hapi.getColumn('lineselection_selected_line',"local_lower_quanta")[0], nu_line[list(coef_line).index(max_coef_line)], max_coef_line, area, integ_line, integ_up/integ_line, 100]
                    else:
                        self.lines_for_overlap[integ_line]=[nu_line, coef_line, coef_line_interp, hapi.getColumn('lineselection_selected_line',"local_lower_quanta")[0], nu_line[list(coef_line).index(max_coef_line)], max_coef_line, 0, integ_line, 100, 100]
                counter+=1
                if counter % 2 == 0:
                    self.emit(self.line_counter_signal, counter, len(line_temp), "line filtering")
            if not self.lines_selected:
                self.emit(self.line_full_overlap_signal)
                self.runs=False
                return
            else:
                self.lines_for_overlap.update(self.lines_selected)
                counter=0
                for i,j in self.lines_selected.items():
                    if not self.runs:
                        return
                    overlap_coef=[]
                    counter_2=0
                    highest_coef=list(zip(*[n[1] for m,n in self.lines_for_overlap.items() if m!=i]))
                    for i in highest_coef:
                        maximum=max(i)
                        if maximum >= j[1][counter_2] and i != 0:
                            overlap_coef.append(j[1][counter_2])
                        else:
                            overlap_coef.append(maximum)
                        counter_2+=1
##                      calculate overlapped area and overlap value
                    overlap_area=np.trapz(overlap_coef, dx=0.01)
                    j[9]=overlap_area/j[7]
                    counter+=1
                    if counter % 2 == 0:
                        self.emit(self.line_counter_signal, counter, len(self.lines_selected), "overlap calculation")
        except Exception as e:
            self.emit(self.error_signal, str(e))
            return

    def stop(self):
        if self.runs:
            self.runs=False
        
        

##HITRAN/HITEMP Tab with functionality to plot stick spectrum, absorption coefficient and absorption/transmittance/radiance spectra
class Spectrum_Tab(QtGui.QWidget):
    def __init__(self, parent=None, tabLayout=None):
        super(Spectrum_Tab, self).__init__()
        uic.loadUi("ui_files/spectrum_buttons.ui", self)
        self.parent = parent
        self.tabLayout = tabLayout
        
        self.profiles    = {"Lorentz":self.radioButton_lorentz, "Gauss":self.radioButton_gauss, "Voigt":self.radioButton_voigt}

        for i,j in self.profiles.items():
            j.toggled.connect(partial(self.set_profile, i, j))
        self.profiles[config["profile_hitran_tab"]].setChecked(1)

        self.spectrum_widget=QtGui.QWidget(self.parent)
        self.hitran_layout=QtGui.QGridLayout(self.spectrum_widget)
##        self.button_layout=QtGui.QGridLayout(self)
        self.tabLayout.addWidget(self.spectrum_widget, 0, 0)
        self.tabLayout.addWidget(self, 1, 0)


        self.statusBar = QtGui.QStatusBar(self)
        self.statusBar.setSizeGripEnabled(False)
        self.progress_label = QtGui.QLabel("")
        self.progressBar = QtGui.QProgressBar(self)
        self.terminate_button = QtGui.QToolButton()
        self.terminate_button.setIcon(QtGui.QIcon("images/cancel.png"))
        self.terminate_button.setAutoRaise(1)
        self.terminate_button.hide()
        self.progressBar.hide()
        self.progress_label.hide()
        self.progressBar.setMaximumSize(QtCore.QSize(200,20))
        self.statusBar.addPermanentWidget(self.progress_label)
        self.statusBar.addPermanentWidget(self.progressBar)
        self.statusBar.addPermanentWidget(self.terminate_button)
        self.tabLayout.addWidget(self.statusBar, 2, 0)
        
        
        self.moltree=MoleculeTree(self, select="hitran")
        self.moltreeWidget=self.moltree.tree
        self.hitran_layout.addWidget(self.moltreeWidget,0,0)

        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self, self.figure)
        self.spectrum_layout = QtGui.QVBoxLayout()
        self.spectrum_layout.addWidget(self.toolbar)
        self.spectrum_layout.addWidget(self.canvas)
        
        self.browser = PointBrowser(self, self.toolbar)
        self.canvas.mpl_connect('pick_event', self.browser.onpick)
        self.canvas.mpl_connect('resize_event', self.browser.resize)

        self.hitran_layout.addLayout(self.spectrum_layout,0,1)


        self.file_loader = QtFileLoader(self)
        
        self.pushButton_linespectra.clicked.connect(partial(self.load_data, "linespectrum", self.plot_linespectrum))
        self.pushButton_absorptioncoef.clicked.connect(partial(self.load_data, "absorptioncoefficient", self.plot_absorptioncoef))

        self.pushButton_absorptionspec.clicked.connect(partial(self.load_data, "absorptionspectrum", self.plot_absorption_spectrum))
        self.pushButton_transmittancespec.clicked.connect(partial(self.load_data, "transmittancespectrum", self.plot_transmittance_spectrum))
        self.pushButton_radiancespec.clicked.connect(partial(self.load_data, "radiancespectrum", self.plot_radiance_spectrum))
        
        self.terminate_button.clicked.connect(self.terminate)

                        
        self.connect(self.file_loader, self.file_loader.line_counter_signal, self.line_counter)

        self.connect(self.file_loader, self.file_loader.error_signal, self.error)
        self.connect(self.file_loader, self.file_loader.log_signal, self.log_print)

##    file loader functions
    def error(self, message):
        self.box("Error", message)
        
    def log_print(self, message, IDS_to_uncheck=[]):
        if config["uncheck_missing"]:
            for i in IDS_to_uncheck:
                self.moltree.selected_ISO_ID[i].setCheckState(0, 0)
        with open("txt/logfile.txt", "a") as f:
            f.write(str(datetime.datetime.now())+" "+message+"\n")
        f.close()
        
    def line_counter(self, RowID, nline):
        self.progress_label.setText("processing data")
        self.progressBar.setRange(0,nline)
        self.progressBar.setValue(RowID)

    def terminate(self):
        self.progress_label.setText("cancel process")
        self.file_loader.stop()

    def closeEvent(self, event):
        json.dump(config, open("txt/config.txt", "w"))
        event.accept()

    def box(self, title, text):
        QtGui.QMessageBox.critical(self, title, text, QtGui.QMessageBox.Ok)

##    load HITRAN/HITEMP data in new QThread
    def load_data(self, spectrum, function):
        if not self.moltree.selected_ISO_ID:
            self.box("Error", "Please select a molecule.")
        elif self.file_loader.isRunning():
            self.box("Warning", "Wait until the calculation has finished.")
        elif config["vmin"][1] >= config["vmax"][1]:
            self.box("Warning", "Please check wave range.")
        else:
            ISO_IDS = list(self.moltree.selected_ISO_ID.keys())
            self.file_loader.update("hitran_tab", ISO_IDS, self.moltree.concentrations, spectrum, config["profile_hitran_tab"])
            if self.file_loader.new_data:
                self.progress_label.setText("loading data")
                self.progressBar.setRange(0,0)
                self.progress_label.show()
                self.progressBar.show()
                self.terminate_button.show()
                self.file_loader.finished.connect(function)
                self.file_loader.start()

        
##    plot stick spectrum
    def plot_linespectrum(self):
        self.file_loader.finished.disconnect()
        self.progress_label.hide()
        self.progressBar.hide()
        self.terminate_button.hide()
        if self.file_loader.processed:
            ax=self.figure.gca()
            ax.plot(self.file_loader.x,self.file_loader.y,label=", ".join([self.moltree.selected_tex_molecules[i] for i in self.file_loader.molecule_names]), picker=5)
            if config["unit_wave"]:
                ax.set_xlabel(r"$Wavenumber\;\nu\;[cm^{-1}]$")
            else:
                ax.set_xlabel(r"$Wavelength\;\lambda\;[nm]$")
            ax.set_ylabel(r"$Intensity\;\;[cm^{^-1}/molecule\;cm^{-2}]$")
            if config["unit_temperature"]:
                ax.set_title(r"$Stick\;spectrum\;296K$")
            else:
                ax.set_title(r"$Stick\;spectrum\;25^\circ C$")
            ax.grid(True)
            ax.legend(framealpha=config["legend_framealpha"], loc=config["legend_loc"])
            self.toolbar._actions["legend_switch"].setChecked(True)
            self.figure.tight_layout()
            ax.autoscale_view(True, True, False)
            ax.get_xaxis().get_major_formatter().set_useOffset(False)
            self.canvas.draw()

##    plot absorption coefficient
    def plot_absorptioncoef(self):
        self.file_loader.finished.disconnect()
        self.progress_label.hide()
        self.progressBar.hide()
        self.terminate_button.hide()
        if self.file_loader.processed:
            ax=self.figure.gca()
            ax.plot(self.file_loader.x, self.file_loader.y, label=", ".join([self.moltree.selected_tex_molecules[i] for i in self.file_loader.molecule_names]), picker=5)
            if config["unit_wave"]:
                ax.set_xlabel(r"$Wavenumber\;\nu\;[cm^{-1}]$")
            else:
                ax.set_xlabel(r"$Wavelength\;\lambda\;[nm]$")
            if config["unit"]:
                ax.set_ylabel(r"$Absorption\;coefficient\;[cm{^2}/molecule]$")
            else:
                ax.set_ylabel(r"$Absorption\;coefficient\;[cm^{-1}]$")
            if config["unit_pressure"]:
                pressure_str=str("%10.2f"%(config["pressure"][1])).strip()+" atm"
            else:
                pressure_str=str("%10.2f"%(config["pressure"][0])).strip() + " bar"
            if config["unit_temperature"]:
                temperature_str=str(config["temperature"][1]) + " K"
            else:
                temperature_str=str("%10.2f"%(config["temperature"][0])).strip() + " °C"
            ax.set_title(r"$Absorption\;coefficient\;"+", ".join([pressure_str, temperature_str]).replace(" ","\;")+"$")
            ax.grid(True)
            ax.legend(framealpha=config["legend_framealpha"], loc=config["legend_loc"])
            self.toolbar._actions["legend_switch"].setChecked(True)
            self.figure.tight_layout()
            ax.autoscale_view(True, True, False)
            ax.get_xaxis().get_major_formatter().set_useOffset(False)
            self.canvas.draw()

##    plot absorption spectrum
    def plot_absorption_spectrum(self):
        self.file_loader.finished.disconnect()
        self.progress_label.hide()
        self.progressBar.hide()
        self.terminate_button.hide()
        if self.file_loader.processed:
            ax=self.figure.gca()
            ax.plot(self.file_loader.x,self.file_loader.y,label=", ".join([self.moltree.selected_tex_molecules[i] for i in self.file_loader.molecule_names]), picker=5)
            if config["unit_wave"]:
                ax.set_xlabel(r"$Wavenumber\;\nu\;[cm^{-1}]$")
            else:
                ax.set_xlabel(r"$Wavelength\;\lambda\;[nm]$")
            ax.set_ylabel(r"$Absorption\;[\%]$")
            if config["unit_pressure"]:
                pressure_str=str("%10.2f"%(config["pressure"][1])).strip()+" atm"
            else:
                pressure_str=str("%10.2f"%(config["pressure"][0])).strip() + " bar"
            if config["unit_temperature"]:
                temperature_str=str(config["temperature"][1]) + " K"
            else:
                temperature_str=str("%10.2f"%(config["temperature"][0])).strip() + " °C"
            ax.set_title(r"$Absorption\;spectrum\;$"+"$"+", ".join([pressure_str, temperature_str]).replace(" ","\;")+"$")
            ax.grid(True)
            ax.legend(framealpha=config["legend_framealpha"], loc=config["legend_loc"])
            self.toolbar._actions["legend_switch"].setChecked(True)
            self.figure.tight_layout()
            ax.autoscale_view(True, True, False)
            ax.get_xaxis().get_major_formatter().set_useOffset(False)
            self.canvas.draw()

##    plot transmittance spectrum    
    def plot_transmittance_spectrum(self):
        self.file_loader.finished.disconnect()
        self.progress_label.hide()
        self.progressBar.hide()
        self.terminate_button.hide()
        if self.file_loader.processed:
            ax=self.figure.gca()
            ax.plot(self.file_loader.x,self.file_loader.y,label=", ".join([self.moltree.selected_tex_molecules[i] for i in self.file_loader.molecule_names]), picker=5)
            if config["unit_wave"]:
                ax.set_xlabel(r"$Wavenumber\;\nu\;[cm^{-1}]$")
            else:
                ax.set_xlabel(r"$Wavelength\;\lambda\;[nm]$")
            ax.set_ylabel(r"$Transmittance\;[\%]$")
            if config["unit_pressure"]:
                pressure_str=str("%10.2f"%(config["pressure"][1])).strip()+" atm"
            else:
                pressure_str=str("%10.2f"%(config["pressure"][0])).strip() + " bar"
            if config["unit_temperature"]:
                temperature_str=str(config["temperature"][1]) + " K"
            else:
                temperature_str=str("%10.2f"%(config["temperature"][0])).strip() + " °C"
            ax.set_title(r"$Transmittance\;spectrum\;$"+"$"+", ".join([pressure_str, temperature_str]).replace(" ","\;")+"$")
            ax.grid(True)
            ax.legend(framealpha=config["legend_framealpha"], loc=config["legend_loc"])
            self.toolbar._actions["legend_switch"].setChecked(True)
            self.figure.tight_layout()
            ax.autoscale_view(True, True, False)
            ax.get_xaxis().get_major_formatter().set_useOffset(False)
            self.canvas.draw()

##    plot radiance spectrum
    def plot_radiance_spectrum(self):
        self.file_loader.finished.disconnect()
        self.progress_label.hide()
        self.progressBar.hide()
        self.terminate_button.hide()
        if self.file_loader.processed:
            ax=self.figure.gca()
            ax.plot(self.file_loader.x,self.file_loader.y,label=", ".join([self.moltree.selected_tex_molecules[i] for i in self.file_loader.molecule_names]), picker=5)
            if config["unit_wave"]:
                ax.set_xlabel(r"$Wavenumber\;\nu\;[cm^{-1}]$")
            else:
                ax.set_xlabel(r"$Wavelength\;\lambda\;[nm]$")
            ax.set_ylabel(r"$Radiance\;[W/sr/cm^{2}/cm^{-1}]$")
            if config["unit_pressure"]:
                pressure_str=str("%10.2f"%(config["pressure"][1])).strip()+" atm"
            else:
                pressure_str=str("%10.2f"%(config["pressure"][0])).strip() + " bar"
            if config["unit_temperature"]:
                temperature_str=str(config["temperature"][1]) + " K"
            else:
                temperature_str=str("%10.2f"%(config["temperature"][0])).strip() + " °C"
            ax.set_title(r"$Radiance\;spectrum\;$"+"$"+", ".join([pressure_str, temperature_str]).replace(" ","\;")+"$")
            ax.grid(True)
            ax.legend(framealpha=config["legend_framealpha"], loc=config["legend_loc"])
            self.toolbar._actions["legend_switch"].setChecked(True)
            self.figure.tight_layout()
            ax.autoscale_view(True, True, False)
            ax.get_xaxis().get_major_formatter().set_useOffset(False)
            self.canvas.draw()
            
    def set_profile(self, i, j, state):
        if state:
            config["profile_hitran_tab"]=i

##file loading thread for HITRAN/HITEMP data
class QtFileLoader(QtCore.QThread):
    def __init__(self,parent=None):
        QtCore.QThread.__init__(self,parent)

        self.vmin=None
        self.vmax=None
        self.unit=None
        self.profile=None
        self.pressure=None
        self.temperature=None
        self.broadening=None
        self.path_length =None
        self.hitemp=None
        self.concentrations={}

        self.spectrum = None
        
        self.table = None
        self.IDS = []
        self.loaded_IDS = []
        self.molecule_names = None
        self.runs = None

        self.reprocess_data = None
        self.reload_data = None

        self.new_data = True
        self.x = None
        self.y = None
        self.processed=False

        self.error_signal = QtCore.SIGNAL("error")
        self.log_signal = QtCore.SIGNAL("log")
        self.line_counter_signal = QtCore.SIGNAL("line_counter")

    def run(self):
        try:
            if self.new_data:
                if not self.runs:
                    try:
                        hapi.dropTable(self.table)
                    except:
                        pass
                self.runs=True
                if self.reload_data:
                    if config["online"]:
                        if config["hitemp"]:
                            hitemp_molecules=[x for x in ["CO","CO2","H2O","NO","OH","SF6","ClONO2","CF4"] if x in self.molecule_names]
##                            print("hitemp")
##                            print(hitemp_molecules)
                            online_molecules=[x for x in self.molecule_names if x not in hitemp_molecules]
                            if hitemp_molecules:
                                hapi.load_data_offline(self.table, self.vmin, self.vmax, hitemp_molecules, self.IDS, config["hitemp"], self)
                                if online_molecules and self.runs:
##                                    print("online molecules")
##                                    print(online_molecules)
                                    hapi.fetch_by_ids(self.table, [x for x in self.IDS if hapi.moleculeName(hapi.ISO_ID[x][0]) in online_molecules], self.vmin, self.vmax, self, wipe=False)
                            else:
                                hapi.fetch_by_ids(self.table, self.IDS, self.vmin, self.vmax, self)
                        else:
                            hapi.fetch_by_ids(self.table, self.IDS, self.vmin, self.vmax, self)
                    else:
                        hapi.load_data_offline(self.table, self.vmin, self.vmax, self.molecule_names, self.IDS, config["hitemp"], self)
                        
                if self.runs:
                    if not hapi.LOCAL_TABLE_CACHE[self.table]['header']['number_of_rows']:
                        raise Exception("No data in given range.")
                    else:
                        ISO_in_table=set([hapi.ISO[i][0] for i in zip(hapi.LOCAL_TABLE_CACHE[self.table]["data"]["molec_id"],hapi.LOCAL_TABLE_CACHE[self.table]["data"]["local_iso_id"])])

                        for i in self.molecule_names.copy():
                            if hapi.moleculeID(i) not in set(hapi.LOCAL_TABLE_CACHE[self.table]["data"]["molec_id"]):
                                self.molecule_names.remove(i)
                            
                        missing_IDS=[i for i in self.IDS if i not in ISO_in_table]
                        if missing_IDS:
                            not_in_table={}
                            for i in missing_IDS:
                                if config["uncheck_missing"]:
                                    self.IDS.remove(i)
                                molecule_name=hapi.ISO_ID[i][5]
                                isoname=hapi.ISO_ID[i][2]
                                if molecule_name not in not_in_table:
                                    not_in_table[molecule_name]=[]
                                    not_in_table[molecule_name].append(isoname)
                                else:
                                    not_in_table[molecule_name].append(isoname)
                            log_text=[]
                            for i,j in not_in_table.items():
                                log_text.append(i+" ("+", ".join(j)+")")
                            log_text="no data available for "+", ".join(log_text)+" in range "+str("%10.2f"%config["vmin"][0]).strip()+"nm - "+str("%10.2f"%config["vmax"][0]).strip()\
                                  +"nm ("+str("%10.2f"%config["vmin"][1]).strip()+"cm^-1 - "+str("%10.2f"%config["vmax"][1]).strip()+"cm^-1)"
                            self.emit(self.log_signal, log_text, missing_IDS)
                    if self.table and self.molecule_names and self.reprocess_data:
                        if self.spectrum is "linespectrum":
                            self.x, self.y = hapi.getStickXY(self.table, self)
                        elif self.spectrum is "absorptioncoefficient":
                            self.x, self.y = self.absorptioncoef(config["unit"])
                        elif self.spectrum is "absorptionspectrum":
                            self.x, self.y = self.absorptioncoef(False)
                            self.x, self.y = hapi.absorptionSpectrum(self.x, self.y, Environment={'l': config["path_length"]})
                            self.y=self.y*100
                        elif self.spectrum is "transmittancespectrum":
                            self.x, self.y = self.absorptioncoef(False)
                            self.x, self.y = hapi.transmittanceSpectrum(self.x, self.y, Environment={'l': config["path_length"]})
                            self.y=self.y*100
                        elif self.spectrum is "radiancespectrum":
                            self.x, self.y = self.absorptioncoef(False)
                            self.x, self.y = hapi.radianceSpectrum(self.x, self.y, Environment={'l': config["path_length"], 'T': config["temperature"][1]})
                        if self.runs:
                            if not config["unit_wave"]:
                                self.x = [1e7/i for i in reversed(self.x)]
                                self.y = [i for i in reversed(self.y)]
                            self.processed=True
                    
                    self.new_data=False
        except Exception as e:
            self.table=None
            self.emit(self.error_signal, str(e))

    def absorptioncoef(self, unit):
        if self.profile=="Voigt":
            return hapi.absorptionCoefficient_Voigt(Concentrations=self.concentrations, IDS=self.IDS, SourceTables=self.table, OmegaRange=(self.vmin,self.vmax), OmegaStep=config["omega_step"], HITRAN_units=unit, GammaL=config["broadening"], Environment={'p':config["pressure"][1],'T':config["temperature"][1]},thread=self)
        elif self.profile=="Lorentz":
            return hapi.absorptionCoefficient_Lorentz(Concentrations=self.concentrations, IDS=self.IDS, SourceTables=self.table, OmegaRange=(self.vmin,self.vmax), OmegaStep=config["omega_step"], HITRAN_units=unit, GammaL=config["broadening"], Environment={'p':config["pressure"][1],'T':config["temperature"][1]},thread=self)
        elif self.profile=="Gauss":
            return hapi.absorptionCoefficient_Doppler(Concentrations=self.concentrations, IDS=self.IDS, SourceTables=self.table, OmegaRange=(self.vmin,self.vmax), OmegaStep=config["omega_step"]/10, HITRAN_units=unit, GammaL=config["broadening"], Environment={'p':config["pressure"][1],'T':config["temperature"][1]},thread=self)


    def stop(self):
        if self.runs:
            self.runs=False
            self.table=None
        
    def update(self, table, ID_List, concentrations, spectrum, profile):
        temp_IDS=ID_List.copy()
        temp_molecules=[]
        temp_concentrations={}
        for i in temp_IDS.copy():
            temp_name=hapi.moleculeName(hapi.ISO_ID[i][0])
            if concentrations[temp_name] == 0:
                temp_IDS.remove(i)
            if temp_name not in temp_molecules:
                if concentrations[temp_name] != 0:
                    temp_molecules.append(temp_name)
                    temp_concentrations[temp_name] = concentrations[temp_name]

        new_molecules = not all([True if i in self.loaded_IDS else False for i in temp_IDS])
        if temp_IDS:
            new_concentrations=bool(self.concentrations != temp_concentrations)
            
            if not self.table or new_molecules or self.vmin != config["vmin"][1] or self.vmax != config["vmax"][1] or self.hitemp != config["hitemp"]:
                self.reload_data=True
                self.reprocess_data=True
                self.loaded_IDS = temp_IDS
                self.new_data=True
            elif set(self.IDS) != set(temp_IDS) or profile != self.profile or spectrum != self.spectrum or self.temperature != config["temperature"][1] or self.pressure != config["pressure"][1]\
                                             or new_concentrations or self.broadening != config["broadening"] or self.path_length != config["path_length"] or config["unit"] != self.unit:
                self.reload_data=False
                self.reprocess_data=True
                self.new_data=True
##                print(config["temperature"][1])
##                print(self.temperature)
##
##                print(config["pressure"][1])
##                print(self.pressure)
##
##                print(new_concentrations)
##
##                print(config["broadening"])
##                print(self.broadening)
##
##                print(config["path_length"])
##                print(self.path_length)
##
##                print(config["unit"])
##                print(self.unit)
                
            if self.new_data:
                self.concentrations=temp_concentrations
                self.processed=False
                self.molecule_names=temp_molecules
                self.IDS = temp_IDS
                self.profile = profile
                self.spectrum = spectrum
                self.unit = config["unit"]
                self.table = table
                self.hitemp = config["hitemp"]
                self.vmin = config["vmin"][1]
                self.vmax = config["vmax"][1]
                self.temperature = config["temperature"][1]
                self.pressure = config["pressure"][1]
                self.broadening = config["broadening"]
                self.path_length = config["path_length"]
        else:
            self.emit(self.error_signal, "All concentrations are zero.")


##Molecule List with concentrations
class MoleculeTree:
    def __init__(self,parent=None, select=None, main= None):
        self.data = None
        self.tree = QtGui.QTreeWidget(parent)
        self.tree.setItemDelegate(HTMLDelegate())
        self.tree.setHeaderLabels(["Molecule","Concentration"])
##        self.tree.header().setResizeMode(1, QtGui.QHeaderView.Stretch)
        self.tree.header().setStretchLastSection(True)
        self.molecule_full_names = json.load(open("txt/molecule_names.txt"))
        self.fill()
        self.selected_ISO_ID = {}
        self.selected_tex_molecules = {}
        self.concentrations_spinboxes = {}
        self.concentrations = {}
        self.select=select
        self.tree.itemChanged.connect(self.refresh)
        
    def fill(self):
        for i in range(1,48):
            if i != 30 and i != 35 and i != 42:
                moleculename=hapi.moleculeName(i)
                molecule = CustomTreeWidgetItem(self.tree)
                molecule.setToolTip(0, self.molecule_full_names[moleculename])
                molecule.setData(0, 32, moleculename)
                if config["sort_molecules_by_names"]:
                    molecule.setSortData(0, self.molecule_full_names[moleculename])
                else:
                    molecule.setSortData(0, moleculename)
                molecule.setFlags(molecule.flags() | QtCore.Qt.ItemIsTristate | QtCore.Qt.ItemIsUserCheckable)
                if config["full_names"]:
                    molecule.setText(0, "<font size=\""+str(config["font-size"])+"\">"+self.transform_html(moleculename)+" ("+self.molecule_full_names[moleculename]+")"+"</font>")
                    self.tree.setMaximumWidth(400)
                    self.tree.setMinimumWidth(400)
                else:
                    molecule.setText(0, "<font size=\""+str(config["font-size"])+"\">"+self.transform_html(moleculename)+"</font>")
                    self.tree.setMaximumWidth(280)
                    self.tree.setMinimumWidth(280)
                if config["sort_molecules_by_names"]:
                    molecule.setSortData(0, self.molecule_full_names[molecule.data(0,32)])
                else:
                    molecule.setSortData(0, molecule.data(0,32))
                temp=[[1,"A"],[2,"B"],[3,"C"],[4,"D"],[5,"E"],[6,"F"],[7,"G"],[8,"H"],[9,"I"],[0,"J"]]
                for j in temp:
                    try:
                        isoname=hapi.isotopologueName(i, j[0])
                        isotope = CustomTreeWidgetItem(molecule)
                        isotope.setSortData(0, j[1])
                        isotope.setData(0, 32, [hapi.ID(i, isoname),isoname])
                        isotope.setFlags(molecule.flags() | QtCore.Qt.ItemIsUserCheckable)
                        isotope.setText(0, "<font size=\""+str(config["font-size"])+"\">"+self.transform_html(isoname)+"</font>")
                        isotope.setCheckState(0, QtCore.Qt.Unchecked)
                    except KeyError:
                        pass
        self.tree.sortItems(0, 0)
##        json.dump(new_m, open("txt/molecule_names.txt", "w"))

    def refresh(self, item):
        if item.parent():
            final_item=item.parent()
            if not item.checkState(0):
                if item.data(0,32)[0] in self.selected_ISO_ID:
                    del self.selected_ISO_ID[item.data(0,32)[0]]
            else:
                self.selected_ISO_ID[item.data(0,32)[0]]=item
        else:
            final_item=item
                
        if final_item.checkState(0):
            if final_item.data(0,32) not in self.selected_tex_molecules:
                self.selected_tex_molecules[final_item.data(0,32)]=self.transform_texstyle(final_item.data(0,32))

                conc_value=100
                objects=[self]
                for obj in gc.get_objects():
                    if isinstance(obj, MoleculeTree):
                        if obj.select is self.select:
                            conc_value-=sum(obj.concentrations.values())
                            objects.append(obj)
                
                if not self.select:
                    conc_value=100-sum(self.concentrations.values())
                self.concentrations[final_item.data(0,32)]=conc_value
                spinBox=QtGui.QDoubleSpinBox()
                self.concentrations_spinboxes[final_item.data(0,32)]=spinBox
                spinBox.setMaximum(conc_value)
                spinBox.setValue(conc_value)
                spinBox.setMinimum(0)
                spinBox.setMaximumWidth(100)
                self.tree.setItemWidget(final_item, 1, spinBox)
        
                def update_conc_value(spinbox, objects, value):
                    temp_sum=sum([self.concentrations[i] for i,j in self.concentrations_spinboxes.items() if j is not spinbox])
                    for i in objects:
                        if i is not self:
                            temp_sum+=sum(i.concentrations.values())
                    for m in objects:
                        for i,j in m.concentrations_spinboxes.items():
                            if j == spinbox:
                                if i in m.concentrations:
                                    m.concentrations[i]=value
                                j.setMaximum(100-temp_sum)
                            else:
                                j.setMaximum(100-value)
                spinBox.valueChanged.connect(partial(update_conc_value, spinBox, objects))

        else:
            if final_item.data(0,32) in self.selected_tex_molecules:
                del self.selected_tex_molecules[final_item.data(0,32)]
                del self.concentrations[final_item.data(0,32)]
                self.concentrations_spinboxes[final_item.data(0,32)].setValue(0)
                del self.concentrations_spinboxes[final_item.data(0,32)]
                self.tree.removeItemWidget(final_item, 1)

##    molecule name in html style
    def transform_html(self, text):
        richtext = ""
        c=0
        while c < len(text):
            if text[c] == "(":
                temp=text[c+1:].split(")",2)
                number = "".join(x for x in temp[0] if x.isdigit())
                richtext+= "<sup>"+number+"</sup>"+temp[0][len(number)]
                c+=len(temp[0])+1
                if c+1 < len(text) and text[c+1].isdigit():
                    richtext+="<sub>"+text[c+1]+"</sub>"
                    c+=1
                c+=1
            else:
                richtext+=text[c]
                if c+1 < len(text):
                    if text[c+1].isdigit():
                        richtext+="<sub>"+text[c+1]+"</sub>"
                        c+=1
                    elif text[c+1]=="p":
                        richtext+="<sup>+</sup>"
                        c+=1
                c+=1
        return richtext

##    molecule name in texstyle
    def transform_texstyle(self, text):
        richtext = r"$"
        c=0
        while c < len(text):
            if text[c] == "(":
                temp=text[c+1:].split(")",2)
                number = "".join(x for x in temp[0] if x.isdigit())
                richtext+= "^{"+number+"}"+temp[0][len(number)]
                c+=len(temp[0])+1
                if c+1 < len(text) and text[c+1].isdigit():
                    richtext+="_{"+text[c+1]+"}"
                    c+=1
                c+=1
            else:
                richtext+=text[c]
                if c+1 < len(text):
                    if text[c+1].isdigit():
                        richtext+="_{"+text[c+1]+"}"
                        c+=1
                    elif text[c+1]=="p":
                        richtext+="^{+}"
                        c+=1
                c+=1
        return richtext+"$"
                

##Custom TableWIdgetItem for numerical sorting
class QCustomTableWidgetItem (QtGui.QTableWidgetItem):
    def __init__ (self, value):
        super(QCustomTableWidgetItem, self).__init__(str('%s' % value))

    def __lt__ (self, other):
        if (isinstance(other, QCustomTableWidgetItem)):
            selfDataValue  = float(self.data(QtCore.Qt.EditRole))
            otherDataValue = float(other.data(QtCore.Qt.EditRole))
            return selfDataValue < otherDataValue
        else:
            return QtGui.QTableWidgetItem.__lt__(self, other)


##Custom Delegate for QTreeWidget to be able to display html style text (for molecule names)
class HTMLDelegate(QtGui.QStyledItemDelegate):    
    def paint(self, painter, option, index):
        options = QtGui.QStyleOptionViewItemV4(option)
        self.initStyleOption(options,index)

        style = QtGui.QApplication.style() if options.widget is None else options.widget.style()

        doc = QtGui.QTextDocument()
        doc.setHtml(options.text)
        doc.setTextWidth(option.rect.width())

        options.text = ""
        style.drawControl(QtGui.QStyle.CE_ItemViewItem, options, painter);

        ctx = QtGui.QAbstractTextDocumentLayout.PaintContext()

        textRect = style.subElementRect(QtGui.QStyle.SE_ItemViewItemText, options)
        painter.save()
        painter.translate(textRect.topLeft())
        painter.setClipRect(textRect.translated(-textRect.topLeft()))
        doc.documentLayout().draw(painter, ctx)

        painter.restore()

    def sizeHint(self, option, index):
        options = QtGui.QStyleOptionViewItemV4(option)
        self.initStyleOption(options,index)

        doc = QtGui.QTextDocument(self)
        doc.setHtml(options.text)
        doc.setTextWidth(options.rect.width())
        return QtCore.QSize(doc.idealWidth(), max(doc.size().height(), options.decorationSize.height()))

##Custom TreeWidgetItem to set custom sort data
class CustomTreeWidgetItem(QtGui.QTreeWidgetItem):
    def __lt__( self, other ):
        if ( not isinstance(other, CustomTreeWidgetItem) ):
            return super(CustomTreeWidgetItem, self).__lt__(other)

        tree = self.treeWidget()
        if ( not tree ):
            column = 0
        else:
            column = tree.sortColumn()

        return self.sortData(column) < other.sortData(column)

    def __init__( self, *args ):
        super(CustomTreeWidgetItem, self).__init__(*args)
        self._sortData = {}

    def sortData( self, column ):
        return self._sortData.get(column, self.text(column))

    def setSortData( self, column, data ):
        self._sortData[column] = data

##Preference window. write values in config and GUI or discard changes if canceled
class prefWindow(QtGui.QWidget):
    def __init__(self, parent=None):
        global config_temp, config
        self.parent=parent
        super(prefWindow, self).__init__()
        uic.loadUi("ui_files/preferences.ui", self)
        config_temp=config.copy()

        self.unit_blocker=False
        
        if not config["unit_wave"]:
            self.label_wave_min.setText("Lambda min [nm]")
            self.label_wave_max.setText("Lambda max [nm]")
        if not config["unit_pressure"]:
            self.label_pressure.setText("Pressure [bar]")
        if not config["unit_temperature"]:
            self.label_temperature.setText("Temperature [°C]")
            self.doubleSpinBox_temperature.setRange(70-273.15, 3000-273.15)
            
        self.doubleSpinBox_wave_min.setValue(self.parent.doubleSpinBox_vmin.value())
        self.doubleSpinBox_wave_max.setValue(self.parent.doubleSpinBox_vmax.value())
        self.doubleSpinBox_path_length.setValue(config["path_length"])
        self.doubleSpinBox_noise_filter.setValue(config["noise_filter"][1])
        self.checkBox_noise_filter.setChecked(config["noise_filter"][0])
        self.doubleSpinBox_overlap_filter.setValue(config["overlap_filter"][1])
        self.checkBox_overlap_filter.setChecked(config["overlap_filter"][0])
        self.doubleSpinBox_pressure.setValue(self.parent.doubleSpinBox_pressure.value())
        self.doubleSpinBox_temperature.setValue(self.parent.doubleSpinBox_temperature.value())
        self.doubleSpinBox_spectrum_filter.setValue(config["spectrum_filter"][1])
        self.spinBox_font_size.setValue(config["font-size"])
        self.checkBox_spectrum_filter.setChecked(config["spectrum_filter"][0])
        self.doubleSpinBox_stretch.setValue(config["measure_range"])
        self.spinBox_frequency.setValue(config["frequency"])
        self.checkBox_hitemp.setChecked(config["hitemp"])
        self.spinBox_savgol_window.setValue(config["savgol_window_length"])
        self.spinBox_savgol_polyorder.setValue(config["savgol_polyorder"])
        self.spinBox_legend_framealpha.setValue(100-config["legend_framealpha"]*100)
        self.comboBox_legend_loc.setCurrentIndex(self.comboBox_legend_loc.findText(config["legend_loc"], QtCore.Qt.MatchFixedString))
        self.spinBox_baseline_polyorder.setValue(config["baseline_polyorder"])
        self.checkBox_automatic_baseline.setChecked(config["automatic_baseline"])


        self.spinBox_baseline_polyorder.valueChanged.connect(self.set_baseline_polyorder)
        self.checkBox_automatic_baseline.stateChanged.connect(self.set_automatic_baseline)
        self.spinBox_legend_framealpha.valueChanged.connect(self.set_legend_framealpha)
        self.spinBox_savgol_polyorder.valueChanged.connect(self.set_savgol_polyorder)
        self.spinBox_savgol_window.valueChanged.connect(self.set_savgol_window)
        self.doubleSpinBox_wave_min.valueChanged.connect(self.set_wave_min)
        self.doubleSpinBox_wave_max.valueChanged.connect(self.set_wave_max)
        self.doubleSpinBox_path_length.valueChanged.connect(self.set_path_length)
        self.doubleSpinBox_noise_filter.valueChanged.connect(self.set_noise_filter_value)
        self.doubleSpinBox_overlap_filter.valueChanged.connect(self.set_overlap_filter_value)
        self.doubleSpinBox_pressure.valueChanged.connect(self.set_pressure_value)
        self.doubleSpinBox_temperature.valueChanged.connect(self.set_temperature_value)
        self.doubleSpinBox_spectrum_filter.valueChanged.connect(self.set_spectrum_filter_value)
        self.doubleSpinBox_stretch.valueChanged.connect(self.set_stretch_factor_value)
        self.spinBox_frequency.valueChanged.connect(self.set_frequency_value)
        self.spinBox_font_size.valueChanged.connect(self.set_font_size_value)

        if config["broadening"]=="gamma_self":
            self.comboBox_broadening.setCurrentIndex(1)
        else:
            self.comboBox_broadening.setCurrentIndex(0)
        self.comboBox_profile.setCurrentIndex(self.comboBox_profile.findText(config["profile"], QtCore.Qt.MatchFixedString))

        self.comboBox_wave.setCurrentIndex(config["unit_wave"])
        self.comboBox_pressure.setCurrentIndex(config["unit_pressure"])
        self.comboBox_temperature.setCurrentIndex(config["unit_temperature"])
        self.comboBox_mode.setCurrentIndex(config["online"])
        self.comboBox_coef.setCurrentIndex(config["unit"])
        self.checkBox_full_names.setChecked(config["full_names"])
        self.checkBox_start_index.setChecked(config["start_index_manually"])
        self.comboBox_sort_molecules.setCurrentIndex(config["sort_molecules_by_names"])

        self.comboBox_legend_loc.currentIndexChanged.connect(self.set_legend_loc)
        self.comboBox_wave.currentIndexChanged.connect(self.unit_wave)
        self.comboBox_pressure.currentIndexChanged.connect(self.unit_pressure)
        self.comboBox_temperature.currentIndexChanged.connect( self.unit_temperature)
        self.comboBox_mode.currentIndexChanged.connect( self.set_mode)
        self.comboBox_broadening.currentIndexChanged.connect( self.set_broadening)
        self.comboBox_profile.currentIndexChanged.connect( self.set_profile)
        self.checkBox_full_names.stateChanged.connect( self.set_full_names)
        self.checkBox_start_index.stateChanged.connect(self.set_start_index_manually)
        self.comboBox_coef.currentIndexChanged.connect( self.set_coef_unit)
        self.comboBox_sort_molecules.currentIndexChanged.connect(self.set_sort_molecule)

        self.checkBox_spectrum_filter.stateChanged.connect(self.set_spectrum_filter)
        self.checkBox_noise_filter.stateChanged.connect(self.set_noise_filter)
        self.checkBox_overlap_filter.stateChanged.connect(self.set_overlap_filter)
        self.checkBox_hitemp.stateChanged.connect(self.set_hitemp)
        
        self.pushButton_save.clicked.connect(self.save)
        self.pushButton_cancel.clicked.connect(self.close)


    def set_sort_molecule(self, index):
        config_temp["sort_molecules_by_names"]=index

    def set_baseline_polyorder(self, value):
        config_temp["baseline_polyorder"] = value
        
    def set_automatic_baseline(self, state):
        config_temp["automatic_baseline"]=state
        
    def set_savgol_polyorder(self, value):
        if value == 1 or value%3 == 0:
            config_temp["savgol_polyorder"]=value
        else:
            self.spinBox_savgol_polyorder.setValue(value+1)
            
    def set_savgol_window(self, value):
        config_temp["savgol_window_length"]=value
        
    def set_legend_framealpha(self, value):
        config_temp["legend_framealpha"]=1-value/100
        
    def set_legend_loc(self):
        config_temp["legend_loc"]=self.comboBox_legend_loc.currentText()
        
    def set_hitemp(self, state):
        config_temp["hitemp"]=state
        
    def set_start_index_manually(self, state):
        config_temp["start_index_manually"]=state
        
    def set_coef_unit(self, state):
        config_temp["unit"]=state
        
    def set_full_names(self, state):
        config_temp["full_names"]=state

    def set_font_size_value(self, value):
        config_temp["font-size"]=value

    def set_frequency_value(self, value):
        config_temp["frequency"]=value
        
    def set_stretch_factor_value(self, value):
        config_temp["measure_range"]=value
        
    def set_spectrum_filter(self, state):
        config_temp["spectrum_filter"][0]=state
        
    def set_spectrum_filter_value(self, value):
        config_temp["spectrum_filter"][1]=value
        
    def set_noise_filter(self, state):
        config_temp["noise_filter"][0]=state
        
    def set_noise_filter_value(self, value):
        config_temp["noise_filter"][1]=value
        
    def set_overlap_filter(self, state):
        config_temp["overlap_filter"][0]=state
        
    def set_overlap_filter_value(self, value):
        config_temp["overlap_filter"][1]=value
        
    def set_path_length(self, value):
        config_temp["path_length"]=value
        
    def set_pressure_value(self, value):
        if config_temp["unit_pressure"]:
            config_temp["pressure"][1] = value
            config_temp["pressure"][0] = value*1.01325
        else:
            config_temp["pressure"][1] = value/1.01325
            config_temp["pressure"][0] = value
        
    def set_temperature_value(self, value):
        if config_temp["unit_temperature"]:
            config_temp["temperature"][1] = value
            config_temp["temperature"][0] = value - 273.15
        else:
            config_temp["temperature"][1] = value + 273.15
            config_temp["temperature"][0] = value
        
    def set_mode(self, index):
        config_temp["online"]=index

    def set_wave_min(self, value):
        if config_temp["unit_wave"]:
            config_temp["vmin"][1]=value
            config_temp["vmax"][0]=1e7/value
        else:
            config_temp["vmin"][0]=value
            config_temp["vmax"][1]=1e7/value
##        if not self.unit_blocker:
##            self.doubleSpinBox_wave_min.setMaximum(self.doubleSpinBox_wave_max.value()-0.01)

    def set_wave_max(self, value):
        if config_temp["unit_wave"]:
            config_temp["vmax"][1]=value
            config_temp["vmin"][0]=1e7/value
        else:
            config_temp["vmax"][0]=value
            config_temp["vmin"][1]=1e7/value
##        if not self.unit_blocker:
##            self.doubleSpinBox_wave_max.setMinimum(self.doubleSpinBox_wave_min.value()+0.01)

    def set_broadening(self, index):
        if index:
            config_temp["broadening"]="gamma_self"
        else:
            config_temp["broadening"]="gamma_air"
        
    def set_profile(self, index):
        if index == 0:
            profile="Voigt"
        if index == 1:
            profile = "Gauss"
        if index == 2:
            profile = "Lorentz"
        config_temp["profile"]=profile
                                                 
    def unit_wave(self, index):
        self.unit_blocker=True
        config_temp["unit_wave"]=index
        if index:
            self.label_wave_min.setText("vmin [cm<sup>-1</sup>]")
            self.label_wave_max.setText("vmax [cm<sup>-1</sup>]")
        else:
            self.label_wave_min.setText("lambda min [nm]")
            self.label_wave_max.setText("lambda max [nm]")
        self.doubleSpinBox_wave_min.setValue(config_temp["vmin"][index])
        self.doubleSpinBox_wave_max.setValue(config_temp["vmax"][index])
        self.unit_blocker=False

    def unit_pressure(self, index):
        config_temp["unit_pressure"]=index
        if index:
            self.label_pressure.setText("Pressure [atm]")
        else:
            self.label_pressure.setText("Pressure [bar]")
        self.doubleSpinBox_pressure.setValue(config_temp["pressure"][index])
        
    def unit_temperature(self, index):
        config_temp["unit_temperature"]=index
        self.doubleSpinBox_temperature.setRange(-300, 3000)
        self.doubleSpinBox_temperature.setValue(config_temp["temperature"][index])
        if index:
            self.label_temperature.setText("Temperature [K]")
            self.doubleSpinBox_temperature.setRange(70, 3000)
        else:
            self.label_temperature.setText("Temperature [°C]")
            self.doubleSpinBox_temperature.setRange(70-273.15, 3000-273.15)
            
    def save(self):
        global font_blocker, config
        font_blocker=1
        if config["full_names"] != config_temp["full_names"] or config["font-size"] != config_temp["font-size"] or config_temp["sort_molecules_by_names"] != config["sort_molecules_by_names"]:
            trees=[]
            for obj in gc.get_objects():
                if isinstance(obj, MoleculeTree):
                    trees.append(obj)
            for j in trees:
                root = j.tree.invisibleRootItem()
                child_count = root.childCount()
                for i in range(child_count):
                    item = root.child(i)
                    if config_temp["full_names"]:
                        item.setText(0, "<font size=\""+str(config_temp["font-size"])+"\">"+j.transform_html(item.data(0,32))+" ("+j.molecule_full_names[item.data(0,32)]+")"+"</font>")
                        j.tree.setMaximumWidth(400)
                        j.tree.setMinimumWidth(400)
                    else:
                        item.setText(0, "<font size=\""+str(config_temp["font-size"])+"\">"+j.transform_html(item.data(0,32))+"</font>")
                        j.tree.setMaximumWidth(280)
                        j.tree.setMinimumWidth(280)
                    if config_temp["sort_molecules_by_names"]:
                        item.setSortData(0, j.molecule_full_names[item.data(0,32)])
                    else:
                        item.setSortData(0, item.data(0,32))
                    if config["font-size"] != config_temp["font-size"]:
                        for k in range(item.childCount()):
                            isotope = item.child(k)
                            for i in range(4):
                                isotope.setText(i, isotope.text(i).replace("font size=\""+str(config["font-size"])+"\"", "font size=\""+str(config_temp["font-size"])+"\""))
                j.tree.sortItems(0,0)
        font_blocker=0
        config=config_temp
        if not config["unit_wave"]:
            self.parent.label_wave.setText("Wave length [nm]")
            self.parent.tableWidget_lineselection.horizontalHeaderItem(3).setText("Wave length")
        else:
            self.parent.label_wave.setText("Wave number [cm<sup>-1</sup>]")
            self.parent.tableWidget_lineselection.horizontalHeaderItem(3).setText("Wave number")
        if not config["unit_pressure"]:
            self.parent.label_pressure.setText("Pressure [bar]")
        else:
            self.parent.label_pressure.setText("Pressure [atm]")
        if config["broadening"] == "gamma_self":
            broadening_text="self"
        elif config["broadening"] == "gamma_air":
            broadening_text="air"
        self.parent.broadening_label.setText("Broadening: "+broadening_text+" ["+config["profile"]+"]")
        self.parent.path_length_label.setText("Path length: "+str(config["path_length"])+" cm")

        self.parent.unit_blocker=True
        self.parent.doubleSpinBox_vmin.setValue(config["vmin"][config["unit_wave"]])
        self.parent.doubleSpinBox_vmax.setValue(config["vmax"][config["unit_wave"]])
        self.parent.unit_blocker=False
        
        self.parent.doubleSpinBox_temperature.setRange(70-273.15, 3000)
        self.parent.doubleSpinBox_temperature.setValue(config["temperature"][config["unit_temperature"]])
        if not config["unit_temperature"]:
            self.parent.label_temperature.setText("Temperature [°C]")
            self.parent.doubleSpinBox_temperature.setRange(70-273.15, 3000-273.15)
        else:
            self.parent.label_temperature.setText("Temperature [K]")
            self.parent.doubleSpinBox_temperature.setRange(70, 3000)
        self.parent.doubleSpinBox_pressure.setValue(config["pressure"][config["unit_pressure"]])

        json.dump(config, open("txt/config.txt", "w"))
        self.close()
        
if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    MainWindow = MainWindow()
    sys.exit(app.exec_())
