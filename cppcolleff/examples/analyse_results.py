#!/usr/bin/env python-sirius

import sys
import numpy as np
import os
from PyQt5.QtWidgets import QApplication, QPushButton, QWidget, QLabel
from PyQt5.QtWidgets import QVBoxLayout, QHBoxLayout, QComboBox
import pyqtgraph as pg
import cppcolleff as coll


class MyWindow(QWidget):

    def __init__(self):
        super().__init__()
        self.results = {}
        self._time = []
        self.planes = ('xx', 'xl', 'de', 'ss')
        self.dwnsmpl = 10
        self.setupUi()
        self._changedir(self.cb_dirs.currentText())

    def setupUi(self):
        self.resize(1200, 900)
        gr = pg.PlotWidget(parent=self, background=[255, 255, 255])
        dtait = pg.PlotDataItem(pen=pg.mkPen('r', width=4))
        dtait2 = pg.PlotDataItem(pen=pg.mkPen('b', width=4))
        dtait.setVisible(True)
        gr.plotItem.addItem(dtait)
        gr.plotItem.addItem(dtait2)
        pl = gr.getPlotItem()
        vl = QVBoxLayout(self)
        vl.addWidget(gr)
        self.graph = gr
        self._plot_rms = dtait
        self._plot_ave = dtait2

        hl = QHBoxLayout()
        vl.addItem(hl)

        self.dirs = sorted([f for f in os.listdir() if os.path.isdir(f)])
        self.cb_dirs = QComboBox(self)
        self.cb_dirs.addItems(self.dirs)
        self.cb_dirs.currentTextChanged.connect(self._changedir)
        hl.addWidget(self.cb_dirs)

        self.cb_planes = QComboBox(self)
        self.cb_planes.addItems(self.planes + ('Wlkick', ))
        self.cb_planes.currentTextChanged.connect(self._changePlane)
        hl.addWidget(self.cb_planes)

    def _changePlane(self, text):
        self._changedir(self.cb_dirs.currentText())

    def _changedir(self, dir_):
        fil = [f for f in os.listdir(dir_) if f == "results.txt"]
        if not fil:
            self.cb_dirs.setCurrentIndex(0)
            return
        fil = dir_ + os.path.sep + fil[0]
        res = self.results.get(dir_, None)
        if res is not None:
            self._update_graphic()
            return
        if not os.path.isfile(fil):
            return
        res = coll.Results_t(1)
        res.from_file(fil)
        dt = 1.726e-6*res.get_calc_every()
        self.results[dir_] = res
        self._time = np.arange(res.ave.size(), step=self.dwnsmpl) * (dt)
        self._update_graphic()

    def _update_graphic(self):
        dir_ = self.cb_dirs.currentText()
        plane = self.cb_planes.currentText()
        res = self.results.get(dir_, None)
        if res is None:
            return
        if plane in self.planes:
            y1 = [getattr(ave, plane) for i, ave in enumerate(res.ave)
                  if not (i % self.dwnsmpl)]
            y2 = [getattr(rms, plane) for i, rms in enumerate(res.std)
                  if not (i % self.dwnsmpl)]
        else:
            y1 = [x for i, x in enumerate(getattr(res, plane))
                  if not (i % self.dwnsmpl)]
            y2 = y1
        self._plot_ave.setData(self._time, y1)
        self._plot_rms.setData(self._time, y2)


if __name__ == "__main__":
    app = QApplication([])
    win = MyWindow()
    win.show()
    sys.exit(app.exec_())
