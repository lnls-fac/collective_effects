#!/usr/bin/env python-sirius

import sys
import numpy as np
import os
from PyQt5.QtWidgets import QApplication, QPushButton, QWidget, QLabel
from PyQt5.QtWidgets import QVBoxLayout, QHBoxLayout, QComboBox, QSpacerItem
from PyQt5.QtCore import QTimer
import pyqtgraph as pg


class MyWindow(QWidget):

    def __init__(self):
        super().__init__()
        self.timer = QTimer()
        self.timer.timeout.connect(self._update_graphic)
        self.timer.setInterval(200)
        self.distrs = {"ss": {}, "de": {}}
        self.turn_numbers = {}
        self.plane = "ss"
        self._ind = 0
        self.pos = None
        self.setupUi()

    def setupUi(self):
        self.resize(1000, 900)
        gr = pg.PlotWidget(parent=self, background=[255, 255, 255])
        dtait = pg.PlotDataItem(pen=pg.mkPen('r', width=4))
        gr.plotItem.addItem(dtait)
        pl = gr.getPlotItem()
        vl = QVBoxLayout(self)
        self.setLayout(vl)
        vl.addWidget(gr)
        self.graph = gr
        self._plot = dtait

        hl = QHBoxLayout()
        vl.addItem(hl)

        hl.addWidget(QLabel("Turn Number", self))
        self._lb_turn = QLabel("None", self)
        hl.addWidget(self._lb_turn)

        # hl.addWidget(QLabel("Current", self))
        # self._lb_curr = QLabel("None", self)
        # hl.addWidget(self._lb_curr)

        pb1 = QPushButton("Start", self)
        pb1.clicked.connect(self.timer.start)
        hl.addWidget(pb1)
        pb2 = QPushButton("Stop", self)
        pb2.clicked.connect(self.timer.stop)
        hl.addWidget(pb2)

        self.dirs = sorted([f for f in os.listdir() if os.path.isdir(f)])
        self.cb_dirs = QComboBox(self)
        self.cb_dirs.addItems(self.dirs)
        self.cb_dirs.currentTextChanged.connect(self._changedir)
        hl.addWidget(self.cb_dirs)

        self.cb_planes = QComboBox(self)
        self.cb_planes.addItems(['ss', 'de'])
        self.cb_planes.currentTextChanged.connect(self._changePlane)
        hl.addWidget(self.cb_planes)

    def _changePlane(self, text):
        self.plane = text
        self._changedir(self.cb_dirs.currentText())

    def _changedir(self, text):
        files = sorted([f for f in os.listdir(text)
                        if f.endswith(self.plane + ".txt")])
        dists = self.distrs[self.plane].get(text, [])
        if not dists:
            turn_numbers = []
            for i, f in enumerate(files):
                fi_ = text + os.path.sep + f
                if i == 0:
                    with open(fi_, 'r') as fp:
                        min_ = float(fp.readline().split()[2])
                        max_ = float(fp.readline().split()[2])
                        bin_ = int(fp.readline().split()[2])
                    self.pos = np.linspace(min_, max_, bin_)
                dists.append(np.loadtxt(fi_, skiprows=4))
                turn_numbers.append(f[4:11])
            self.distrs[self.plane][text] = dists
            self.turn_numbers[text] = turn_numbers

    def _update_graphic(self):
        text = self.cb_dirs.currentText()
        dists = self.distrs[self.plane].get(text, [])
        if dists:
            self._plot.setData(self.pos, dists[self._ind])
            self._lb_turn.setText(self.turn_numbers[text][self._ind])
            self._ind = (self._ind + 1) % len(dists)
            if not self._ind:
                self.timer.stop()


if __name__ == "__main__":
    app = QApplication([])
    win = MyWindow()
    win.show()
    sys.exit(app.exec_())
