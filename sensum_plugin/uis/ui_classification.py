# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'uis/ui_classification.ui'
#
# Created: Wed Jun 11 12:14:53 2014
#      by: PyQt4 UI code generator 4.10.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_Classification(object):
    def setupUi(self, Classification):
        Classification.setObjectName(_fromUtf8("Classification"))
        Classification.resize(319, 415)
        Classification.setMinimumSize(QtCore.QSize(319, 415))
        Classification.setMaximumSize(QtCore.QSize(319, 415))
        self.widget = QtGui.QWidget(Classification)
        self.widget.setGeometry(QtCore.QRect(30, 50, 261, 151))
        self.widget.setObjectName(_fromUtf8("widget"))
        self.gridLayout = QtGui.QGridLayout(self.widget)
        self.gridLayout.setMargin(0)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.label_input = QtGui.QLabel(self.widget)
        self.label_input.setObjectName(_fromUtf8("label_input"))
        self.gridLayout.addWidget(self.label_input, 0, 0, 1, 1)
        self.pushButton_input = QtGui.QPushButton(self.widget)
        self.pushButton_input.setObjectName(_fromUtf8("pushButton_input"))
        self.gridLayout.addWidget(self.pushButton_input, 1, 1, 1, 1)
        self.pushButton_output = QtGui.QPushButton(self.widget)
        self.pushButton_output.setObjectName(_fromUtf8("pushButton_output"))
        self.gridLayout.addWidget(self.pushButton_output, 3, 1, 1, 1)
        self.lineEdit_input = QtGui.QLineEdit(self.widget)
        self.lineEdit_input.setObjectName(_fromUtf8("lineEdit_input"))
        self.gridLayout.addWidget(self.lineEdit_input, 1, 0, 1, 1)
        self.label_output = QtGui.QLabel(self.widget)
        self.label_output.setObjectName(_fromUtf8("label_output"))
        self.gridLayout.addWidget(self.label_output, 2, 0, 1, 1)
        self.lineEdit_output = QtGui.QLineEdit(self.widget)
        self.lineEdit_output.setObjectName(_fromUtf8("lineEdit_output"))
        self.gridLayout.addWidget(self.lineEdit_output, 3, 0, 1, 1)
        self.comboBox_supervised = QtGui.QComboBox(self.widget)
        self.comboBox_supervised.setObjectName(_fromUtf8("comboBox_supervised"))
        self.comboBox_supervised.addItem(_fromUtf8(""))
        self.comboBox_supervised.addItem(_fromUtf8(""))
        self.gridLayout.addWidget(self.comboBox_supervised, 4, 0, 1, 1)
        self.groupBox = QtGui.QGroupBox(Classification)
        self.groupBox.setGeometry(QtCore.QRect(10, 210, 301, 151))
        self.groupBox.setObjectName(_fromUtf8("groupBox"))
        self.frame_supervised = QtGui.QFrame(self.groupBox)
        self.frame_supervised.setGeometry(QtCore.QRect(11, 26, 279, 114))
        self.frame_supervised.setFrameShape(QtGui.QFrame.StyledPanel)
        self.frame_supervised.setFrameShadow(QtGui.QFrame.Raised)
        self.frame_supervised.setObjectName(_fromUtf8("frame_supervised"))
        self.formLayout_2 = QtGui.QFormLayout(self.frame_supervised)
        self.formLayout_2.setFieldGrowthPolicy(QtGui.QFormLayout.AllNonFixedFieldsGrow)
        self.formLayout_2.setObjectName(_fromUtf8("formLayout_2"))
        self.label_supervised_type = QtGui.QLabel(self.frame_supervised)
        self.label_supervised_type.setObjectName(_fromUtf8("label_supervised_type"))
        self.formLayout_2.setWidget(0, QtGui.QFormLayout.LabelRole, self.label_supervised_type)
        self.label_training_type = QtGui.QLabel(self.frame_supervised)
        self.label_training_type.setObjectName(_fromUtf8("label_training_type"))
        self.formLayout_2.setWidget(0, QtGui.QFormLayout.FieldRole, self.label_training_type)
        self.comboBox_supervised_type = QtGui.QComboBox(self.frame_supervised)
        self.comboBox_supervised_type.setObjectName(_fromUtf8("comboBox_supervised_type"))
        self.comboBox_supervised_type.addItem(_fromUtf8(""))
        self.comboBox_supervised_type.addItem(_fromUtf8(""))
        self.comboBox_supervised_type.addItem(_fromUtf8(""))
        self.comboBox_supervised_type.addItem(_fromUtf8(""))
        self.comboBox_supervised_type.addItem(_fromUtf8(""))
        self.comboBox_supervised_type.addItem(_fromUtf8(""))
        self.comboBox_supervised_type.addItem(_fromUtf8(""))
        self.formLayout_2.setWidget(1, QtGui.QFormLayout.LabelRole, self.comboBox_supervised_type)
        self.lineEdit_training_field = QtGui.QLineEdit(self.frame_supervised)
        self.lineEdit_training_field.setObjectName(_fromUtf8("lineEdit_training_field"))
        self.formLayout_2.setWidget(1, QtGui.QFormLayout.FieldRole, self.lineEdit_training_field)
        self.label_training = QtGui.QLabel(self.frame_supervised)
        self.label_training.setObjectName(_fromUtf8("label_training"))
        self.formLayout_2.setWidget(2, QtGui.QFormLayout.LabelRole, self.label_training)
        self.lineEdit_training = QtGui.QLineEdit(self.frame_supervised)
        self.lineEdit_training.setObjectName(_fromUtf8("lineEdit_training"))
        self.formLayout_2.setWidget(3, QtGui.QFormLayout.LabelRole, self.lineEdit_training)
        self.pushButton_training = QtGui.QPushButton(self.frame_supervised)
        self.pushButton_training.setObjectName(_fromUtf8("pushButton_training"))
        self.formLayout_2.setWidget(3, QtGui.QFormLayout.FieldRole, self.pushButton_training)
        self.frame_unsupervised = QtGui.QFrame(self.groupBox)
        self.frame_unsupervised.setGeometry(QtCore.QRect(20, 50, 266, 64))
        self.frame_unsupervised.setFrameShape(QtGui.QFrame.StyledPanel)
        self.frame_unsupervised.setFrameShadow(QtGui.QFrame.Raised)
        self.frame_unsupervised.setObjectName(_fromUtf8("frame_unsupervised"))
        self.formLayout = QtGui.QFormLayout(self.frame_unsupervised)
        self.formLayout.setObjectName(_fromUtf8("formLayout"))
        self.label_nclasses = QtGui.QLabel(self.frame_unsupervised)
        self.label_nclasses.setObjectName(_fromUtf8("label_nclasses"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.LabelRole, self.label_nclasses)
        self.label_niteration = QtGui.QLabel(self.frame_unsupervised)
        self.label_niteration.setObjectName(_fromUtf8("label_niteration"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.FieldRole, self.label_niteration)
        self.spinBox_nclasses = QtGui.QSpinBox(self.frame_unsupervised)
        self.spinBox_nclasses.setObjectName(_fromUtf8("spinBox_nclasses"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.LabelRole, self.spinBox_nclasses)
        self.spinBox_niteration = QtGui.QSpinBox(self.frame_unsupervised)
        self.spinBox_niteration.setObjectName(_fromUtf8("spinBox_niteration"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.FieldRole, self.spinBox_niteration)
        self.buttonBox = QtGui.QDialogButtonBox(Classification)
        self.buttonBox.setGeometry(QtCore.QRect(80, 370, 160, 27))
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.label_title = QtGui.QLabel(Classification)
        self.label_title.setGeometry(QtCore.QRect(70, 20, 201, 21))
        self.label_title.setObjectName(_fromUtf8("label_title"))

        self.retranslateUi(Classification)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("accepted()")), Classification.accept)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("rejected()")), Classification.reject)
        QtCore.QMetaObject.connectSlotsByName(Classification)

    def retranslateUi(self, Classification):
        Classification.setWindowTitle(_translate("Classification", "Sensum", None))
        self.label_input.setText(_translate("Classification", "Input File:", None))
        self.pushButton_input.setText(_translate("Classification", "...", None))
        self.pushButton_output.setText(_translate("Classification", "...", None))
        self.label_output.setText(_translate("Classification", "Output File:", None))
        self.comboBox_supervised.setItemText(0, _translate("Classification", "Supervised", None))
        self.comboBox_supervised.setItemText(1, _translate("Classification", "Unsupervised", None))
        self.groupBox.setTitle(_translate("Classification", "Options", None))
        self.label_supervised_type.setText(_translate("Classification", "Supervised Type", None))
        self.label_training_type.setText(_translate("Classification", "Discriminant Field", None))
        self.comboBox_supervised_type.setItemText(0, _translate("Classification", "libsvm", None))
        self.comboBox_supervised_type.setItemText(1, _translate("Classification", "svm", None))
        self.comboBox_supervised_type.setItemText(2, _translate("Classification", "dt", None))
        self.comboBox_supervised_type.setItemText(3, _translate("Classification", "gbt", None))
        self.comboBox_supervised_type.setItemText(4, _translate("Classification", "bayes", None))
        self.comboBox_supervised_type.setItemText(5, _translate("Classification", "rf", None))
        self.comboBox_supervised_type.setItemText(6, _translate("Classification", "knn", None))
        self.label_training.setText(_translate("Classification", "Input Training File:", None))
        self.pushButton_training.setText(_translate("Classification", "...", None))
        self.label_nclasses.setText(_translate("Classification", "Number of classes", None))
        self.label_niteration.setText(_translate("Classification", "Number of Iterations", None))
        self.label_title.setText(_translate("Classification", "<html><head/><body><p><span style=\" font-size:16pt;\">CLASSIFICATION</span></p></body></html>", None))

