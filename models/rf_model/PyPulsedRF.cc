#include "PulsedRF.hh"

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <string>
#include <iostream>

struct PyPulsedRF {
    PyObject_HEAD
    PulsedRF structure;
};

std::vector<double> ConvertList(PyObject* pylist) {
    if (!PyList_Check(pylist)) {
        throw std::string("Could not parse argument as a list");
    }
    std::vector<double> list(PyList_Size(pylist));
    for (Py_ssize_t i = 0; i < PyList_Size(pylist); ++i) {
        PyObject* pyfloat = PyList_GetItem(pylist, i);
        double cfloat = PyFloat_AsDouble(pyfloat);
        if (PyErr_Occurred()) {
            throw std::string("Could not parse item as a float");
        }
        list[i] = cfloat;
    }
    return list;
}

static PyObject* SetTanhFieldModel(PyObject* self, PyObject* args) {
    double b0, centreLength, endLength = 0;
    if (!PyArg_ParseTuple(args, "ddd", &b0, &centreLength, &endLength)) {
        PyErr_SetString(PyExc_ValueError, "Failed to parse arguments");
        return NULL;
    }
    PyPulsedRF* pyrf = (PyPulsedRF*)self;
    PulsedRF& sol = pyrf->structure;
    TanhFieldModel* tfm = new TanhFieldModel();
    tfm->_b0 = b0;
    tfm->_x0 = centreLength/2;
    tfm->_lambda = endLength;
    sol.fieldModel = std::unique_ptr<OnAxisFieldModel>(tfm);
    tfm->Initialise(20);
    Py_RETURN_NONE;
}

static PyObject* SetMaxDerivative(PyObject* self, PyObject* args) {
    int maxDerivative = 0;
    if (!PyArg_ParseTuple(args, "i", &maxDerivative)) {
        PyErr_SetString(PyExc_ValueError, "Failed to parse arguments");
        return NULL;
    }
    if (maxDerivative < 0) {
        PyErr_SetString(PyExc_ValueError, "Maximum derivative should be >= 0");
        return NULL;
    }
    PyPulsedRF* pyrf = (PyPulsedRF*)self;
    pyrf->structure.maxDerivative = maxDerivative;
    Py_RETURN_NONE;
}

static PyObject* GetMaxDerivative(PyObject* self, PyObject* args) {
    PyPulsedRF* pyrf = (PyPulsedRF*)self;
    unsigned int maxDerivative = pyrf->structure.maxDerivative;
    PyObject* pyDeriv = PyLong_FromSize_t(maxDerivative);
    Py_INCREF(pyDeriv);
    return pyDeriv;
}

static PyObject* SetLength(PyObject* self, PyObject* args) {
    double length = 0;
    if (!PyArg_ParseTuple(args, "d", &length)) {
        PyErr_SetString(PyExc_ValueError, "Failed to parse arguments");
        return NULL;
    }
    PyPulsedRF* pyrf = (PyPulsedRF*)self;
    pyrf->structure.length = length;
    Py_RETURN_NONE;
}

static PyObject* GetLength(PyObject* self, PyObject* args) {
    PyPulsedRF* pyrf = (PyPulsedRF*)self;
    double length = pyrf->structure.length;
    PyObject* pyLength = PyFloat_FromDouble(length);
    Py_INCREF(pyLength);
    return pyLength;
}

static PyObject* SetMaxR(PyObject* self, PyObject* args) {
    double maxR = 0;
    if (!PyArg_ParseTuple(args, "d", &maxR)) {
        PyErr_SetString(PyExc_ValueError, "Failed to parse arguments");
        return NULL;
    }
    PyPulsedRF* pyrf = (PyPulsedRF*)self;
    pyrf->structure.maxR = maxR;
    Py_RETURN_NONE;
}

static PyObject* GetMaxR(PyObject* self, PyObject* args) {
    PyPulsedRF* pyrf = (PyPulsedRF*)self;
    double maxR = pyrf->structure.maxR;
    PyObject* pyMaxR = PyFloat_FromDouble(maxR);
    Py_INCREF(pyMaxR);
    return pyMaxR;
}

static PyObject* SetV0(PyObject* self, PyObject* args) {
    double v_0 = 0;
    if (!PyArg_ParseTuple(args, "d", &v_0)) {
        PyErr_SetString(PyExc_ValueError, "Failed to parse arguments");
        return NULL;
    }
    PyPulsedRF* pyrf = (PyPulsedRF*)self;
    pyrf->structure.v_0 = v_0*pyrf->structure.c_light;
    Py_RETURN_NONE;
}

static PyObject* GetV0(PyObject* self, PyObject* args) {
    PyPulsedRF* pyrf = (PyPulsedRF*)self;
    double v_0 = pyrf->structure.v_0/pyrf->structure.c_light;
    PyObject* pyV0 = PyFloat_FromDouble(v_0);
    Py_INCREF(pyV0);
    return pyV0;
}

static PyObject* SetZOffset(PyObject* self, PyObject* args) {
    double zOffset = 0;
    if (!PyArg_ParseTuple(args, "d", &zOffset)) {
        PyErr_SetString(PyExc_ValueError, "Failed to parse arguments");
        return NULL;
    }
    PyPulsedRF* pyrf = (PyPulsedRF*)self;
    pyrf->structure.zOffset = zOffset;
    Py_RETURN_NONE;
}

static PyObject* GetZOffset(PyObject* self, PyObject* args) {
    PyPulsedRF* pyrf = (PyPulsedRF*)self;
    double zOffset = pyrf->structure.zOffset;
    PyObject* pyZOffset = PyFloat_FromDouble(zOffset);
    Py_INCREF(pyZOffset);
    return pyZOffset;
}


static PyObject* GetFieldValue(PyObject* self, PyObject* args) {
    double x=0, y=0, z=0, time=0;
    if (!PyArg_ParseTuple(args, "dddd", &x, &y, &z, &time)) {
        PyErr_SetString(PyExc_ValueError, "Failed to parse arguments");
        return NULL;
    }
    PyPulsedRF* pyrf = (PyPulsedRF*)self;
    PulsedRF& rf = pyrf->structure;
    std::vector<double> pos = {x, y, z};
    std::vector<double> bfield = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    try {
        rf.GetFieldValue(pos, time, bfield);
    } catch (std::string& str) {
        PyErr_SetString(PyExc_RuntimeError, str.c_str());
        return nullptr;
    }
    PyObject* pyfield = Py_BuildValue("(dddddd)", bfield[0], bfield[1], bfield[2], bfield[3], bfield[4], bfield[5]);
    Py_INCREF(pyfield);
    return pyfield;
}

static PyObject* GetOnAxisField(PyObject* self, PyObject* args) {
    double z = 0;
    int derivativeOrder = 0;
    if (!PyArg_ParseTuple(args, "di", &z, &derivativeOrder)) {
        PyErr_SetString(PyExc_ValueError, "Failed to parse arguments");
        return NULL;
    }
    PyPulsedRF* pyrf = (PyPulsedRF*)self;
    PulsedRF& rf = pyrf->structure;
    if (rf.fieldModel.get() == nullptr) {
        PyErr_SetString(PyExc_RuntimeError, "On axis field was not set");
        return nullptr;
    }
    double efield = 0;
    try {
        rf.fieldModel->GetField(z, derivativeOrder, efield);
    } catch (std::string& str) {
        PyErr_SetString(PyExc_RuntimeError, str.c_str());
        return nullptr;
    }
    PyObject* pyefield = Py_BuildValue("d", efield);
    Py_INCREF(pyefield);
    return pyefield;

}


static PyMethodDef pulsed_rf_module_methods[] = {
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static PyMethodDef pulsed_rf_methods[] = {
    {"get_field_value",  GetFieldValue, METH_VARARGS,      "Get the field value on the axis, or its derivative."},
    {"get_on_axis_field",  GetOnAxisField, METH_VARARGS, "Get the field value."},
    {"set_tanh_field_model", SetTanhFieldModel, METH_VARARGS, "Set the field model - using a tanh model. 3 arguments are e0, centre length, end length"},
    {"set_max_derivative", SetMaxDerivative, METH_VARARGS, "Set the maximum derivative in the expansion"},
    {"get_max_derivative", GetMaxDerivative, METH_VARARGS, "Get the maximum derivative in the expansion"},
    {"set_length", SetLength, METH_VARARGS, "Set the length of the bounding box"},
    {"get_length", GetLength, METH_VARARGS, "Get the length of the bounding box"},
    {"set_max_r", SetMaxR, METH_VARARGS, "Set the maximum radius of the bounding box"},
    {"get_max_r", GetMaxR, METH_VARARGS, "Get the maximum radius of the bounding box"},
    {"set_v0", SetV0, METH_VARARGS, "Set the pulse velocity as a fraction of speed of light"},
    {"get_v0", GetV0, METH_VARARGS, "Get the pulse velocity as a fraction of speed of light"},
    {"set_z_offset", SetZOffset, METH_VARARGS, "Set the z position of the pulse at t = 0"},
    {"get_z_offset", GetZOffset, METH_VARARGS, "Get the z position of the pulse at t = 0"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef pulsed_rf_module = {
    PyModuleDef_HEAD_INIT,
    "pulsed_rf",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    pulsed_rf_module_methods
};

static PyTypeObject PyPulsedRFType = {
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "models.rf_model.pulsed_rf.PulsedRF",
    .tp_basicsize = sizeof(PyPulsedRF),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_doc = PyDoc_STR("Pulsed RF"),
    .tp_methods = pulsed_rf_methods,
    .tp_new = PyType_GenericNew
};

PyMODINIT_FUNC
PyInit_pulsed_rf(void)
{
     if (PyType_Ready(&PyPulsedRFType) < 0)
        return NULL;

    PyObject* rf_module = PyModule_Create(&pulsed_rf_module);
    if (rf_module == nullptr) {
        return NULL;
    }

    Py_INCREF(&PyPulsedRFType);
    if (PyModule_AddObject(rf_module, "PulsedRF", (PyObject *) &PyPulsedRFType) < 0) {
        Py_DECREF(&PyPulsedRFType);
        Py_DECREF(rf_module);
        return NULL;
    }
    return rf_module;
}
