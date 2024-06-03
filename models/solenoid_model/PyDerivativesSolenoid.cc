#include "DerivativesSolenoid.hh"

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <string>

struct PyDerivativesSolenoid {
    PyObject_HEAD
    DerivativesSolenoid solenoid;
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

static PyObject* SetFourierFieldModel(PyObject* self, PyObject* args) {
    double length = 0;    
    PyObject* pyHarmonicList;
    if (!PyArg_ParseTuple(args, "dO", &length, &pyHarmonicList)) {
        PyErr_SetString(PyExc_ValueError, "Failed to parse arguments");
        return NULL;
    }
    std::vector<double> harmonicList;
    try {
        harmonicList = ConvertList(pyHarmonicList);
    } catch (std::string str) {
        PyErr_SetString(PyExc_ValueError, str.c_str());
        return NULL;
    }
    PyDerivativesSolenoid* pysol = (PyDerivativesSolenoid*)self;
    DerivativesSolenoid& sol = pysol->solenoid;
    FourierFieldModel* ffm = new FourierFieldModel();
    ffm->cellLength = length;
    ffm->harmonicList = harmonicList;
    sol.fieldModel = std::unique_ptr<OnAxisFieldModel>(ffm);
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
    PyDerivativesSolenoid* pysol = (PyDerivativesSolenoid*)self;
    pysol->solenoid.maxDerivative = maxDerivative;
    Py_RETURN_NONE;    
}

static PyObject* GetMaxDerivative(PyObject* self, PyObject* args) {
    PyDerivativesSolenoid* pysol = (PyDerivativesSolenoid*)self;
    unsigned int maxDerivative = pysol->solenoid.maxDerivative;
    PyObject* pyDeriv = PyLong_FromSize_t(maxDerivative);
    Py_INCREF(pyDeriv);
    return pyDeriv;
}

static PyObject* SetMaxZ(PyObject* self, PyObject* args) {
    double maxZ = 0;
    if (!PyArg_ParseTuple(args, "d", &maxZ)) {
        PyErr_SetString(PyExc_ValueError, "Failed to parse arguments");
        return NULL;
    }
    PyDerivativesSolenoid* pysol = (PyDerivativesSolenoid*)self;
    pysol->solenoid.maxZ = maxZ;
    Py_RETURN_NONE;
}

static PyObject* GetMaxZ(PyObject* self, PyObject* args) {
    PyDerivativesSolenoid* pysol = (PyDerivativesSolenoid*)self;
    unsigned int maxZ = pysol->solenoid.maxZ;
    PyObject* pyMaxZ = PyFloat_FromDouble(maxZ);
    Py_INCREF(pyMaxZ);
    return pyMaxZ;
}

static PyObject* SetMinZ(PyObject* self, PyObject* args) {
    double minZ = 0;
    if (!PyArg_ParseTuple(args, "d", &minZ)) {
        PyErr_SetString(PyExc_ValueError, "Failed to parse arguments");
        return NULL;
    }
    PyDerivativesSolenoid* pysol = (PyDerivativesSolenoid*)self;
    pysol->solenoid.minZ = minZ;
    Py_RETURN_NONE;
}

static PyObject* GetMinZ(PyObject* self, PyObject* args) {
    PyDerivativesSolenoid* pysol = (PyDerivativesSolenoid*)self;
    unsigned int minZ = pysol->solenoid.minZ;
    PyObject* pyMinZ = PyFloat_FromDouble(minZ);
    Py_INCREF(pyMinZ);
    return pyMinZ;
}

static PyObject* SetMaxR(PyObject* self, PyObject* args) {
    double maxR = 0;
    if (!PyArg_ParseTuple(args, "d", &maxR)) {
        PyErr_SetString(PyExc_ValueError, "Failed to parse arguments");
        return NULL;
    }
    PyDerivativesSolenoid* pysol = (PyDerivativesSolenoid*)self;
    pysol->solenoid.maxR = maxR;
    Py_RETURN_NONE;
}

static PyObject* GetMaxR(PyObject* self, PyObject* args) {
    PyDerivativesSolenoid* pysol = (PyDerivativesSolenoid*)self;
    unsigned int maxR = pysol->solenoid.maxR;
    PyObject* pyMaxR = PyFloat_FromDouble(maxR);
    Py_INCREF(pyMaxR);
    return pyMaxR;
}

static PyObject* GetFieldValue(PyObject* self, PyObject* args) {
    double x=0, y=0, z=0, time=0;
    if (!PyArg_ParseTuple(args, "dddd", &x, &y, &z, &time)) {
        PyErr_SetString(PyExc_ValueError, "Failed to parse arguments");
        return NULL;
    }
    PyDerivativesSolenoid* pysol = (PyDerivativesSolenoid*)self;
    DerivativesSolenoid& sol = pysol->solenoid;
    std::vector<double> pos = {x, y, z};
    std::vector<double> bfield = {0.0, 0.0, 0.0};
    try {
        sol.GetFieldValue(pos, time, bfield);
    } catch (std::string& str) {
        PyErr_SetString(PyExc_RuntimeError, str.c_str());
        return nullptr;
    }
    PyObject* pybfield = Py_BuildValue("(ddd)", bfield[0], bfield[1], bfield[2]);
    Py_INCREF(pybfield);
    return pybfield;
}

static PyObject* GetOnAxisField(PyObject* self, PyObject* args) {
    double z = 0;
    int derivativeOrder = 0;
    if (!PyArg_ParseTuple(args, "d|i", &z, &derivativeOrder)) {
        PyErr_SetString(PyExc_ValueError, "Failed to parse arguments");
        return NULL;
    }
    PyDerivativesSolenoid* pysol = (PyDerivativesSolenoid*)self;
    DerivativesSolenoid& sol = pysol->solenoid;
    if (sol.fieldModel.get() == nullptr) {
        PyErr_SetString(PyExc_RuntimeError, "On axis field was not set");
        return nullptr;
    }
    double bfield = 0;
    try {
        sol.fieldModel->GetField(z, derivativeOrder, bfield);
    } catch (std::string& str) {
        PyErr_SetString(PyExc_RuntimeError, str.c_str());
        return nullptr;
    }
    PyObject* pybfield = Py_BuildValue("d", bfield);
    Py_INCREF(pybfield);
    return pybfield;

}

static PyMethodDef derivatives_solenoid_module_methods[] = {
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static PyMethodDef derivatives_solenoid_methods[] = {
    {"get_field_value",  GetFieldValue, METH_VARARGS,
     "Get the field value."},
    {"get_on_axis_field",  GetOnAxisField, METH_VARARGS,
     "Get the field value on the axis, or its derivative."},
    {"set_fourier_field_model", SetFourierFieldModel, METH_VARARGS, "Set the field model"},
    {"set_max_derivative", SetMaxDerivative, METH_VARARGS, "Set the maximum derivative in the expansion"},
    {"get_max_derivative", GetMaxDerivative, METH_VARARGS, "Get the maximum derivative in the expansion"},
    {"set_max_z", SetMaxZ, METH_VARARGS, "Set the maximum z of the bounding box"},
    {"get_max_z", GetMaxZ, METH_VARARGS, "Get the maximum z of the bounding box"},
    {"set_min_z", SetMinZ, METH_VARARGS, "Set the minimum z of the bounding box"},
    {"get_min_z", GetMinZ, METH_VARARGS, "Get the minimum z of the bounding box"},
    {"set_max_r", SetMaxR, METH_VARARGS, "Set the maximum radius of the bounding box"},
    {"get_max_r", GetMaxR, METH_VARARGS, "Get the maximum radius of the bounding box"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef derivatives_solenoid_module = {
    PyModuleDef_HEAD_INIT,
    "derivatives_solenoid",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    derivatives_solenoid_module_methods
};

static PyTypeObject PyDerivativesSolenoidType = {
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "derivatives_solenoid.DerivativesSolenoid",
    .tp_basicsize = sizeof(PyDerivativesSolenoid),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_doc = PyDoc_STR("Derivatives Solenoid"),
    .tp_methods = derivatives_solenoid_methods,
    .tp_new = PyType_GenericNew
};

PyMODINIT_FUNC
PyInit_derivatives_solenoid(void)
{
     if (PyType_Ready(&PyDerivativesSolenoidType) < 0)
        return NULL;

    PyObject* deriv_module = PyModule_Create(&derivatives_solenoid_module);
    if (deriv_module == nullptr) {
        return NULL;
    }

    Py_INCREF(&PyDerivativesSolenoidType);
    if (PyModule_AddObject(deriv_module, "DerivativesSolenoid", (PyObject *) &PyDerivativesSolenoidType) < 0) {
        Py_DECREF(&PyDerivativesSolenoidType);
        Py_DECREF(deriv_module);
        return NULL;
    }
    return deriv_module;
}