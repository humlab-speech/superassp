/*
 * _snack.c - Python C extension for the Snack Sound Toolkit
 *
 * Provides a direct Python API over the Snack C library.  No Tcl or Tk
 * APIs are exposed to Python callers; an embedded Tcl interpreter is
 * used internally for resource management and complex operations.
 *
 * The Sound_Init() entry point (from generic/sound.c) initialises Snack
 * without Tk, so this extension has no GUI dependency.
 *
 * Copyright (C) 1997-2005 Kare Sjolander
 * Python bindings copyright (C) 2024 contributors
 *
 * This file is part of the Snack Sound Toolkit and follows the
 * repository license terms in LICENSE.txt.
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "tcl.h"
#include "snack.h"

/* Single embedded Tcl interpreter, initialised once at module load. */
static Tcl_Interp *g_interp = NULL;
static int g_sound_counter  = 0;
static int g_filter_counter = 0;

/* Forward declarations for types referenced by method implementations. */
static PyTypeObject PySoundType;
static PyTypeObject PyFilterType;

/* -----------------------------------------------------------------------
 * Helpers
 * ----------------------------------------------------------------------- */

/* Convert the current Tcl result into a Python exception and return NULL. */
static PyObject *
raise_tcl_error(const char *prefix)
{
    const char *msg = g_interp ? Tcl_GetStringResult(g_interp) : "(no interp)";
    if (prefix)
        PyErr_Format(PyExc_RuntimeError, "%s: %s", prefix, msg);
    else
        PyErr_SetString(PyExc_RuntimeError, msg);
    return NULL;
}

/* Convert a Tcl list object into a Python list of floats. */
static PyObject *
tcl_list_to_pylist_doubles(Tcl_Obj *obj)
{
    int n;
    Tcl_Obj **elems;
    if (Tcl_ListObjGetElements(g_interp, obj, &n, &elems) != TCL_OK)
        return raise_tcl_error("ListObjGetElements");

    PyObject *list = PyList_New(n);
    if (!list) return NULL;

    for (int i = 0; i < n; i++) {
        double v = 0.0;
        Tcl_GetDoubleFromObj(g_interp, elems[i], &v);
        PyObject *f = PyFloat_FromDouble(v);
        if (!f) { Py_DECREF(list); return NULL; }
        PyList_SET_ITEM(list, i, f);
    }
    return list;
}

/*
 * Convert a Tcl list of sub-lists into a Python list of lists of floats.
 * Used by formant() and pitch() with non-AMDF methods.
 */
static PyObject *
tcl_list_to_pylist_of_double_lists(Tcl_Obj *obj)
{
    int n;
    Tcl_Obj **rows;
    if (Tcl_ListObjGetElements(g_interp, obj, &n, &rows) != TCL_OK)
        return raise_tcl_error("ListObjGetElements (outer)");

    PyObject *outer = PyList_New(n);
    if (!outer) return NULL;

    for (int i = 0; i < n; i++) {
        PyObject *inner = tcl_list_to_pylist_doubles(rows[i]);
        if (!inner) { Py_DECREF(outer); return NULL; }
        PyList_SET_ITEM(outer, i, inner);
    }
    return outer;
}

/*
 * Convert a Tcl string list result into a Python list of str.
 */
static PyObject *
tcl_list_result_to_pylist_str(const char *context)
{
    Tcl_Obj *result = Tcl_GetObjResult(g_interp);
    int n;
    Tcl_Obj **elems;
    if (Tcl_ListObjGetElements(g_interp, result, &n, &elems) != TCL_OK)
        return raise_tcl_error(context);

    PyObject *list = PyList_New(n);
    if (!list) return NULL;
    for (int i = 0; i < n; i++) {
        PyObject *s = PyUnicode_FromString(Tcl_GetString(elems[i]));
        if (!s) { Py_DECREF(list); return NULL; }
        PyList_SET_ITEM(list, i, s);
    }
    return list;
}

/*
 * Append keyword arguments to a Tcl command string as -key value pairs.
 * Returns the new offset into cmd, or -1 on truncation/error.
 */
static int
append_kwargs_as_tcl_opts(char *cmd, int off, int bufsz, PyObject *kwds)
{
    if (!kwds) return off;
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    while (PyDict_Next(kwds, &pos, &key, &value)) {
        const char *k = PyUnicode_AsUTF8(key);
        if (!k) return -1;
        PyObject *vs = PyObject_Str(value);
        if (!vs) return -1;
        const char *v = PyUnicode_AsUTF8(vs);
        if (!v) { Py_DECREF(vs); return -1; }
        int n = snprintf(cmd + off, (size_t)(bufsz - off), " -%s %s", k, v);
        Py_DECREF(vs);
        if (n < 0 || n >= bufsz - off) return -1;
        off += n;
    }
    return off;
}

/* -----------------------------------------------------------------------
 * Filter Python type
 * ----------------------------------------------------------------------- */

typedef struct {
    PyObject_HEAD
    char name[32];
} PyFilterObject;

static void
PyFilter_dealloc(PyFilterObject *self)
{
    if (self->name[0] && g_interp) {
        char cmd[64];
        snprintf(cmd, sizeof(cmd), "%s destroy", self->name);
        Tcl_Eval(g_interp, cmd);
        self->name[0] = '\0';
    }
    Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject *
PyFilter_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    const char *filttype;
    if (!PyArg_ParseTuple(args, "s", &filttype)) return NULL;

    PyFilterObject *self = (PyFilterObject *)type->tp_alloc(type, 0);
    if (!self) return NULL;

    snprintf(self->name, sizeof(self->name), "_pyflt%d", ++g_filter_counter);

    /* Build: snack::filter <type> <options...> */
    char cmd[512];
    int off = snprintf(cmd, sizeof(cmd), "snack::filter %s", filttype);
    if (kwds) {
        off = append_kwargs_as_tcl_opts(cmd, off, (int)sizeof(cmd), kwds);
        if (off < 0) {
            Py_DECREF(self);
            PyErr_SetString(PyExc_RuntimeError, "Filter options too long");
            return NULL;
        }
    }

    if (Tcl_Eval(g_interp, cmd) != TCL_OK) {
        Py_DECREF(self);
        return raise_tcl_error("Filter creation");
    }

    /* Tcl returns the filter object name; store it. */
    const char *returned_name = Tcl_GetStringResult(g_interp);
    if (returned_name && returned_name[0])
        snprintf(self->name, sizeof(self->name), "%s", returned_name);

    return (PyObject *)self;
}

/* configure(*args) */
static PyObject *
PyFilter_configure(PyFilterObject *self, PyObject *args)
{
    /* Build objv manually to handle arbitrary positional args. */
    Py_ssize_t nargs = PyTuple_GET_SIZE(args);
    Tcl_Obj **objv   = (Tcl_Obj **)PyMem_Malloc((size_t)(nargs + 2) * sizeof(Tcl_Obj *));
    if (!objv) return PyErr_NoMemory();

    objv[0] = Tcl_NewStringObj(self->name, -1);
    Tcl_IncrRefCount(objv[0]);
    objv[1] = Tcl_NewStringObj("configure", -1);
    Tcl_IncrRefCount(objv[1]);
    for (Py_ssize_t i = 0; i < nargs; i++) {
        PyObject *s = PyObject_Str(PyTuple_GET_ITEM(args, i));
        if (!s) {
            for (Py_ssize_t j = 0; j < i + 2; j++) Tcl_DecrRefCount(objv[j]);
            PyMem_Free(objv);
            return NULL;
        }
        objv[i + 2] = Tcl_NewStringObj(PyUnicode_AsUTF8(s), -1);
        Tcl_IncrRefCount(objv[i + 2]);
        Py_DECREF(s);
    }

    int rc = Tcl_EvalObjv(g_interp, (int)(nargs + 2), objv, 0);
    for (Py_ssize_t i = 0; i < nargs + 2; i++) Tcl_DecrRefCount(objv[i]);
    PyMem_Free(objv);

    if (rc != TCL_OK)
        return raise_tcl_error("Filter.configure");
    Py_RETURN_NONE;
}

/* destroy() */
static PyObject *
PyFilter_destroy(PyFilterObject *self, PyObject *Py_UNUSED(ignored))
{
    if (self->name[0] && g_interp) {
        char cmd[64];
        snprintf(cmd, sizeof(cmd), "%s destroy", self->name);
        Tcl_Eval(g_interp, cmd);
        self->name[0] = '\0';
    }
    Py_RETURN_NONE;
}

static PyObject *
PyFilter_repr(PyFilterObject *self)
{
    return PyUnicode_FromFormat("Filter(%s)", self->name);
}

static PyMethodDef PyFilter_methods[] = {
    {"configure", (PyCFunction)PyFilter_configure, METH_VARARGS,
     "configure(*args)\n\nReconfigure the filter coefficients or parameters."},
    {"destroy",   (PyCFunction)PyFilter_destroy,   METH_NOARGS,
     "destroy()\n\nFree the Snack filter object."},
    {NULL, NULL, 0, NULL}
};

static PyTypeObject PyFilterType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name      = "_snack.Filter",
    .tp_basicsize = sizeof(PyFilterObject),
    .tp_dealloc   = (destructor)PyFilter_dealloc,
    .tp_repr      = (reprfunc)PyFilter_repr,
    .tp_flags     = Py_TPFLAGS_DEFAULT,
    .tp_doc       =
        "Filter(type, **options)\n\n"
        "Snack digital filter object.  ``type`` is a Snack filter type string\n"
        "such as 'echo', 'iir', 'map', or any type registered by an extension.\n"
        "Use Sound.apply_filter() to apply the filter to a sound.",
    .tp_methods   = PyFilter_methods,
    .tp_new       = PyFilter_new,
};

/* -----------------------------------------------------------------------
 * Sound Python type
 * ----------------------------------------------------------------------- */

typedef struct {
    PyObject_HEAD
    Sound *snd;
    char   name[32];
} PySoundObject;

static void
PySound_dealloc(PySoundObject *self)
{
    if (self->snd && g_interp) {
        char cmd[64];
        snprintf(cmd, sizeof(cmd), "%s destroy", self->name);
        Tcl_Eval(g_interp, cmd);
        self->snd = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject *
PySound_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    int         rate     = 16000;
    int         channels = 1;
    const char *encoding = "Lin16";

    static char *kwlist[] = {"sample_rate", "channels", "encoding", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|iis", kwlist,
                                     &rate, &channels, &encoding))
        return NULL;

    PySoundObject *self = (PySoundObject *)type->tp_alloc(type, 0);
    if (!self) return NULL;

    snprintf(self->name, sizeof(self->name), "_pysnd%d", ++g_sound_counter);

    char cmd[256];
    snprintf(cmd, sizeof(cmd),
             "sound %s -rate %d -channels %d -encoding %s",
             self->name, rate, channels, encoding);

    if (Tcl_Eval(g_interp, cmd) != TCL_OK) {
        Py_DECREF(self);
        return raise_tcl_error("Sound creation");
    }

    self->snd = Snack_GetSound(g_interp, self->name);
    if (!self->snd) {
        char destroy[64];
        snprintf(destroy, sizeof(destroy), "%s destroy", self->name);
        Tcl_Eval(g_interp, destroy);
        Py_DECREF(self);
        PyErr_SetString(PyExc_RuntimeError, "Failed to retrieve sound object");
        return NULL;
    }

    return (PyObject *)self;
}

/* read(filename, start=0, end=-1) */
static PyObject *
PySound_read(PySoundObject *self, PyObject *args, PyObject *kwds)
{
    const char *filename;
    int         start = 0, end = -1;
    static char *kwlist[] = {"filename", "start", "end", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|ii", kwlist,
                                     &filename, &start, &end))
        return NULL;

    char cmd[2048];
    if (end >= 0)
        snprintf(cmd, sizeof(cmd),
                 "%s read {%s} -start %d -end %d", self->name, filename, start, end);
    else
        snprintf(cmd, sizeof(cmd), "%s read {%s}", self->name, filename);

    if (Tcl_Eval(g_interp, cmd) != TCL_OK)
        return raise_tcl_error("Sound.read");
    Py_RETURN_NONE;
}

/* write(filename, format=None, start=0, end=-1) */
static PyObject *
PySound_write(PySoundObject *self, PyObject *args, PyObject *kwds)
{
    const char *filename;
    const char *fmt   = NULL;
    int         start = 0, end = -1;
    static char *kwlist[] = {"filename", "format", "start", "end", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|sii", kwlist,
                                     &filename, &fmt, &start, &end))
        return NULL;

    char cmd[2048];
    int  off = snprintf(cmd, sizeof(cmd), "%s write {%s}", self->name, filename);
    if (fmt)
        off += snprintf(cmd + off, sizeof(cmd) - (size_t)off,
                        " -fileformat %s", fmt);
    if (end >= 0)
        snprintf(cmd + off, sizeof(cmd) - (size_t)off,
                 " -start %d -end %d", start, end);

    if (Tcl_Eval(g_interp, cmd) != TCL_OK)
        return raise_tcl_error("Sound.write");
    Py_RETURN_NONE;
}

/* play(blocking=True) */
static PyObject *
PySound_play(PySoundObject *self, PyObject *args, PyObject *kwds)
{
    int blocking = 1;
    static char *kwlist[] = {"blocking", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|p", kwlist, &blocking))
        return NULL;

    char cmd[64];
    snprintf(cmd, sizeof(cmd), "%s play%s",
             self->name, blocking ? " -blocking 1" : "");

    if (Tcl_Eval(g_interp, cmd) != TCL_OK)
        return raise_tcl_error("Sound.play");
    Py_RETURN_NONE;
}

/* record(**kw) */
static PyObject *
PySound_record(PySoundObject *self, PyObject *args, PyObject *kwds)
{
    char cmd[256];
    int  off = snprintf(cmd, sizeof(cmd), "%s record", self->name);
    off = append_kwargs_as_tcl_opts(cmd, off, (int)sizeof(cmd), kwds);
    if (off < 0) {
        PyErr_SetString(PyExc_RuntimeError, "record options too long");
        return NULL;
    }
    if (Tcl_Eval(g_interp, cmd) != TCL_OK)
        return raise_tcl_error("Sound.record");
    Py_RETURN_NONE;
}

/* stop() */
static PyObject *
PySound_stop(PySoundObject *self, PyObject *Py_UNUSED(ignored))
{
    char cmd[64];
    snprintf(cmd, sizeof(cmd), "%s stop", self->name);
    Tcl_Eval(g_interp, cmd);
    Py_RETURN_NONE;
}

/* pause() */
static PyObject *
PySound_pause(PySoundObject *self, PyObject *Py_UNUSED(ignored))
{
    char cmd[64];
    snprintf(cmd, sizeof(cmd), "%s pause", self->name);
    Tcl_Eval(g_interp, cmd);
    Py_RETURN_NONE;
}

/* crop(start, end) */
static PyObject *
PySound_crop(PySoundObject *self, PyObject *args)
{
    int start, end;
    if (!PyArg_ParseTuple(args, "ii", &start, &end)) return NULL;
    char cmd[128];
    snprintf(cmd, sizeof(cmd), "%s crop %d %d", self->name, start, end);
    if (Tcl_Eval(g_interp, cmd) != TCL_OK)
        return raise_tcl_error("Sound.crop");
    Py_RETURN_NONE;
}

/* cut(start, end) */
static PyObject *
PySound_cut(PySoundObject *self, PyObject *args)
{
    int start, end;
    if (!PyArg_ParseTuple(args, "ii", &start, &end)) return NULL;
    char cmd[128];
    snprintf(cmd, sizeof(cmd), "%s cut %d %d", self->name, start, end);
    if (Tcl_Eval(g_interp, cmd) != TCL_OK)
        return raise_tcl_error("Sound.cut");
    Py_RETURN_NONE;
}

/* reverse() */
static PyObject *
PySound_reverse(PySoundObject *self, PyObject *Py_UNUSED(ignored))
{
    char cmd[64];
    snprintf(cmd, sizeof(cmd), "%s reverse", self->name);
    if (Tcl_Eval(g_interp, cmd) != TCL_OK)
        return raise_tcl_error("Sound.reverse");
    Py_RETURN_NONE;
}

/* flush() */
static PyObject *
PySound_flush(PySoundObject *self, PyObject *Py_UNUSED(ignored))
{
    char cmd[64];
    snprintf(cmd, sizeof(cmd), "%s flush", self->name);
    if (Tcl_Eval(g_interp, cmd) != TCL_OK)
        return raise_tcl_error("Sound.flush");
    Py_RETURN_NONE;
}

/* destroy() — explicit release without waiting for GC */
static PyObject *
PySound_destroy(PySoundObject *self, PyObject *Py_UNUSED(ignored))
{
    if (self->snd && g_interp) {
        char cmd[64];
        snprintf(cmd, sizeof(cmd), "%s destroy", self->name);
        Tcl_Eval(g_interp, cmd);
        self->snd = NULL;
    }
    Py_RETURN_NONE;
}

/* changed(flag) — notify Snack that external code modified the buffer */
static PyObject *
PySound_changed(PySoundObject *self, PyObject *args)
{
    const char *flag;
    if (!PyArg_ParseTuple(args, "s", &flag)) return NULL;
    char cmd[128];
    snprintf(cmd, sizeof(cmd), "%s changed %s", self->name, flag);
    if (Tcl_Eval(g_interp, cmd) != TCL_OK)
        return raise_tcl_error("Sound.changed");
    Py_RETURN_NONE;
}

/* swap() — swap the byte order of all samples */
static PyObject *
PySound_swap(PySoundObject *self, PyObject *Py_UNUSED(ignored))
{
    char cmd[64];
    snprintf(cmd, sizeof(cmd), "%s swap", self->name);
    if (Tcl_Eval(g_interp, cmd) != TCL_OK)
        return raise_tcl_error("Sound.swap");
    Py_RETURN_NONE;
}

/* convert(**kw) — change sample rate, encoding, or channel count in-place */
static PyObject *
PySound_convert(PySoundObject *self, PyObject *args, PyObject *kwds)
{
    char cmd[256];
    int  off = snprintf(cmd, sizeof(cmd), "%s convert", self->name);
    off = append_kwargs_as_tcl_opts(cmd, off, (int)sizeof(cmd), kwds);
    if (off < 0) {
        PyErr_SetString(PyExc_RuntimeError, "convert options too long");
        return NULL;
    }
    if (Tcl_Eval(g_interp, cmd) != TCL_OK)
        return raise_tcl_error("Sound.convert");
    Py_RETURN_NONE;
}

/* info() -> dict */
static PyObject *
PySound_info(PySoundObject *self, PyObject *Py_UNUSED(ignored))
{
    char cmd[64];
    snprintf(cmd, sizeof(cmd), "%s info", self->name);
    if (Tcl_Eval(g_interp, cmd) != TCL_OK)
        return raise_tcl_error("Sound.info");

    Tcl_Obj *result = Tcl_GetObjResult(g_interp);
    int n;
    Tcl_Obj **elems;
    if (Tcl_ListObjGetElements(g_interp, result, &n, &elems) != TCL_OK)
        return raise_tcl_error("Sound.info (parse)");
    if (n < 8) {
        PyErr_SetString(PyExc_RuntimeError, "Sound.info: unexpected result length");
        return NULL;
    }

    /* Fields: length rate max min encoding channels fileFormat headerSize */
    int    length   = 0, rate = 0, channels = 0, header_size = 0;
    double max_val  = 0.0, min_val = 0.0;
    Tcl_GetIntFromObj(g_interp, elems[0], &length);
    Tcl_GetIntFromObj(g_interp, elems[1], &rate);
    Tcl_GetDoubleFromObj(g_interp, elems[2], &max_val);
    Tcl_GetDoubleFromObj(g_interp, elems[3], &min_val);
    const char *encoding    = Tcl_GetString(elems[4]);
    Tcl_GetIntFromObj(g_interp, elems[5], &channels);
    const char *file_format = Tcl_GetString(elems[6]);
    Tcl_GetIntFromObj(g_interp, elems[7], &header_size);

    return Py_BuildValue(
        "{s:i, s:i, s:d, s:d, s:s, s:i, s:s, s:i}",
        "length",      length,
        "sample_rate", rate,
        "max",         max_val,
        "min",         min_val,
        "encoding",    encoding,
        "channels",    channels,
        "file_format", file_format,
        "header_size", header_size);
}

/* length(n=None) — get or set the number of sample frames */
static PyObject *
PySound_length(PySoundObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *n_obj = NULL;
    static char *kwlist[] = {"n", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &n_obj))
        return NULL;

    if (n_obj == NULL || n_obj == Py_None) {
        return PyLong_FromLong(Snack_GetLength(self->snd));
    } else {
        long n = PyLong_AsLong(n_obj);
        if (n == -1 && PyErr_Occurred()) return NULL;
        char cmd[128];
        snprintf(cmd, sizeof(cmd), "%s length %ld", self->name, n);
        if (Tcl_Eval(g_interp, cmd) != TCL_OK)
            return raise_tcl_error("Sound.length");
        Py_RETURN_NONE;
    }
}

/* max_sample(**kw) -> float */
static PyObject *
PySound_max_sample(PySoundObject *self, PyObject *args, PyObject *kwds)
{
    char cmd[256];
    int  off = snprintf(cmd, sizeof(cmd), "%s max", self->name);
    off = append_kwargs_as_tcl_opts(cmd, off, (int)sizeof(cmd), kwds);
    if (off < 0) {
        PyErr_SetString(PyExc_RuntimeError, "max options too long");
        return NULL;
    }
    if (Tcl_Eval(g_interp, cmd) != TCL_OK)
        return raise_tcl_error("Sound.max_sample");
    double v = 0.0;
    Tcl_GetDoubleFromObj(g_interp, Tcl_GetObjResult(g_interp), &v);
    return PyFloat_FromDouble(v);
}

/* min_sample(**kw) -> float */
static PyObject *
PySound_min_sample(PySoundObject *self, PyObject *args, PyObject *kwds)
{
    char cmd[256];
    int  off = snprintf(cmd, sizeof(cmd), "%s min", self->name);
    off = append_kwargs_as_tcl_opts(cmd, off, (int)sizeof(cmd), kwds);
    if (off < 0) {
        PyErr_SetString(PyExc_RuntimeError, "min options too long");
        return NULL;
    }
    if (Tcl_Eval(g_interp, cmd) != TCL_OK)
        return raise_tcl_error("Sound.min_sample");
    double v = 0.0;
    Tcl_GetDoubleFromObj(g_interp, Tcl_GetObjResult(g_interp), &v);
    return PyFloat_FromDouble(v);
}

/*
 * data(data=None) — get or set raw PCM bytes.
 * With no argument: returns the sound's raw sample data as bytes.
 * With a bytes argument: loads that raw data into the sound.
 */
static PyObject *
PySound_data(PySoundObject *self, PyObject *args, PyObject *kwds)
{
    Py_buffer view = {0};
    static char *kwlist[] = {"data", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|y*", kwlist, &view))
        return NULL;

    if (view.buf != NULL) {
        /* Set: pass bytes as a Tcl byte array to avoid string escaping. */
        Tcl_Obj *objv[3];
        objv[0] = Tcl_NewStringObj(self->name, -1);
        Tcl_IncrRefCount(objv[0]);
        objv[1] = Tcl_NewStringObj("data", -1);
        Tcl_IncrRefCount(objv[1]);
        objv[2] = Tcl_NewByteArrayObj((unsigned char *)view.buf, (int)view.len);
        Tcl_IncrRefCount(objv[2]);

        int rc = Tcl_EvalObjv(g_interp, 3, objv, 0);
        Tcl_DecrRefCount(objv[0]);
        Tcl_DecrRefCount(objv[1]);
        Tcl_DecrRefCount(objv[2]);
        PyBuffer_Release(&view);

        if (rc != TCL_OK)
            return raise_tcl_error("Sound.data (set)");
        Py_RETURN_NONE;
    } else {
        /* Get: retrieve raw sample bytes. */
        char cmd[64];
        snprintf(cmd, sizeof(cmd), "%s data", self->name);
        if (Tcl_Eval(g_interp, cmd) != TCL_OK)
            return raise_tcl_error("Sound.data (get)");

        Tcl_Obj *result = Tcl_GetObjResult(g_interp);
        int blen = 0;
        unsigned char *bytes = Tcl_GetByteArrayFromObj(result, &blen);
        return PyBytes_FromStringAndSize((char *)bytes, blen);
    }
}

/*
 * append_data(data) — append raw PCM bytes to the end of the sound.
 */
static PyObject *
PySound_append_data(PySoundObject *self, PyObject *args)
{
    Py_buffer view = {0};
    if (!PyArg_ParseTuple(args, "y*", &view)) return NULL;

    Tcl_Obj *objv[3];
    objv[0] = Tcl_NewStringObj(self->name, -1);
    Tcl_IncrRefCount(objv[0]);
    objv[1] = Tcl_NewStringObj("append", -1);
    Tcl_IncrRefCount(objv[1]);
    objv[2] = Tcl_NewByteArrayObj((unsigned char *)view.buf, (int)view.len);
    Tcl_IncrRefCount(objv[2]);

    int rc = Tcl_EvalObjv(g_interp, 3, objv, 0);
    Tcl_DecrRefCount(objv[0]);
    Tcl_DecrRefCount(objv[1]);
    Tcl_DecrRefCount(objv[2]);
    PyBuffer_Release(&view);

    if (rc != TCL_OK)
        return raise_tcl_error("Sound.append_data");
    Py_RETURN_NONE;
}

/* concatenate(other) — concatenate another sound onto the end of self */
static PyObject *
PySound_concatenate(PySoundObject *self, PyObject *args)
{
    PySoundObject *other;
    if (!PyArg_ParseTuple(args, "O!", &PySoundType, &other)) return NULL;
    char cmd[128];
    snprintf(cmd, sizeof(cmd), "%s concatenate %s", self->name, other->name);
    if (Tcl_Eval(g_interp, cmd) != TCL_OK)
        return raise_tcl_error("Sound.concatenate");
    Py_RETURN_NONE;
}

/* insert(other, position) — insert other sound at sample position */
static PyObject *
PySound_insert(PySoundObject *self, PyObject *args)
{
    PySoundObject *other;
    int position;
    if (!PyArg_ParseTuple(args, "O!i", &PySoundType, &other, &position))
        return NULL;
    char cmd[128];
    snprintf(cmd, sizeof(cmd), "%s insert %s %d", self->name, other->name, position);
    if (Tcl_Eval(g_interp, cmd) != TCL_OK)
        return raise_tcl_error("Sound.insert");
    Py_RETURN_NONE;
}

/* copy(other, **kw) — copy sample data from other into self */
static PyObject *
PySound_copy(PySoundObject *self, PyObject *args, PyObject *kwds)
{
    PySoundObject *other;
    if (!PyArg_ParseTuple(args, "O!", &PySoundType, &other)) return NULL;
    char cmd[256];
    int  off = snprintf(cmd, sizeof(cmd), "%s copy %s", self->name, other->name);
    off = append_kwargs_as_tcl_opts(cmd, off, (int)sizeof(cmd), kwds);
    if (off < 0) {
        PyErr_SetString(PyExc_RuntimeError, "copy options too long");
        return NULL;
    }
    if (Tcl_Eval(g_interp, cmd) != TCL_OK)
        return raise_tcl_error("Sound.copy");
    Py_RETURN_NONE;
}

/* mix(other, **kw) — mix other into self */
static PyObject *
PySound_mix(PySoundObject *self, PyObject *args, PyObject *kwds)
{
    PySoundObject *other;
    if (!PyArg_ParseTuple(args, "O!", &PySoundType, &other)) return NULL;
    char cmd[256];
    int  off = snprintf(cmd, sizeof(cmd), "%s mix %s", self->name, other->name);
    off = append_kwargs_as_tcl_opts(cmd, off, (int)sizeof(cmd), kwds);
    if (off < 0) {
        PyErr_SetString(PyExc_RuntimeError, "mix options too long");
        return NULL;
    }
    if (Tcl_Eval(g_interp, cmd) != TCL_OK)
        return raise_tcl_error("Sound.mix");
    Py_RETURN_NONE;
}

/* power(**kw) -> list of float */
static PyObject *
PySound_power(PySoundObject *self, PyObject *args, PyObject *kwds)
{
    char cmd[256];
    int  off = snprintf(cmd, sizeof(cmd), "%s power", self->name);
    off = append_kwargs_as_tcl_opts(cmd, off, (int)sizeof(cmd), kwds);
    if (off < 0) {
        PyErr_SetString(PyExc_RuntimeError, "power options too long");
        return NULL;
    }
    if (Tcl_Eval(g_interp, cmd) != TCL_OK)
        return raise_tcl_error("Sound.power");
    return tcl_list_to_pylist_doubles(Tcl_GetObjResult(g_interp));
}

/* stretch(**kw) — time-scale modification */
static PyObject *
PySound_stretch(PySoundObject *self, PyObject *args, PyObject *kwds)
{
    char cmd[256];
    int  off = snprintf(cmd, sizeof(cmd), "%s stretch", self->name);
    off = append_kwargs_as_tcl_opts(cmd, off, (int)sizeof(cmd), kwds);
    if (off < 0) {
        PyErr_SetString(PyExc_RuntimeError, "stretch options too long");
        return NULL;
    }
    if (Tcl_Eval(g_interp, cmd) != TCL_OK)
        return raise_tcl_error("Sound.stretch");
    Py_RETURN_NONE;
}

/* shape(**kw) — apply amplitude envelope */
static PyObject *
PySound_shape(PySoundObject *self, PyObject *args, PyObject *kwds)
{
    char cmd[512];
    int  off = snprintf(cmd, sizeof(cmd), "%s shape", self->name);
    off = append_kwargs_as_tcl_opts(cmd, off, (int)sizeof(cmd), kwds);
    if (off < 0) {
        PyErr_SetString(PyExc_RuntimeError, "shape options too long");
        return NULL;
    }
    if (Tcl_Eval(g_interp, cmd) != TCL_OK)
        return raise_tcl_error("Sound.shape");
    Py_RETURN_NONE;
}

/* apply_filter(filt, **kw) — apply a Filter object in-place */
static PyObject *
PySound_apply_filter(PySoundObject *self, PyObject *args, PyObject *kwds)
{
    PyFilterObject *filt;
    if (!PyArg_ParseTuple(args, "O!", &PyFilterType, &filt)) return NULL;
    char cmd[256];
    int  off = snprintf(cmd, sizeof(cmd), "%s filter %s", self->name, filt->name);
    off = append_kwargs_as_tcl_opts(cmd, off, (int)sizeof(cmd), kwds);
    if (off < 0) {
        PyErr_SetString(PyExc_RuntimeError, "filter options too long");
        return NULL;
    }
    if (Tcl_Eval(g_interp, cmd) != TCL_OK)
        return raise_tcl_error("Sound.apply_filter");
    Py_RETURN_NONE;
}

/* get_sample(channel, index) -> float */
static PyObject *
PySound_get_sample(PySoundObject *self, PyObject *args)
{
    int channel, index;
    if (!PyArg_ParseTuple(args, "ii", &channel, &index)) return NULL;

    int nlen  = Snack_GetLength(self->snd);
    int nchan = Snack_GetNumChannels(self->snd);
    if (index < 0 || index >= nlen) {
        PyErr_SetString(PyExc_IndexError, "sample index out of range");
        return NULL;
    }
    if (channel < 0 || channel >= nchan) {
        PyErr_SetString(PyExc_IndexError, "channel index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(Snack_GetSample(self->snd, channel, index));
}

/* set_sample(channel, index, value) */
static PyObject *
PySound_set_sample(PySoundObject *self, PyObject *args)
{
    int    channel, index;
    double value;
    if (!PyArg_ParseTuple(args, "iid", &channel, &index, &value)) return NULL;

    int nlen  = Snack_GetLength(self->snd);
    int nchan = Snack_GetNumChannels(self->snd);
    if (index < 0 || index >= nlen) {
        PyErr_SetString(PyExc_IndexError, "sample index out of range");
        return NULL;
    }
    if (channel < 0 || channel >= nchan) {
        PyErr_SetString(PyExc_IndexError, "channel index out of range");
        return NULL;
    }
    Snack_SetSample(self->snd, channel, index, value);
    Py_RETURN_NONE;
}

/* get_samples() -> flat list [ch0s0, ch1s0, ch0s1, ...] */
static PyObject *
PySound_get_samples(PySoundObject *self, PyObject *Py_UNUSED(ignored))
{
    int nframes = Snack_GetLength(self->snd);
    int nchan   = Snack_GetNumChannels(self->snd);

    PyObject *list = PyList_New((Py_ssize_t)nframes * nchan);
    if (!list) return NULL;

    for (int i = 0; i < nframes; i++) {
        for (int c = 0; c < nchan; c++) {
            PyObject *v = PyFloat_FromDouble(Snack_GetSample(self->snd, c, i));
            if (!v) { Py_DECREF(list); return NULL; }
            PyList_SET_ITEM(list, (Py_ssize_t)i * nchan + c, v);
        }
    }
    return list;
}

/*
 * pitch(**kwargs) — compute F0 contour.
 *
 * With the default AMDF method returns a flat list of floats.
 * With any other method (e.g. method='ESPS') returns a list of
 * [pitch_value, voiced_flag] sub-lists per frame.
 */
static PyObject *
PySound_pitch(PySoundObject *self, PyObject *args, PyObject *kwds)
{
    char cmd[512];
    int  off = snprintf(cmd, sizeof(cmd), "%s pitch", self->name);
    int  nested = 0;   /* becomes 1 for non-AMDF methods */

    if (kwds) {
        PyObject *meth_obj = PyDict_GetItemString(kwds, "method");
        if (meth_obj) {
            const char *meth = PyUnicode_AsUTF8(meth_obj);
            if (meth && strcasecmp(meth, "amdf") != 0)
                nested = 1;
        }
        off = append_kwargs_as_tcl_opts(cmd, off, (int)sizeof(cmd), kwds);
        if (off < 0) {
            PyErr_SetString(PyExc_RuntimeError, "pitch options too long");
            return NULL;
        }
    }

    if (Tcl_Eval(g_interp, cmd) != TCL_OK)
        return raise_tcl_error("Sound.pitch");

    Tcl_Obj *result = Tcl_GetObjResult(g_interp);
    return nested
        ? tcl_list_to_pylist_of_double_lists(result)
        : tcl_list_to_pylist_doubles(result);
}

/*
 * formant(**kwargs) -> list of lists of floats.
 * Each inner list contains the formant frequencies for one analysis frame.
 */
static PyObject *
PySound_formant(PySoundObject *self, PyObject *args, PyObject *kwds)
{
    char cmd[512];
    int  off = snprintf(cmd, sizeof(cmd), "%s formant", self->name);
    off = append_kwargs_as_tcl_opts(cmd, off, (int)sizeof(cmd), kwds);
    if (off < 0) {
        PyErr_SetString(PyExc_RuntimeError, "formant options too long");
        return NULL;
    }
    if (Tcl_Eval(g_interp, cmd) != TCL_OK)
        return raise_tcl_error("Sound.formant");
    return tcl_list_to_pylist_of_double_lists(Tcl_GetObjResult(g_interp));
}

/* power_spectrum(fft_length=512, start=0, end=-1) -> list of dB floats */
static PyObject *
PySound_power_spectrum(PySoundObject *self, PyObject *args, PyObject *kwds)
{
    int fft_length = 512, start = 0, end = -1;
    static char *kwlist[] = {"fft_length", "start", "end", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|iii", kwlist,
                                     &fft_length, &start, &end))
        return NULL;

    char cmd[256];
    int  off = snprintf(cmd, sizeof(cmd),
                        "%s dBPowerSpectrum -fftlength %d -start %d",
                        self->name, fft_length, start);
    if (end >= 0)
        snprintf(cmd + off, sizeof(cmd) - (size_t)off, " -end %d", end);

    if (Tcl_Eval(g_interp, cmd) != TCL_OK)
        return raise_tcl_error("Sound.power_spectrum");
    return tcl_list_to_pylist_doubles(Tcl_GetObjResult(g_interp));
}

/* Properties (read-only): length, sample_rate, channels, encoding_id */
static PyObject *
PySound_get_length(PySoundObject *self, void *Py_UNUSED(closure))
{
    return PyLong_FromLong(Snack_GetLength(self->snd));
}

static PyObject *
PySound_get_sample_rate(PySoundObject *self, void *Py_UNUSED(closure))
{
    return PyLong_FromLong(Snack_GetSampleRate(self->snd));
}

static PyObject *
PySound_get_channels(PySoundObject *self, void *Py_UNUSED(closure))
{
    return PyLong_FromLong(Snack_GetNumChannels(self->snd));
}

static PyObject *
PySound_get_encoding(PySoundObject *self, void *Py_UNUSED(closure))
{
    return PyLong_FromLong(Snack_GetSampleEncoding(self->snd));
}

static PyObject *
PySound_repr(PySoundObject *self)
{
    return PyUnicode_FromFormat(
        "Sound(sample_rate=%d, channels=%d, length=%d)",
        Snack_GetSampleRate(self->snd),
        Snack_GetNumChannels(self->snd),
        Snack_GetLength(self->snd));
}

static PyMethodDef PySound_methods[] = {
    /* I/O */
    {"read",           (PyCFunction)PySound_read,           METH_VARARGS | METH_KEYWORDS,
     "read(filename, start=0, end=-1)\n\nLoad sound data from a file."},
    {"write",          (PyCFunction)PySound_write,          METH_VARARGS | METH_KEYWORDS,
     "write(filename, format=None, start=0, end=-1)\n\nSave sound data to a file."},
    {"data",           (PyCFunction)PySound_data,           METH_VARARGS | METH_KEYWORDS,
     "data(data=None) -> bytes or None\n\n"
     "With no argument, return the raw PCM bytes of the sound.\n"
     "With a bytes argument, load that raw data into the sound."},
    {"append_data",    (PyCFunction)PySound_append_data,    METH_VARARGS,
     "append_data(data)\n\nAppend raw PCM bytes to the end of the sound."},
    /* Playback / recording */
    {"play",           (PyCFunction)PySound_play,           METH_VARARGS | METH_KEYWORDS,
     "play(blocking=True)\n\nPlay the sound."},
    {"record",         (PyCFunction)PySound_record,         METH_VARARGS | METH_KEYWORDS,
     "record(**options)\n\nStart recording from the audio device into this sound."},
    {"stop",           (PyCFunction)PySound_stop,           METH_NOARGS,
     "stop()\n\nStop playback or recording."},
    {"pause",          (PyCFunction)PySound_pause,          METH_NOARGS,
     "pause()\n\nToggle pause for the current play or record operation."},
    /* Editing */
    {"crop",           (PyCFunction)PySound_crop,           METH_VARARGS,
     "crop(start, end)\n\nRemove all samples outside [start, end]."},
    {"cut",            (PyCFunction)PySound_cut,            METH_VARARGS,
     "cut(start, end)\n\nRemove all samples in the range [start, end]."},
    {"reverse",        (PyCFunction)PySound_reverse,        METH_NOARGS,
     "reverse()\n\nReverse the sample order."},
    {"flush",          (PyCFunction)PySound_flush,          METH_NOARGS,
     "flush()\n\nRemove all audio data from the sound (set length to zero)."},
    {"destroy",        (PyCFunction)PySound_destroy,        METH_NOARGS,
     "destroy()\n\nExplicitly free the Snack sound object without waiting for GC."},
    {"changed",        (PyCFunction)PySound_changed,        METH_VARARGS,
     "changed(flag)\n\n"
     "Notify Snack that the sound buffer has been modified externally.\n"
     "flag should be 'more' (data was appended) or 'new' (data was replaced)."},
    {"concatenate",    (PyCFunction)PySound_concatenate,    METH_VARARGS,
     "concatenate(other)\n\nAppend all samples from other onto the end of self."},
    {"insert",         (PyCFunction)PySound_insert,         METH_VARARGS,
     "insert(other, position)\n\nInsert other at the given sample position."},
    {"copy",           (PyCFunction)PySound_copy,           METH_VARARGS | METH_KEYWORDS,
     "copy(other, **options)\n\nCopy sample data from other into self."},
    /* Processing */
    {"convert",        (PyCFunction)PySound_convert,        METH_VARARGS | METH_KEYWORDS,
     "convert(**options)\n\n"
     "Convert the sound in-place.  Common options:\n"
     "  rate=<int>, channels=<int>, encoding=<str>"},
    {"swap",           (PyCFunction)PySound_swap,           METH_NOARGS,
     "swap()\n\nSwap the byte order of all samples."},
    {"mix",            (PyCFunction)PySound_mix,            METH_VARARGS | METH_KEYWORDS,
     "mix(other, **options)\n\nMix sample data from other into self."},
    {"stretch",        (PyCFunction)PySound_stretch,        METH_VARARGS | METH_KEYWORDS,
     "stretch(**options)\n\nTime-scale the sound (tempo change without pitch shift)."},
    {"shape",          (PyCFunction)PySound_shape,          METH_VARARGS | METH_KEYWORDS,
     "shape(**options)\n\nApply an amplitude envelope to the sound."},
    {"apply_filter",   (PyCFunction)PySound_apply_filter,   METH_VARARGS | METH_KEYWORDS,
     "apply_filter(filt, **options)\n\nApply a Filter object to the sound in-place."},
    /* Sample access */
    {"get_sample",     (PyCFunction)PySound_get_sample,     METH_VARARGS,
     "get_sample(channel, index) -> float\n\nReturn the value of one sample."},
    {"set_sample",     (PyCFunction)PySound_set_sample,     METH_VARARGS,
     "set_sample(channel, index, value)\n\nSet the value of one sample."},
    {"get_samples",    (PyCFunction)PySound_get_samples,    METH_NOARGS,
     "get_samples() -> list\n\n"
     "Return all samples as a flat list [ch0s0, ch1s0, ch0s1, ...]."},
    /* Metadata */
    {"info",           (PyCFunction)PySound_info,           METH_NOARGS,
     "info() -> dict\n\n"
     "Return a dict with keys: length, sample_rate, max, min, encoding,\n"
     "channels, file_format, header_size."},
    {"length",         (PyCFunction)PySound_length,         METH_VARARGS | METH_KEYWORDS,
     "length(n=None) -> int or None\n\n"
     "Without argument: return the number of sample frames (same as .length property).\n"
     "With an integer argument: resize the sound to n frames."},
    {"max_sample",     (PyCFunction)PySound_max_sample,     METH_VARARGS | METH_KEYWORDS,
     "max_sample(**options) -> float\n\nReturn the largest positive sample value."},
    {"min_sample",     (PyCFunction)PySound_min_sample,     METH_VARARGS | METH_KEYWORDS,
     "min_sample(**options) -> float\n\nReturn the largest negative sample value."},
    /* Analysis */
    {"pitch",          (PyCFunction)PySound_pitch,          METH_VARARGS | METH_KEYWORDS,
     "pitch(**options) -> list\n\n"
     "Compute the F0 (pitch) contour.\n"
     "Default method (AMDF) returns a flat list of floats.\n"
     "Other methods (e.g. method='ESPS') return a list of\n"
     "[pitch_value, voiced_flag] sub-lists."},
    {"formant",        (PyCFunction)PySound_formant,        METH_VARARGS | METH_KEYWORDS,
     "formant(**options) -> list of lists\n\n"
     "Compute formant trajectories.  Returns a list of frames;\n"
     "each frame is a list of formant frequency values."},
    {"power_spectrum", (PyCFunction)PySound_power_spectrum, METH_VARARGS | METH_KEYWORDS,
     "power_spectrum(fft_length=512, start=0, end=-1) -> list of floats\n\n"
     "Compute the log FFT power spectrum in dB."},
    {"power",          (PyCFunction)PySound_power,          METH_VARARGS | METH_KEYWORDS,
     "power(**options) -> list of floats\n\nCompute RMS power per frame."},
    {NULL, NULL, 0, NULL}
};

static PyGetSetDef PySound_getsetters[] = {
    {"length",      (getter)PySound_get_length,      NULL,
     "Number of sample frames (read-only; use length() method to resize)", NULL},
    {"sample_rate", (getter)PySound_get_sample_rate, NULL, "Sample rate in Hz",        NULL},
    {"channels",    (getter)PySound_get_channels,    NULL, "Number of channels",       NULL},
    {"encoding",    (getter)PySound_get_encoding,    NULL,
     "Sample encoding as an integer constant (LIN16, ALAW, etc.)",         NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyTypeObject PySoundType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name      = "_snack.Sound",
    .tp_basicsize = sizeof(PySoundObject),
    .tp_dealloc   = (destructor)PySound_dealloc,
    .tp_repr      = (reprfunc)PySound_repr,
    .tp_flags     = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_doc       =
        "Sound(sample_rate=16000, channels=1, encoding='Lin16')\n\n"
        "Snack sound object.  Backed directly by the Snack C library;\n"
        "no Tcl or Tk APIs are visible to callers of this type.\n\n"
        "encoding can be any string accepted by Snack: 'Lin16', 'Lin8',\n"
        "'Alaw', 'Mulaw', 'Lin32', 'Float', 'Double', etc.",
    .tp_methods   = PySound_methods,
    .tp_getset    = PySound_getsetters,
    .tp_new       = PySound_new,
};

/* -----------------------------------------------------------------------
 * Module-level audio / mixer functions
 * ----------------------------------------------------------------------- */

static PyObject *
mod_get_output_devices(PyObject *Py_UNUSED(self), PyObject *Py_UNUSED(args))
{
    if (Tcl_Eval(g_interp, "snack::audio outputDevices") != TCL_OK)
        return raise_tcl_error("get_output_devices");
    return tcl_list_result_to_pylist_str("get_output_devices (parse)");
}

static PyObject *
mod_get_input_devices(PyObject *Py_UNUSED(self), PyObject *Py_UNUSED(args))
{
    if (Tcl_Eval(g_interp, "snack::audio inputDevices") != TCL_OK)
        return raise_tcl_error("get_input_devices");
    return tcl_list_result_to_pylist_str("get_input_devices (parse)");
}

/* audio_frequencies() — list of supported sample rates */
static PyObject *
mod_audio_frequencies(PyObject *Py_UNUSED(self), PyObject *Py_UNUSED(args))
{
    if (Tcl_Eval(g_interp, "snack::audio frequencies") != TCL_OK)
        return raise_tcl_error("audio_frequencies");

    Tcl_Obj *result = Tcl_GetObjResult(g_interp);
    int n;
    Tcl_Obj **elems;
    if (Tcl_ListObjGetElements(g_interp, result, &n, &elems) != TCL_OK)
        return raise_tcl_error("audio_frequencies (parse)");

    PyObject *list = PyList_New(n);
    if (!list) return NULL;
    for (int i = 0; i < n; i++) {
        int rate = 0;
        Tcl_GetIntFromObj(g_interp, elems[i], &rate);
        PyObject *r = PyLong_FromLong(rate);
        if (!r) { Py_DECREF(list); return NULL; }
        PyList_SET_ITEM(list, i, r);
    }
    return list;
}

/* audio_encodings() — list of supported encoding format names */
static PyObject *
mod_audio_encodings(PyObject *Py_UNUSED(self), PyObject *Py_UNUSED(args))
{
    if (Tcl_Eval(g_interp, "snack::audio encodings") != TCL_OK)
        return raise_tcl_error("audio_encodings");
    return tcl_list_result_to_pylist_str("audio_encodings (parse)");
}

/* audio_elapsed_time() — seconds since playback/recording started */
static PyObject *
mod_audio_elapsed_time(PyObject *Py_UNUSED(self), PyObject *Py_UNUSED(args))
{
    if (Tcl_Eval(g_interp, "snack::audio elapsedTime") != TCL_OK)
        return raise_tcl_error("audio_elapsed_time");
    double t = 0.0;
    Tcl_GetDoubleFromObj(g_interp, Tcl_GetObjResult(g_interp), &t);
    return PyFloat_FromDouble(t);
}

/* audio_play_gain(gain=None) — get or set the playback gain (0-100) */
static PyObject *
mod_audio_play_gain(PyObject *Py_UNUSED(self), PyObject *args)
{
    PyObject *gain_obj = NULL;
    if (!PyArg_ParseTuple(args, "|O", &gain_obj)) return NULL;

    char cmd[64];
    if (gain_obj && gain_obj != Py_None) {
        double gain = PyFloat_AsDouble(gain_obj);
        if (gain == -1.0 && PyErr_Occurred()) return NULL;
        snprintf(cmd, sizeof(cmd), "snack::audio play_gain %g", gain);
    } else {
        snprintf(cmd, sizeof(cmd), "snack::audio play_gain");
    }
    if (Tcl_Eval(g_interp, cmd) != TCL_OK)
        return raise_tcl_error("audio_play_gain");
    double v = 0.0;
    Tcl_GetDoubleFromObj(g_interp, Tcl_GetObjResult(g_interp), &v);
    return PyFloat_FromDouble(v);
}

/* audio_record_gain(gain=None) — get or set the recording gain (0-100) */
static PyObject *
mod_audio_record_gain(PyObject *Py_UNUSED(self), PyObject *args)
{
    PyObject *gain_obj = NULL;
    if (!PyArg_ParseTuple(args, "|O", &gain_obj)) return NULL;

    char cmd[64];
    if (gain_obj && gain_obj != Py_None) {
        double gain = PyFloat_AsDouble(gain_obj);
        if (gain == -1.0 && PyErr_Occurred()) return NULL;
        snprintf(cmd, sizeof(cmd), "snack::audio record_gain %g", gain);
    } else {
        snprintf(cmd, sizeof(cmd), "snack::audio record_gain");
    }
    if (Tcl_Eval(g_interp, cmd) != TCL_OK)
        return raise_tcl_error("audio_record_gain");
    double v = 0.0;
    Tcl_GetDoubleFromObj(g_interp, Tcl_GetObjResult(g_interp), &v);
    return PyFloat_FromDouble(v);
}

/* audio_play_latency(ms=None) — get or set the play latency in milliseconds */
static PyObject *
mod_audio_play_latency(PyObject *Py_UNUSED(self), PyObject *args)
{
    PyObject *ms_obj = NULL;
    if (!PyArg_ParseTuple(args, "|O", &ms_obj)) return NULL;

    char cmd[64];
    if (ms_obj && ms_obj != Py_None) {
        long ms = PyLong_AsLong(ms_obj);
        if (ms == -1 && PyErr_Occurred()) return NULL;
        snprintf(cmd, sizeof(cmd), "snack::audio playLatency %ld", ms);
    } else {
        snprintf(cmd, sizeof(cmd), "snack::audio playLatency");
    }
    if (Tcl_Eval(g_interp, cmd) != TCL_OK)
        return raise_tcl_error("audio_play_latency");
    int v = 0;
    Tcl_GetIntFromObj(g_interp, Tcl_GetObjResult(g_interp), &v);
    return PyLong_FromLong(v);
}

/* audio_select_output(device) */
static PyObject *
mod_audio_select_output(PyObject *Py_UNUSED(self), PyObject *args)
{
    const char *device;
    if (!PyArg_ParseTuple(args, "s", &device)) return NULL;
    char cmd[256];
    snprintf(cmd, sizeof(cmd), "snack::audio selectOutput {%s}", device);
    if (Tcl_Eval(g_interp, cmd) != TCL_OK)
        return raise_tcl_error("audio_select_output");
    Py_RETURN_NONE;
}

/* audio_select_input(device) */
static PyObject *
mod_audio_select_input(PyObject *Py_UNUSED(self), PyObject *args)
{
    const char *device;
    if (!PyArg_ParseTuple(args, "s", &device)) return NULL;
    char cmd[256];
    snprintf(cmd, sizeof(cmd), "snack::audio selectInput {%s}", device);
    if (Tcl_Eval(g_interp, cmd) != TCL_OK)
        return raise_tcl_error("audio_select_input");
    Py_RETURN_NONE;
}

/* audio_play() — resume paused audio */
static PyObject *
mod_audio_play(PyObject *Py_UNUSED(self), PyObject *Py_UNUSED(args))
{
    if (Tcl_Eval(g_interp, "snack::audio play") != TCL_OK)
        return raise_tcl_error("audio_play");
    Py_RETURN_NONE;
}

/* audio_stop() — stop all audio playback on the device */
static PyObject *
mod_audio_stop(PyObject *Py_UNUSED(self), PyObject *Py_UNUSED(args))
{
    if (Tcl_Eval(g_interp, "snack::audio stop") != TCL_OK)
        return raise_tcl_error("audio_stop");
    Py_RETURN_NONE;
}

/* audio_pause() — toggle pause on all audio playback */
static PyObject *
mod_audio_pause(PyObject *Py_UNUSED(self), PyObject *Py_UNUSED(args))
{
    if (Tcl_Eval(g_interp, "snack::audio pause") != TCL_OK)
        return raise_tcl_error("audio_pause");
    Py_RETURN_NONE;
}

static PyMethodDef snack_methods[] = {
    /* Device enumeration */
    {"get_output_devices", mod_get_output_devices, METH_NOARGS,
     "get_output_devices() -> list of str\n\nReturn available audio output device names."},
    {"get_input_devices",  mod_get_input_devices,  METH_NOARGS,
     "get_input_devices() -> list of str\n\nReturn available audio input device names."},
    {"audio_select_output", mod_audio_select_output, METH_VARARGS,
     "audio_select_output(device)\n\nSelect a default audio output device by name."},
    {"audio_select_input",  mod_audio_select_input,  METH_VARARGS,
     "audio_select_input(device)\n\nSelect a default audio input device by name."},
    /* Device capabilities */
    {"audio_frequencies",  mod_audio_frequencies,  METH_NOARGS,
     "audio_frequencies() -> list of int\n\nReturn supported sample rates for the current device."},
    {"audio_encodings",    mod_audio_encodings,    METH_NOARGS,
     "audio_encodings() -> list of str\n\nReturn supported sample encodings for the current device."},
    /* Playback / recording control */
    {"audio_play",         mod_audio_play,         METH_NOARGS,
     "audio_play()\n\nResume paused audio playback on the device."},
    {"audio_stop",         mod_audio_stop,         METH_NOARGS,
     "audio_stop()\n\nStop all audio playback on the device."},
    {"audio_pause",        mod_audio_pause,        METH_NOARGS,
     "audio_pause()\n\nToggle pause for all audio playback on the device."},
    {"audio_elapsed_time", mod_audio_elapsed_time, METH_NOARGS,
     "audio_elapsed_time() -> float\n\nReturn seconds elapsed since playback/recording started."},
    /* Gain / latency */
    {"audio_play_gain",    mod_audio_play_gain,    METH_VARARGS,
     "audio_play_gain(gain=None) -> float\n\n"
     "Get or set the playback gain (0–100).  Returns the current gain."},
    {"audio_record_gain",  mod_audio_record_gain,  METH_VARARGS,
     "audio_record_gain(gain=None) -> float\n\n"
     "Get or set the recording gain (0–100).  Returns the current gain."},
    {"audio_play_latency", mod_audio_play_latency, METH_VARARGS,
     "audio_play_latency(ms=None) -> int\n\n"
     "Get or set how much audio (in milliseconds) is buffered to the device.\n"
     "Returns the current latency."},
    {NULL, NULL, 0, NULL}
};

/* -----------------------------------------------------------------------
 * Module initialisation
 * ----------------------------------------------------------------------- */

static struct PyModuleDef snackmodule = {
    PyModuleDef_HEAD_INIT,
    "_snack",
    "Snack sound toolkit – Python C extension.\n\n"
    "Uses an embedded Tcl interpreter for resource management; Tk is NOT\n"
    "required.  The Sound type provides direct C-level access to audio\n"
    "I/O, file formats, DSP analysis, and sample data.\n\n"
    "Types:\n"
    "  Sound  — core audio object\n"
    "  Filter — digital filter (apply with Sound.apply_filter())\n\n"
    "Module functions cover audio device control (get_output_devices,\n"
    "audio_frequencies, audio_play_gain, etc.).",
    -1,
    snack_methods
};

PyMODINIT_FUNC
PyInit__snack(void)
{
    /* Initialise the embedded Tcl interpreter. */
    g_interp = Tcl_CreateInterp();
    if (!g_interp) {
        PyErr_SetString(PyExc_RuntimeError, "Tcl_CreateInterp failed");
        return NULL;
    }
    if (Tcl_Init(g_interp) != TCL_OK) {
        PyErr_Format(PyExc_RuntimeError,
                     "Tcl_Init failed: %s", Tcl_GetStringResult(g_interp));
        return NULL;
    }

    /* Sound_Init is the Tcl-only (no-Tk) entry point from generic/sound.c. */
    extern int Sound_Init(Tcl_Interp *interp);
    if (Sound_Init(g_interp) != TCL_OK) {
        PyErr_Format(PyExc_RuntimeError,
                     "Sound_Init failed: %s", Tcl_GetStringResult(g_interp));
        return NULL;
    }

    if (PyType_Ready(&PySoundType) < 0) return NULL;
    if (PyType_Ready(&PyFilterType) < 0) return NULL;

    PyObject *m = PyModule_Create(&snackmodule);
    if (!m) return NULL;

    Py_INCREF(&PySoundType);
    if (PyModule_AddObject(m, "Sound", (PyObject *)&PySoundType) < 0) {
        Py_DECREF(&PySoundType);
        Py_DECREF(m);
        return NULL;
    }

    Py_INCREF(&PyFilterType);
    if (PyModule_AddObject(m, "Filter", (PyObject *)&PyFilterType) < 0) {
        Py_DECREF(&PyFilterType);
        Py_DECREF(m);
        return NULL;
    }

    /* Encoding integer constants (from jkAudIO.h) */
    PyModule_AddIntConstant(m, "LIN16",        LIN16);
    PyModule_AddIntConstant(m, "ALAW",         ALAW);
    PyModule_AddIntConstant(m, "MULAW",        MULAW);
    PyModule_AddIntConstant(m, "LIN8OFFSET",   LIN8OFFSET);
    PyModule_AddIntConstant(m, "LIN8",         LIN8);
    PyModule_AddIntConstant(m, "LIN24",        LIN24);
    PyModule_AddIntConstant(m, "LIN32",        LIN32);
    PyModule_AddIntConstant(m, "SNACK_FLOAT",  SNACK_FLOAT);
    PyModule_AddIntConstant(m, "SNACK_DOUBLE", SNACK_DOUBLE);

    return m;
}
