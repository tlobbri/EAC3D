
%module EAC3D_module

%typemap(in) double (*func)(double) {
    $1 = (double (*)(double xx))PyLong_AsVoidPtr($input);;
}

%typemap(in) double (*func_J_relax)(double, double) {
    $1 = (double (*)(double t, double tprime))PyLong_AsVoidPtr($input);;
}
%typemap(in) double (*bnd)(double, double , double, double) {
    $1 = (double (*)(double y, double z, double t,  double val))PyLong_AsVoidPtr($input);;
}
%typemap(in) double (*func_restrain)(double, double,double) {
    $1 = (double (*)(double x, double y, double z))PyLong_AsVoidPtr($input);;
}

%typemap(in) double (*init)(double, double,double) {
    $1 = (double (*)(double x, double y, double z))PyLong_AsVoidPtr($input);;
}
%typemap(in) double (*DiffusionCoef)(double) {
    $1 = (double (*)(double x))PyLong_AsVoidPtr($input);;
}

%{

#include "EAC3D.h"

%}

%include "std_string.i" 

%include "EAC3D.h"

%extend EAC3D {
    
%pythoncode %{
    def set_T_BNDXS_py(self, bnd, bnd_type):
        import ctypes
        py_callback_type = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double)
        bnd_py = py_callback_type(bnd)
        bnd_ptr = ctypes.cast(bnd_py, ctypes.c_void_p).value
        return EAC3D.set_T_BNDXS(self, bnd_ptr, bnd_type)

    def set_T_BNDXF_py(self, bnd, bnd_type):
        import ctypes
        py_callback_type = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double)
        bnd_py = py_callback_type(bnd)
        bnd_ptr = ctypes.cast(bnd_py, ctypes.c_void_p).value
        return EAC3D.set_T_BNDXF(self, bnd_ptr, bnd_type)

    def set_T_BNDYS_py(self, bnd, bnd_type):
        import ctypes
        py_callback_type = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double)
        bnd_py = py_callback_type(bnd)
        bnd_ptr = ctypes.cast(bnd_py, ctypes.c_void_p).value
        return EAC3D.set_T_BNDYS(self, bnd_ptr, bnd_type)

    def set_T_BNDYF_py(self, bnd, bnd_type):
        import ctypes
        py_callback_type = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double)
        bnd_py = py_callback_type(bnd)
        bnd_ptr = ctypes.cast(bnd_py, ctypes.c_void_p).value
        return EAC3D.set_T_BNDYF(self, bnd_ptr, bnd_type)

    def set_T_BNDZS_py(self, bnd, bnd_type):
        import ctypes
        py_callback_type = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double)
        bnd_py = py_callback_type(bnd)
        bnd_ptr = ctypes.cast(bnd_py, ctypes.c_void_p).value
        return EAC3D.set_T_BNDZS(self, bnd_ptr, bnd_type)

    def set_T_BNDZF_py(self, bnd, bnd_type):
        import ctypes
        py_callback_type = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double)
        bnd_py = py_callback_type(bnd)
        bnd_ptr = ctypes.cast(bnd_py, ctypes.c_void_p).value
        return EAC3D.set_T_BNDZF(self, bnd_ptr, bnd_type)

    def set_T_lambda_py(self, DCoef):
        import ctypes
        py_callback_type = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
        DCoef_py = py_callback_type(DCoef)
        DCoef_ptr = ctypes.cast(DCoef_py, ctypes.c_void_p).value
        return EAC3D.set_T_lambda(self, DCoef_ptr)

    def set_H_BNDXS_py(self, bnd, bnd_type):
        import ctypes
        py_callback_type = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double)
        bnd_py = py_callback_type(bnd)
        bnd_ptr = ctypes.cast(bnd_py, ctypes.c_void_p).value
        return EAC3D.set_H_BNDXS(self, bnd_ptr, bnd_type)

    def set_H_BNDXF_py(self, bnd, bnd_type):
        import ctypes
        py_callback_type = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double)
        bnd_py = py_callback_type(bnd)
        bnd_ptr = ctypes.cast(bnd_py, ctypes.c_void_p).value
        return EAC3D.set_H_BNDXF(self, bnd_ptr, bnd_type)

    def set_H_BNDYS_py(self, bnd, bnd_type):
        import ctypes
        py_callback_type = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double)
        bnd_py = py_callback_type(bnd)
        bnd_ptr = ctypes.cast(bnd_py, ctypes.c_void_p).value
        return EAC3D.set_H_BNDYS(self, bnd_ptr, bnd_type)

    def set_H_BNDYF_py(self, bnd, bnd_type):
        import ctypes
        py_callback_type = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double)
        bnd_py = py_callback_type(bnd)
        bnd_ptr = ctypes.cast(bnd_py, ctypes.c_void_p).value
        return EAC3D.set_H_BNDYF(self, bnd_ptr, bnd_type)

    def set_H_BNDZS_py(self, bnd, bnd_type):
        import ctypes
        py_callback_type = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double)
        bnd_py = py_callback_type(bnd)
        bnd_ptr = ctypes.cast(bnd_py, ctypes.c_void_p).value
        return EAC3D.set_H_BNDZS(self, bnd_ptr, bnd_type)

    def set_H_BNDZF_py(self, bnd, bnd_type):
        import ctypes
        py_callback_type = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double)
        bnd_py = py_callback_type(bnd)
        bnd_ptr = ctypes.cast(bnd_py, ctypes.c_void_p).value
        return EAC3D.set_H_BNDZF(self, bnd_ptr, bnd_type)

    def set_H_lambda_py(self, DCoef):
        import ctypes
        py_callback_type = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
        DCoef_py = py_callback_type(DCoef)
        DCoef_ptr = ctypes.cast(DCoef_py, ctypes.c_void_p).value
        return EAC3D.set_H_lambda(self, DCoef_ptr)
    
    def init_field_py(self, init, bnd_id):
        import ctypes
        py_callback_type = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double,ctypes.c_double,ctypes.c_double)
        init_py = py_callback_type(init)
        init_ptr = ctypes.cast(init_py, ctypes.c_void_p).value
        return EAC3D.init_field(self, init_ptr, bnd_id)
    def set_DStrainODH_py(self, func):
        import ctypes
        py_callback_type = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
        func_py = py_callback_type(func)
        func_ptr = ctypes.cast(func_py, ctypes.c_void_p).value
        return EAC3D.set_DStrainODH(self, func_ptr)
    def set_DStrainODT_py(self, func):
        import ctypes
        py_callback_type = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
        func_py = py_callback_type(func)
        func_ptr = ctypes.cast(func_py, ctypes.c_void_p).value
        return EAC3D.set_DStrainODT(self, func_ptr)
    
    def set_Restrain_py(self, func):
        import ctypes
        py_callback_type = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double,ctypes.c_double,ctypes.c_double)
        func_py = py_callback_type(func)
        func_ptr = ctypes.cast(func_py, ctypes.c_void_p).value
        return EAC3D.set_Restrain(self, func_ptr)
    def set_Relaxation_py(self, func):
        import ctypes
        py_callback_type = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double, ctypes.c_double)
        func_py = py_callback_type(func)
        func_ptr = ctypes.cast(func_py, ctypes.c_void_p).value
        return EAC3D.set_Relaxation(self, func_ptr)
%}

};





