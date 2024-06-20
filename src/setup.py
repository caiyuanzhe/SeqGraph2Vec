from setuptools import setup
from Cython.Build import cythonize

from setuptools import setup, Extension
from Cython.Build import cythonize

ext_modules = [
    Extension(
        "cython_function",
        sources=["cython_function.pyx"],
    )
]

setup(
    ext_modules=cythonize(ext_modules, annotate=True),
)