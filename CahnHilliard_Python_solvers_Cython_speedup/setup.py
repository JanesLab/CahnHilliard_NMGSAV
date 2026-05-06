from setuptools import setup
from Cython.Build import cythonize
import numpy as np

setup(
    name="CHsolvers",
    ext_modules=cythonize(
        [
            "CHsolvers/CahnHilliard_NMG.pyx",
            "CHsolvers/aux_functions_NMG.pyx",
            "CHsolvers/NMG_solver.pyx",
            "CHsolvers/relax.pyx",
            "CHsolvers/laplacian.pyx",
        ],
        language_level=3,
    ),
    include_dirs=[np.get_include()],
)
