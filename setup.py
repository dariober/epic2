# require pysam is pre-installed
try:
    import pysam
except ImportError:
    raise Exception('pysam not found; please install pysam first')

import os
import sys

from distutils.core import setup

from setuptools import find_packages, Extension, Command
from Cython.Build import cythonize
macros = []

install_requires = ["scipy", "numpy", "natsort", "cython", "pysam"]


if os.getenv("TRAVIS"):
    install_requires.append("coveralls")


if sys.version_info[0] == 2:
    install_requires.append("functools32")

compile_options = ["-Ofast", "-Wall", "-std=c++11"] #, "-frename-registers", "-funroll-loops"] # , "-lgzstream", "-lz"

# bamtools_args = ["-I ./lib/bamtools-master/include/ -L ./lib/bamtools-master/src/"]
                   # -fprofile-generate
                   #"-fopenmp", "-D_GLIBCXX_PARALLEL"]

from subprocess import check_output

conda_path = check_output("which conda", shell=True).decode().strip()

conda_include = []
conda_lib = []
if conda_path:
    conda_base = conda_path.replace("bin/conda", "")
    conda_include.append(conda_base + "include/")
    conda_lib.append(conda_base + "lib/")



extensions = [Extension("epic2.src.reads_to_bins",
                        ["epic2/src/reads_to_bins.pyx", "epic2/src/gzstream.cpp"], language="c++",
                        include_dirs=conda_include,
                        library_dirs=conda_lib,
                        extra_compile_args=compile_options,
                        libraries=["z"]),
              Extension("epic2.src.SICER_stats",
                        ["epic2/src/SICER_stats.pyx"], language="c++",
                        extra_compile_args=compile_options),
              Extension("epic2.src.statistics",
                        ["epic2/src/statistics.pyx"], language="c++",
                        extra_compile_args=compile_options),
              Extension("epic2.src.find_islands",
                        ["epic2/src/find_islands.pyx"], language="c++",
                        extra_compile_args=compile_options),
              Extension("epic2.src.read_bam",
                        ["epic2/src/read_bam.pyx"], language="c++",
                        extra_compile_args=compile_options,
                        include_dirs=pysam.get_include(), #  + conda_include
                        library_dirs=conda_lib,
                        define_macros=pysam.get_defines()),
              Extension("epic2.src.genome_info",
                        ["epic2/src/genome_info.pyx"], language="c++",
                        extra_compile_args=compile_options)]

setup(
    name="epic2",
    packages=find_packages(),
    ext_modules = cythonize(extensions, annotate=True),
    scripts=["bin/epic2"],
    package_data={'epic2': ['effective_sizes/*.txt', 'chromsizes/*chromsizes'],
                  '': ['*.pyx', '*.pxd', '*.h', '*.c']},
    version="0.0.2",
    description="Ultraperformant ChIP-Seq broad peak/domain finder.",
    author="Endre Bakken Stovner",
    author_email="endrebak85@gmail.com",
    url="http://github.com/endrebak/epic2",
    keywords=["ChIP-Seq"],
    license=["MIT"],
    install_requires=install_requires,
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Development Status :: 4 - Beta",
        "Environment :: Other Environment", "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
        "Topic :: Scientific/Engineering"
    ],
    include_dirs=["."],
    long_description=
    ("Chip-Seq broad peak/domain finder based on SICER. See the url for more info."
     ))
