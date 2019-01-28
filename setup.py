from skbuild import setup

setup(
    name="fastscapelib-fortran",
    version="2.6.0",
    description=("A library of efficient algorithms"
                 "for landscape evolution modeling"),
    author='Jean Braun',
    license="GPLv3",
    packages=['fastscapelib_fortran'],
    cmake_args=['-DBUILD_PYTHON_MODULE=ON',
                '-DBUILD_FASTSCAPELIB_STATIC=OFF',
                '-DUSE_FLEXURE=ON'],
    cmake_languages=('C', 'Fortran'),
    cmake_minimum_required_version='3.5'
)
