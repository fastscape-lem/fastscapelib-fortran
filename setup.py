from skbuild import setup

setup(
    name="fastscapelib-fortran",
    version="2.9.0dev",
    description=("A library of efficient algorithms"
                 "for landscape evolution modeling"),
    author='Jean Braun',
    license="GPLv3",
    packages=['fastscapelib_fortran'],
    package_dir={"": "src_python"},
    cmake_args=['-DBUILD_FASTSCAPELIB_STATIC=OFF',
                '-DUSE_FLEXURE=ON'],
    cmake_languages=('C', 'Fortran'),
    cmake_minimum_required_version='3.5'
)
