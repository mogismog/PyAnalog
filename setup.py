def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('analog', parent_package, top_path)
    config.add_extension('fortran_analog',
                         sources=['./analog/fortran_routines.f90','./analog/fortran_analog.pyf',
                                  ])
    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
