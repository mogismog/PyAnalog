#! /usr/bin/env python
#
# Copyright (C) 2013-2014 Francisco Alvarez <francisco.m.alvarez@gmail.com>

DESCRIPTION = "PyAnalog - Weather forecasting with analogous past dates."
DISTNAME = 'pyanalog'
MAINTAINER = 'Francisco M. Alvarez'
MAINTAINER_EMAIL = 'francisco.m.alvarez@gmail.com'
URL = 'https://github.com/mogismog/PyAnalog'
LICENSE = 'Apache v2 License'
DOWNLOAD_URL = 'https://github.com/mogismog/PyAnalog'
VERSION = '0.1 dev'


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('pyanalog', parent_package, top_path)
    config.add_extension('fortran_analog',
                         sources=['./pyanalog/analog/fortran_routines.f90', './pyanalog/analog/fortran_analog.pyf',
                         ])
    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup

    setup(maintainer=MAINTAINER,
          maintainer_email=MAINTAINER_EMAIL,
          description=DESCRIPTION,
          license=LICENSE,
          url=URL,
          version=VERSION,
          download_url=DOWNLOAD_URL,
          packages = ['pyanalog','pyanalog.analog'],
          author='Francisco Alvarez',
          author_email='francisco.m.alvarez@gmail.com',
          install_requires=["numpy"],
          **configuration(top_path='').todict())
