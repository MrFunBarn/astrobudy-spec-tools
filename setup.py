#! /usr/bin/python

from distutils.core import setup

setup(name='astrobudy-spec-tools',
      version='1.0',
      description='Some tools to view and normalize specral data',
      author='Brandon Bell',
      author_email='bmbellat@gmail.com',
      url='git://github.com/MrFunBarn/astrobudy-spec-tools.git',
      py_modules=['interactivenorm'],
      scripts=['median-chiron','inspect-chiron','spec-normalize']
      )
