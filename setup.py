#! /usr/bin/python

from distutils.core import setup

setup(name='chiron-spec-tools',
      version='1.0',
      description='Some tools to view and normalize SMARTS 1.5m CHIRON data',
      author='Brandon Bell',
      author_email='bmbellat@gmail.com',
      url='git://github.com/MrFunPy/chiron-spec-tools.git',
      py_modules=['InteractiveNorm'],
      scripts=['Median-CHIRON','Inspect-CHIRON','Normalize-CHIRON']
      )
