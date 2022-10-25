from setuptools import setup

setup(name='GIL',
      version='1.0.0',
      description='GIL: A python package for designing custom indexing primers',
      author='de Boer Lab',
      python_requires='>=3.6',
      license='MIT',
      packages=['GIL'],
      entry_points={
        'console_scripts': ['GIL = GIL.GIL:main']
      },
      include_package_data=True,
)