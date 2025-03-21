from io import open

from setuptools import find_packages, setup

with open('pemtk/__init__.py', 'r') as f:
    for line in f:
        if line.startswith('__version__'):
            version = line.strip().split('=')[1].strip(' \'"')
            break
    else:
        version = '0.0.1'

with open('README.rst', 'r', encoding='utf-8') as f:
    readme = f.read()

REQUIRES = []

setup(
    name='PEMtk',
    version=version,
    description='Quantum Metrology with Photoelectrons platform data & analysis layer.',
    long_description=readme,
    long_description_content_type = 'text/x-rst',
    author='Paul Hockett',
    author_email='pemtoolkit@gmail.com',
    maintainer='Paul Hockett',
    maintainer_email='pemtoolkit@gmail.com',
    url='https://github.com/phockett/PEMtk',
    license='GNU3',

    keywords=[
        '',
    ],

    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python :: Implementation :: PyPy',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Chemistry'
    ],

    install_requires=REQUIRES,
    tests_require=['coverage', 'pytest'],

    packages=find_packages(),
)
