from distutils.core import setup
import alnclst

setup(name='alnclst',
        version=alnclst.__version__,
        description=alnclst.__doc__,
        author="Christopher Small",
        author_email="csmall@fhcrc.org",
        scripts=['alnclst.py'],
        py_modules=['alnclst'],
        requires=['scipy', "numpy", "biopython"])

