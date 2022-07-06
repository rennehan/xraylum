from setuptools import setup, find_packages


setup(name = 'xraylum',
      version = '0.1',
      description = 'xraylum: Compute xray luminosity using AtomDB tables.',
      url = 'https://github.com/rennehan/xraylum',
      author = 'Doug Rennehan',
      author_email = 'douglas.rennehan@gmail.com',
      license = 'GPLv3',
      packages = find_packages(),
      zip_safe = False, install_requires=['numpy', 'scipy', 'h5py', 'pyatomdb', 'astropy'])
