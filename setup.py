import setuptools

VERSION_MAJOR = 0
VERSION_MINOR = 1
VERSION_MICRO = 1

with open("README.txt", "r") as fh:
    long_description = fh.read()

setuptools.setup(name='Pycroscopy3D',
                 version='{}.{}.{}'.format(
                     VERSION_MAJOR, VERSION_MINOR, VERSION_MICRO),
                 description='Python tools for microscopy',
                 author='Marco Pascucci',
                 long_description=long_description,
                 long_description_content_type="text/markdown",
                 author_email='marpas.paris@gmail.com',
                 url='https://github.com/mpascucci/pycroscopy',
                 packages=setuptools.find_packages(),
                 include_package_data=True,
                 install_requires=['multipagetiff',
                                   'numpy', 'matplotlib', 'tqdm', 'connected-components-3d', 'ants'],
                 classifiers=[
                     "Programming Language :: Python",
                     "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
                     "Operating System :: OS Independent",
                 ]
                 )
