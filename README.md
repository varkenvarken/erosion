erosion
=======

python scripts to simulate erosion.

In its final form it will be styled as a Blender add-on, currently it's a commandline script.

More information on http://blenderthings.blogspot.com  (articles will be tagged with an 'erosion'tag )

erosion/            is the directory containing the files necessary for the Blender addon

erosion.zip         is the zipped version of the directory mentioned above and can be used to install the addon
                    (File->User Preferences->Addons->Install from file)

Dependencies:
=============

Blender 2.69

Python 3.3.x        (Blender 2.69 comes with Python 3.3)

Numpy               (See installation instructions for Windows below)
                    http://www.numpy.org/
                    
Numexpr             (Optional but *highly* recommended. See installation instructions for Windows below)
                    http://code.google.com/p/numexpr/

psutil              (Optional. See installation instructions for Windows below)
                    http://code.google.com/p/psutil/
                    
Installation of dependencies on Windows
=======================================

Blender comes bundled with its own Python distribution and this means that you have to install any packages this addon depends upon in this distribution.

For unix version this might be a simple as downloading the package from its source an running 'python setup --install'  from the python directory with the Blender installation, but for windows this is a bit more involved.

1. make sure that you have installed Pyhon 3.3.x outside Blender as well and that the platform macthes (i.e. if you use 64bit Blender, install a 64bit Python as well). The version of Python that comes with Blender is not added to the registry so any windows installers looking for a Python 3.3 installation won't find it)

2. download the windows installer suitable for your platform (i.e. again 64 bits if that's your Blender and Python version). Using a binary installer saves you the hassle of compiling packages yourself and many friendly people already did the hard work and provided suitable version to save you the trouble. Christoph Gohlke for example is an excellent source of installers: check http://www.lfd.uci.edu/~gohlke/pythonlibs/

3. run the installer and follow instructions. The installer should find the Python version you installed outside Blender in the registry

4. now copy the installed package from the external Python's lib to Blender's lib. For example on my system I have the site packages for the external Python in E:\Python33\Lib\site-packages and the site packages of my latest Blender version live in E:\blender-2.69-fd0b104-win64\blender-2.69-fd0b104-win64\2.69\python\lib\site-packages so after installation I simply copy the whole directory of a newly installed package across.

5. Test if it works: Open a Python console in Blender and type for example 'import numpy'

6. Done, after installing the erosion addon it should work. 
