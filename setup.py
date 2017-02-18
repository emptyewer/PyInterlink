"""
This is a setup.py script generated by py2applet

Usage:
    python setup.py py2app --includes sip
    python setup.py py2exe --includes sip

"""
import os
import sys
import glob
from setuptools import setup
if sys.platform == 'darwin':
    import py2app
elif sys.platform == 'win32':
    import py2exe
    import py2exe.build_exe
    import py2exe.build_exe.isSystemDLL
sys.setrecursionlimit(100000)

def find_data_files(sources, targets, patterns):
    """Locates the specified data-files and returns the matches
    in a data_files compatible format.

    source is the root of the source data tree.
        Use '' or '.' for current directory.
    target is the root of the target data tree.
        Use '' or '.' for the distribution directory.
    patterns is a sequence of glob-patterns for the
        files you want to copy.
    """

    ret = {}
    for i, source in enumerate(sources):
        target = targets[i]
        if glob.has_magic(source) or glob.has_magic(target):
            raise ValueError("Magic not allowed in src, target")
        pattern = os.path.join(source, patterns[i])
        for filename in glob.glob(pattern):
            if os.path.isfile(filename):
                targetpath = os.path.join(target, os.path.relpath(filename, source))
                path = os.path.dirname(targetpath)
                ret.setdefault(path, []).append(filename)
    return sorted(ret.items())


APP = ['interlink.py']
INCLUDES = ['sip', 'PyQt4', 'operator', 'sys', 'os']
OPTIONS = {'argv_emulation': True,
           'includes': INCLUDES,
           'iconfile' : 'icons/Icon.icns',
           'plist': {'CFBundleGetInfoString': 'Interlink',
                     'CFBundleIdentifier': 'edu.uiowa.ahern.interlink',
                     'CFBundleShortVersionString': '0.1',
                     'CFBundleName': 'Interlink',
                     'CFBundleVersion': '01',
                     'NSHumanReadableCopyright': '(c) 2016 Venkatramanan Krishnamani'},
           'excludes': ['PyQt4.QtDesigner', 'PyQt4.QtNetwork', 'PyQt4.QtOpenGL', 'PyQt4.QtScript', 'PyQt4.QtSql', 'PyQt4.QtTest', 'PyQt4.QtWebKit', 'PyQt4.QtXml', 'PyQt4.phonon'],}

DATA_FILES_MAC = find_data_files(['structs', 'structs/pyteomics', 'ui/macos', 'icons/48/User_Interface'],
                                 ['structs', 'structs/pyteomics', 'ui/macos', 'icons/48/User_Interface'],
                                 ['*.py', '*.py', '*.ui', '*'])

DATA_FILES_WIN = find_data_files(['structs', 'structs/pyteomics', 'ui/win', 'icons/48/User_Interface'],
                             ['structs', 'structs/pyteomics', 'ui/win', 'icons/48/User_Interface'],
                             ['*.py', '*.py', '*.ui', '*'])
if sys.platform == 'darwin':
    setup(
        app=APP,
        name='Interlink',
        options={'py2app': OPTIONS},
        setup_requires=['py2app'],
        author='Venkatramanan Krishnamani',
        data_files=DATA_FILES_MAC,
    )
elif sys.platform == 'win32':
    origIsSystemDLL = py2exe.build_exe.isSystemDLL
    def isSystemDLL(pathname):
            if os.path.basename(pathname).lower() in ("msvcp71.dll", "dwmapi.dll", "'msvcp90.dll'"):
                    return 0
            return origIsSystemDLL(pathname)
    py2exe.build_exe.isSystemDLL = isSystemDLL
    setup(
        version='0.1',
        description='Interlink',
        author='Venkatramanan Krishnamani',
        windows=[{"script":'interlink.py',
                   # "icon_resources": [(1, "icon/Icon.ico")],
                   "dest_base":"Interlink"
                }],
        data_files=DATA_FILES_WIN,
        options={"py2exe": {'includes': INCLUDES,
                            "optimize": 2,
                            "compressed": 2,
                            "bundle_files": 1,
                            }}
    )

