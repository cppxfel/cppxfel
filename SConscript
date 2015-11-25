import os
import libtbx.load_env
Import("env_etc")

env_etc.cppxfel_dist = libtbx.env.dist_path("cppxfel")
env_etc.cppxfel_include = os.path.dirname(env_etc.cppxfel_dist)
if (not env_etc.no_boost_python and hasattr(env_etc, "boost_adaptbx_include")):
    Import("env_no_includes_boost_python_ext")
    env = env_no_includes_boost_python_ext.Clone()
    env_etc.enable_more_warnings(env=env)
    env_etc.include_registry.append(
        env=env,
        paths=[
            env_etc.libtbx_include,
            env_etc.scitbx_include,
            env_etc.cctbx_include,
            env_etc.cctbx_include,
            env_etc.ccp4io_include,
            env_etc.boost_include,
            env_etc.boost_adaptbx_include,
            env_etc.python_include,
            env_etc.dxtbx_include,
            env_etc.cppxfel_include])
    print env_etc.cppxfel_dist
    env.Append(
		LIBS=env_etc.libm + ["scitbx_boost_python", "boost_thread-mt", "boost_system-mt",
		"boost_python",
		"cctbx",
		"ccp4io"])
    env.Replace(SHCCFLAGS=['-std=c++0x', '-fPIC'])
    
if env_etc.clang_version:
  wd = ["-Wno-unused-variable"]
  env.Append(CCFLAGS=wd)

if 'BOOST_LOCATION' in os.environ:
	env.Append(LIBPATH = os.environ['BOOST_LOCATION'])
	print ("Appending directory containing boost libraries: " + os.environ['BOOST_LOCATION'])

source = [
'boost_python/cppxfel_ext.cc',
'source/AmbiguityBreaker.cpp',
'source/FileParser.cpp',
'source/FileReader.cpp',
'source/GraphDrawer.cpp',
'source/Holder.cpp',
'source/Image.cpp',
'source/IOMRefiner.cpp',
'source/InputFileParser.cpp',
'source/Logger.cpp',
'source/LoggableObject.cpp',
'source/Matrix.cpp',
'source/Miller.cpp',
'source/MtzGrouper.cpp',
'source/MtzManager.cpp',
'source/MtzManagerMinimize.cpp',
'source/MtzManagerRefine.cpp',
'source/MtzRefiner.cpp',
'source/NelderMead.cpp',
'source/Panel.cpp',
'source/PanelParser.cpp',
'source/PythonExt.cpp',
'source/ReflectionManager.cpp',
'source/Scaler.cpp',
'source/ScalingManager.cpp',
'source/Shoebox.cpp',
'source/Spot.cpp',
'source/SpotVector.cpp',
'source/StatisticsManager.cpp',
'source/Vector.cpp',
'source/Wiki.cpp',
'source/XManager.cpp',
'source/gaussianfit.cpp',
'source/lbfgs_scaling.cpp',
'source/main.cpp',
'source/misc.cpp']

env.SharedLibrary(
    target='#/lib/cppxfel_ext', 
    source=source,
    LIBS=env["LIBS"])
