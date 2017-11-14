/*
 * flex_ext.cc
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../source/MtzRefiner.h"
#include "../source/parameters.h"
#include <iostream>
#include <string>
#include "../source/PythonExt.h"
#include "../source/InputFileParser.h"
#include <boost/container/vector.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "../source/MtzManager.h"
#include <dials/model/data/shoebox.h>
#include <dials/algorithms/centroid/centroid.h>

namespace cppxfel { namespace boost_python {

        using namespace boost::python;

        template<class T>
        boost::python::list vector_to_py_list(const vector<T>& v)
        {
                boost::python::object get_iter = boost::python::iterator<vector<T> >();
                boost::python::object iter = get_iter(v);
                boost::python::list l(iter);
                return l;
        }

        template<class T>
        std::vector<T> py_list_to_vector(boost::python::list theList)
        {
                std::vector<T> vec;

                for (int i = 0; i < len(theList); ++i)
    {
        vec.push_back(boost::python::extract<T>(theList[i]));
    }

                return vec;
        }

        boost::python::list getMtzs(InputFileParser parser)
        {
                vector<MtzPtr> mtzs = parser.getRefiner()->getMtzManagers();
                boost::python::list list = vector_to_py_list(mtzs);
                return list;
        }

        void cppxfelScript(std::string scriptFile)
        {
                runScriptFromPython(scriptFile);
        }

        BOOST_PYTHON_MODULE(cppxfel_ext)
        {
                def ("run", &cppxfelScript);
                def ("runCommandLineArgs", &runCommandLine);

                class_<MtzManager, MtzPtr, boost::noncopyable>("MtzManager", no_init)
                        .def("gridSearch", &MtzManager::gridSearch)
                        .def("refCorrelation", &MtzManager::getRefCorrelation)
                        .def("applyUnrefinedPartiality", &MtzManager::applyUnrefinedPartiality)
                ;

                class_<dials::model::Shoebox<float> >("DialsShoebox", no_init);

                class_<std::vector<MtzPtr> >("MtzArray")
                        .def(vector_indexing_suite<vector<MtzPtr> >());

                class_<InputFileParser>("cppParser", init<std::string>())
                        .def("parse", &InputFileParser::parseFromPython)
                        .def("refine", &InputFileParser::refine)
                        .def("mtzs", &getMtzs)
//                      .def("setPythonSelf", &InputFileParser::setPythonSelf)
                        .def("loadImage", &InputFileParser::loadDxtbxImage)
                        .def("addMatrixToLastImage", &InputFileParser::addMatrixToLastImage)
                        .def("integrate", &InputFileParser::integrate)
                ;

        }
}

}
