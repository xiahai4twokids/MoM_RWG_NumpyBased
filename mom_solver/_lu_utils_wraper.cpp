# include <boost/python.hpp>
# include "_lu_utils.h"
using namespace boost::python;

BOOST_PYTHON_MODULE(_lu_utils) {
       	// 导出普通函数 def("fun_name_in_python", &fun_name_in_c);
	//  导出类及部分成员 
	class_<RankRevealingMethod_tradition_python>("RankRevealingMethod_tradition_python") 
	//  类名，默认构造函数 
		//.def("LUDec", &RankRevealingMethod_tradition_python::LUDec, "args_list"); 
		.def("LUDec", &RankRevealingMethod_tradition_python::LUDec, args("field","source","wavenumber", "rows","columns", "threshold_remain", "threshold_rank")); 
}


