//-------------------------------------------------------------------------
//
// Description:
//      make progress_display from BOOST library available
//
//-------------------------------------------------------------------------


#ifndef PROGRESS_DISPLAY_HPP
#define PROGRESS_DISPLAY_HPP


#include <boost/version.hpp>
// progress_display was moved to different header file and name space in version 1.72.0
#if (BOOST_VERSION >= 107200)
#include <boost/timer/progress_display.hpp>
using boost::timer::progress_display;
#else
#include <boost/progress.hpp>
using boost::progress_display;
#endif


#endif  // PROGRESS_DISPLAY_HPP
