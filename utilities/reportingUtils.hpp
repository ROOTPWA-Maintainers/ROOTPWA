///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of rootpwa
//
//    rootpwa is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rootpwa is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rootpwa. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
//
// Description:
//      some primitive standardized streams for reporting plus some
//      stream operators for common STL classes
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef REPORTINGUTILS_HPP
#define REPORTINGUTILS_HPP


#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <limits>
#include <vector>
#include <cmath>
#include <unistd.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "RVersion.h"


namespace rpwa {


	//////////////////////////////////////////////////////////////////////////////
	// general output stream manipulator with one argument

	template <typename T>
	class omanip {

	private:

		typedef std::ostream& (*funcPointer)(std::ostream&, const T&);


	public:

		omanip(funcPointer func,
		       const T&    val)
			: _func(func),
			  _val (val )
		{	}

		friend
		std::ostream&
		operator << (std::ostream& out,
		             const omanip& manip)
		{	return manip._func(out, manip._val); }

	private:

		funcPointer _func;
		T           _val;
	};


	//////////////////////////////////////////////////////////////////////////////
	// functions for colored output on interactive terminals

	// based on code taken from cmake (www.cmake.org)
	// Source/kwsys/Terminal.c
	inline
	bool
	terminalSupportsColor()
	{
		// terminals running inside emacs do not support escape sequences
		const char* envEmacs = getenv("EMACS");
		if (envEmacs and (*envEmacs == 't'))
			return false;

		// check terminal name
		const std::string colorTermNames[] = {
	    "Eterm",
	    "ansi",
	    "color-xterm",
	    "con132x25", "con132x30", "con132x43", "con132x60",
	    "con80x25", "con80x28", "con80x30", "con80x43", "con80x50", "con80x60",
	    "cons25",
	    "console",
	    "cygwin",
	    "dtterm",
	    "eterm-color",
	    "gnome", "gnome-256color",
	    "konsole", "konsole-256color",
	    "kterm",
	    "linux",
	    "msys",
	    "linux-c",
	    "mach-color",
	    "mlterm",
	    "putty",
	    "rxvt", "rxvt-256color", "rxvt-cygwin", "rxvt-cygwin-native",
	    "rxvt-unicode", "rxvt-unicode-256color",
	    "screen", "screen-256color", "screen-256color-bce", "screen-bce", "screen-w", "screen.linux",
	    "vt100",
	    "xterm", "xterm-16color", "xterm-256color", "xterm-88color", "xterm-color", "xterm-debian"
    };
    const char* envTerm = getenv("TERM");
    if (envTerm) {
	    const std::string t = envTerm;
	    for (unsigned i = 0; i < sizeof(colorTermNames) / sizeof(std::string) ; ++i)
		    if (colorTermNames[i] == t)
			    return true;
    }
		return false;
	}
#ifndef __CINT__
	const bool isColorTerminal = terminalSupportsColor();
#endif

	inline
	bool
	streamIsNotInteractive(const int fileDescriptor)
	{
		struct stat streamStatus;
		if (fstat(fileDescriptor, &streamStatus) == 0)
			// check wether stream is a regular file
			if(streamStatus.st_mode & S_IFREG)
				return true;
		return false;
	}

	inline
	bool
	stdoutIsColorTerminal()
	{
		return isColorTerminal and (isatty(STDOUT_FILENO) == 1)
			and not streamIsNotInteractive(STDOUT_FILENO);
	}

	inline
	bool
	stderrIsColorTerminal()
	{
		return isColorTerminal and (isatty(STDERR_FILENO) == 1)
			and not streamIsNotInteractive(STDERR_FILENO);
	}


	// VT100 escape sequences
	enum vt100EscapeCodesEnum {
		NORMAL     = 0,
		BOLD       = 1,
		UNDERLINE  = 4,
		BLINK      = 5,
		INVERSE    = 7,
		FG_BLACK   = 30,
		FG_RED     = 31,
		FG_GREEN   = 32,
		FG_YELLOW  = 33,
		FG_BLUE    = 34,
		FG_MAGENTA = 35,
		FG_CYAN    = 36,
		FG_WHITE   = 37,
		BG_BLACK   = 40,
		BG_RED     = 41,
		BG_GREEN   = 42,
		BG_YELLOW  = 43,
		BG_BLUE    = 44,
		BG_MAGENTA = 45,
		BG_CYAN    = 46,
		BG_WHITE   = 47
	};

	// ostream manipulator function that inserts VT100 escape sequence into stream
	inline
	std::ostream&
	vt100SequenceFor(std::ostream&               out,
	                 const vt100EscapeCodesEnum& vt100Code)
	{
		bool isVt100 = false;
		if (out == std::cout)
			isVt100 = stdoutIsColorTerminal();
		else if (out == std::cerr)
			isVt100 = stderrIsColorTerminal();
		if (isVt100)
			out << "\33[" << vt100Code << "m";
		return out;
	}

	// ostream manipulator for VT100 codes
	inline
	omanip<vt100EscapeCodesEnum>
	setStreamTo(const vt100EscapeCodesEnum vt100Code)
	{ return omanip<vt100EscapeCodesEnum>(&vt100SequenceFor, vt100Code); }


	//////////////////////////////////////////////////////////////////////////////
	// macros for printing errors, warnings, and infos

	// cuts out block "className::methodName" from __PRETTY_FUNCTION__ output
	inline
	std::string
	getClassMethod__(std::string prettyFunction)
	{
		size_t pos = prettyFunction.find("(");
		if (pos == std::string::npos)
			return prettyFunction;           // something is not right
		prettyFunction.erase(pos);         // cut away signature
		pos = prettyFunction.rfind(" ");
		if (pos == std::string::npos)
			return prettyFunction;           // something is not right
		prettyFunction.erase(0, pos + 1);  // cut away return type
		return prettyFunction;
	}

#ifndef __CINT__
#define printErr   std::cerr << rpwa::setStreamTo(rpwa::FG_RED    ) << "!!! " << __PRETTY_FUNCTION__ << " [" << __FILE__ << ":" << __LINE__ << "]: error: "   << rpwa::setStreamTo(rpwa::NORMAL) << std::flush
#define printWarn  std::cerr << rpwa::setStreamTo(rpwa::FG_YELLOW ) << "??? " << __PRETTY_FUNCTION__ << " [" << __FILE__ << ":" << __LINE__ << "]: warning: " << rpwa::setStreamTo(rpwa::NORMAL) << std::flush
#define printSucc  std::cout << rpwa::setStreamTo(rpwa::FG_GREEN  ) << "*** " << rpwa::getClassMethod__(__PRETTY_FUNCTION__) << "(): success: " << rpwa::setStreamTo(rpwa::NORMAL) << std::flush
#define printInfo  std::cout << rpwa::setStreamTo(rpwa::BOLD      ) << ">>> " << rpwa::getClassMethod__(__PRETTY_FUNCTION__) << "(): info: "    << rpwa::setStreamTo(rpwa::NORMAL) << std::flush
#define printDebug std::cout << rpwa::setStreamTo(rpwa::FG_MAGENTA) << "+++ " << rpwa::getClassMethod__(__PRETTY_FUNCTION__) << "(): debug: "   << rpwa::setStreamTo(rpwa::NORMAL) << std::flush
#else
// rootcint crashes in dictionary generation (sigh)
#define printErr   std::cerr << std::flush
#define printWarn  std::cerr << std::flush
#define printSucc  std::cout << std::flush
#define printInfo  std::cout << std::flush
#define printDebug std::cout << std::flush
#endif


	//////////////////////////////////////////////////////////////////////////////
	// functions to print version and compilation info

	// check macro variables set by Makefile
#ifndef CMAKE_HOST_SYSTEM_NAME
#define CMAKE_HOST_SYSTEM_NAME "undefined"
#endif
#ifndef CMAKE_HOST_SYSTEM_PROCESSOR
#define CMAKE_HOST_SYSTEM_PROCESSOR "undefined"
#endif
#ifndef CMAKE_HOST_SYSTEM_VERSION
#define CMAKE_HOST_SYSTEM_VERSION "undefined"
#endif
#ifndef NMB_CPU_CORES
#define NMB_CPU_CORES "undefined"
#endif
#ifndef HOSTNAME
#define HOSTNAME "undefined"
#endif
#ifndef USER
#define USER "undefined"
#endif
#ifndef CMAKE_SOURCE_DIR
#define CMAKE_SOURCE_DIR "undefined"
#endif
#ifndef CMAKE_BUILD_TYPE
#define CMAKE_BUILD_TYPE "undefined"
#endif
#ifndef GIT_HASH
#define GIT_HASH "undefined"
#endif
#ifndef Boost_LIBRARY_VERSION
#define Boost_LIBRARY_VERSION "undefined"
#endif
#ifndef Boost_INCLUDE_DIRS
#define Boost_INCLUDE_DIRS "undefined"
#endif
#ifndef Libconfig_VERSION
#define Libconfig_VERSION "undefined"
#endif
#ifndef Libconfig_DIR
#define Libconfig_DIR "undefined"
#endif
#ifndef ROOTSYS
#define ROOTSYS "undefined"
#endif
#ifndef CUDA_VERSION
#define CUDA_VERSION "undefined"
#endif
#ifndef CUDA_LIB_DIRS
#define CUDA_LIB_DIRS "undefined"
#endif
#ifndef PYTHONLIBS_VERSION_STRING
#define PYTHONLIBS_VERSION_STRING "undefined"
#endif
#ifndef PYTHON_INCLUDE_DIRS
#define PYTHON_INCLUDE_DIRS "undefined"
#endif


	inline
	void
	printGitHash()
	{
		printInfo << "git repository hash at compile time was "
		          << "'" << GIT_HASH << "'" << std::endl;
	}


	inline
	void
	printCompilerInfo()
	{
		printInfo << "this executable was compiled" << std::endl
		          << "    on a '" << CMAKE_HOST_SYSTEM_NAME << " " << CMAKE_HOST_SYSTEM_VERSION << " "
		          << CMAKE_HOST_SYSTEM_PROCESSOR << "' machine" << std::endl
		          << "    in directory '" << USER << "@" << HOSTNAME ":"
		          << CMAKE_SOURCE_DIR << "'" << std::endl
		          << "    on " << __DATE__ << " " << __TIME__
		          << " by compiler version " << __VERSION__ << ","
		          << " build type is '" << CMAKE_BUILD_TYPE << "'" << std::endl;
	}


	inline
	void
	printLibraryInfo()
	{
		printInfo << "this project was linked against" << std::endl
		          << "    Boost version " << Boost_LIBRARY_VERSION << " in '" << Boost_INCLUDE_DIRS << "'"  << std::endl
		          << "    libConfig version " << Libconfig_VERSION << " in '" << Libconfig_DIR << "'" << std::endl
		          << "    ROOT version " << ROOT_RELEASE << " in '" << ROOTSYS << "'"  << std::endl;
#ifdef USE_CUDA
		std::cout << "    CUDA version " << CUDA_VERSION << " in '" << CUDA_LIB_DIRS << "'" << std::endl;
#endif
#ifdef USE_MPI
		std::cout << "    MPI libraries; using Boost.MPI" << std::endl;
#endif
#ifdef USE_PYTHON
		std::cout << "    Python version " << PYTHONLIBS_VERSION_STRING << " in '" << PYTHON_INCLUDE_DIRS << "'; using Boost.Python" << std::endl;
#endif
	}


	//////////////////////////////////////////////////////////////////////////////
	// output stream manipulator that prints a value with its maximum precision

	template<typename T> class maxPrecisionValue__;

	template<typename T>
	inline
	maxPrecisionValue__<T>
	maxPrecision(const T& value)
	{ return maxPrecisionValue__<T>(value); }

	// output stream manipulator that prints a value with its maximum precision
	// in addition manipulator reserves space so that values will align
	template<typename T>
	inline
	maxPrecisionValue__<T>
	maxPrecisionAlign(const T& value)
	{ return maxPrecisionValue__<T>(value, maxPrecisionValue__<T>::ALIGN); }

	// output stream manipulator that prints a value with maximum precision for double
	template<typename T>
	inline
	maxPrecisionValue__<T>
	maxPrecisionDouble(const T& value)
	{ return maxPrecisionValue__<T>(value, maxPrecisionValue__<T>::DOUBLE); }

	// general helper class that encapsulates a value of type T
	template<typename T>
	class maxPrecisionValue__ {
	public:
		enum modeEnum { PLAIN,
		                ALIGN,
		                DOUBLE};  // forces precision for double
		maxPrecisionValue__(const T&       value,
		                    const modeEnum mode = PLAIN)
			: _value(value),
			  _mode (mode)
		{ }
		std::ostream& print(std::ostream& out) const
		{
			const int nmbDigits = (_mode != DOUBLE) ? std::numeric_limits<T>::digits10 + 1
				: std::numeric_limits<double>::digits10 + 1;
			std::ostringstream s;
			s.precision(nmbDigits);
			s.setf(std::ios_base::scientific, std::ios_base::floatfield);
			s << _value;
			switch (_mode) {
			case ALIGN:
				return out << std::setw(nmbDigits + 7) << s.str();  // make space for sign, dot, and exponent
			case PLAIN: case DOUBLE: default:
				return out << s.str();
			}
		}
	private:
		const T& _value;
		modeEnum _mode;
	};

	template<typename T>
	inline
	std::ostream& operator << (std::ostream&                 out,
	                           const maxPrecisionValue__<T>& value)
	{ return value.print(out); }


	//////////////////////////////////////////////////////////////////////////////
	// indenting

	inline
	void
	indent(std::ostream&      out,
	       const unsigned int offset)
	{
		for (unsigned int i = 0; i < offset; ++i)
			out << " ";
	}


	//////////////////////////////////////////////////////////////////////////////
	// simple stream operators for some STL classes

	template<typename T1, typename T2>
	inline
	std::ostream&
	operator << (std::ostream&            out,
	             const std::pair<T1, T2>& pair)
	{
		return out << "(" << pair.first << ", " << pair.second << ")";
	}


	template<typename T>
	inline
	std::ostream&
	operator << (std::ostream&         out,
	             const std::vector<T>& vec)
	{
		out << "{";
		for (unsigned int i = 0; i < (vec.size() - 1); ++i)
			out << "[" << i << "] = " << vec[i] << ", ";
		return out << "[" << vec.size() - 1 << "] = " << vec[vec.size() - 1] << "}";
	}


	//////////////////////////////////////////////////////////////////////////////
	// various stuff

	template<typename T>
	inline
	unsigned int
	nmbOfDigits(const T& val)
	{
		double logVal = 0;
		if (val > 0)
			logVal = log(val);
		if (val < 0)
			logVal = log(-val);
		return (unsigned int)(logVal / log(10)) + 1;
	}


}  // namespace rpwa


#endif  // REPORTINGUTILS_HPP
