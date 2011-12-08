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
// File and Version Information:
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
//
// Description:
//      Boost.Iostream filter for output indentation
//      taken from Boost mailing list:
//      http://lists.boost.org/Archives/boost/2008/02/133679.php
//      copyright Roland Schwarz
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef INDENTSTREAM_H
#define INDENTSTREAM_H


#ifndef ROOT_CINT
#include <ios>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/concepts.hpp>
#include <boost/iostreams/operations.hpp>


namespace rpwa {


	// a small helper class to get global variable behaviour for the
	// stream pword index for header only implementation the expression
	// xdent__<int>()() will always get the same unique index
	template <typename T>
	struct pwordIndex__ {
		int operator() ()
		{
			if (!initialized_) {
				index_       = std::ios::xalloc();
				initialized_ = true;
			}
			return index_;
		}

	private:

		static T    index_;
		static bool initialized_;
	};

	template <typename T> T    pwordIndex__<T>::index_;
	template <typename T> bool pwordIndex__<T>::initialized_;


	// the ctor is private, so the filter can only be created from its
	// static push function
	// this is to assert registration to the stream
	class indentFilter : public boost::iostreams::output_filter {

	public:

		template<typename sink>
		bool put(sink& dest,
		         int   c)
		{
			if (c == '\n')
				lineStart_ = true;
			else
				if (lineStart_) {
					for(int i = 0; i < indent_; ++i)
						for(int j = 0; j < width_; ++j)
							boost::iostreams::put(dest, ' ');
					lineStart_ = false;
				}
			return boost::iostreams::put(dest, c);
		}
    
		template<typename sink>
		void close(sink&)
		{
			indent_    = 0;
			lineStart_ = true;
		}
    
		void indentIn()
		{
			++indent_;
		}
    
		void indentOut()
		{
			if (indent_ > 0)
				--indent_;
		}
    
		// of course it would be more elegant to modify the
		// filtering_ostream push function instead, ...
		static void push(boost::iostreams::filtering_ostream& out,
		                 int                                  width = 4)
		{
			out.push(indentFilter(width));
			indentFilter* filter =	out.component<indentFilter>(out.size() - 1);
			out.pword(pwordIndex__<int>()()) = filter;
			filter->out_ = &out;
		}
    
		// when the filter is destroyed, it must be deregistered.
		~indentFilter()
		{
			if (out_)
				out_->pword(pwordIndex__<int>()()) = 0;
		}
    
	private:

		explicit indentFilter(int width)
			: indent_   (0),
			  lineStart_(true),
			  width_    (width),
			  out_      (0)
		{ }

		int           indent_;
		bool          lineStart_;
		int           width_;
		std::ostream* out_;

	};


	// these are the manipulators that change indentation.
	// note that this will even work when the filter_stream
	// is accessed through its basic_ostream.
	// uniqueness of pwordIndex__<int>()() guarantees correct cast
	// from void* to indentFilter* .
	template<class charT, class traits>
	inline std::basic_ostream<charT, traits>&
	indentIn(std::basic_ostream<charT, traits>& out)
	{
		indentFilter* filter = (indentFilter*)out.pword(pwordIndex__<int>()());
		if (filter) { 
			out.flush(); 
			filter->indentIn(); 
		}
		return out;
	}


	template<class charT, class traits>
	inline std::basic_ostream<charT, traits>&
	indentOut(std::basic_ostream<charT, traits>& out)
	{
		indentFilter* filter = (indentFilter*)out.pword(pwordIndex__<int>()());
		if (out) { 
			out.flush(); 
			filter->indentOut(); 
		}
		return out;
	}


	// define indentable versions of cout end cerr
	class indentStream : public boost::iostreams::filtering_ostream {
		
	public:
		
		indentStream(std::ostream& out)
			: boost::iostreams::filtering_ostream()
		{
			indentFilter::push(*this, 4);
			this->push(out);
		}
		
	};
	extern indentStream cout;
	extern indentStream cerr;


}  // namespace rpwa


#else


#include <iostream>


namespace rpwa {

	extern std::ostream cout;
	extern std::ostream cerr;

}


#endif  // ROOT_CINT


// override STL streams
using rpwa::cout;
using rpwa::cerr;


#endif  // INDENTSTREAM_H
